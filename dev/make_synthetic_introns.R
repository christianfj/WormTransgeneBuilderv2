# Makes the pool of synthetic introns, inst/extdata/Synthetic_Introns.csv (column Sequence). The app picks
# introns from this pool at random, and swaps one out when it would create a recognition site of an enzyme
# the user wants removed (insert_introns(avoid_sites =, pool =)).
# Run from the project root:   Rscript dev/make_synthetic_introns.R [number_of_introns] [seed]
# It replaces the python script used for the first pool. The same seed gives the same pool.
#
# Where the design comes from: Blumenthal & Steward, "RNA processing and gene structure", in C. elegans II
# (1997), section "Cis-splicing in worms" (https://www.ncbi.nlm.nih.gov/books/NBK20075/, see also
# https://pubmed.ncbi.nlm.nih.gov/21278743/):
#  * more than half of the C. elegans introns are shorter than 60 nt; the very short ones peak at 48 nt, with more
#    introns than expected at 44 and 52 nt; only 4 of 659 are shorter than 40 nt.
#    -> lengths 44 to 52, weighted to peak at 48 with side peaks at 44 and 52 (the weights below are an
#       approximation of that description, not counts read from the figure);
#  * 5' splice site GU..., 3' splice site consensus UUUUCAG (the U at -5 is almost invariant), no branch-site
#    consensus, but A is more frequent at -16, -17 and -18 from the 3' splice site;
#    -> every intron starts gtaagttt and ends ttttcag, and has an A at -16, -17 and -18 with probability 0.7 each;
#  * introns are A+U-rich (about 70%) compared with exons (about 54%); the short ones are richer in U.
#    -> the part between the fixed ends has 72 to 78% A+T (as in the first pool), with a little more T than A.
# Added in this version, to keep introns from being (mis)spliced or silenced and from interfering with cloning:
#  * no cryptic splice sites: gtaag only at the 5' end, tcag (so also ttcag, ttttcag) only at the 3' end;
#  * no recognition site of the commonly used and Golden Gate enzymes (Enzymes.csv, columns Common and GG), on either
#    strand, inside the intron. (Sites that appear only where an intron meets its exon are handled in the app, at
#    the time of insertion, by choosing another intron from this pool.)
#  * no 20-base window within 3 mismatches of a known piRNA (minimum Hamming distance >= 4, the app's default
#    threshold, Priyadarshini et al. 2022) of any of the four species' databases, type 2 piRNAs included.
suppressMessages(devtools::load_all(quiet = TRUE))

args <- commandArgs(trailingOnly = TRUE)
n_wanted <- if (length(args) >= 1) as.integer(args[1]) else 1000L
seed <- if (length(args) >= 2) as.integer(args[2]) else 2026L
set.seed(seed)

five_prime <- "gtaagttt"
three_prime <- "ttttcag"
lengths_allowed <- 44:52
length_weights <- c(`44` = 3, `45` = 2, `46` = 3, `47` = 4, `48` = 6, `49` = 4, `50` = 3, `51` = 2, `52` = 3)
min_piRNA_distance <- 4

# Sites to keep out: both strands of every enzyme that is "commonly used" or Golden Gate.
enzymes <- utils::read.csv(file.path("inst", "extdata", "Enzymes.csv"), stringsAsFactors = FALSE)
avoid_enzymes <- enzymes$Enzyme[enzymes$Common | enzymes$GG]
avoid_regex <- paste(vapply(enzyme_search_sites(avoid_enzymes), enzyme_site_regex, ""), collapse = "|")

# Every piRNA of every species, as encoded 20-mers (the app's engine), to compare the 20-base windows with.
all_pirnas <- unique(select_pirna_database(c("elegans", "briggsae", "remanei", "brenneri"), include_type2 = TRUE)$Target_20mer)
reference <- make_pirna_reference(all_pirnas)

draw_one <- function() {
  len <- sample(lengths_allowed, 1, prob = length_weights[as.character(lengths_allowed)])
  mid_len <- len - nchar(five_prime) - nchar(three_prime)
  at_fraction <- runif(1, 0.72, 0.78)
  at <- round(mid_len * at_fraction)
  # Within A+T a little more T than A (short C. elegans introns are U-rich); within G+C, equal.
  p <- c(a = 0.45 * at, t = 0.55 * at, c = 0.5 * (mid_len - at), g = 0.5 * (mid_len - at))
  mid <- sample(c("a", "t", "c", "g"), mid_len, replace = TRUE, prob = p)
  intron <- strsplit(paste0(five_prime, paste(mid, collapse = ""), three_prime), "")[[1]]
  # A at -16, -17, -18 from the 3' splice site (the last base of the intron is -1).
  for (k in 16:18) if (runif(1) < 0.7) intron[len - k + 1] <- "a"
  paste(intron, collapse = "")
}

passes_sequence_rules <- function(x) {
  inner <- substr(x, 2, nchar(x))                                  # everything after the first base
  !grepl("gtaag", inner, fixed = TRUE) &&                           # no second 5' splice site
    !grepl("tcag", substr(x, 1, nchar(x) - 4), fixed = TRUE) &&     # no 3' splice site before the real one
    !grepl("ttttcag", substr(x, 1, nchar(x) - 1), fixed = TRUE) &&
    !grepl(avoid_regex, toupper(x), perl = TRUE)
}

pool <- character(0)
rejected <- c(rules = 0L, piRNA = 0L)
while (length(pool) < n_wanted) {
  # A batch is drawn, screened for the sequence rules, then all windows of the survivors are compared with the piRNAs at once.
  batch <- unique(replicate(400, draw_one()))
  batch <- setdiff(batch, pool)
  ok <- vapply(batch, passes_sequence_rules, TRUE)
  rejected["rules"] <- rejected["rules"] + sum(!ok)
  batch <- batch[ok]
  if (length(batch) == 0) next
  windows <- lapply(batch, function(x) substring(toupper(x), 1:(nchar(x) - 19), 20:nchar(x)))
  row_min <- scan_kmers_against_pirnas(unlist(windows), reference)$row_min
  intron_min <- vapply(split(row_min, rep(seq_along(batch), lengths(windows))), min, 0)
  keep <- intron_min >= min_piRNA_distance
  rejected["piRNA"] <- rejected["piRNA"] + sum(!keep)
  pool <- c(pool, batch[keep])
}
pool <- pool[seq_len(n_wanted)]

stopifnot(!anyDuplicated(pool), all(grepl("^gtaagttt[acgt]+ttttcag$", pool)))
out <- file.path("inst", "extdata", "Synthetic_Introns.csv")
utils::write.csv(data.frame(Sequence = pool), out, row.names = FALSE, quote = FALSE)

at_pct <- 100 * vapply(strsplit(pool, ""), function(x) mean(x %in% c("a", "t")), 0)
cat(sprintf("Wrote %d introns to %s (seed %d)\n", length(pool), out, seed))
cat("Rejected while drawing: ", rejected[["rules"]], " by the sequence rules, ", rejected[["piRNA"]], " by the piRNA screen\n", sep = "")
cat("Lengths:\n"); print(table(nchar(pool)))
cat("A+T of the whole intron (%): "); cat(paste(names(summary(at_pct)), round(summary(at_pct), 1), collapse = "  "), "\n")
