# Phases 3 and 5 of the piRNA rebuild. Usage, from the project root:
#   Rscript data-raw/validate_new21U_and_type2.R genomic.fa piRNA_sequence.csv
# 1. Finds where the old "type 2" sequences (Mello lab, Gu et al. 2012) lie in the WS298 genome.
# 2. Tests which reading of the Claycomb Supp 6 `location` column (chrom:pos:strand) gives T-starting 21-mers.
# 3. Matches the type 2 sequences and the annotated piRNAs to those 466 PRG-1-IP-enriched new 21U loci.
args <- commandArgs(trailingOnly = TRUE)
source("data-raw/fasta_helpers.R")
fa <- open_fasta(args[1])
old <- read.csv(args[2], header = FALSE, col.names = c("Name", "Target_20mer", "Full_21mer"), stringsAsFactors = FALSE)
new <- read.csv("data-raw/staging/pirna_extracted.csv", stringsAsFactors = FALSE)
s6 <- read.csv("data-raw/inputs/claycomb_supp6_466_new_21U.csv", stringsAsFactors = FALSE)
say <- function(...) cat(sprintf(...), "\n")

# ---- type 2 sequences in the genome ----------------------------------------------------------------------
t2 <- old[grepl(" \\(type 2\\)$", old$Name), ]
h <- find_in_genome(fa, t2$Full_21mer, chrom_lengths_ws298)
h$Name <- setNames(t2$Name, t2$Full_21mer)[ifelse(h$strand == "+", h$seq, revcomp(h$seq))]
h$five_prime <- ifelse(h$strand == "+", h$start, h$end)           # position of the piRNA's first base
write.csv(h[, c("Name", "seq", "chrom", "start", "end", "strand", "five_prime")], "data-raw/staging/type2_genome_hits.csv",
          row.names = FALSE, quote = FALSE)
say("type 2: %d old entries; found in the genome: %d (%d hits); not found: %d", nrow(t2), length(unique(h$Name)), nrow(h),
    sum(!t2$Name %in% h$Name))

# ---- which reading of the Supp 6 location gives T-starting 21-mers? ------------------------------------------
# `location` is chrom:pos:strand. Candidate reading: pos (shifted by an offset) is the LEFTMOST genomic base of a 21-base
# interval; the piRNA is the plus-strand bases for "+", and their reverse complement for "-" (so its 5' end is the
# right end of the interval). The other candidate, "pos is the 5' end of the piRNA", is scored for comparison.
p <- do.call(rbind, strsplit(s6$location, ":"))
loc <- data.frame(chrom = p[, 1], pos = as.integer(p[, 2]), strand = p[, 3], stringsAsFactors = FALSE)
piece <- function(left, ...) {
  vapply(seq_len(nrow(loc)), function(i) {
    bases <- toupper(read_region(fa, loc$chrom[i], left[i], left[i] + 20))
    if (loc$strand[i] == "+") bases else revcomp(bases)
  }, "")
}
left_of_5prime <- ifelse(loc$strand == "+", loc$pos, loc$pos - 20)   # reading B: pos is the piRNA's 5' end
cat("\nFraction of the 466 21-mers that start with T:\n")
scores <- sapply(-2:2, function(o) mean(substr(piece(loc$pos + o), 1, 1) == "T"))
print(data.frame(`offset (pos is the leftmost genomic base)` = -2:2, fraction_T = round(scores, 3), check.names = FALSE))
say("reading B (pos is the piRNA's own 5' end): %.3f", mean(substr(piece(left_of_5prime), 1, 1) == "T"))
best <- (-2:2)[which.max(scores)]
loc$left <- loc$pos + best
loc$Full_21mer <- piece(loc$left)
loc$five_prime_chrom_pos <- ifelse(loc$strand == "+", loc$left, loc$left + 20)
say("chosen: pos + %d is the leftmost base; 21-mers starting with T: %d of %d; ACGT only: %d", best,
    sum(substr(loc$Full_21mer, 1, 1) == "T"), nrow(loc), sum(grepl("^[ACGT]{21}$", loc$Full_21mer)))

# ---- compare with annotated piRNAs and with the type 2 hits ----------------------------------------------------
new$five_prime <- ifelse(new$strand == "+", new$start, new$end)
key <- function(chrom, pos, strand) paste(chrom, pos, strand)
loc$key <- key(loc$chrom, loc$five_prime_chrom_pos, loc$strand)
loc$same_5prime_as_annotated <- loc$key %in% key(new$chrom, new$five_prime, new$strand)
loc$identical_21mer_annotated <- loc$Full_21mer %in% new$Full_21mer
overlap_annot <- function(i) {
  s <- new[new$chrom == loc$chrom[i], ]
  lo <- loc$left[i]
  hi <- lo + 20
  any(s$start <= hi & s$end >= lo)
}
loc$overlaps_annotated <- vapply(seq_len(nrow(loc)), overlap_annot, NA)
loc$same_5prime_as_type2 <- loc$key %in% key(h$chrom, h$five_prime, h$strand)
loc$identical_21mer_type2 <- loc$Full_21mer %in% t2$Full_21mer
say("\nSupp 6 loci (466): identical 21-mer to an annotated WS298 piRNA: %d; same 5' end as an annotated piRNA: %d; overlap an annotated piRNA at all: %d",
    sum(loc$identical_21mer_annotated), sum(loc$same_5prime_as_annotated), sum(loc$overlaps_annotated))
say("                    identical 21-mer to an old type 2 sequence: %d; same 5' end as a type 2 hit: %d",
    sum(loc$identical_21mer_type2), sum(loc$same_5prime_as_type2))
say("                    neither annotated-overlapping nor type 2: %d", sum(!loc$overlaps_annotated & !loc$same_5prime_as_type2))
out <- cbind(s6, loc[, c("Full_21mer", "five_prime_chrom_pos", "overlaps_annotated", "same_5prime_as_annotated", "same_5prime_as_type2", "identical_21mer_type2")])
write.csv(out, "data-raw/staging/claycomb_466_with_sequences.csv", row.names = FALSE, quote = FALSE)

# which type 2 sequences are supported by the PRG-1 IP data
h$in_supp6 <- key(h$chrom, h$five_prime, h$strand) %in% loc$key
say("\ntype 2 old entries whose genomic locus is one of the 466 PRG-1-IP-enriched loci (same 5' end): %d (of %d found in the genome)",
    length(unique(h$Name[h$in_supp6])), length(unique(h$Name)))
