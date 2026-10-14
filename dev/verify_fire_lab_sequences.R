# Compares the Fire lab sequences of the app (the "Fire_lab_introns" in Introns.csv and the "Fire lab synthetic spliced" 5' UTR in
# 5pUTRs.csv) with pPD95.75 (Addgene plasmid 1494, record "pPD95_75" of the Fire lab kit file in data-raw/fire_lab_vectors/), and with the
# pPD95.75 sequence file supplied by the authors (ApE, 14-NOV-2018; excerpt of positions 841 to 2160 in data-raw/fire_lab/pPD95.75_excerpt.txt).
# Run from the project root:   Rscript dev/verify_fire_lab_sequences.R
# Results are written up in docs/FIRE_LAB_SEQUENCES.md. (dev/fire_lab_vector_search.R searches all 286 vectors of the kit.)
suppressMessages(devtools::load_all(quiet = TRUE))
source("dev/fire_lab_parse.R")

addgene <- read_fire_lab_vectors()[["pPD95_75"]]$seq
cat("Addgene pPD95_75:", nchar(addgene), "bp (circular; the numbering below is the numbering of that record)\n\n")
gc <- function(x) 100 * mean(strsplit(toupper(x), "")[[1]] %in% c("G", "C"))
mismatches <- function(a, b) if (nchar(a) == nchar(b)) sum(strsplit(a, "")[[1]] != strsplit(b, "")[[1]]) else NA

introns <- utils::read.csv("inst/extdata/Introns.csv", stringsAsFactors = FALSE)
fire <- toupper(unlist(introns[introns$Name == "Fire_lab_introns", c("Sequence_1", "Sequence_2", "Sequence_3")]))
names(fire) <- c("Sequence_1", "Sequence_2", "Sequence_3")
utr <- toupper(utils::read.csv("inst/extdata/5pUTRs.csv", stringsAsFactors = FALSE)$Sequence[2])
bases_before_gfp_intron <- NULL

cat("Fire lab introns of the app in pPD95.75 (Addgene):\n")
gfp_start <- find_in_circular("ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTC", addgene)       # start of the gfp coding region
cat("  gfp coding region starts at", gfp_start$start, "on the", gfp_start$strand, "strand\n")
for (i in 1:3) {
  hit <- find_in_circular(fire[i], addgene)
  # phase: bases of gfp before the intron (modulo 3); the gfp lies on the + strand if gfp_start$strand is "+"
  before <- if (!is.na(hit$start) && gfp_start$strand == "+") (hit$start - gfp_start$start) %% nchar(addgene) else NA
  cat(sprintf("  %s: %d bases, GC %.1f%%; found at %s (%s strand); gfp bases before it %s (phase %s)\n", names(fire)[i], nchar(fire[i]), gc(fire[i]),
              hit$start, hit$strand, before, if (is.na(before)) NA else before %% 3))
}
cat(sprintf("  Mean GC of the three: %.1f%% (the menu says 24%%)\n\n", mean(vapply(fire, gc, 0))))

hit <- find_in_circular(utr, addgene)
cat(sprintf("5' UTR of the app: %d bases, GC %.1f%%; found in pPD95.75 at %s..%s (%s strand)\n", nchar(utr), gc(utr), hit$start, hit$start + nchar(utr) - 1, hit$strand))
intron_a <- substr(utr, 26, 67)
cat("  exon part:", substr(utr, 1, 25), "| intron A (", nchar(intron_a), "bases):", intron_a, "| exon part:", substr(utr, 68, 92), "\n\n")

# The authors' own pPD95.75 file (ApE, 2018): same sequence except for the third gfp intron?
lines <- grep("^ +[0-9]+ ", readLines("data-raw/fire_lab/pPD95.75_excerpt.txt"), value = TRUE)
excerpt <- toupper(paste(gsub("[^A-Za-z]", "", sub("^ +[0-9]+ ", "", lines)), collapse = ""))
anchor <- find_in_circular(substr(excerpt, 1, 40), addgene)
same_stretch <- substr(paste0(addgene, addgene), anchor$start, anchor$start + nchar(excerpt) - 1)
d <- which(strsplit(excerpt, "")[[1]] != strsplit(same_stretch, "")[[1]])
cat(sprintf("Excerpt of the authors' pPD95.75 file (%d bases, positions 841 to 2160 of that file) against Addgene's: %d differing bases, all inside the third gfp intron (positions %s of the file)\n",
            nchar(excerpt), length(d), paste(d + 840, collapse = " ")))
cat("  Addgene (= the app):", substr(same_stretch, min(d) - 14, max(d) + 14), "\n  authors' file       :", substr(excerpt, min(d) - 14, max(d) + 14), "\n\n")

enzymes <- utils::read.csv("inst/extdata/Enzymes.csv", stringsAsFactors = FALSE)
has_site <- function(x, enzyme) any(vapply(enzyme_search_sites(enzyme, warn = FALSE), function(s) nrow(find_enzyme_site_ranges(toupper(x), s)) > 0, TRUE))
listed <- function(x) { keep <- enzymes[enzymes$Common | enzymes$GG, "Enzyme"]; keep[vapply(keep, function(e) has_site(x, e), TRUE)] }
cat("Sites of the commonly used or Golden Gate enzymes inside:\n  5' UTR:", paste(listed(utr), collapse = " "), "\n")
for (i in 1:3) cat("  ", names(fire)[i], ":", paste(listed(fire[i]), collapse = " "), "\n")
