# Phase 2 of the piRNA database rebuild (see data-raw/pirna_rebuild_brief.md).
#
# Reads the piRNA-only subset of the WormBase GFF3 (data-raw/inputs/gff3_piRNA_subset.gff3, made by streaming
# the 7 GB annotation file with grep) and cuts each piRNA out of the genomic FASTA. The FASTA (102 MB) is
# never loaded: it is read through a small index of byte positions, one 21-base piece at a time.
#
# Usage, from the project root:
#   Rscript data-raw/extract_pirna_from_gff3.R path/to/genomic.fa
#
# Output: data-raw/staging/pirna_extracted.csv with one row per piRNA:
#   Name (the 21ur name), WBGene, sequence_name, chrom, start, end, strand (1-based, inclusive),
#   Full_21mer (the piRNA as an RNA read 5' to 3', written as DNA: the genomic bases for a "+" feature, their
#   reverse complement for a "-" feature), Target_20mer (the reverse complement of Full_21mer without its first
#   base: the stretch of the sense transcript the piRNA pairs with; this is what the app compares).

fasta <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(!is.na(fasta), file.exists(fasta))
gff <- "data-raw/inputs/gff3_piRNA_subset.gff3"
out <- "data-raw/staging/pirna_extracted.csv"

revcomp <- function(x) {
  vapply(strsplit(x, ""), function(b) paste(rev(chartr("ACGT", "TGCA", b)), collapse = ""), "")
}

# ---- the annotation ----------------------------------------------------------------------------------
g <- utils::read.delim(gff, header = FALSE, comment.char = "#", stringsAsFactors = FALSE, quote = "",
                       col.names = c("chrom", "source", "type", "start", "end", "score", "strand", "phase", "attr"))
g <- g[g$type == "piRNA", ]
attr_value <- function(attr, key) {
  m <- regmatches(attr, regexpr(paste0("(^|;)", key, "=[^;]*"), attr))
  sub(paste0("^;?", key, "="), "", m)
}
stopifnot(all(g$end - g$start + 1 == 21), all(g$strand %in% c("+", "-")))
pi <- data.frame(Name = attr_value(g$attr, "locus"),
                 WBGene = sub("^Gene:", "", attr_value(g$attr, "Parent")),
                 sequence_name = attr_value(g$attr, "Name"),
                 chrom = g$chrom, start = g$start, end = g$end, strand = g$strand,
                 stringsAsFactors = FALSE)

# ---- an index of the FASTA: where the bases of each chromosome begin -----------------------------------
hdr <- system2("grep", c("-b", "-a", shQuote("^>"), shQuote(fasta)), stdout = TRUE, env = "LC_ALL=C")
byte <- as.numeric(sub(":.*", "", hdr))
name <- sub("^[0-9]+:>", "", hdr)
# (offset of the first base = byte offset of the header line + length of the header line + 1 newline)
seq_offset <- setNames(byte + nchar(sub("^[0-9]+:", "", hdr)) + 1, name)
width <- nchar(readLines(fasta, n = 2)[2])      # bases per line (the first sequence line)
stopifnot(width == 60)

con <- file(fasta, "rb")
on.exit(close(con), add = TRUE)
read_region <- function(chrom, start, end) {
  first <- start - 1
  last <- end - 1
  seek(con, seq_offset[[chrom]] + first + first %/% width)
  n_bytes <- (end - start + 1) + (last %/% width - first %/% width)
  gsub("\n", "", rawToChar(readBin(con, "raw", n_bytes)), fixed = TRUE)
}

plus_strand <- toupper(mapply(read_region, pi$chrom, pi$start, pi$end))
stopifnot(all(nchar(plus_strand) == 21))
pi$Full_21mer <- ifelse(pi$strand == "+", plus_strand, revcomp(plus_strand))
pi$Target_20mer <- revcomp(substr(pi$Full_21mer, 2, 21))
pi <- pi[order(as.integer(sub("^21ur-", "", pi$Name))), ]

utils::write.csv(pi, out, row.names = FALSE, quote = FALSE)
cat("Wrote", nrow(pi), "piRNAs to", out, "\n")
cat("ACGT only:", sum(grepl("^[ACGT]{21}$", pi$Full_21mer)), "of", nrow(pi), "\n")
cat("Start with T:", sum(substr(pi$Full_21mer, 1, 1) == "T"), "of", nrow(pi), "\n")
