# Small helpers to read pieces of the genomic FASTA without loading it (a byte index, then seek).
# Source this file, then call open_fasta(path); read_region(chrom, start, end) gives the plus-strand bases
# (1-based, inclusive). Assumes 60 bases per line, which open_fasta() checks.
revcomp <- function(x) {
  vapply(strsplit(x, ""), function(b) paste(rev(chartr("ACGT", "TGCA", b)), collapse = ""), "")
}
open_fasta <- function(path) {
  hdr <- system2("grep", c("-b", "-a", shQuote("^>"), shQuote(path)), stdout = TRUE, env = "LC_ALL=C")
  byte <- as.numeric(sub(":.*", "", hdr))
  text <- sub("^[0-9]+:", "", hdr)
  fa <- list(path = path, con = file(path, "rb"), offset = setNames(byte + nchar(text) + 1, sub("^>", "", text)))
  stopifnot(nchar(readLines(path, n = 2)[2]) == 60)
  fa
}
read_region <- function(fa, chrom, start, end) {
  first <- start - 1
  last <- end - 1
  seek(fa$con, fa$offset[[chrom]] + first + first %/% 60)
  gsub("\n", "", rawToChar(readBin(fa$con, "raw", (end - start + 1) + (last %/% 60 - first %/% 60))), fixed = TRUE)
}
# Positions of every exact occurrence (either strand) of many equally long sequences in the whole genome. The genome is
# read in blocks of `block` bases (never all at once) and every window of each block is looked up in a hash table of the
# queries, so no outside program is needed. Returns chrom, start, end (1-based) and the strand on which the query reads
# as written ("+" if the genome has the query itself, "-" if it has its reverse complement).
find_in_genome <- function(fa, queries, chrom_lengths, block = 2e6) {
  k <- unique(nchar(queries))
  stopifnot(length(k) == 1)
  patterns <- unique(c(queries, revcomp(queries)))
  found <- list()
  for (chrom in names(chrom_lengths)) {
    n <- chrom_lengths[[chrom]]
    for (s in seq(1, n, by = block)) {
      e <- min(n, s + block - 1 + (k - 1))
      bases <- toupper(read_region(fa, chrom, s, e))
      starts <- seq_len(min(block, nchar(bases) - k + 1))
      windows <- substring(bases, starts, starts + k - 1)
      i <- which(windows %in% patterns)
      if (length(i)) found[[length(found) + 1]] <- data.frame(chrom = chrom, start = s + starts[i] - 1, seq = windows[i])
    }
  }
  h <- do.call(rbind, found)
  h$end <- h$start + k - 1
  h$strand <- ifelse(h$seq %in% queries, "+", "-")
  h
}
chrom_lengths_ws298 <- c(I = 15072434, II = 15279421, III = 13783801, IV = 17493829, V = 20924180, X = 17718942, MtDNA = 13794)
