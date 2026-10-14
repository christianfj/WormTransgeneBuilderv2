# Shared helpers for the other-species piRNA databases (see data-raw/pirna_other_species_brief.md).
revcomp <- function(x) {
  vapply(strsplit(x, ""), function(b) paste(rev(chartr("ACGT", "TGCA", b)), collapse = ""), "")
}

# Reads a whole (small) genome FASTA into a named character vector, one string per sequence (upper case).
read_fasta <- function(path) {
  lines <- readLines(path)
  header <- grep("^>", lines)
  names_ <- sub("^>([^ \t]+).*", "\\1", lines[header])
  ends <- c(header[-1] - 1, length(lines))
  seqs <- mapply(function(a, b) toupper(paste(lines[(a + 1):b], collapse = "")), header, ends)
  setNames(unname(seqs), names_)
}

# The piRNA's own 5' end and its upstream sequence, read in the piRNA's direction.
# `start`, `end` are 1-based, inclusive genome coordinates of the 21 bases; `strand` is "+" or "-".
# Returns the `n` bases upstream of the first base (5' to 3', so the last character is position -1).
upstream_of <- function(genome, chrom, start, end, strand, n = 60) {
  mapply(function(ch, s, e, sd) {
    g <- genome[[ch]]
    if (sd == "+") {
      lo <- s - n
      if (lo < 1) return(NA_character_)
      substr(g, lo, s - 1)
    } else {
      hi <- e + n
      if (hi > nchar(g)) return(NA_character_)
      revcomp(substr(g, e + 1, hi))
    }
  }, chrom, start, end, strand, USE.NAMES = FALSE)
}

# Ruby upstream motif flag: the core GTTTC lies entirely within positions -50 to -35 relative to the first base of the
# piRNA (position +1 is the first base; the motif core is about 35 to 42 bases upstream of it, Billi et al. 2013).
# `upstream` must be the 60 bases before the piRNA as returned by upstream_of().
has_ruby_core <- function(upstream) {
  window <- ifelse(is.na(upstream), NA_character_, substr(upstream, 11, 26))   # positions -50 .. -35 of a 60-base string
  grepl("GTTTC", window, fixed = TRUE) & !is.na(window)
}
