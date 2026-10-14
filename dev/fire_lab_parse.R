# Reads the Addgene GenBank file of the whole Fire lab vector kit (data-raw/fire_lab_vectors/, 286 records; sequences as deposited at
# Addgene) into a list with the elements id (Addgene record), name (vector) and seq (upper case). Used by dev/fire_lab_vector_search.R,
# dev/verify_fire_lab_sequences.R and tests/testthat/test-fire_lab_sequences.R.
read_fire_lab_vectors <- function(path = "data-raw/fire_lab_vectors/fire-lab-c-elegans-vector-kit-genbank-sequences.txt") {
  x <- readLines(path, warn = FALSE)
  starts <- grep("^LOCUS", x)
  ends <- c(starts[-1] - 1, length(x))
  one <- function(i) {
    r <- x[starts[i]:ends[i]]
    o <- grep("^ORIGIN", r); e <- grep("^//", r)
    e <- if (length(e)) e[1] else length(r) + 1
    list(id = sub("^LOCUS +(\\S+).*", "\\1", r[1]), name = sub("^SOURCE +", "", grep("^SOURCE", r, value = TRUE)[1]),
         seq = toupper(paste(gsub("[^A-Za-z]", "", r[(o + 1):(e - 1)]), collapse = "")))
  }
  v <- lapply(seq_along(starts), one)
  stats::setNames(v, vapply(v, `[[`, "", "name"))
}
reverse_complement_dna <- function(s) chartr("ACGT", "TGCA", paste(rev(strsplit(s, "")[[1]]), collapse = ""))
# Where `pattern` occurs in a circular sequence: start position (1-based, in the numbering of the record), or NA; the strand it was found on.
find_in_circular <- function(pattern, seq) {
  twice <- paste0(seq, seq)
  p <- regexpr(pattern, twice, fixed = TRUE)[1]
  if (p > 0) return(list(start = (p - 1) %% nchar(seq) + 1, strand = "+"))
  p <- regexpr(reverse_complement_dna(pattern), twice, fixed = TRUE)[1]
  if (p > 0) return(list(start = (p - 1) %% nchar(seq) + 1, strand = "-"))
  list(start = NA_integer_, strand = NA_character_)
}
