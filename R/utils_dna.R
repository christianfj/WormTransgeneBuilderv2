#' Small DNA Helpers in Base R
#'
#' The app needs only three small things that used to come from the Biostrings (Bioconductor)
#' and stringr packages: the reverse complement of a restriction site, approximate matching of a
#' piRNA target for the viewer, and counting a base. They are short, so they are written in base
#' R. That removes Bioconductor, whose releases are tied to one R version each and are a common
#' cause of an app failing to reinstall years later, from the list of dependencies.
#'
#' Letters are the IUPAC DNA codes: A C G T and the ambiguity codes R (A or G), Y (C or T),
#' S (C or G), W (A or T), K (G or T), M (A or C), B (C, G or T), D (A, G or T), H (A, C or T),
#' V (A, C or G) and N (any base).
#' @name dna_utils
#' @noRd
NULL

# Each IUPAC letter as the set of bases it stands for, written as a 4-bit number
# (A = 1, C = 2, G = 4, T = 8). Two letters are compatible when their sets overlap.
iupac_masks <- c(A = 1L, C = 2L, G = 4L, T = 8L,
                 R = 5L, Y = 10L, S = 6L, W = 9L, K = 12L, M = 3L,
                 B = 14L, D = 13L, H = 11L, V = 7L, N = 15L)

# The complement of each IUPAC letter (the letter for the complementary bases).
iupac_complement <- c(A = "T", C = "G", G = "C", T = "A",
                      R = "Y", Y = "R", S = "S", W = "W", K = "M", M = "K",
                      B = "V", D = "H", H = "D", V = "B", N = "N")

#' Reverse Complement
#'
#' @param seq Character vector of DNA sequences (IUPAC letters, any case).
#' @return The reverse complement of each, in upper case. Stops with an error if a sequence
#'   contains a letter that is not an IUPAC DNA code.
#' @noRd
reverse_complement <- function(seq) {
  vapply(seq, function(x) {
    comp <- unname(iupac_complement[rev(strsplit(toupper(x), "")[[1]])])
    if (is.na(x) || anyNA(comp)) {
      stop("Not a DNA sequence; only the IUPAC letters A C G T R Y S W K M B D H V N are allowed: ", x, call. = FALSE)
    }
    paste(comp, collapse = "")
  }, "", USE.NAMES = FALSE)
}

#' Count One Base
#'
#' @param dna Character string.
#' @param base A single letter.
#' @return How many times the letter occurs, ignoring case.
#' @noRd
count_base <- function(dna, base) {
  nchar(dna) - nchar(gsub(base, "", toupper(dna), fixed = TRUE))
}

#' Windows of a Sequence That Match a Pattern, Allowing Mismatches
#'
#' Slides the pattern along the sequence one base at a time and reports every position where the
#' pattern and the window of the sequence differ at no more than `max_mismatch` places. A pair of
#' letters counts as matching when their IUPAC sets overlap, so an `N` in either matches anything.
#' A character that is not an IUPAC letter matches nothing. Only whole windows inside the
#' sequence are considered; a pattern hanging over an end is not a match.
#'
#' @param pattern Character string; the pattern (for example a 20-base piRNA target).
#' @param subject Character string; the sequence to search.
#' @param max_mismatch Integer; the most differing positions allowed.
#' @return An integer vector of 1-based start positions (empty if there are none).
#' @noRd
approximate_match_starts <- function(pattern, subject, max_mismatch = 0) {
  p <- unname(iupac_masks[strsplit(toupper(pattern), "")[[1]]])
  s <- unname(iupac_masks[strsplit(toupper(subject), "")[[1]]])
  p[is.na(p)] <- 0L
  s[is.na(s)] <- 0L
  len <- length(p)
  n <- length(s)
  if (len == 0 || len > n) return(integer(0))

  windows <- n - len + 1L
  mismatches <- integer(windows)
  for (j in seq_len(len)) {
    mismatches <- mismatches + (bitwAnd(p[j], s[j:(j + windows - 1L)]) == 0L)
  }
  as.integer(which(mismatches <= max_mismatch))
}
