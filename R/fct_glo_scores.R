#' GLO Score Table and Lookup
#'
#' The Germline Optimization (GLO) score of a coding sequence is the sum, over its 12-base
#' windows, of a score that depends on the 12-mer (Dickinson et al.). There is a score for
#' every possible 12-mer of A/C/G/T: 4^12 = 16,777,216 of them.
#'
#' The scores are kept in one plain numeric vector, in alphabetical order of the 12-mers
#' (A < C < G < T). A 12-mer's position in that vector is therefore just its bases read as a
#' base-4 number (A = 0, C = 1, G = 2, T = 3), plus one: `AAAAAAAAAAAA` is entry 1,
#' `AAAAAAAAAAAC` is entry 2, and `TTTTTTTTTTTT` is the last. A lookup is a single
#' calculation and a single index, with no names to store and nothing to search.
#'
#' (The original dictionary, `glo_full`, stored the 12-mer text as the names of the vector,
#' which took about 90% of its memory and most of its load time, and needed a binary search
#' to find a window. It has been removed from `data/`; its last version is in git history,
#' commit bdee5b0. See `dev/make_glo_scores.R` for how the table was derived from it.)
#'
#' The table is read from `inst/extdata/glo_scores.rds` the first time it is needed, not when
#' the package is loaded, and then kept.
#' @noRd
NULL

# Weight of each of the 12 positions in a window's base-4 number.
glo_window_weights <- as.integer(4^(11:0))

#' The GLO Score Table
#'
#' @return A numeric vector of 4^12 scores; entry k is the score of the 12-mer whose base-4
#'   value is k - 1.
#' @noRd
glo_score_table <- local({
  table <- NULL
  function() {
    if (is.null(table)) {
      path <- system.file("extdata", "glo_scores.rds", package = "WormTransgeneBuilder")
      if (!nzchar(path)) stop("The GLO score table (inst/extdata/glo_scores.rds) was not found.", call. = FALSE)
      table <<- readRDS(path)
    }
    table
  }
})

#' Numeric Code of Each Base
#'
#' @param chars Character vector of single letters.
#' @return An integer vector: A = 0, C = 1, G = 2, T = 3, and `NA` for anything else (such as `N`).
#' @noRd
glo_base_codes <- function(chars) match(chars, c("A", "C", "G", "T")) - 1L

#' GLO Scores of Windows
#'
#' @param codes Integer vector of base codes from `glo_base_codes()`.
#' @param starts Integer vector; the position (1-based) of the first base of each 12-base window.
#' @return A numeric vector, one score per window. A window containing a base that is not
#'   A/C/G/T has no entry in the table and gets the median score, as it always has.
#' @noRd
glo_window_scores <- function(codes, starts) {
  k <- length(starts)
  if (k == 0) return(numeric(0))
  positions <- outer(starts, 0:11, "+")
  window_code <- rowSums(matrix(codes[positions], nrow = k) * rep(glo_window_weights, each = k))
  scores <- glo_score_table()[window_code + 1]
  scores[is.na(window_code)] <- WormTransgeneBuilder::glo_median_score
  scores
}
