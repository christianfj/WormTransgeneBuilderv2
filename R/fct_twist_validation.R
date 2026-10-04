#' Safely Calculate GC Content
#'
#' Internal helper function to calculate the GC content of a DNA sequence
#' robustly, handling edge cases like empty strings or NULL inputs.
#'
#' @param dna_seq Character string; the DNA sequence to analyze.
#' @return A numeric value representing the GC fraction, or NA if invalid.
#' @noRd
safe_calc_gc <- function(dna_seq) {
  if (is.null(dna_seq) || nchar(dna_seq) == 0) return(NA_real_)

  # Count G and C (case-insensitive)
  g_count <- count_base(dna_seq, "G")
  c_count <- count_base(dna_seq, "C")

  return((g_count + c_count) / nchar(dna_seq))
}

#' Validate Sequence for Commercial Synthesis
#'
#' Evaluates a DNA sequence against standard commercial synthesis constraints
#' (specifically mimicking Twist Bioscience guidelines). Checks for homopolymers,
#' extreme global GC content, and extreme local GC content windows.
#'
#' @param dna_seq Character string; the DNA sequence to validate.
#' @return A list containing:
#'   \item{pass}{Logical; TRUE if all synthesis rules are satisfied.}
#'   \item{warnings}{Character vector; specific violation messages if the sequence fails.}
#' @export
validate_for_synthesis <- function(dna_seq) {
  if (is.null(dna_seq) || length(dna_seq) == 0 || dna_seq == "") {
    return(list(pass = FALSE, warnings = "Error: Input sequence is empty."))
  }

  dna_seq <- toupper(gsub("[^ATGCN]", "", dna_seq))
  report <- list(pass = TRUE, warnings = c())

  # 1. Check Homopolymers (Limit 14bp)
  for (base in c("A", "T", "G", "C")) {
    pattern <- paste0(rep(base, 14), collapse = "")
    if (isTRUE(grepl(pattern, dna_seq))) {
      report$pass <- FALSE
      report$warnings <- c(report$warnings, paste("Critical: Homopolymer of", base, ">= 14bp detected."))
    }
  }

  # 2. Global GC Content Check (Limit 25-75%)
  gc <- safe_calc_gc(dna_seq)

  if (is.na(gc)) {
    report$pass <- FALSE
    report$warnings <- c(report$warnings, "Error: Could not calculate GC content.")
  } else if (gc < 0.25 || gc > 0.75) {
    report$pass <- FALSE
    report$warnings <- c(report$warnings, paste0("Critical: Global GC content (", round(gc * 100, 1), "%) is outside the 25-75% range."))
  }

  # 3. Local GC Windows Check (50bp window, limit 80%)
  if (!is.na(gc) && nchar(dna_seq) >= 50) {
    for (i in 1:(nchar(dna_seq) - 49)) {
      window <- substr(dna_seq, i, i + 49)
      window_gc <- safe_calc_gc(window)

      if (!is.na(window_gc) && window_gc > 0.80) {
        report$pass <- FALSE
        report$warnings <- c(report$warnings, paste("Warning: High local GC (>80%) detected at position", i))
        break
      }
    }
  }

  return(report)
}
