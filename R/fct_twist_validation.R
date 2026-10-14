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

# The protein encoded by ccdB (the toxin of the type II toxin-antitoxin system of the F plasmid; the full coding sequence of 101 amino acids, supplied by the authors). Twist asks that it is not included
# in a design. Because the app recodes sequences, the check is made on the protein in all six reading frames, not on the DNA.
ccdb_protein <- function() {
  dna <- paste0("atgcagtttaaggtttacacctataaaagagagagccgttatcgtctgtttgtggatgtacagagtgatattattgacacgcccgggcgacggatggtgatccccctggccagtgcacgtctgctgtcagataaagtctcccgtgaactttacccggtggtgcatatcggggatgaaagctggcgcatgatgaccaccgatatggccagtgtgccggtctccgttatcggggaagaagtggctgatctcagccaccgcgaaaatgacatcaaaaacgccattaacctgatgttctggggaatataa")
  translate_frame(toupper(dna))
}

# Translation of one reading frame (starting at base 1) with the app's codon table; a stop codon is "X", an unknown codon "?".
translate_frame <- function(dna) {
  n <- nchar(dna) %/% 3
  if (n == 0) return("")
  cf <- WormTransgeneBuilder::codon_freqs
  aa <- stats::setNames(cf$Amino, toupper(cf$Codon))[substring(dna, seq(1, 3 * n - 2, 3), seq(3, 3 * n, 3))]
  aa[is.na(aa)] <- "?"
  paste(aa, collapse = "")
}

# TRUE if a stretch of 30 amino acids of ccdB is encoded anywhere in the sequence (either strand, any frame).
contains_ccdb <- function(dna_seq, window = 30L) {
  target <- ccdb_protein()
  target <- substring(target, 1:(nchar(target) - window + 1), window:nchar(target))
  for (strand in list(dna_seq, reverse_complement(dna_seq))) {
    for (frame in 1:3) {
      if (nchar(strand) - frame + 1 < 3 * window) next
      protein <- translate_frame(substring(strand, frame))
      if (any(vapply(target, function(w) grepl(w, protein, fixed = TRUE), TRUE))) return(TRUE)
    }
  }
  FALSE
}

# Positions (1-based, of the second copy) at which a stretch of `len` bases repeats an earlier, non-overlapping stretch on the same strand.
direct_repeat_starts <- function(dna_seq, len) {
  n <- nchar(dna_seq)
  if (n < 2 * len) return(integer(0))
  kmers <- substring(dna_seq, 1:(n - len + 1), len:n)
  first <- match(kmers, kmers)
  which(seq_along(kmers) - first >= len)
}

#' Validate Sequence for Commercial Synthesis
#'
#' Checks a DNA sequence against Twist Bioscience's published gene synthesis guidelines (read on 2026-10-14; see `docs/METHODS.md` for the
#' pages). Twist scores the complexity of a design (Low, Medium or High) and decides itself; the result here separates two levels:
#' \itemize{
#'   \item **standard guidelines**: homopolymers of 14 bases or more; global GC content outside 25 to 75 percent; a 50-base window with more than
#'     80 percent GC (a softer warning than Twist's own limit of 90); a length below 300 or above 7,000 bases (and above 5,000, too long for gene
#'     fragments); a ccdB coding region (translated in all six frames, so a recoded copy is found too).
#'   \item **high complexity**: any homopolymer longer than 30 bases; a 50-base window below 10 or above 90 percent GC; a direct repeat of more than
#'     200 bases. Such a design is not eligible for the Express service and is made as a Complex Gene, which takes longer.
#' }
#' Not checked: short tandem repeats longer than 100 bases, hairpin-forming repeats longer than 100 bases, and Twist's machine-learning complexity score.
#'
#' @param dna_seq Character string; the DNA sequence to validate.
#' @return A list containing:
#'   \item{pass}{Logical; TRUE if there is no issue at either level.}
#'   \item{high_complexity}{Logical; TRUE if a high-complexity criterion is met.}
#'   \item{standard_issues}{Character vector; messages about the standard guidelines.}
#'   \item{high_complexity_issues}{Character vector; messages about high-complexity criteria.}
#'   \item{warnings}{Character vector; all the messages (standard issues first), for display.}
#' @export
validate_for_synthesis <- function(dna_seq) {
  if (is.null(dna_seq) || length(dna_seq) == 0 || dna_seq == "") {
    return(list(pass = FALSE, high_complexity = FALSE, standard_issues = "Error: Input sequence is empty.",
                high_complexity_issues = character(0), warnings = "Error: Input sequence is empty."))
  }

  dna_seq <- gsub("[^ATGCN]", "", toupper(dna_seq))   # upper case first: introns, UTRs and promoters of the final sequence are lower case
  standard <- character(0)
  high <- character(0)
  n <- nchar(dna_seq)

  # 1. Homopolymers: standard guideline >= 14 bases; more than 30 bases is high complexity.
  for (base in c("A", "T", "G", "C")) {
    runs <- gregexpr(paste0(base, "{14,}"), dna_seq)[[1]]
    if (runs[1] > 0) {
      longest <- max(attr(runs, "match.length"))
      if (longest > 30) {
        high <- c(high, paste0("Homopolymer of ", base, " of ", longest, " bp (longer than 30 bp)."))
      } else {
        standard <- c(standard, paste0("Homopolymer of ", base, " of ", longest, " bp (14 bp or more)."))
      }
    }
  }

  # 2. Global GC content (25 to 75 percent).
  gc <- safe_calc_gc(dna_seq)
  if (is.na(gc)) {
    standard <- c(standard, "Error: Could not calculate GC content.")
  } else if (gc < 0.25 || gc > 0.75) {
    standard <- c(standard, paste0("Global GC content (", round(gc * 100, 1), "%) is outside the 25-75% range."))
  }

  # 3. Local GC in 50-base windows: above 80 percent is a warning, below 10 or above 90 percent is high complexity.
  if (n >= 50) {
    is_gc <- c(0, cumsum(strsplit(dna_seq, "")[[1]] %in% c("G", "C")))
    window_gc <- (is_gc[51:(n + 1)] - is_gc[1:(n - 49)]) / 50
    low <- which(window_gc < 0.10)
    top <- which(window_gc > 0.90)
    mid <- which(window_gc > 0.80 & window_gc <= 0.90)
    if (length(low) > 0) high <- c(high, paste0("A 50 bp window with less than 10% GC (first at position ", low[1], ")."))
    if (length(top) > 0) high <- c(high, paste0("A 50 bp window with more than 90% GC (first at position ", top[1], ")."))
    if (length(mid) > 0 && length(top) == 0) standard <- c(standard, paste0("A 50 bp window with more than 80% GC (first at position ", mid[1], ")."))
  }

  # 4. Direct repeats of more than 200 bases (high complexity).
  rep_start <- direct_repeat_starts(dna_seq, 201L)
  if (length(rep_start) > 0) high <- c(high, paste0("A direct repeat of more than 200 bp (a second copy starts at position ", rep_start[1], ")."))

  # 5. Length (gene fragments 300 to 5,000 bp; clonal genes up to 7,000 bp).
  if (n < 300) {
    standard <- c(standard, paste0("Length (", n, " bp) is below the 300 bp minimum."))
  } else if (n > 7000) {
    standard <- c(standard, paste0("Length (", n, " bp) is above the 7,000 bp maximum for clonal genes."))
  } else if (n > 5000) {
    standard <- c(standard, paste0("Length (", n, " bp) is above the 5,000 bp maximum for gene fragments (clonal genes only)."))
  }

  # 6. ccdB must not be included.
  if (contains_ccdb(dna_seq)) standard <- c(standard, "Contains a ccdB coding region (Twist: do not include ccdB).")

  list(pass = length(standard) == 0 && length(high) == 0,
       high_complexity = length(high) > 0,
       standard_issues = standard,
       high_complexity_issues = high,
       warnings = c(standard, high))
}
