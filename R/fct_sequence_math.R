#' Calculate GC Content
#'
#' Computes the global GC percentage of a DNA sequence.
#'
#' @param dna_seq Character string; the input DNA sequence.
#' @return A numeric value representing the GC fraction (e.g., 0.45 for 45%).
#' @export
calc_gc_content <- function(dna_seq) {
  if (!is.character(dna_seq) || nchar(dna_seq) == 0) return(NA_real_)
  dna_seq <- toupper(gsub("[^ATGCNatgcn]", "", dna_seq))

  g_count <- count_base(dna_seq, "G")
  c_count <- count_base(dna_seq, "C")

  return((g_count + c_count) / nchar(dna_seq))
}

#' Calculate Codon Adaptation Index (CAI)
#'
#' Calculates the CAI as the geometric mean of the relative adaptiveness
#' of all codons in a sequence. To prevent arithmetic underflow from multiplying
#' many small probabilities, this uses the identity: exp(mean(log(weights))).
#'
#' @param dna_seq Character string; the input DNA sequence.
#' @param usage_column Character string; the column name in `codon_freqs` to use for relative adaptiveness weights.
#' @return A numeric value representing the CAI score.
#' @export
calc_cai <- function(dna_seq, usage_column = "Ubiq") {
  dna_seq <- toupper(gsub("[^ATGCNatgcn]", "", dna_seq))
  if (nchar(dna_seq) < 3) return(NA_real_)

  # Drop dangling non-codon bases at the end of the sequence
  if (nchar(dna_seq) %% 3 != 0) {
    dna_seq <- substr(dna_seq, 1, nchar(dna_seq) - (nchar(dna_seq) %% 3))
  }

  codons <- substring(dna_seq, seq(1, nchar(dna_seq) - 2, 3), seq(3, nchar(dna_seq), 3))
  idx <- match(codons, toupper(WormTransgeneBuilder::codon_freqs$Codon))

  weights <- WormTransgeneBuilder::codon_freqs[[usage_column]][idx]
  weights <- weights[!is.na(weights) & weights > 0]

  if (length(weights) == 0) return(0)

  return(exp(mean(log(weights))))
}

#' Calculate GLO Fitness Score
#'
#' Evaluates a sequence against the Germline Optimized (GLO) 12-mer score table.
#' The sequence is divided into overlapping 12-base pair windows (stepping one codon at
#' a time), and the score is the sum of the empirically derived fitness values for each
#' window. A window's score is found from its bases (see `glo_window_scores()`).
#'
#' @param dna_seq Character string; the input DNA sequence.
#' @return A numeric value representing the cumulative GLO fitness score.
#' @export
calc_glo_score <- function(dna_seq) {
  dna_seq <- toupper(gsub("[^ATGCNatgcn]", "", dna_seq))
  cds_vec <- unlist(strsplit(dna_seq, ""))
  aa_len <- length(cds_vec) / 3

  if (aa_len < 4) return(0)

  num_windows <- aa_len - 3
  windows <- seq_len(floor(num_windows))
  scores <- glo_window_scores(glo_base_codes(cds_vec), (windows - 1) * 3 + 1)

  # Added one at a time, left to right, exactly as before, so the total is bit-for-bit the same.
  total_score <- 0
  for (score in scores) total_score <- total_score + score

  return(total_score)
}
