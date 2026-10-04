#' Germline Optimize (GLO) a DNA Sequence
#'
#' Implementation of GLO algorithm developed by Daniel Dickinson (https://pubmed.ncbi.nlm.nih.gov/30109984/).
#'
#' Employs a greedy, localized perturbation algorithm to optimize a coding sequence
#' for germline expression based on a fitness score for every possible 12-mer. A 12-mer's
#' score is found by working out its position in the score table from its bases (see
#' `glo_window_scores()`), so a lookup is one calculation and one index.
#'
#' Furthermore, when evaluating a codon mutation, it does not recalculate the fitness
#' of the entire gene. It maps and recalculates only the 4 overlapping 12-mer windows
#' directly impacted by the substitution, radically improving iteration speed.
#'
#' @param dna_seq Character string; the input DNA sequence to optimize.
#' @return A character string containing the GLO-optimized DNA sequence.
#' @export
optimize_glo <- function(dna_seq) {
  current_seq_str <- toupper(dna_seq)
  cds_vec <- unlist(strsplit(current_seq_str, ""))
  aa_len <- length(cds_vec) / 3

  if (aa_len < 4) return(dna_seq)

  upper_codons <- toupper(WormTransgeneBuilder::codon_freqs$Codon)
  codon_to_aa <- setNames(WormTransgeneBuilder::codon_freqs$Amino, upper_codons)
  aa_to_syn <- split(upper_codons, WormTransgeneBuilder::codon_freqs$Amino)

  num_windows <- aa_len - 3
  window_scores <- numeric(num_windows)
  loop_vec <- cds_vec
  # The same sequence as numbers (A = 0 ... T = 3), which is what the scoring works on.
  codes <- glo_base_codes(loop_vec)

  # Establish baseline fitness
  baseline_windows <- seq_len(floor(num_windows))
  window_scores[baseline_windows] <- glo_window_scores(codes, (baseline_windows - 1) * 3 + 1)

  iteration <- 1
  max_iterations <- 5
  changes_made <- TRUE

  while (iteration <= max_iterations && changes_made) {
    changes_made <- FALSE
    iteration <- iteration + 1
    loop_window_scores <- window_scores

    codon_indices <- sample(1:aa_len)
    counter <- 1

    while (counter <= length(codon_indices)) {
      codon_idx <- codon_indices[counter]
      pos <- (codon_idx - 1) * 3 + 1
      old_codon <- paste(loop_vec[pos:(pos + 2)], collapse = "")

      aa <- unname(codon_to_aa[old_codon])
      if (is.na(aa)) { counter <- counter + 1; next }

      synonymous <- aa_to_syn[[aa]]
      if (length(synonymous) < 2) { counter <- counter + 1; next }

      syn_options <- synonymous[synonymous != old_codon]
      new_codon <- if(length(syn_options) == 1) syn_options else sample(syn_options, 1)
      new_codon_vec <- unlist(strsplit(new_codon, ""))

      # Identify ONLY the 12-mer windows touching this mutated codon
      affected_windows <- max(1, codon_idx - 3):min(num_windows, codon_idx)
      old_local_sum <- sum(loop_window_scores[affected_windows])

      # Try the new codon in place; it is put back below if it does not help.
      old_codes <- codes[pos:(pos + 2)]
      codes[pos:(pos + 2)] <- glo_base_codes(new_codon_vec)

      # Calculate fitness strictly for the overlapping mutated zone
      temp_new_scores <- glo_window_scores(codes, (affected_windows - 1) * 3 + 1)
      new_local_sum <- 0
      for (score in temp_new_scores) new_local_sum <- new_local_sum + score

      # Greedy acceptance: if it improves the local zone, keep it
      if (new_local_sum > old_local_sum) {
        loop_vec[pos:(pos + 2)] <- new_codon_vec
        loop_window_scores[affected_windows] <- temp_new_scores
        changes_made <- TRUE
      } else {
        codes[pos:(pos + 2)] <- old_codes
      }
      counter <- counter + 1
    }
    window_scores <- loop_window_scores
  }

  return(paste(loop_vec, collapse = ""))
}
