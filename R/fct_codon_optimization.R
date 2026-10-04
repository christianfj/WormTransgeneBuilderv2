#' Check Biological Validity of CDS
#'
#' Implements strict biological rules to validate a coding sequence (CDS):
#' requires a valid start codon (ATG), a valid stop codon (TAA, TAG, TGA),
#' a length that is a multiple of 3, and the absence of internal stop codons.
#'
#' @param dna_seq Character string; the input DNA sequence to validate.
#' @return A list containing:
#'   \item{valid}{Logical; TRUE if the sequence passes all checks, FALSE otherwise.}
#'   \item{warnings}{Character vector; specific error messages if validation fails.}
#' @export
check_cds_validity <- function(dna_seq) {
  dna_seq <- toupper(gsub("[^ATGCNatgcn]", "", dna_seq))
  len <- nchar(dna_seq)
  warnings <- c()

  if (len < 6) return(list(valid = FALSE, warnings = "Sequence must be at least 6 bp."))

  if (len %% 3 != 0) {
    warnings <- c(warnings, paste0("Sequence length (", len, ") is not a multiple of 3."))
  }

  start_codon <- substr(dna_seq, 1, 3)
  if (start_codon != "ATG") {
    warnings <- c(warnings, paste("Sequence does not begin with a start codon (ATG). Found:", start_codon))
  }

  stop_codon <- substr(dna_seq, len - 2, len)
  if (!stop_codon %in% c("TAA", "TAG", "TGA")) {
    warnings <- c(warnings, paste("Sequence does not end with a valid stop codon. Found:", stop_codon))
  }

  if (len >= 6) {
    codons <- substring(dna_seq, seq(1, len - 5, 3), seq(3, len - 3, 3))
    internal_stops <- codons %in% c("TAA", "TAG", "TGA")
    if (any(internal_stops)) {
      warnings <- c(warnings, paste("Internal stop codon found at position:", paste(which(internal_stops) * 3 - 2, collapse = ", ")))
    }
  }

  return(list(valid = length(warnings) == 0, warnings = warnings))
}

#' Sample a Synonymous Codon
#'
#' Selects a different synonymous codon for a given amino acid based on the
#' specified usage frequencies. Used primarily to mutate sequences when
#' destroying restriction enzyme sites or piRNA targets.
#'
#' @param aa Character string; the single-letter amino acid code.
#' @param old_codon Character string; the 3-nucleotide codon to be replaced.
#' @param usage_column_name Character string; the column name in the internal codon frequency table to use for weighted sampling.
#' @return A character string of the newly sampled 3-nucleotide codon.
#' @export
sample_new_codon <- function(aa, old_codon, usage_column_name) {
  df <- WormTransgeneBuilder::codon_freqs[WormTransgeneBuilder::codon_freqs$Amino == aa, ]
  if (nrow(df) < 2) return(toupper(old_codon))

  df_subset <- df[toupper(df$Codon) != toupper(old_codon), ]
  if (nrow(df_subset) == 0) return(toupper(old_codon))

  new_codon <- sample(df_subset$Codon, size = 1, prob = df_subset[[usage_column_name]])

  return(toupper(as.character(new_codon)))
}
