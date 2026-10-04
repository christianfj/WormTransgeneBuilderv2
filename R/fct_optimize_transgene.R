#' Biologically-Aware Intron Inserter
#'
#' Inserts a vector of intron sequences into a coding sequence (CDS) while
#' attempting to place them at optimal boundaries (e.g., matching the AAG/CAG consensus).
#'
#' @param cds_seq Character string; the coding DNA sequence to modify.
#' @param introns Character vector; the intron sequences to insert.
#' @param spacing_type Character string; either "early" (front-loaded) or "equidistant".
#' @param in_frame Logical; if TRUE (the default) an intron is always inserted between two codons
#'   (phase 0). If FALSE it may be inserted after any base, so it can fall inside a codon.
#' @details The best place for each intron is looked for within 15 codons (45 bases) either side of
#'   its target. A place scores 3 if the three bases before the intron are AAG, 2 if they are CAG and 1
#'   if they end in AG; the highest score wins and the place nearest the target wins a tie. With
#'   `in_frame = TRUE` the three bases before the intron are one codon.
#' @return A list containing the modified sequence (`seq`) and a numeric vector
#'   of the insertion base pairs (`insertions`), the number of bases of `cds_seq` before each intron.
#' @export
insert_introns <- function(cds_seq, introns, spacing_type = "early", in_frame = TRUE) {
  cds_len <- nchar(cds_seq)
  num_introns <- length(introns)
  if (num_introns == 0) return(list(seq = cds_seq, insertions = numeric(0)))

  total_codons <- floor(cds_len / 3)
  target_codons <- numeric(num_introns)

  if (spacing_type == "equidistant") {
    step <- total_codons / (num_introns + 1)
    for (i in 1:num_introns) target_codons[i] <- round(i * step)
  } else {
    current_target <- 50
    for (i in 1:num_introns) {
      if (current_target > total_codons - 15) current_target <- total_codons - 15
      target_codons[i] <- current_target
      current_target <- current_target + 50
    }
    target_codons <- sort(unique(target_codons))
    if (length(target_codons) < num_introns) {
      step <- total_codons / (num_introns + 1)
      for (i in 1:num_introns) target_codons[i] <- round(i * step)
    }
  }

  actual_insertion_bps <- numeric(num_introns)

  for (i in 1:num_introns) {
    target_c <- target_codons[i]
    search_radius <- 15

    if (!in_frame) {
      # Any base may precede the intron: score every base position near the target, not just codon ends.
      target_b <- target_c * 3
      prev_b <- if (i == 1) 15 else actual_insertion_bps[i - 1] + 30
      start_b <- max(prev_b, target_b - search_radius * 3)
      end_b <- min(cds_len - 15, target_b + search_radius * 3)
      if (start_b > end_b) {
        actual_insertion_bps[i] <- target_b
        next
      }
      positions <- start_b:end_b
      last3 <- toupper(substring(cds_seq, positions - 2, positions))
      scores <- ifelse(last3 == "AAG", 3, ifelse(last3 == "CAG", 2, ifelse(substr(last3, 2, 3) == "AG", 1, 0)))
      actual_insertion_bps[i] <- positions[order(-scores, abs(positions - target_b))[1]]
      next
    }
    prev_boundary <- if(i == 1) 5 else (actual_insertion_bps[i-1]/3 + 10)

    possible_codons <- list()
    start_c <- max(prev_boundary, target_c - search_radius)
    end_c <- min(total_codons - 5, target_c + search_radius)

    if (start_c > end_c) {
      actual_insertion_bps[i] <- target_c * 3
      next
    }

    for (c in start_c:end_c) {
      prec_codon <- substr(cds_seq, (c - 1) * 3 + 1, c * 3)
      possible_codons[[length(possible_codons) + 1]] <- list(codon = toupper(prec_codon), pos = c)
    }

    if (length(possible_codons) == 0) {
      actual_insertion_bps[i] <- target_c * 3
      next
    }

    for (j in seq_along(possible_codons)) {
      seq <- possible_codons[[j]]$codon
      score <- 0
      if (seq == "AAG") { score <- 3 }
      else if (seq == "CAG") { score <- 2 }
      else if (substr(seq, 2, 3) == "AG") { score <- 1 }
      possible_codons[[j]]$score <- score
      possible_codons[[j]]$distance <- abs(possible_codons[[j]]$pos - target_c)
    }

    sorted_matches <- possible_codons[order(-sapply(possible_codons, function(x) x$score),
                                            sapply(possible_codons, function(x) x$distance))]
    actual_insertion_bps[i] <- sorted_matches[[1]]$pos * 3
  }

  pieces <- list()
  last_pos <- 1
  for (i in 1:num_introns) {
    insert_pos <- actual_insertion_bps[i]
    pieces[[length(pieces) + 1]] <- substr(cds_seq, last_pos, insert_pos)
    pieces[[length(pieces) + 1]] <- tolower(introns[i])
    last_pos <- insert_pos + 1
  }
  pieces[[length(pieces) + 1]] <- substr(cds_seq, last_pos, cds_len)

  return(list(seq = paste(pieces, collapse = ""), insertions = actual_insertion_bps))
}

#' Probabilistic DNA Recoder
#'
#' Recodes a DNA sequence by substituting synonymous codons according to a
#' specified usage frequency column in the internal codon table.
#'
#' @param dna_seq Character string; the DNA sequence to recode.
#' @param usage_column_name Character string; the name of the column in the internal codon frequency table to use.
#' @return A character string of the recoded DNA sequence.
#' @export
probabilistic_recode <- function(dna_seq, usage_column_name) {
  current_seq <- toupper(gsub("[^ATGCNatgcn]", "", dna_seq))
  cds_vec <- unlist(strsplit(current_seq, ""))
  aa_len <- length(cds_vec) / 3

  # FIXED: Added package namespace to codon_freqs
  upper_codons <- toupper(WormTransgeneBuilder::codon_freqs$Codon)
  codon_to_aa <- setNames(WormTransgeneBuilder::codon_freqs$Amino, upper_codons)
  weights_vec <- setNames(WormTransgeneBuilder::codon_freqs[[usage_column_name]], upper_codons)
  aa_to_syn <- split(upper_codons, WormTransgeneBuilder::codon_freqs$Amino)

  for (i in 1:aa_len) {
    pos <- (i - 1) * 3 + 1
    old_codon <- paste(cds_vec[pos:(pos + 2)], collapse = "")
    aa <- unname(codon_to_aa[old_codon])

    # In the CAI1 column all three stop codons have the same top weight, so the
    # choice among them would be random. Keep the stop codon that was supplied so
    # that the "High expression" recoding is fully deterministic.
    if (identical(usage_column_name, "CAI1") && !is.na(aa) && aa == "X") next

    if (!is.na(aa)) {
      synonymous <- aa_to_syn[[aa]]
      if (length(synonymous) > 1) {
        syn_weights <- weights_vec[synonymous]
        if (sum(syn_weights, na.rm = TRUE) > 0) new_codon <- sample(synonymous, 1, prob = syn_weights)
        else new_codon <- sample(synonymous, 1)
        cds_vec[pos:(pos+2)] <- unlist(strsplit(new_codon, ""))
      }
    }
  }
  return(paste(cds_vec, collapse = ""))
}

#' Find the RNAfold Program
#'
#' Looks for ViennaRNA's `RNAfold` first on the PATH and then in a few usual
#' install locations (a `~` means the current user's home directory).
#'
#' @param on_path Character string; what `Sys.which("RNAfold")` found ("" if nothing).
#' @param fallbacks Character vector; locations to try, in order, when it is not on the PATH.
#' @return The path to `RNAfold` as a plain string, or `""` if it was not found.
#' @noRd
find_rnafold <- function(on_path = Sys.which("RNAfold"),
                         fallbacks = c("bin/RNAfold", "inst/bin/RNAfold", "~/miniconda3/bin/RNAfold",
                                       "~/anaconda3/bin/RNAfold", "/opt/homebrew/bin/RNAfold",
                                       "/usr/local/bin/RNAfold")) {
  if (any(nzchar(on_path))) return(unname(on_path[nzchar(on_path)][1]))
  for (p in fallbacks) {
    expanded <- path.expand(p)
    if (file.exists(expanded)) return(expanded)
  }
  ""
}

#' Pure-R Ribosomal Binding Site Optimization
#'
#' Generates candidate start sequences and uses ViennaRNA's RNAfold to select the
#' variant with the least stable secondary structure (highest MFE). Any rare piRNA
#' homologies introduced here will be natively handled by the downstream pipeline.
#'
#' @param cds_start_seq Character string; the first 39bp of the CDS.
#' @param usage_column_name Character string; the frequency column used for generating candidate sequences.
#' @param num_candidates Integer; the number of candidate sequences to generate and test.
#' @return A character string containing the optimized start sequence.
#' @export
optimize_rbs <- function(cds_start_seq, usage_column_name, num_candidates = 100) {
  rnafold_bin <- find_rnafold()

  if (rnafold_bin == "") {
    warning("RNAfold binary not found in path. Skipping RBS optimization.")
    return(cds_start_seq)
  }

  candidates <- c(cds_start_seq)
  for (i in 1:num_candidates) candidates <- c(candidates, probabilistic_recode(cds_start_seq, usage_column_name))
  candidates <- unique(candidates)
  fold_targets <- paste0("AAAA", candidates)

  rnafold_out <- tryCatch({
    system2(rnafold_bin, args = c("--noPS"), input = fold_targets, stdout = TRUE, stderr = FALSE)
  }, error = function(e) return(NULL))

  if (is.null(rnafold_out) || length(rnafold_out) == 0) return(cds_start_seq)

  struct_lines <- grep("\\(", rnafold_out, value = TRUE)
  energies <- numeric(length(candidates))
  for (i in seq_along(candidates)) {
    if (i <= length(struct_lines)) {
      line <- struct_lines[i]
      mfe_str <- gsub(".*\\(([^)]+)\\).*", "\\1", line)
      mfe_str <- gsub(" ", "", mfe_str)
      mfe_val <- suppressWarnings(as.numeric(mfe_str))
      if (is.na(mfe_val)) mfe_val <- -999
      energies[i] <- mfe_val
    } else energies[i] <- -999
  }

  # Simply return the candidate with the highest MFE
  best_idx <- which.max(energies)
  return(candidates[best_idx])
}

#' Naive Back-Translator
#'
#' Converts an amino acid sequence into a dummy DNA sequence using common C. elegans codons.
#' This allows the optimization engine to process protein inputs natively.
#'
#' @param prot_seq Character string; the input amino acid sequence.
#' @return A character string of the translated dummy DNA sequence.
#' @noRd
back_translate <- function(prot_seq) {
  bt_map <- c(
    'A'='GCT', 'C'='TGT', 'D'='GAT', 'E'='GAA', 'F'='TTT', 'G'='GGT',
    'H'='CAT', 'I'='ATT', 'K'='AAG', 'L'='CTT', 'M'='ATG', 'N'='AAT',
    'P'='CCT', 'Q'='CAA', 'R'='CGT', 'S'='TCT', 'T'='ACT', 'V'='GTT',
    'W'='TGG', 'Y'='TAT', '*'='TAA'
  )
  prot_vec <- unlist(strsplit(toupper(prot_seq), ""))
  dna_vec <- bt_map[prot_vec]
  dna_vec[is.na(dna_vec)] <- "NNN"
  return(paste(dna_vec, collapse = ""))
}

#' Master Transgene Optimization Assembly Line
#'
#' Coordinates the full transgene assembly pipeline, including molecular tag attachment,
#' codon optimization, piRNA shielding, restriction site removal, and intron insertion.
#'
#' @param dna_seq Character string; the raw input sequence (DNA or Protein).
#' @param opt_method Character string; specific codon recoding method to deploy.
#' @param enzymes Character vector; restriction enzymes to screen and remove.
#' @param remove_pirnas Logical; whether to actively shield against piRNA targets.
#' @param target_min_hamming Numeric; threshold for piRNA exclusion.
#' @param introns Character vector; sequence strings for intron insertion.
#' @param spacing_type Character string; intron distribution logic ("early" or "equidistant").
#' @param force_reading_frame Logical; if TRUE introns are inserted only between codons (in the reading
#'   frame, phase 0); if FALSE they may be inserted after any base.
#' @param add_consensus_start Logical; whether to pad the 5' end with "aaaa".
#' @param rbs_opt Logical; whether to invoke the ViennaRNA optimization step.
#' @param promoter_seq Character string; sequence to prepend as a promoter.
#' @param utr5_seq Character string; sequence to prepend as a 5' UTR.
#' @param utr3_seq Character string; sequence to append as a 3' UTR.
#' @param is_protein Logical; signals the engine to apply naive back-translation to the input.
#' @param tag_n_seqs Character vector; sequences for N-terminal molecular tags.
#' @param tag_c_seqs Character vector; sequences for C-terminal molecular tags.
#' @return A complex list containing the final assembled sequence, structural coordinates,
#'   optimization logs, and piRNA distance metrics.
#' @export
optimize_transgene <- function(dna_seq, opt_method = "Ubiq", enzymes = c("BsaI"), remove_pirnas = TRUE,
                               target_min_hamming = 2, introns = NULL, spacing_type = "early",
                               force_reading_frame = TRUE, add_consensus_start = FALSE, rbs_opt = FALSE,
                               promoter_seq = NULL, utr5_seq = NULL, utr3_seq = NULL, is_protein = FALSE,
                               tag_n_seqs = NULL, tag_c_seqs = NULL) {

  if (is_protein) {
    current_seq <- back_translate(dna_seq)
    if (!grepl("(TAA|TAG|TGA)$", toupper(current_seq))) current_seq <- paste0(current_seq, "TAA")
  } else {
    val_init <- check_cds_validity(dna_seq)
    if (!val_init$valid) stop(paste("Invalid CDS:", paste(val_init$warnings, collapse = " | ")))
    current_seq <- toupper(gsub("[^ATGCNatgcn]", "", dna_seq))
  }

  report_log <- c(paste("Original Input Type:", ifelse(is_protein, "Protein", "DNA")))

  clean_n_tags <- c()
  clean_c_tags <- c()
  if (!is.null(tag_n_seqs) && length(tag_n_seqs) > 0) clean_n_tags <- toupper(gsub("[^ATGCNatgcn]", "", tag_n_seqs))
  if (!is.null(tag_c_seqs) && length(tag_c_seqs) > 0) clean_c_tags <- toupper(gsub("[^ATGCNatgcn]", "", tag_c_seqs))

  tag_n_starts_base <- numeric(0)
  tag_n_ends_base <- numeric(0)
  if (length(clean_n_tags) > 0) {
    curr_pos <- 4
    for (t in clean_n_tags) {
      tag_n_starts_base <- c(tag_n_starts_base, curr_pos)
      curr_pos <- curr_pos + nchar(t)
      tag_n_ends_base <- c(tag_n_ends_base, curr_pos - 1)
    }
    current_seq <- paste0(
      substr(current_seq, 1, 3),
      paste(clean_n_tags, collapse=""),
      substr(current_seq, 4, nchar(current_seq))
    )
  }

  tag_c_starts_base <- numeric(0)
  tag_c_ends_base <- numeric(0)
  if (length(clean_c_tags) > 0) {
    curr_pos <- nchar(current_seq) - 2
    for (t in clean_c_tags) {
      tag_c_starts_base <- c(tag_c_starts_base, curr_pos)
      curr_pos <- curr_pos + nchar(t)
      tag_c_ends_base <- c(tag_c_ends_base, curr_pos - 1)
    }
    current_seq <- paste0(
      substr(current_seq, 1, nchar(current_seq) - 3),
      paste(clean_c_tags, collapse=""),
      substr(current_seq, nchar(current_seq) - 2, nchar(current_seq))
    )
  }

  true_initial_dists <- if (!is_protein) get_pirna_distances(current_seq) else NULL

  # 2. Base Recoding & RBS Optimization
  if (rbs_opt && nchar(current_seq) >= 39) {
    cds_start <- substr(current_seq, 1, 39)
    cds_rest  <- substr(current_seq, 40, nchar(current_seq))
    shield_col <- if(opt_method %in% c("None", "GLO", "Chance")) "Ubiq" else opt_method

    # RBS Optimization executes instantly based solely on RNAfold MFE
    best_start <- optimize_rbs(cds_start, shield_col, num_candidates = 100)

    if (opt_method == "GLO") best_rest <- optimize_glo(cds_rest)
    else if (opt_method != "None") best_rest <- probabilistic_recode(cds_rest, usage_column_name = opt_method)
    else best_rest <- cds_rest
    current_seq <- paste0(best_start, best_rest)

  } else {
    if (opt_method == "GLO") current_seq <- optimize_glo(current_seq)
    else if (opt_method != "None") current_seq <- probabilistic_recode(current_seq, usage_column_name = opt_method)
  }

  final_dists <- NULL

  # 3. Global piRNA Shielding (Will naturally fix any issues from RBS generation)
  if (remove_pirnas) {
    shield_col <- if(opt_method %in% c("None", "GLO", "Chance")) "Ubiq" else opt_method
    pirna_out <- remove_pirna_sites(current_seq, usage_column_name = shield_col, target_min_hamming = target_min_hamming)
    current_seq <- pirna_out$sequence
    final_dists <- pirna_out$final_dists
  } else {
    final_dists <- get_pirna_distances(current_seq)
  }

  if (length(enzymes) > 0) {
    shield_col <- if(opt_method %in% c("None", "GLO", "Chance")) "Ubiq" else opt_method
    current_seq <- remove_enzyme_sites(current_seq, enzymes, usage_column_name = shield_col)
  }

  optimized_cds <- current_seq

  intron_insertions <- numeric(0)  # bases of optimized_cds after which each intron was inserted
  intron_lengths <- numeric(0)     # length of each inserted intron

  if (!is.null(introns)) {
    introns <- tolower(gsub("[^ATGCatgc]", "", introns[!is.na(introns) & nchar(introns) > 0]))
    if (length(introns) > 0) {
      intron_res <- insert_introns(current_seq, introns, spacing_type = spacing_type, in_frame = force_reading_frame)
      current_seq <- intron_res$seq
      intron_insertions <- intron_res$insertions
      intron_lengths <- nchar(introns)

      if (length(tag_n_starts_base) > 0) {
        for (j in seq_along(tag_n_starts_base)) {
          shift_start <- sum(nchar(introns)[intron_res$insertions < tag_n_starts_base[j]])
          shift_end <- sum(nchar(introns)[intron_res$insertions < tag_n_ends_base[j]])
          tag_n_starts_base[j] <- tag_n_starts_base[j] + shift_start
          tag_n_ends_base[j] <- tag_n_ends_base[j] + shift_end
        }
      }

      if (length(tag_c_starts_base) > 0) {
        for (j in seq_along(tag_c_starts_base)) {
          shift_start <- sum(nchar(introns)[intron_res$insertions < tag_c_starts_base[j]])
          shift_end <- sum(nchar(introns)[intron_res$insertions < tag_c_ends_base[j]])
          tag_c_starts_base[j] <- tag_c_starts_base[j] + shift_start
          tag_c_ends_base[j] <- tag_c_ends_base[j] + shift_end
        }
      }
    }
  }

  cds_start_idx <- 1
  cds_end_idx <- nchar(current_seq)

  if (add_consensus_start || rbs_opt) {
    if (substr(current_seq, 1, 3) == "ATG") current_seq <- substr(current_seq, 4, nchar(current_seq))
    current_seq <- paste0("aaaaATG", current_seq)
    cds_start_idx <- 5
    cds_end_idx <- nchar(current_seq)
  }

  if (!is.null(utr5_seq)) {
    clean_utr5 <- tolower(gsub("[^ATGCatgc]", "", utr5_seq))
    current_seq <- paste0(clean_utr5, current_seq)
    cds_start_idx <- cds_start_idx + nchar(clean_utr5)
    cds_end_idx <- cds_end_idx + nchar(clean_utr5)
  }

  if (!is.null(promoter_seq)) {
    clean_promoter <- tolower(gsub("[^ATGCatgc]", "", promoter_seq))
    current_seq <- paste0(clean_promoter, current_seq)
    cds_start_idx <- cds_start_idx + nchar(clean_promoter)
    cds_end_idx <- cds_end_idx + nchar(clean_promoter)
  }

  final_tag_n_starts <- if(length(tag_n_starts_base)>0) tag_n_starts_base + cds_start_idx - 1 else NULL
  final_tag_n_ends <- if(length(tag_n_ends_base)>0) tag_n_ends_base + cds_start_idx - 1 else NULL
  final_tag_c_starts <- if(length(tag_c_starts_base)>0) tag_c_starts_base + cds_start_idx - 1 else NULL
  final_tag_c_ends <- if(length(tag_c_ends_base)>0) tag_c_ends_base + cds_start_idx - 1 else NULL

  if (!is.null(utr3_seq)) current_seq <- paste0(current_seq, tolower(gsub("[^ATGCatgc]", "", utr3_seq)))

  return(list(
    original_seq = dna_seq,
    optimized_cds = optimized_cds,
    final_seq = current_seq,
    log = report_log,
    initial_pirna_dists = true_initial_dists,
    final_pirna_dists = final_dists,
    cds_start = cds_start_idx,
    cds_end = cds_end_idx,
    tag_n_starts = final_tag_n_starts,
    tag_n_ends = final_tag_n_ends,
    tag_c_starts = final_tag_c_starts,
    tag_c_ends = final_tag_c_ends,
    intron_insertions = intron_insertions,
    intron_lengths = intron_lengths
  ))
}
