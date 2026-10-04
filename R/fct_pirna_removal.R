#' Get Exact Hamming Distances for all piRNAs
#'
#' Slices the input sequence into 20-mer sliding windows and computes the
#' minimum Hamming distance between any 20-mer in the sequence and the
#' established C. elegans piRNA database. The distances are exact; they are computed
#' with whole-number arithmetic on a 2-bit encoding of the bases (see
#' `scan_kmers_against_pirnas()`), which is much faster than comparing characters.
#'
#' @param dna_seq Character string; the input DNA sequence.
#' @return A numeric vector of the minimum Hamming distances for each known piRNA.
#' @export
get_pirna_distances <- function(dna_seq) {
  current_seq <- toupper(gsub("[^ATGCNatgcn]", "", dna_seq))
  seq_len <- nchar(current_seq)
  if (seq_len < 20) return(rep(NA, length(WormTransgeneBuilder::pirna_sequences)))

  starts <- 1:(seq_len - 19)
  kmers <- substring(current_seq, starts, starts + 19)

  # The minimum distance for each piRNA across the entire sequence (NA if it cannot be compared).
  scan_kmers_against_pirnas(kmers)$col_min
}

#' Shield piRNA Sites with the Smallest CAI Loss (deterministic)
#'
#' Used for the "High expression" (CAI1) method, where every codon has already been
#' set to the single best codon and there is nothing to sample from. Works one
#' change at a time: take the first 20-mer window that is closer to a piRNA than
#' the threshold, try every synonymous change in the codons under it, keep only
#' the changes that (a) increase that window's distance to the piRNA database and
#' (b) do not bring any other window closer to a piRNA than it already was (never
#' below the threshold), and apply the one that costs the least Codon Adaptation
#' Index. The cost of replacing codon `old` by `new` is `log(w_new) - log(w_old)`.
#' Ties go to the larger gain in distance, then the earlier codon, then the
#' alphabetically first codon, so the same input always gives the same output.
#' Because every change strictly improves a bounded score, the search always ends.
#'
#' @param current_seq Character string; upper-case coding sequence (A, C, G, T, N).
#' @param target_min_hamming Numeric; windows closer than this to a piRNA are fixed.
#' @param score_column Character string; graded codon-weight column in `codon_freqs`.
#' @param max_steps Integer; safety limit on the number of codon changes.
#' @return The shielded sequence. A warning is given if windows remain below the
#'   threshold because no synonymous change could improve them.
#' @noRd
remove_pirna_sites_max_cai <- function(current_seq, target_min_hamming,
                                       score_column = "BringMansWeights", max_steps = 2000,
                                       initial_scan = NULL) {
  cf <- WormTransgeneBuilder::codon_freqs
  codons_up <- toupper(cf$Codon)
  codon_to_aa <- setNames(cf$Amino, codons_up)
  weights <- setNames(cf[[score_column]], codons_up)
  aa_to_syn <- split(codons_up, cf$Amino)
  log_w <- function(codon) log(pmax(weights[codon], 1e-6))

  pirnas <- WormTransgeneBuilder::pirna_sequences

  # 20-mer windows of `dna`.
  windows_of <- function(dna) {
    starts <- seq_len(nchar(dna) - 19)
    substring(dna, starts, starts + 19)
  }
  # Minimum distance of every window to any piRNA in `reference` (encoded with
  # make_pirna_reference()), capped at the threshold: only values below the threshold
  # matter, so anything larger is stored as the threshold.
  window_min_dists <- function(dna, reference) {
    if (length(reference$lo) == 0) return(rep(target_min_hamming, nchar(dna) - 19))
    pmin(scan_kmers_against_pirnas(windows_of(dna), reference)$row_min, target_min_hamming)
  }

  seq_vec <- unlist(strsplit(current_seq, ""))
  n <- length(seq_vec)
  # The caller has usually scanned the whole sequence already; reuse that instead of
  # scanning it again.
  wmin <- if (is.null(initial_scan)) window_min_dists(current_seq, pirna_reference())
          else pmin(initial_scan$row_min, target_min_hamming)
  unfixable <- integer(0)

  for (step in seq_len(max_steps)) {
    violating <- setdiff(which(wmin < target_min_hamming), unfixable)
    if (length(violating) == 0) break
    i <- violating[1]  # first offending window

    codon_range <- ((i - 1) %/% 3 + 1):((i + 18) %/% 3 + 1)
    codon_range <- codon_range[(codon_range - 1) * 3 + 3 <= n]

    # One codon change alters at most 3 bases, so it can move a window's distance
    # to a piRNA by at most 3. Only piRNAs within (threshold + 2) of the region can
    # therefore fall below the threshold after this step: compare candidates to
    # those alone. This gives exactly the same answer as using the full database.
    region_lo <- max(1, (codon_range[1] - 1) * 3 + 1 - 19)
    region_hi <- min(length(wmin), (codon_range[length(codon_range)] - 1) * 3 + 3)
    region_dna <- paste(seq_vec[region_lo:(region_hi + 19)], collapse = "")
    region_min <- scan_kmers_against_pirnas(windows_of(region_dna), pirna_reference())$col_min
    relevant <- make_pirna_reference(pirnas[which(region_min <= target_min_hamming + 2)])

    best <- NULL  # list(score, gain, lo, new_local, codon_index, new_codon)
    for (ci in codon_range) {
      p <- (ci - 1) * 3 + 1
      old_codon <- paste(seq_vec[p:(p + 2)], collapse = "")
      aa <- unname(codon_to_aa[old_codon])
      if (is.na(aa)) next
      alternatives <- sort(setdiff(aa_to_syn[[aa]], old_codon))
      if (length(alternatives) == 0) next

      lo <- max(1, p - 19)                     # windows that overlap this codon
      hi <- min(length(wmin), p + 2)
      old_local <- wmin[lo:hi]
      region <- lo:(hi + 19)

      for (alt in alternatives) {
        trial <- seq_vec[region]
        trial[(p - lo + 1):(p - lo + 3)] <- unlist(strsplit(alt, ""))
        new_local <- window_min_dists(paste(trial, collapse = ""), relevant)

        # No window may get worse (below the threshold), and window i must improve.
        if (any(new_local < old_local)) next
        gain <- new_local[i - lo + 1] - old_local[i - lo + 1]
        if (gain <= 0) next

        score <- unname(log_w(alt) - log_w(old_codon))
        better <- is.null(best) || score > best$score + 1e-12 ||
          (abs(score - best$score) <= 1e-12 && gain > best$gain)
        if (better) best <- list(score = score, gain = gain, lo = lo, new_local = new_local,
                                 codon_index = ci, new_codon = alt)
      }
    }

    if (is.null(best)) {
      unfixable <- c(unfixable, i)
      next
    }
    p <- (best$codon_index - 1) * 3 + 1
    seq_vec[p:(p + 2)] <- unlist(strsplit(best$new_codon, ""))
    wmin[best$lo:(best$lo + length(best$new_local) - 1)] <- best$new_local
  }

  if (any(wmin < target_min_hamming)) {
    warning(sum(wmin < target_min_hamming), " 20-mer window(s) remain closer than ", target_min_hamming,
            " mismatches to a piRNA because no synonymous change could improve them.", call. = FALSE)
  }
  paste(seq_vec, collapse = "")
}

#' Find and Remove piRNA Binding Sites
#'
#' Identifies 20-mer regions within the sequence that are dangerously homologous
#' to known piRNAs (falling below the target Hamming distance). Mutates codons
#' strictly within those vulnerable windows using synonymous alternatives to
#' push the Hamming distance into the safe threshold.
#'
#' @param dna_seq Character string; the coding DNA sequence to shield.
#' @param usage_column_name Character string; the codon table column used to select synonymous mutations.
#' @param target_min_hamming Numeric; the safety threshold (default is usually 2 or 3 mismatches).
#' @details With `usage_column_name = "CAI1"` (the "High expression" method) the codon changes are chosen
#'   deterministically, one at a time, as the change with the smallest loss of Codon Adaptation Index.
#' @return A list containing the shielded `sequence`, `initial_dists`, and `final_dists`.
#' @export
remove_pirna_sites <- function(dna_seq, usage_column_name = "Ubiq", target_min_hamming = 2) {
  current_seq <- toupper(gsub("[^ATGCNatgcn]", "", dna_seq))

  if (nchar(current_seq) < 20) {
    initial_dists <- get_pirna_distances(current_seq)  # NA for every piRNA
    return(list(sequence = current_seq, initial_dists = initial_dists, final_dists = initial_dists))
  }

  seq_len <- nchar(current_seq)
  starts <- 1:(seq_len - 19)
  kmers <- substring(current_seq, starts, starts + 19)

  # One scan of the sequence gives both the distance to each piRNA and the windows that
  # are closer than the threshold (this used to be done in two separate full passes).
  scan <- scan_kmers_against_pirnas(kmers, threshold = target_min_hamming)
  initial_dists <- scan$col_min
  if (all(is.na(initial_dists))) {
    return(list(sequence = current_seq, initial_dists = initial_dists, final_dists = initial_dists))
  }

  # "High expression" (CAI1) has nothing to sample from (all other codon weights are 0),
  # so use the deterministic, smallest-CAI-loss search instead of random re-rolling.
  if (identical(usage_column_name, "CAI1")) {
    final_seq <- remove_pirna_sites_max_cai(current_seq, target_min_hamming, "BringMansWeights",
                                            initial_scan = scan)
    final_dists <- if (identical(final_seq, current_seq)) initial_dists else get_pirna_distances(final_seq)
    return(list(sequence = final_seq, initial_dists = initial_dists, final_dists = final_dists))
  }

  # The windows that are too close to a piRNA, in the order the random search visits them.
  hits <- scan$pairs

  if (nrow(hits) == 0) {
    return(list(sequence = current_seq, initial_dists = initial_dists, final_dists = initial_dists))
  }

  # Snap hit bases to their respective codon starting boundaries
  stpos <- unique(hits$window)
  stpos <- unique(c(stpos, stpos + 3, stpos + 6, stpos + 9, stpos + 12, stpos + 15))
  stpos <- unique(stpos - (stpos %% 3) + 1)

  seq_vec <- unlist(strsplit(current_seq, ""))
  stpos <- stpos[stpos <= (length(seq_vec) - 2)]

  upper_codons <- toupper(WormTransgeneBuilder::codon_freqs$Codon)
  codon_to_aa <- setNames(WormTransgeneBuilder::codon_freqs$Amino, upper_codons)
  weights_vec <- setNames(WormTransgeneBuilder::codon_freqs[[usage_column_name]], upper_codons)
  aa_to_syn <- split(upper_codons, WormTransgeneBuilder::codon_freqs$Amino)

  for (pos in stpos) {
    old_codon <- paste(seq_vec[pos:(pos + 2)], collapse = "")
    aa <- unname(codon_to_aa[old_codon])

    if (!is.na(aa)) {
      synonymous <- aa_to_syn[[aa]]
      synonymous <- synonymous[synonymous != old_codon]

      if (length(synonymous) > 0) {
        syn_weights <- weights_vec[synonymous]
        if (sum(syn_weights, na.rm = TRUE) > 0) new_codon <- sample(synonymous, 1, prob = syn_weights)
        else new_codon <- sample(synonymous, 1)
        seq_vec[pos:(pos + 2)] <- unlist(strsplit(new_codon, ""))
      }
    }
  }

  final_seq <- paste(seq_vec, collapse = "")
  final_dists <- get_pirna_distances(final_seq)

  return(list(sequence = final_seq, initial_dists = initial_dists, final_dists = final_dists))
}
