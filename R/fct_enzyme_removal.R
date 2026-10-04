#' Find Restriction Site Positions (IUPAC-aware)
#'
#' Searches a DNA sequence for one or more recognition sites. Sites may contain
#' IUPAC ambiguity codes (for example `Y` = C or T, `R` = A or G, `N` = any base),
#' so a site such as BanI's `GGYRCC` matches `GGCGCC`, `GGTACC`, and so on.
#' Ambiguity codes are interpreted in the *site* only: an `N` in the sequence being
#' searched is treated as an unknown base and does not match a specific base.
#'
#' The caller is responsible for including the reverse-complement of each site.
#'
#' @param dna_seq Character string; the (upper-case) sequence to search. Any
#'   character that is not an IUPAC DNA code (for example `*` or stray letters
#'   in a pasted sequence) is treated as an unknown base, which keeps every
#'   reported position aligned with the original string.
#' @param sites Character vector; recognition sites, which may contain IUPAC codes.
#' @return `find_enzyme_site_ranges()` returns a data frame with 1-based `start`
#'   and `end` columns, one row per match (unsorted, may contain duplicates; zero
#'   rows if nothing matched). `find_enzyme_site_starts()` returns just the start
#'   positions as an integer vector.
#' @noRd
find_enzyme_site_ranges <- function(dna_seq, sites) {
  subject <- gsub("[^ACGTRYSWKMBDHVN]", "N", dna_seq)
  starts <- integer(0)
  ends <- integer(0)
  for (pat in sites) {
    found <- as.integer(gregexpr(enzyme_site_regex(pat), subject, perl = TRUE)[[1]])
    found <- found[found > 0L]
    starts <- c(starts, found)
    ends <- c(ends, found + nchar(pat) - 1L)
  }
  data.frame(start = starts, end = ends)
}

# Which sequence letters each letter of a site matches. A site letter matches the letters whose
# bases it includes: R = A or G matches A, G and also a sequence letter R; N matches every
# letter; a plain A matches only A (a sequence letter R or N is not a plain A). This is the
# rule Biostrings::matchPattern(fixed = "subject") followed, measured letter by letter, which
# the previous implementation used.
iupac_match_sets <- c(A = "A", C = "C", G = "G", T = "T",
                      R = "AGR", Y = "CTY", S = "CGS", W = "ATW", K = "GTK", M = "ACM",
                      B = "CGTYSKB", D = "AGTRWKD", H = "ACTYWMH", V = "ACGRSMV",
                      N = "ACGTRYSWKMBDHVN")

# The regular expression that finds every place a site occurs, including overlapping ones (a
# look-ahead matches without using up the letters, so the next start can be found). Built
# once per site and remembered.
enzyme_site_regex <- local({
  cache <- new.env(parent = emptyenv())
  function(site) {
    hit <- cache[[site]]
    if (!is.null(hit)) return(hit)
    letters <- strsplit(toupper(site), "")[[1]]
    sets <- iupac_match_sets[letters]
    if (anyNA(sets)) {
      stop("A restriction site may only contain IUPAC DNA letters (A C G T R Y S W K M B D H V N): ", site, call. = FALSE)
    }
    pieces <- ifelse(nchar(sets) == 1, sets, paste0("[", sets, "]"))
    regex <- paste0("(?=", paste(pieces, collapse = ""), ")")
    assign(site, regex, envir = cache)
    regex
  }
})

#' @rdname find_enzyme_site_ranges
#' @noRd
find_enzyme_site_starts <- function(dna_seq, sites) {
  find_enzyme_site_ranges(dna_seq, sites)$start
}

#' List the Synonymous Changes That Reduce the Number of Sites
#'
#' For the codons from `first_codon` to `last_codon`, tries every synonymous
#' alternative and keeps those that leave fewer restriction sites in the
#' neighbourhood than there were before (so a change that only swaps one site for
#' another is never listed). Any site that overlaps a codon lies inside the window
#' that is checked.
#'
#' @param seq_vec Character vector; the sequence, one base per element.
#' @param first_codon,last_codon Integer; the range of codon numbers to try.
#' @param sites Character vector; recognition sites and reverse complements (IUPAC allowed).
#' @param codon_to_aa Named character vector; amino acid for each codon.
#' @param aa_to_syn List; the codons for each amino acid.
#' @return A list with one entry per usable change, in codon order and then
#'   alphabetical order of the new codon, each a list with `codon_index`,
#'   `old_codon` and `new_codon`.
#' @noRd
find_site_breaking_changes <- function(seq_vec, first_codon, last_codon, sites, codon_to_aa, aa_to_syn) {
  n <- length(seq_vec)
  max_site_len <- max(nchar(sites))
  changes <- list()
  for (ci in first_codon:last_codon) {
    p <- (ci - 1) * 3 + 1
    if (p + 2 > n) next
    old_codon <- paste(seq_vec[p:(p + 2)], collapse = "")
    aa <- unname(codon_to_aa[old_codon])
    if (is.na(aa)) next
    alternatives <- sort(setdiff(aa_to_syn[[aa]], old_codon))
    if (length(alternatives) == 0) next

    win_lo <- max(1, p - max_site_len + 1)
    win_hi <- min(n, p + 2 + max_site_len - 1)
    sites_before <- nrow(find_enzyme_site_ranges(paste(seq_vec[win_lo:win_hi], collapse = ""), sites))

    for (alt in alternatives) {
      trial <- seq_vec[win_lo:win_hi]
      trial[(p - win_lo + 1):(p - win_lo + 3)] <- unlist(strsplit(alt, ""))
      if (nrow(find_enzyme_site_ranges(paste(trial, collapse = ""), sites)) < sites_before) {
        changes[[length(changes) + 1]] <- list(codon_index = ci, old_codon = old_codon, new_codon = alt)
      }
    }
  }
  changes
}

#' Warn About Restriction Sites That Could Not Be Removed
#'
#' @param remaining Data frame from `find_enzyme_site_ranges()` listing the sites
#'   still present in the final sequence.
#' @return Nothing; gives one warning naming how many sites remain and where they start.
#' @noRd
warn_unremoved_sites <- function(remaining) {
  warning(nrow(remaining), " restriction site(s) could not be removed by synonymous codon changes ",
          "(at position(s) ", paste(sort(unique(remaining$start)), collapse = ", "), "). ",
          "Please check the sequence before ordering it.", call. = FALSE)
}

#' Remove Restriction Sites with the Smallest CAI Loss (deterministic)
#'
#' Used for the "High expression" (CAI1) method. That method recodes every codon
#' to the best codon for its amino acid, so there is no random choice to make:
#' for each restriction site left in the sequence this function looks at every
#' codon that overlaps the site and every synonymous alternative, keeps only the
#' changes that reduce the number of sites in that neighbourhood, and applies the
#' one that costs the least Codon Adaptation Index. The cost of replacing codon
#' `old` by `new` is `log(w_new) - log(w_old)`, so the best change has the highest
#' value (a codon with weight 0 is treated as a very small weight). Ties go to
#' the earliest codon, then to the alphabetically first codon, so the same input
#' always gives the same output.
#'
#' @param dna_seq Character string; upper-case coding sequence (A, C, G, T, N).
#' @param sites Character vector; recognition sites and their reverse
#'   complements, which may contain IUPAC codes.
#' @param score_column Character string; the graded codon-weight column in
#'   `codon_freqs` (best codon = 1) used to score changes.
#' @param max_fixes Integer; safety limit on the number of codon changes.
#' @return The sequence with removable sites removed. A warning is given if any
#'   site could not be removed.
#' @noRd
remove_enzyme_sites_max_cai <- function(dna_seq, sites, score_column = "BringMansWeights",
                                        max_fixes = 1000) {
  cf <- WormTransgeneBuilder::codon_freqs
  codons_up <- toupper(cf$Codon)
  codon_to_aa <- setNames(cf$Amino, codons_up)
  weights <- setNames(cf[[score_column]], codons_up)
  aa_to_syn <- split(codons_up, cf$Amino)
  log_w <- function(codon) log(pmax(weights[codon], 1e-6))

  seq_vec <- unlist(strsplit(dna_seq, ""))
  unfixable <- character(0)  # "start-end" keys of sites no change could break

  for (fix in seq_len(max_fixes)) {
    current <- paste(seq_vec, collapse = "")
    hits <- find_enzyme_site_ranges(current, sites)
    if (nrow(hits) == 0) break
    hits <- hits[order(hits$start, hits$end), , drop = FALSE]
    keys <- paste(hits$start, hits$end, sep = "-")
    todo <- which(!keys %in% unfixable)
    if (length(todo) == 0) break

    h <- todo[1]
    first_codon <- (hits$start[h] - 1) %/% 3 + 1
    last_codon <- (hits$end[h] - 1) %/% 3 + 1

    best <- NULL  # list(score, codon_index, new_codon)
    for (ch in find_site_breaking_changes(seq_vec, first_codon, last_codon, sites, codon_to_aa, aa_to_syn)) {
      score <- unname(log_w(ch$new_codon) - log_w(ch$old_codon))
      better <- is.null(best) || score > best$score + 1e-12
      if (better) best <- list(score = score, codon_index = ch$codon_index, new_codon = ch$new_codon)
    }

    if (is.null(best)) {
      unfixable <- c(unfixable, keys[h])
      next
    }
    p <- (best$codon_index - 1) * 3 + 1
    seq_vec[p:(p + 2)] <- unlist(strsplit(best$new_codon, ""))
  }

  remaining <- find_enzyme_site_ranges(paste(seq_vec, collapse = ""), sites)
  if (nrow(remaining) > 0) warn_unremoved_sites(remaining)
  paste(seq_vec, collapse = "")
}

#' Remove Restriction Enzyme Sites
#'
#' Iteratively scans the DNA sequence for specified restriction enzyme recognition
#' sites (checking both forward and reverse complement strands for non-palindromic sites).
#' Recognition sites may contain IUPAC ambiguity codes (for example `GGYRCC` for BanI).
#' Mutates identified sites using synonymous codons to preserve the amino acid sequence.
#'
#' @param dna_seq Character string; the coding DNA sequence to process.
#' @param enzymes_to_remove Character vector; names of the restriction enzymes to screen for (e.g., "BsaI").
#' @param usage_column_name Character string; the column name in the internal codon table used to weight synonymous mutations.
#'   The special value `"CAI1"` (the "High expression" method) switches to a deterministic search that makes
#'   the change with the smallest loss of Codon Adaptation Index.
#' @return A character string of the sequence with targeted restriction sites removed.
#' @export
remove_enzyme_sites <- function(dna_seq, enzymes_to_remove, usage_column_name) {
  dna_seq <- gsub("[^ATGCNatgcn]", "", dna_seq)

  if (length(enzymes_to_remove) == 0 || is.null(dna_seq) || dna_seq == "") return(dna_seq)

  unknown <- setdiff(enzymes_to_remove, rownames(WormTransgeneBuilder::enzy))
  if (length(unknown) > 0) {
    warning("Ignoring restriction enzyme(s) that are not in the enzyme table: ",
            paste(unknown, collapse = ", "), ".", call. = FALSE)
  }

  sites_to_check <- c()
  for (enz in enzymes_to_remove) {
    site <- toupper(WormTransgeneBuilder::enzy[enz, "Site"])
    if (!is.na(site)) {
      if (!grepl("^[ACGTRYSWKMBDHVN]+$", site)) {
        stop("Recognition site for ", enz, " contains characters that are not valid IUPAC DNA codes: ", site)
      }
      sites_to_check <- c(sites_to_check, site)
      rc_site <- reverse_complement(site)
      sites_to_check <- c(sites_to_check, rc_site)
    }
  }
  sites_to_check <- unique(sites_to_check)

  # "High expression" (CAI1) recodes to the single best codon, so a random
  # weighted draw is impossible (all other weights are 0). Use the deterministic,
  # smallest-CAI-loss search instead.
  if (identical(usage_column_name, "CAI1") && length(sites_to_check) > 0) {
    return(remove_enzyme_sites_max_cai(toupper(dna_seq), sites_to_check, "BringMansWeights"))
  }

  remove_enzyme_sites_weighted(toupper(dna_seq), sites_to_check, usage_column_name)
}

#' Remove Restriction Sites by Weighted Random Synonymous Changes
#'
#' The method used for every codon-usage choice except High expression. It scans for
#' sites, and for each one replaces the first codon the site touches with a synonymous
#' codon drawn at random, weighted by the usage column, then scans again (up to 100
#' passes). If that first codon cannot change (Met and Trp have no synonyms), it falls
#' back to the other codons the site covers: among the changes that leave fewer sites
#' in the neighbourhood, one is drawn at random with the same weights. The fallback is
#' used only when the first codon cannot change, so every case that worked before gives
#' exactly the same sequence as before. A warning is given if any site remains.
#'
#' @param dna_seq Character string; upper-case coding sequence (A, C, G, T, N).
#' @param sites_to_check Character vector; recognition sites and reverse complements.
#' @param usage_column_name Character string; the codon-table column used as weights.
#' @return The sequence with removable sites removed.
#' @noRd
remove_enzyme_sites_weighted <- function(dna_seq, sites_to_check, usage_column_name) {
  cf <- WormTransgeneBuilder::codon_freqs
  codons_up <- toupper(cf$Codon)
  codon_to_aa <- setNames(cf$Amino, codons_up)
  aa_to_syn <- split(codons_up, cf$Amino)
  usage_weights <- setNames(cf[[usage_column_name]], codons_up)

  iteration <- 0
  max_iterations <- 100
  seq_vec <- unlist(strsplit(dna_seq, ""))

  while (iteration < max_iterations) {
    iteration <- iteration + 1
    current_seq_str <- paste(seq_vec, collapse = "")

    hits <- find_enzyme_site_ranges(current_seq_str, sites_to_check)
    if (nrow(hits) == 0) break
    hits$snapped <- hits$start - (hits$start %% 3) + 1
    stpos <- unique(hits$snapped)

    for (pos in stpos) {
      if (pos > (length(seq_vec) - 2)) next
      old_codon <- paste(seq_vec[pos:(pos + 2)], collapse = "")

      aa <- WormTransgeneBuilder::codon_freqs$Amino[toupper(WormTransgeneBuilder::codon_freqs$Codon) == old_codon][1]

      changed <- FALSE
      if (!is.na(aa)) {
        new_codon <- sample_new_codon(aa, old_codon, usage_column_name)
        seq_vec[pos:(pos + 2)] <- unlist(strsplit(new_codon, ""))
        changed <- new_codon != old_codon
      }

      if (!changed) {
        # The first codon of the site cannot change, so try the codons the site covers.
        covered <- hits[hits$snapped == pos, , drop = FALSE]
        options <- find_site_breaking_changes(
          seq_vec, min((covered$start - 1) %/% 3 + 1), max((covered$end - 1) %/% 3 + 1),
          sites_to_check, codon_to_aa, aa_to_syn
        )
        if (length(options) > 0) {
          w <- unname(usage_weights[vapply(options, function(o) o$new_codon, "")])
          w[is.na(w)] <- 0
          pick <- if (sum(w) > 0) sample.int(length(options), 1, prob = w) else sample.int(length(options), 1)
          p <- (options[[pick]]$codon_index - 1) * 3 + 1
          seq_vec[p:(p + 2)] <- unlist(strsplit(options[[pick]]$new_codon, ""))
        }
      }
    }

    # If a whole pass changed nothing (every codon that had to change is Met, Trp or
    # unknown and no neighbouring change helps), further passes would be identical, so stop.
    if (identical(paste(seq_vec, collapse = ""), current_seq_str)) break
  }

  final_seq <- paste(seq_vec, collapse = "")
  remaining <- find_enzyme_site_ranges(final_seq, sites_to_check)
  if (nrow(remaining) > 0) warn_unremoved_sites(remaining)
  return(final_seq)
}
