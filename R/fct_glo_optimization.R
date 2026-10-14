#' Germline Optimize (GLO) a DNA Sequence
#'
#' The GLO (GermLine Optimized) method of Daniel Dickinson, first described and used in Fielmich et al. (2018) eLife 7:e38198
#' (https://doi.org/10.7554/eLife.38198; PubMed 30109984); the original Perl code is at https://github.com/dannyhmg/germline.
#' The score of a sequence is the sum of the score of every 12-base window that starts at a codon boundary (see `glo_window_scores()`);
#' the score of a 12-mer reflects how it is used in native germline-expressed genes. The amino-acid sequence never changes.
#'
#' Two methods are available.
#'
#' **`"exact"` (the default)** finds the coding with the highest possible score. The score is a sum of terms that each involve four
#' consecutive codons, so the best coding can be found exactly by dynamic programming: go along the sequence and, for every choice of the last
#' three codons, keep the best score reached so far, then trace back. This takes a fraction of a second even for 20 kb. Because 82% of the 12-mers
#' score zero, very many different codings share the best score; they are told apart by a second criterion that can never override the
#' score: among the best-scoring codings, the one with the highest Codon Adaptation Index (computed with the `Ubiq` codon table, as the app reports it) wins.
#' (The CAI is the geometric mean of the codon weights, so maximizing the sum of the log weights is maximizing the CAI.) Optionally
#' a small random jitter (`diversity`) lets repeated runs differ a little among near-equal choices.
#'
#' **`"heuristic"`** is a port of the original algorithm (`OptimizerTools::optimize()`), which works with whole four-codon words:
#' \enumerate{
#'   \item every synonymous coding of each group of four amino acids is scored;
#'   \item a first sequence is assembled by placing the highest-scoring words first, then the best that fit;
#'   \item it is improved by taking positions in random order, scoring every coding word there in the context of the neighbouring codons
#'     (own window counted fully, the six overlapping windows 0.75, 0.5 and 0.25 on each side) and keeping the best word if the score of the
#'     whole sequence rises;
#'   \item this is done in `rounds` (10) rounds from the same first sequence; the best result is kept.
#' }
#' It differs from the original only in its stopping rule: a position is revisited only when a neighbouring codon has changed, so the search
#' stops when no position can be improved. It reaches the scores of the unchanged original code (`docs/GLO_VERIFICATION.md`) and falls 0.1 to 2%
#' short of the exact optimum.
#'
#' Codons that cannot be recoded stay as they are: stop codons, and anything that is not a codon of A, C, G and T. A sequence of
#' fewer than four codons is returned unchanged. Bases after the last complete codon are kept.
#'
#' @param dna_seq Character string; the input DNA sequence to optimize.
#' @param method `"exact"` (default) or `"heuristic"`, see above.
#' @param rounds Integer; the number of improvement rounds of the heuristic (the original uses 10). Ignored by the exact method.
#' @param diversity Numeric; size of the random jitter added to the CAI tie-break of the exact method (the default 0 is deterministic: the
#'   highest-CAI coding among the best-scoring ones). It is on the scale of the natural logarithm of the codon usage weights and never changes the GLO score.
#' @param avoid_sites Optional character vector of DNA recognition sites (IUPAC codes allowed, both strands as separate entries, for example from
#'   `enzyme_search_sites()`). The exact method then finds the best coding **that contains none of them**: a site of up to 10 bases always lies inside
#'   one window of four codons, and any such window containing a site is ruled out (with a penalty so large that it is only accepted when no
#'   alternative exists). Longer sites, and sites that cannot be avoided, are left for the later restriction-site step. Ignored by the heuristic.
#' @return A character string containing the GLO-optimized DNA sequence (upper case).
#' @export
optimize_glo <- function(dna_seq, method = c("exact", "heuristic"), rounds = 10L, diversity = 0, avoid_sites = NULL) {
  method <- match.arg(method)
  seq_str <- toupper(dna_seq)
  n <- nchar(seq_str) %/% 3
  if (n < 4) return(dna_seq)
  leftover <- substring(seq_str, 3 * n + 1)
  codon_chr <- substring(seq_str, seq(1, 3 * n - 2, 3), seq(3, 3 * n, 3))
  cod0 <- glo_codon_index(codon_chr)

  # What each position may be: the codons of its amino acid, or itself if it is a stop codon or not a plain codon.
  amino <- glo_codon_amino()[codon_chr]
  variable <- !is.na(cod0) & !is.na(amino) & amino != "X"
  synonyms <- glo_synonym_indices()
  options <- lapply(seq_len(n), function(i) if (variable[i]) synonyms[[amino[i]]] else cod0[i])
  signature <- ifelse(variable, amino, paste0("=", cod0))

  if (method == "exact") {
    cod <- glo_exact(options, n, diversity, avoid_sites, codon_chr)
    out <- ifelse(is.na(cod), codon_chr, glo_codon_names()[cod + 1])
    return(paste0(paste(out, collapse = ""), leftover))
  }

  P <- n - 3L                                         # the number of four-codon words
  words <- vector("list", P)
  cache <- new.env(parent = emptyenv())
  for (p in seq_len(P)) {
    key <- paste(signature[p:(p + 3)], collapse = "|")
    w <- cache[[key]]
    if (is.null(w)) { w <- glo_word_table(options[p:(p + 3)]); cache[[key]] <- w }
    words[[p]] <- w
  }

  # ---- first guess: the highest-scoring words first, then the best that fit ----------------------------------------
  coding <- rep(NA_real_, n); filled <- rep(FALSE, n)
  assigned <- rep(1L, P)
  assigned_score <- vapply(words, function(w) w$s[1], 0)
  # The positions still to place are kept in a binary heap (smallest key first); the key is exact in doubles because scores are whole numbers.
  sort_key <- function(score, position) -score * 1e6 + position
  heap <- sort_key(assigned_score, seq_len(P))
  heap_size <- P
  heap <- heap[order(heap)]                           # a sorted array is a valid heap
  heap_pop <- function() {
    top <- heap[1]; heap[1] <<- heap[heap_size]; heap_size <<- heap_size - 1L
    i <- 1L
    while (TRUE) {
      l <- 2L * i; r <- l + 1L; m <- i
      if (l <= heap_size && heap[l] < heap[m]) m <- l
      if (r <= heap_size && heap[r] < heap[m]) m <- r
      if (m == i) break
      tmp <- heap[i]; heap[i] <<- heap[m]; heap[m] <<- tmp; i <- m
    }
    top
  }
  heap_push <- function(key) {
    heap_size <<- heap_size + 1L; i <- heap_size; heap[i] <<- key
    while (i > 1L && heap[i %/% 2L] > heap[i]) { j <- i %/% 2L; tmp <- heap[i]; heap[i] <<- heap[j]; heap[j] <<- tmp; i <- j }
  }
  while (heap_size > 0) {
    key <- heap_pop()
    e <- as.integer(key %% 1e6)                          # the position is the part below 1e6 (%% is never negative here)
    w <- words[[e]]$W[assigned[e], ]
    seg <- e:(e + 3)
    same <- (coding[seg] == w) | (is.na(coding[seg]) & is.na(w)); same[is.na(same)] <- FALSE
    if (all(!filled[seg] | same)) {
      coding[seg] <- w; filled[seg] <- TRUE
    } else {
      assigned[e] <- assigned[e] + 1L
      if (assigned[e] > nrow(words[[e]]$W)) stop("GLO: no coding word fits position ", e, call. = FALSE)
      assigned_score[e] <- words[[e]]$s[assigned[e]]
      heap_push(sort_key(assigned_score[e], e))
    }
  }

  # ---- improvement rounds -----------------------------------------------------------------------------------------------
  best <- NULL
  for (round in seq_len(rounds)) {
    cur <- coding
    dirty <- rep(TRUE, P)
    while (any(dirty)) {
      for (pos in sample(which(dirty))) {              # pos is 1-based here; the codon position is pos - 1
        if (!dirty[pos]) next
        dirty[pos] <- FALSE
        top <- glo_best_word(cur, pos - 1L, words[[pos]], n)
        now <- cur[pos:(pos + 3)]
        if (all((top == now) | (is.na(top) & is.na(now)), na.rm = FALSE) %in% TRUE) next
        trial <- cur; trial[pos:(pos + 3)] <- top
        starts <- max(0L, pos - 4L):min(n - 4L, pos + 2L)        # windows (0-based starts) that contain a changed codon
        if (sum(glo_window_s100(trial, starts)) - sum(glo_window_s100(cur, starts)) > 0) {
          cur <- trial
          dirty[max(1L, pos - 6L):min(P, pos + 6L)] <- TRUE
        }
      }
    }
    total <- sum(glo_window_s100(cur, 0:(n - 4L)))
    if (is.null(best) || total > best$total) best <- list(cod = cur, total = total)
  }

  out <- ifelse(is.na(best$cod), codon_chr, glo_codon_names()[best$cod + 1])
  paste0(paste(out, collapse = ""), leftover)
}

# ---- helpers ------------------------------------------------------------------------------------------------------------

# The 64 codons in the order of their base-4 value (AAA = 0, AAC = 1, ..., TTT = 63): the same order as the digits of a 12-mer's position in the
# score table, so four codon indices combine into a window's table position by reading them as base-64 digits.
glo_codon_names <- local({
  names_ <- NULL
  function() {
    if (is.null(names_)) {
      b <- c("A", "C", "G", "T")
      g <- expand.grid(d3 = 0:3, d2 = 0:3, d1 = 0:3)
      names_ <<- paste0(b[g$d1 + 1], b[g$d2 + 1], b[g$d3 + 1])
    }
    names_
  }
})

glo_codon_index <- function(chr) match(chr, glo_codon_names()) - 1

# Amino acid (one letter, "X" for stop) of each codon, named by the upper-case codon.
glo_codon_amino <- function() {
  cf <- WormTransgeneBuilder::codon_freqs
  stats::setNames(cf$Amino, toupper(cf$Codon))
}

# Codon indices of each amino acid's synonymous codons.
glo_synonym_indices <- function() {
  a <- glo_codon_amino()
  lapply(split(names(a), unname(a)), glo_codon_index)
}

# Scores as whole numbers (hundredths): every score in the table is a multiple of 0.01, so sums and comparisons are exact.
# `idx` is a window's position in the table (0-based); NA (a window with a base that is not A/C/G/T) gets the median score.
glo_s100 <- function(idx) {
  v <- round(glo_score_table()[idx + 1] * 100)
  v[is.na(v)] <- round(WormTransgeneBuilder::glo_median_score * 100)
  v
}

# Scores of the 12-base windows that start at codons `starts` (0-based) of the codon-index vector `cod`.
glo_window_s100 <- function(cod, starts) {
  glo_s100(cod[starts + 1] * 262144 + cod[starts + 2] * 4096 + cod[starts + 3] * 64 + cod[starts + 4])
}

# Every coding word for four positions, given the codon indices allowed at each, best first. Returns the words (W, one row each), their
# scores (s, hundredths) and the pieces needed to score the neighbouring windows quickly: p1..p3 are the first 1..3 codons read as base-64
# digits, s1..s3 the last 1..3.
glo_word_table <- function(opts) {
  n <- lengths(opts); total <- prod(n); k <- seq_len(total) - 1
  c1 <- opts[[1]][(k %/% (n[2] * n[3] * n[4])) %% n[1] + 1]
  c2 <- opts[[2]][(k %/% (n[3] * n[4])) %% n[2] + 1]
  c3 <- opts[[3]][(k %/% n[4]) %% n[3] + 1]
  c4 <- opts[[4]][k %% n[4] + 1]
  s <- glo_s100(c1 * 262144 + c2 * 4096 + c3 * 64 + c4)
  o <- order(-s, seq_len(total))
  c1 <- c1[o]; c2 <- c2[o]; c3 <- c3[o]; c4 <- c4[o]
  list(W = cbind(c1, c2, c3, c4), s = s[o],
       p1 = c1, p2 = c1 * 64 + c2, p3 = (c1 * 64 + c2) * 64 + c3,
       s1 = c4, s2 = c3 * 64 + c4, s3 = (c2 * 64 + c3) * 64 + c4)
}

# The best word at 0-based position `pos` given the codons around it: each word is scored as its own window (weight 4) plus the six
# windows that overlap it (weights 3, 2, 1 on each side; the original's 0.75, 0.5, 0.25 scaled by 4 to stay whole numbers).
# Windows that would run off either end of the sequence do not exist and count nothing.
glo_best_word <- function(cur, pos, wt, n) {
  sc <- 4 * wt$s
  for (m in 1:3) {                                     # windows that start m codons before the word
    if (pos - m >= 0) {
      left <- Reduce(function(a, b) a * 64 + b, cur[(pos - m + 1):pos])
      sc <- sc + (4 - m) * glo_s100(left * 64^(4 - m) + wt[[paste0("p", 4 - m)]])
    }
  }
  for (m in 1:3) {                                     # windows that start m codons after the word's first codon
    if (pos + m + 3 <= n - 1) {
      right <- Reduce(function(a, b) a * 64 + b, cur[(pos + 5):(pos + 4 + m)])
      sc <- sc + (4 - m) * glo_s100(wt[[paste0("s", 4 - m)]] * 64^m + right)
    }
  }
  wt$W[which.max(sc), ]
}

# The exact optimum by dynamic programming. `options[[i]]` are the codon indices allowed at position i. Returns the chosen codon index at every
# position (NA where the codon is not a plain codon). The score is the sum, over i, of the score of the window of codons i..i+3, so after choosing
# codon i+3 only the three codons before it matter for what comes next: V(b, c, d) = max over a of V(a, b, c) + score(a, b, c, d).
# A second, much smaller term breaks ties (never larger than a quarter of one hundredth of a score point in total): the log of the usage weight of
# each chosen codon (the sum of these is n times the log of the CAI), plus optional random jitter, so among equally scoring codings the one with the
# highest CAI wins.
glo_exact <- function(options, n, diversity, avoid_sites = NULL, codon_chr = NULL) {
  cf <- WormTransgeneBuilder::codon_freqs
  usage <- numeric(64); usage[glo_codon_index(toupper(cf$Codon)) + 1] <- log(pmax(cf$Ubiq, 1e-6))
  tie <- lapply(options, function(o) if (length(o) > 1) usage[o + 1] + stats::rnorm(length(o), 0, diversity) else rep(0, length(o)))
  span <- sum(vapply(tie, function(x) max(abs(x)), 0))
  eps <- if (span > 0) 0.25 / span else 0
  tie <- lapply(tie, function(x) x * eps)

  # Windows of four codons that contain a recognition site are ruled out. The strings of a window depend only on which amino acids (or fixed codons)
  # it covers, so the answer is remembered by that signature.
  avoid_regex <- if (length(avoid_sites) > 0) paste(vapply(unique(avoid_sites), enzyme_site_regex, ""), collapse = "|") else NULL
  names_of <- function(i) if (anyNA(options[[i]])) codon_chr[i] else glo_codon_names()[options[[i]] + 1]
  opt_names <- if (is.null(avoid_regex)) NULL else lapply(seq_len(n), names_of)
  signature <- if (is.null(avoid_regex)) NULL else vapply(seq_len(n), function(i) paste(opt_names[[i]], collapse = ","), "")
  forbidden_cache <- new.env(parent = emptyenv())
  forbidden <- function(i) {                                          # logical array (a, b, c, d) for the window of codons i-3..i
    key <- paste(signature[(i - 3):i], collapse = "|")
    f <- forbidden_cache[[key]]
    if (is.null(f)) {
      nm <- opt_names[(i - 3):i]
      strings <- outer(outer(outer(nm[[1]], nm[[2]], paste0), nm[[3]], paste0), nm[[4]], paste0)
      f <- array(grepl(avoid_regex, strings, perl = TRUE), dim = lengths(nm))
      forbidden_cache[[key]] <- f
    }
    f
  }
  penalty <- 1e9                                                      # in hundredths of a score point: far more than any window can score

  V <- array(outer(outer(tie[[1]], tie[[2]], "+"), tie[[3]], "+"), dim = lengths(options[1:3]))
  back <- vector("list", n)
  for (i in 4:n) {
    o <- options[(i - 3):i]; dm <- lengths(o)
    g <- outer(outer(outer(o[[1]] * 262144, o[[2]] * 4096, "+"), o[[3]] * 64, "+"), o[[4]], "+")
    G <- array(glo_s100(g), dim = dm)                                  # score of the window a, b, c, d
    G <- G + rep(tie[[i]], each = prod(dm[1:3]))                       # tie-break term of the new codon d
    if (!is.null(avoid_regex)) G <- G - penalty * forbidden(i)
    best <- array(-Inf, dim = dm[2:4]); arg <- array(1L, dim = dm[2:4])
    for (a in seq_len(dm[1])) {
      cand <- array(G[a, , , ], dim = dm[2:4]) + array(V[a, , ], dim = dm[2:4])[, , rep(1L, dm[4]), drop = FALSE]
      better <- cand > best
      best[better] <- cand[better]; arg[better] <- a
    }
    V <- best; back[[i]] <- arg
  }
  k <- arrayInd(which.max(V), dim(V))
  pick <- integer(n); pick[n] <- k[3]; pick[n - 1] <- k[2]; pick[n - 2] <- k[1]
  for (i in n:4) pick[i - 3] <- back[[i]][pick[i - 2], pick[i - 1], pick[i]]
  vapply(seq_len(n), function(i) options[[i]][pick[i]], 0)
}

# The change in GLO score (in hundredths) if codon `ci` (1-based) of the sequence `seq_vec` (a vector of single bases) is replaced by `new_codon`:
# only the up to four 12-base windows that contain that codon change. A tiny extra term (never more than a thousandth of a score point) prefers
# the more common codon when two changes score the same. Used to choose the least harmful synonymous change when a restriction site or a piRNA match
# has to be removed from a GLO-optimized sequence.
glo_change_gain <- function(seq_vec, ci, new_codon) {
  n <- length(seq_vec) %/% 3
  lo <- max(1L, ci - 3L); hi <- min(n, ci + 3L)
  old_chr <- vapply(lo:hi, function(k) paste(seq_vec[((k - 1) * 3 + 1):(k * 3)], collapse = ""), "")
  cod <- glo_codon_index(old_chr)
  new <- cod; new[ci - lo + 1] <- glo_codon_index(new_codon)
  starts <- (max(lo, ci - 3L):min(ci, hi - 3L)) - lo               # 0-based starts of the windows that contain codon ci
  gain <- if (length(starts) == 0) 0 else sum(glo_window_s100(new, starts)) - sum(glo_window_s100(cod, starts))
  cf <- WormTransgeneBuilder::codon_freqs
  w <- stats::setNames(cf$Ubiq, toupper(cf$Codon))
  gain + 1e-4 * unname(log(pmax(w[new_codon], 1e-6)) - log(pmax(w[old_chr[ci - lo + 1]], 1e-6)))
}

