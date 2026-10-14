# R port of the ORIGINAL GLO code (Seqscore.pm and OptimizerTools.pm of https://github.com/dannyhmg/germline = dickinson-lab/germline), used to
# check the app (dev/verify_glo_against_original.R) and to study how the result depends on the number of rounds and iterations.
# Needs `orig` (the 4^12 original scores as a numeric vector, alphabetical order) to exist before it is sourced.
# ---- the original code, ported -----------------------------------------------------------------------------------------
codon_names <- c(outer(c("A","C","G","T"), outer(c("A","C","G","T"), c("A","C","G","T"), paste0), paste0))   # not the base-4 order; see below
codon_chr <- function(i) { b <- c("A","C","G","T"); paste0(b[i %/% 16 + 1], b[(i %/% 4) %% 4 + 1], b[i %% 4 + 1]) }
codon_idx <- function(chr) vapply(strsplit(chr, ""), function(x) sum((match(x, c("A","C","G","T")) - 1) * c(16, 4, 1)), 0)   # base-4 index 0..63
orig_codons <- list(A = c("GCA","GCC","GCG","GCT"), C = c("TGC","TGT"), D = c("GAC","GAT"), E = c("GAA","GAG"), F = c("TTC","TTT"),
  G = c("GGA","GGC","GGG","GGT"), H = c("CAC","CAT"), I = c("ATA","ATC","ATT"), K = c("AAA","AAG"), L = c("CTA","CTC","CTG","CTT","TTA","TTG"),
  M = "ATG", N = c("AAC","AAT"), P = c("CCA","CCC","CCG","CCT"), Q = c("CAA","CAG"), R = c("AGA","AGG","CGA","CGC","CGG","CGT"),
  S = c("AGC","AGT","TCA","TCC","TCG","TCT"), T = c("ACA","ACC","ACG","ACT"), V = c("GTA","GTC","GTG","GTT"), W = "TGG", Y = c("TAC","TAT"))
orig_codons <- lapply(orig_codons, codon_idx)
win <- function(cod, start) orig[1 + cod[start] * 262144 + cod[start + 1] * 4096 + cod[start + 2] * 64 + cod[start + 3]]   # start: 1-based vector, in range

# Seqscore::score_word(position, sequence): `position` is 0-based; weights 0.25, 0.5, 0.75, 1, 0.75, 0.5, 0.25 for words starting at position-3 ... position+3
score_word <- function(position, cod) {
  n <- length(cod); total <- 0
  for (k in -3:3) { s <- position + k; if (s >= 0 && s + 3 <= n - 1) total <- total + (1 - abs(k) / 4) * win(cod, s + 1) }
  total
}
score_sequence <- function(cod) { n <- length(cod); sum(vapply(seq(0, n, by = 4), score_word, 0, cod = cod)) / n }
app_score <- function(cod) { n <- length(cod); sum(win(cod, seq_len(n - 3))) / n }

# Optimizer::optimize(sequence_lib, AA)
optimize_original <- function(aa, rounds, iter_scale = 1) {
  aa <- strsplit(aa, "")[[1]]; L <- length(aa)
  words <- lapply(seq_len(L - 3), function(c0) { g <- expand.grid(rev(lapply(aa[c0:(c0 + 3)], function(a) orig_codons[[a]])), KEEP.OUT.ATTRS = FALSE); as.matrix(g[, 4:1]) })   # rows = coding 12-mers (4 codon indices)
  wscore <- lapply(words, function(w) orig[1 + w[, 1] * 262144 + w[, 2] * 4096 + w[, 3] * 64 + w[, 4]])
  ord <- lapply(wscore, function(s) order(-s, seq_along(s)))
  assigned <- rep(1L, L - 3); assigned_score <- vapply(seq_along(ord), function(i) wscore[[i]][ord[[i]][1]], 0)
  top_score <- assigned_score
  coding <- rep(NA_real_, L)
  pending <- order(-top_score, seq_along(top_score))
  while (length(pending) > 0) {
    e <- pending[1]; pending <- pending[-1]
    w <- words[[e]][ord[[e]][assigned[e]], ]
    cur <- coding[e:(e + 3)]
    if (all(is.na(cur)) || all(is.na(cur) | cur == w)) coding[e:(e + 3)] <- w
    else {
      assigned[e] <- assigned[e] + 1L
      assigned_score[e] <- wscore[[e]][ord[[e]][assigned[e]]]
      pending <- c(pending, e); pending <- pending[order(-assigned_score[pending], pending)]
    }
  }
  first <- score_sequence(coding)
  maxiter <- iter_scale * -4.6052 / log(1 - 1 / (L - 3))
  best <- NULL
  for (h in seq_len(rounds)) {
    cur <- coding; cur_score <- first; counter <- 0
    while (counter <= maxiter) {
      pos <- sample(0:(L - 4), 1)                                 # 0-based, as random_uniform_integer(1, 0, L-4)
      lo <- if (pos >= 3) pos - 3 else 0; hi <- if (pos <= L - 7) pos + 6 else L - 1
      rel <- if (pos >= 3) 3 else pos
      W <- words[[pos + 1]][ord[[pos + 1]], , drop = FALSE]
      ctx <- cur[(lo + 1):(hi + 1)]; n <- length(ctx)
      sc <- numeric(nrow(W))
      for (k in -3:3) { s <- rel + k; if (s >= 0 && s + 3 <= n - 1) {
        wt <- 1 - abs(k) / 4
        # window starting at codon s of ctx; the new word occupies codons rel..rel+3 of ctx
        cols <- lapply(0:3, function(j) { p <- s + j; if (p >= rel && p <= rel + 3) W[, p - rel + 1] else rep(ctx[p + 1], nrow(W)) })
        sc <- sc + wt * orig[1 + cols[[1]] * 262144 + cols[[2]] * 4096 + cols[[3]] * 64 + cols[[4]]] } }
      top <- which.max(sc)                                         # first maximum, like a stable sort descending
      topw <- W[top, ]
      if (all(topw == cur[(pos + 1):(pos + 4)])) { counter <- counter + 1; next }
      test <- cur; test[(pos + 1):(pos + 4)] <- topw
      ts <- score_sequence(test)
      if (ts > cur_score) { cur <- test; cur_score <- ts; counter <- 0 } else counter <- counter + 1
    }
    if (is.null(best) || cur_score > best$score) best <- list(cod = cur, score = cur_score)
  }
  list(cod = best$cod, score = best$score, first = first)
}


protein_of_codons <- function(cod) { map <- setNames(rep(names(orig_codons), lengths(orig_codons)), unlist(orig_codons)); paste(map[as.character(cod)], collapse = "") }
