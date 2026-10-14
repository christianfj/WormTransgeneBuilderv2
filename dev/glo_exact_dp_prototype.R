# PROTOTYPE (2026-10-07, not used by the app): the exact optimum of the GLO score by dynamic programming, compared with optimize_glo().
# The score is a sum of terms that each involve 4 consecutive codons, so the best coding can be found exactly, position by position, keeping the best score for
# every choice of the last three codons. See docs/GLO_VERIFICATION.md. Run from the project root after setting the paths below.
suppressMessages(devtools::load_all(quiet = TRUE)); source("dev/benchmarks/benchmark_helpers.R")
glo_dp <- function(dna) {
  n <- nchar(dna) %/% 3; codon_chr <- substring(dna, seq(1, 3*n-2, 3), seq(3, 3*n, 3))
  cod0 <- glo_codon_index(codon_chr); amino <- glo_codon_amino()[codon_chr]
  variable <- !is.na(cod0) & !is.na(amino) & amino != "X"; syn <- glo_synonym_indices()
  opts <- lapply(seq_len(n), function(i) if (variable[i]) syn[[amino[i]]] else cod0[i])
  # V over (c1,c2,c3); then for i = 4..n: V'(b,c,d) = max_a V(a,b,c) + s(a,b,c,d)
  V <- array(0, dim = lengths(opts[1:3])); back <- vector("list", n)
  for (i in 4:n) {
    o <- opts[(i-3):i]; dm <- lengths(o)
    G <- array(glo_s100(outer(outer(outer(o[[1]]*262144, o[[2]]*4096, "+"), o[[3]]*64, "+"), o[[4]], "+")), dim = dm)   # a,b,c,d
    best <- array(-Inf, dim = dm[2:4]); arg <- array(1L, dim = dm[2:4])
    for (a in seq_len(dm[1])) {
      cand <- G[a, , , , drop = TRUE]; dim(cand) <- dm[2:4]
      cand <- cand + array(V[a, , ], dim = dm[2:4])  # V[a,b,c] broadcast over d
      better <- cand > best; best[better] <- cand[better]; arg[better] <- a
    }
    V <- best; back[[i]] <- arg
  }
  # traceback
  k <- arrayInd(which.max(V), dim(V)); b <- k[1]; c <- k[2]; d <- k[3]
  pick <- integer(n); pick[n] <- d; pick[n-1] <- c; pick[n-2] <- b
  for (i in n:4) { a <- back[[i]][pick[i-2], pick[i-1], pick[i]]; pick[i-3] <- a }
  cod <- vapply(seq_len(n), function(i) opts[[i]][pick[i]], 0)
  list(seq = paste(ifelse(is.na(cod), codon_chr, glo_codon_names()[cod+1]), collapse = ""), total = max(V))
}
a <- glo_codon_amino(); prot <- function(s) paste(a[substring(s, seq(1,nchar(s)-2,3), seq(3,nchar(s),3))], collapse = "")
gfp <- toupper(paste(readLines("/Volumes/CFJ_backup1S1TB/pirna_other_species/codon_check/germline/GFP_S65C.fasta", warn=FALSE), collapse=""))
cases <- list(`GFP 238 aa` = gfp, `benchmark 300 bp` = make_benchmark_cds(300), `benchmark 3 kb` = make_benchmark_cds(3000), `benchmark 20 kb` = make_benchmark_cds(20000))
for (nm in names(cases)) { g <- cases[[nm]]; k <- nchar(g)/3
  t1 <- system.time(d <- glo_dp(g))[["elapsed"]]; stopifnot(identical(prot(d$seq), prot(g)))
  set.seed(1); t2 <- system.time(h <- optimize_glo(g))[["elapsed"]]
  cat(sprintf("%-17s exact optimum (dynamic programming) %9.3f per codon in %5.2f s | heuristic %9.3f in %5.2f s | difference %.4f%%\n", nm, calc_glo_score(d$seq)/k, t1, calc_glo_score(h)/k, t2, 100*(calc_glo_score(d$seq)-calc_glo_score(h))/calc_glo_score(d$seq))) }
# ---- brute-force check on small proteins: DP total equals the best of ALL codings, and equals the score of the sequence it returns
set.seed(5); syn <- glo_synonym_indices(); aas <- names(syn)[lengths(syn) > 1 & lengths(syn) <= 4]
ok <- TRUE
for (rep in 1:12) {
  p <- sample(aas, 7, TRUE); dna <- paste(vapply(p, function(x) glo_codon_names()[syn[[x]][1] + 1], ""), collapse = "")
  all <- expand.grid(lapply(p, function(x) syn[[x]]))
  tot <- apply(all, 1, function(cd) sum(glo_window_s100(cd, 0:3)))
  d <- glo_dp(dna); ds <- sum(glo_window_s100(glo_codon_index(substring(d$seq, seq(1, 19, 3), seq(3, 21, 3))), 0:3))
  ok <- ok && d$total == max(tot) && ds == max(tot)
}
cat("brute force (12 random 7-codon proteins, all codings enumerated): DP is exactly optimal:", ok, "\n")
