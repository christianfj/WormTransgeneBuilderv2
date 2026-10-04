# Compare the GLO functions before and after the score lookup was rewritten, back to back.
#
#   Rscript dev/benchmarks/compare_glo.R
#
# The previous versions of calc_glo_score() and optimize_glo() (which searched the named
# 54 MB dictionary with a binary search) are taken from git history, commit 3bf0298. Each
# size is run with the old and the new code alternately, several times, and the best time
# of each is reported, which is much less affected by other programs using the CPU than a
# single run. The results are also checked to be exactly identical (same seed for the random
# optimization). Run it from the project root.

suppressMessages(devtools::load_all(quiet = TRUE))
source("dev/benchmarks/benchmark_helpers.R")

# The old dictionary, data/glo_full.rda, was removed from the repository afterwards; its last
# version is in git history (commit bdee5b0). The old functions refer to it as
# WormTransgeneBuilder::glo_full, which is redirected here to the copy loaded from history.
glo_rda <- tempfile(fileext = ".rda")
system2("git", c("show", "bdee5b0:data/glo_full.rda"), stdout = glo_rda)
old <- new.env(parent = asNamespace("WormTransgeneBuilder"))
load(glo_rda, envir = old)
for (f in c("R/fct_glo_optimization.R", "R/fct_sequence_math.R")) {
  code <- system2("git", c("show", paste0("3bf0298:", f)), stdout = TRUE)
  code <- gsub("WormTransgeneBuilder::glo_full", "glo_full", code, fixed = TRUE)
  eval(parse(text = code), envir = old)
}

t1 <- function(f) system.time(f())[["elapsed"]]
best_of <- function(reps, fs) {
  times <- lapply(fs, function(f) numeric(0))
  for (r in seq_len(reps)) for (nm in names(fs)) times[[nm]] <- c(times[[nm]], t1(fs[[nm]]))
  vapply(times, min, 0)
}

cat("Loading both score tables once before timing (the old dictionary takes several seconds) ...\n")
invisible(old$calc_glo_score(make_benchmark_cds(300)))
invisible(WormTransgeneBuilder:::glo_score_table())
cat("\n")

# `--quick` skips the 20 kb gene, which takes the old code several seconds per call. (The old
# dictionary needs 1.3 GB of memory, so on a computer that is short of memory the old code is
# far slower than it should be; the new code is not affected.)
sizes <- if ("--quick" %in% commandArgs(TRUE)) c(300, 3000) else c(300, 3000, 20000)
for (n in sizes) {
  g <- make_benchmark_cds(n)
  reps <- if (n >= 20000) 2 else 5
  # calc_glo_score is fast, so time many calls at once and divide (a single call is below
  # the resolution of R's timer).
  k <- if (n >= 20000) 5 else 50
  s <- best_of(reps, list(old = function() for (i in seq_len(k)) old$calc_glo_score(g),
                          new = function() for (i in seq_len(k)) calc_glo_score(g))) / k
  cat(sprintf("%5d bp  calc_glo_score: old %8.2f ms | new %7.2f ms (%.1fx)\n", n, s[["old"]] * 1000, s[["new"]] * 1000, s[["old"]] / s[["new"]]))
  o <- best_of(reps, list(old = function() { set.seed(1); old$optimize_glo(g) }, new = function() { set.seed(1); optimize_glo(g) }))
  cat(sprintf("%5d bp  optimize_glo:   old %8.2f s  | new %7.2f s  (%.1fx)\n", n, o[["old"]], o[["new"]], o[["old"]] / o[["new"]]))
  set.seed(1); a <- old$optimize_glo(g); set.seed(1); b <- optimize_glo(g)
  stopifnot(identical(a, b), identical(old$calc_glo_score(g), calc_glo_score(g)), identical(old$calc_glo_score(a), calc_glo_score(b)))
}
cat("\nResults identical (optimized sequence and scores, same seed) at all sizes: TRUE\n")
cat("\nMemory used by the score data:\n")
cat(sprintf("  old dictionary (glo_full, with names): %.0f MB\n", as.numeric(object.size(old$glo_full)) / 1024^2))
cat(sprintf("  new score table:                       %.0f MB\n", as.numeric(object.size(WormTransgeneBuilder:::glo_score_table())) / 1024^2))
