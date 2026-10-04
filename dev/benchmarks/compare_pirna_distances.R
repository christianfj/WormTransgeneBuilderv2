# Compare the old and the new piRNA distance computation, back to back.
#
#   Rscript dev/benchmarks/compare_pirna_distances.R
#
# Whole-suite benchmark timings are noisy whenever something else is using the CPU. This
# script alternates the three methods (the original stringdist comparison, the new engine
# on one core, and the new engine on its default number of cores), repeats each, and reports
# the best time of each, which is much less affected by background load. It also checks that
# all three give exactly the same distances. Run it from the project root.
suppressMessages(devtools::load_all(quiet=TRUE)); source("dev/benchmarks/benchmark_helpers.R")
old_distances <- function(seq) { s <- toupper(gsub("[^ATGCNatgcn]", "", seq)); st <- seq_len(nchar(s)-19); k <- substring(s, st, st+19)
  m <- stringdist::stringdistmatrix(k, WormTransgeneBuilder::pirna_sequences, method="hamming", nthread = max(1, parallel::detectCores()-1))
  r <- suppressWarnings(apply(m, 2, min, na.rm=TRUE)); r[is.infinite(r)] <- NA; r }
t1 <- function(f) system.time(f())[["elapsed"]]
cat("load average before:", system("uptime | sed 's/.*load averages*: //'", intern=TRUE), "\n\n")
for (n in c(300, 3000, 20000)) {
  g <- make_benchmark_cds(n); starts <- seq_len(nchar(g)-19); kmers <- substring(g, starts, starts+19)
  reps <- if (n >= 20000) 2 else 5
  res <- list(old = numeric(0), one = numeric(0), multi = numeric(0))
  for (r in seq_len(reps)) {   # interleaved so any background load hits all three equally
    res$old   <- c(res$old,   t1(function() old_distances(g)))
    res$one   <- c(res$one,   t1(function() scan_kmers_against_pirnas(kmers, cores = 1)))
    res$multi <- c(res$multi, t1(function() get_pirna_distances(g)))
  }
  cat(sprintf("%5d bp (%d windows): OLD stringdist %.2f s | NEW 1 core %.2f s (%.1fx) | NEW multi-core (default: %d cores) %.2f s (%.1fx)   [best of %d, interleaved]\n",
      n, length(kmers), min(res$old), min(res$one), min(res$old)/min(res$one), pirna_scan_cores(length(kmers)), min(res$multi), min(res$old)/min(res$multi), reps))
  stopifnot(identical(old_distances(g), get_pirna_distances(g)))
}
cat("\nidentical results old vs new on all three genes: TRUE\nload average after:", system("uptime | sed 's/.*load averages*: //'", intern=TRUE), "\n")
