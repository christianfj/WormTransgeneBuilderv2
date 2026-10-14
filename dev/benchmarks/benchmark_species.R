# How long does the piRNA screen take for each choice of species? (2026-10, after the species selection was added.)
# Run from the project root:   Rscript dev/benchmarks/benchmark_species.R
# For each gene size and each selection: get_pirna_distances(), the random ("Ubiq") shield, the "High expression" (CAI1) shield,
# and the whole pipeline with the shield on. Best of 3 runs, in seconds. The genes are the benchmark genes (they contain planted
# C. elegans piRNA targets taken from the database as it was before the rebuild), so the shields have work to do.
suppressMessages(devtools::load_all(quiet = TRUE))
source("dev/benchmarks/benchmark_helpers.R")
best <- function(expr, n = 3) { expr <- substitute(expr); min(vapply(seq_len(n), function(i) system.time(eval(expr, parent.frame()))[["elapsed"]], 1)) }
selections <- list(
  "old database (17,849, as before)" = benchmark_pirna_targets(),
  "C. elegans (default)" = select_pirna_database("elegans")$Target_20mer,
  "C. elegans + C. briggsae" = select_pirna_database(c("elegans", "briggsae"))$Target_20mer,
  "all four species" = select_pirna_database(c("elegans", "briggsae", "remanei", "brenneri"))$Target_20mer)
rows <- list()
for (n in c(3000, 20000)) {
  g <- make_benchmark_cds(n); bg <- make_best_codon_benchmark_cds(n)
  for (nm in names(selections)) {
    p <- selections[[nm]]
    for (thr in c(3, 4)) {
      if (nm == "old database (17,849, as before)" && thr == 4) next
      t_dist <- best(get_pirna_distances(g, p))
      t_ubiq <- best({ set.seed(1); remove_pirna_sites(g, "Ubiq", thr, pirnas = p) }, 2)
      t_cai1 <- best(remove_pirna_sites(bg, "CAI1", thr, pirnas = p), 2)
      t_all <- best({ set.seed(1); optimize_transgene(g, opt_method = "Ubiq", enzymes = NULL, remove_pirnas = TRUE, target_min_hamming = thr, pirna_targets = p) }, 2)
      rows[[length(rows) + 1]] <- data.frame(bp = n, selection = nm, piRNAs = length(p), threshold = thr, distances_s = round(t_dist, 2),
        shield_random_s = round(t_ubiq, 2), shield_cai1_s = round(t_cai1, 2), pipeline_s = round(t_all, 2))
      print(rows[[length(rows)]], row.names = FALSE)
    }
  }
}
res <- do.call(rbind, rows)
utils::write.csv(res, "dev/benchmarks/results/species_timing.csv", row.names = FALSE)
