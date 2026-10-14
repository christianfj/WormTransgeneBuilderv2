# How the GLO score depends on the number of rounds and on the stopping count, for the port of the original code and for the app (2026-10-07).
# Run from the project root; needs the original scores as orig_scores.rds (see dev/verify_glo_against_original.R); edit the two paths below.
suppressMessages(devtools::load_all(quiet = TRUE))
orig <- readRDS("/Volumes/CFJ_backup1S1TB/pirna_other_species/codon_check/orig_scores.rds")
source("dev/glo_original_port.R")
source("dev/benchmarks/benchmark_helpers.R")
a <- glo_codon_amino()
prot <- function(s) paste(a[substring(s, seq(1, nchar(s) - 2, 3), seq(3, nchar(s), 3))], collapse = "")
gfp <- toupper(paste(readLines("/Volumes/CFJ_backup1S1TB/pirna_other_species/codon_check/germline/GFP_S65C.fasta", warn = FALSE), collapse = ""))
tests <- list(GFP_238aa = prot(gfp), gene_1kb = sub("X$", "", prot(make_benchmark_cds(3000))))
cat("score per codon (mean of 4 seeds); first-guess-only value in the first column\n")
for (nm in names(tests)) {
  aa <- tests[[nm]]; n <- nchar(aa)
  first <- optimize_original(aa, 0)$first
  res <- list()
  for (rounds in c(1, 2, 4, 10)) {
    sc <- vapply(1:4, function(s) { set.seed(s); optimize_original(aa, rounds)$score }, 0)
    res[[paste0("original, ", rounds, " round(s)")]] <- mean(sc)
  }
  for (sc_ in c(0.25, 0.5)) { s_ <- vapply(1:4, function(s) { set.seed(s); optimize_original(aa, 4, iter_scale = sc_)$score }, 0); res[[sprintf("original, 4 rounds, %.2fx the stopping count", sc_)]] <- mean(s_) }
  for (rounds in c(1, 2, 4, 10)) {
    sc <- vapply(1:4, function(s) { set.seed(s); x <- optimize_glo(back_translate(aa), method = "heuristic", rounds = rounds); calc_glo_score(x) / n }, 0)
    res[[paste0("app, ", rounds, " round(s)")]] <- mean(sc)
  }
  cat(sprintf("\n%s (%d aa): first guess only %.2f\n", nm, n, first)); for (k in names(res)) cat(sprintf("  %-48s %9.3f\n", k, res[[k]]))
}
