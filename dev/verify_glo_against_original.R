# Compares this app's GLO (GermLine Optimized) code with the original implementation by D. J. Dickinson
# (Fielmich et al. 2018, eLife 7:e38198; source code https://github.com/dannyhmg/germline, Perl, MIT licence).
# Run from the project root (needs the original score database; see below):
#   Rscript dev/verify_glo_against_original.R <orig_scores.rds> <GFP_S65C.fasta of the original repository> [n_rounds]
#
# `orig_scores.rds` is the original score database (data/sequence_lib_scores.db in that repository, a 737 MB Git LFS file; download it from
# https://media.githubusercontent.com/media/elifesciences-publications/germline/master/data/sequence_lib_scores.db), dumped with `db_dump -p`
# and saved as the vector of the 4^12 scores in alphabetical order of the 12-mers (as inst/extdata/glo_scores.rds). The checks:
#  1. The scores in the app's table are identical to the original database (all 16,777,216 values).
#  2. The original scoring function (Seqscore.pm: sum over words of a 7-word weighted window, divided by the number of codons) equals the app's
#     sum over all 12-mer windows, divided by the number of codons.
#  3. The original optimization (OptimizerTools.pm: best-word first guess, then random-position refinement, 10 rounds) is ported here and
#     compared with the app's optimize_glo() on the same proteins. The port was checked against the unchanged Perl code (dev/glo_original_perl/):
#     both give a score of 703.8 for the 238-residue GFP protein of the original repository.
suppressMessages(devtools::load_all(quiet = TRUE))
args <- commandArgs(trailingOnly = TRUE)
orig <- readRDS(args[1])
gfp_file <- args[2]
rounds <- if (length(args) > 2) as.integer(args[3]) else 10L
tab <- readRDS("inst/extdata/glo_scores.rds")
say <- function(...) cat(sprintf(...), "\n")
say("1. app table vs original database: %d entries, max |difference| %g", length(tab), max(abs(tab - orig)))

source("dev/glo_original_port.R")

set.seed(1)
random_cds <- function(n) sample(0:63, n, TRUE)
d <- vapply(c(5, 8, 13, 40, 101, 238), function(n) { cod <- random_cds(n); score_sequence(cod) - app_score(cod) }, 0)
say("2. original scoring minus the app's all-windows scoring, for random sequences of 5 to 238 codons: max |difference| %g", max(abs(d)))
# and with the app's own function, on a real DNA sequence
dna <- paste(codon_chr(random_cds(120)), collapse = "")
say("   app's calc_glo_score() / number of codons vs the original score of the same DNA: %.6f vs %.6f", calc_glo_score(dna) / 120, score_sequence(codon_idx(substring(dna, seq(1, 358, 3), seq(3, 360, 3)))))

protein <- function(cod) { map <- setNames(rep(names(orig_codons), lengths(orig_codons)), unlist(orig_codons)); paste(map[as.character(cod)], collapse = "") }
gfp <- toupper(paste(readLines(gfp_file, warn = FALSE), collapse = ""))
gfp_dna <- gsub("[^ACGT]", "", sub("^>.*", "", gfp)); gfp_aa <- protein(codon_idx(substring(gfp_dna, seq(1, nchar(gfp_dna) - 2, 3), seq(3, nchar(gfp_dna), 3))))
tests <- list(`GFP S65C (239 aa)` = gfp_aa, `first 60 aa of GFP` = substr(gfp_aa, 1, 60), `gld-1-like random 120 aa` = paste(sample(names(orig_codons), 120, TRUE), collapse = ""))
for (nm in names(tests)) {
  aa <- tests[[nm]]; n <- nchar(aa)
  t0 <- Sys.time(); o <- optimize_original(aa, rounds); t1 <- as.numeric(Sys.time() - t0, units = "secs")
  stopifnot(identical(protein(o$cod), aa))
  # the app: back-translate, then optimize_glo (random)
  ours <- vapply(1:5, function(i) { set.seed(100 + i); s <- optimize_glo(back_translate(aa), method = "heuristic"); codon <- codon_idx(substring(s, seq(1, nchar(s) - 2, 3), seq(3, nchar(s), 3))); stopifnot(identical(protein(codon), aa)); score_sequence(codon) }, 0)
  say("3. %-26s original algorithm (%d rounds, %.0f s): first guess %.1f, final %.1f | app's optimize_glo, 5 runs: %.1f to %.1f (mean %.1f) | gap to original: %.1f%%",
      nm, rounds, t1, o$first, o$score, min(ours), max(ours), mean(ours), 100 * (o$score - mean(ours)) / o$score)
}
