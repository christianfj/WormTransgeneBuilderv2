# optimize_rbs() makes 100 random synonymous versions of the first 39 bp, folds
# each with ViennaRNA's RNAfold, and keeps the one with the highest minimum free
# energy (the least stable structure). The original is always among the candidates.
# These tests need RNAfold, and are skipped without it.

skip_without_rnafold <- function() {
  skip_if(Sys.which("RNAfold") == "", "RNAfold (ViennaRNA) is not on the PATH")
}

# Minimum free energy of a start sequence, folded the same way the function does it
# (with "AAAA" in front, standing for the upstream consensus start).
fold_mfe <- function(seq) {
  out <- system2(Sys.which("RNAfold"), args = "--noPS", input = paste0("AAAA", seq),
                 stdout = TRUE, stderr = FALSE)
  line <- grep("\\(", out, value = TRUE)[1]
  as.numeric(gsub(" ", "", sub(".*\\(([^)]+)\\)\\s*$", "\\1", line)))
}

rbs_start <- function() substr(glo_cds(), 1, 39)  # helper from test-glo_optimization.R

test_that("the RBS fixture is 39 bp and encodes 13 codons", {
  expect_equal(nchar(rbs_start()), 39)
})

test_that("RNAfold gives an energy for the fixture", {
  skip_without_rnafold()
  expect_true(is.finite(fold_mfe(rbs_start())))
})

test_that("the RBS start keeps its length, protein and alphabet", {
  skip_without_rnafold()
  for (s in 1:5) {
    set.seed(s)
    out <- optimize_rbs(rbs_start(), "Ubiq")
    expect_equal(nchar(out), 39, info = paste("seed", s))
    expect_equal(protein_of(out), protein_of(rbs_start()), info = paste("seed", s))
    expect_match(out, "^[ACGT]+$", info = paste("seed", s))
  }
})

test_that("the RBS start is never more stable than the input", {
  skip_without_rnafold()
  before <- fold_mfe(rbs_start())
  for (s in 1:5) {
    set.seed(s)
    out <- optimize_rbs(rbs_start(), "Ubiq")
    # Higher (less negative) energy = weaker structure; the input is a candidate.
    expect_gte(fold_mfe(out), before - 1e-9)
  }
})

test_that("the RBS start is reproducible for a fixed seed", {
  skip_without_rnafold()
  set.seed(3); a <- optimize_rbs(rbs_start(), "Ubiq")
  set.seed(3); b <- optimize_rbs(rbs_start(), "Ubiq")
  expect_identical(a, b)
})

test_that("in CAI1 mode the result is the input or the all-best-codon version", {
  # CAI1 recoding is deterministic, so every candidate is one of these two.
  skip_without_rnafold()
  set.seed(1)
  out <- optimize_rbs(rbs_start(), "CAI1")
  expect_true(out %in% c(rbs_start(), probabilistic_recode(rbs_start(), "CAI1")))
})

test_that("a different number of candidates still returns a valid start", {
  skip_without_rnafold()
  set.seed(1)
  out <- optimize_rbs(rbs_start(), "Ubiq", num_candidates = 10)
  expect_equal(protein_of(out), protein_of(rbs_start()))
})
