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

test_that("with include_input = FALSE the input start is never returned unchanged", {
  skip_without_rnafold()
  # A start in which no codon is the best codon of its amino acid, so any recoding by the graded weights
  # that picks a best codon somewhere differs from it; over several seeds it must differ every time.
  worst <- paste0("ATG", "AAA", "AAT", "GCA", "CTA", "GAA", "ACA", "TCA", "CCA", "GTA", "GGA", "ATA", "CAA")
  expect_equal(nchar(worst), 39)
  for (s in 1:5) {
    set.seed(s)
    out <- optimize_rbs(worst, "BringMansWeights", include_input = FALSE)
    expect_false(identical(out, worst), info = paste("seed", s))
    expect_equal(protein_of(out), protein_of(worst), info = paste("seed", s))
  }
})

test_that("graded weights give RNAfold many different High-expression candidates, not just one", {
  set.seed(1)
  candidates <- unique(replicate(100, probabilistic_recode(rbs_start(), "BringMansWeights")))
  expect_gt(length(candidates), 50)
  expect_length(unique(replicate(20, probabilistic_recode(rbs_start(), "CAI1"))), 1)  # the old behaviour
})

test_that("with RNAfold unavailable the fallback is returned (and a warning given)", {
  withr::local_envvar(PATH = "/nonexistent")
  testthat::local_mocked_bindings(find_rnafold = function(...) "")
  expect_warning(out <- optimize_rbs(rbs_start(), "BringMansWeights", include_input = FALSE,
                                     fallback_seq = "FALLBACK"), "RNAfold")
  expect_identical(out, "FALLBACK")
  expect_warning(out2 <- optimize_rbs(rbs_start(), "Ubiq"), "RNAfold")
  expect_identical(out2, rbs_start())   # other methods: unchanged behaviour, the input start
})

test_that("a different number of candidates still returns a valid start", {
  skip_without_rnafold()
  set.seed(1)
  out <- optimize_rbs(rbs_start(), "Ubiq", num_candidates = 10)
  expect_equal(protein_of(out), protein_of(rbs_start()))
})
