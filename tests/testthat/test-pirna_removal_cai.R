# "High expression" (CAI1) mode for the piRNA shield. As with restriction sites,
# the CAI1 column has nothing to sample from, so the shield makes as few changes
# as it can, choosing at each step the synonymous change that costs the least CAI
# (scored with the graded BringMansWeights column). The result is deterministic.

# A CDS that contains a real piRNA target exactly (Hamming distance 0), read in
# frame as MDTSTKNT*.
# This one is made only of best (CAI1) codons, so it survives the "High
# expression" recoding with the piRNA target still intact, as in the real pipeline.
pirna_cds <- function() "ATGAAGGTCATCGTCATCTTCGACTAA"

cods <- function(dna) substring(dna, seq(1, nchar(dna), 3), seq(3, nchar(dna), 3))
w_of <- function(codon) {
  cf <- WormTransgeneBuilder::codon_freqs
  cf$BringMansWeights[match(toupper(codon), toupper(cf$Codon))]
}

test_that("the fixture really contains a piRNA target", {
  expect_equal(min(get_pirna_distances(pirna_cds()), na.rm = TRUE), 0)
  expect_equal(protein_of(pirna_cds()), "MKVIVIFD*")
  expect_identical(probabilistic_recode(pirna_cds(), "CAI1"), pirna_cds())
})

test_that("CAI1 shielding lifts every piRNA above the threshold and keeps the protein", {
  for (target in c(2, 3)) {
    out <- remove_pirna_sites(pirna_cds(), "CAI1", target)
    expect_true(all(out$final_dists >= target, na.rm = TRUE), info = paste("target", target))
    expect_equal(protein_of(out$sequence), protein_of(pirna_cds()), info = paste("target", target))
    expect_equal(nchar(out$sequence), nchar(pirna_cds()), info = paste("target", target))
  }
})

test_that("the CAI1 piRNA shield is deterministic", {
  set.seed(1); a <- remove_pirna_sites(pirna_cds(), "CAI1", 3)$sequence
  set.seed(2); b <- remove_pirna_sites(pirna_cds(), "CAI1", 3)$sequence
  set.seed(3); d <- remove_pirna_sites(pirna_cds(), "CAI1", 3)$sequence
  expect_identical(a, b)
  expect_identical(a, d)
})

test_that("the CAI1 shield changes fewer codons than the random shield does", {
  det <- remove_pirna_sites(pirna_cds(), "CAI1", 3)$sequence
  n_det <- sum(cods(det) != cods(pirna_cds()))
  n_rand <- vapply(1:10, function(s) {
    set.seed(s)
    sum(cods(remove_pirna_sites(pirna_cds(), "Ubiq", 3)$sequence) != cods(pirna_cds()))
  }, numeric(1))
  expect_lt(n_det, median(n_rand))
})

test_that("the CAI1 shield keeps CAI higher than the random shield does", {
  det <- remove_pirna_sites(pirna_cds(), "CAI1", 3)$sequence
  cai_det <- calc_cai(det, "BringMansWeights")
  cai_rand <- vapply(1:10, function(s) {
    set.seed(s)
    calc_cai(remove_pirna_sites(pirna_cds(), "Ubiq", 3)$sequence, "BringMansWeights")
  }, numeric(1))
  expect_gt(cai_det, median(cai_rand))
})

test_that("when one codon change is enough, the cheapest such change is made", {
  # Brute force every single-codon synonymous change; if any reaches the
  # threshold, the shield must have made exactly one change and the best one.
  cf <- WormTransgeneBuilder::codon_freqs
  start <- pirna_cds()
  target <- 1  # one changed base is enough to leave distance 0
  cd <- cods(start)
  best_score <- -Inf; best_seq <- NA
  for (i in seq_along(cd)) {
    aa <- cf$Amino[match(cd[i], toupper(cf$Codon))]
    for (alt in setdiff(toupper(cf$Codon[cf$Amino == aa]), cd[i])) {
      trial <- cd; trial[i] <- alt; s <- paste(trial, collapse = "")
      if (all(get_pirna_distances(s) >= target, na.rm = TRUE)) {
        sc <- log(max(w_of(alt), 1e-6)) - log(max(w_of(cd[i]), 1e-6))
        if (sc > best_score + 1e-12) { best_score <- sc; best_seq <- s }
      }
    }
  }
  skip_if(is.na(best_seq), "No single-codon change reaches the threshold for this fixture")
  out <- remove_pirna_sites(start, "CAI1", target)$sequence
  expect_equal(sum(cods(out) != cd), 1)
  expect_identical(out, best_seq)
})

test_that("a sequence without piRNA matches is returned unchanged", {
  cds <- "ATGGCTAAAGGTCTTTCTGAACCTACTGTTAATCAAATTTATTGGTAA"
  out <- remove_pirna_sites(cds, "CAI1", 2)
  expect_identical(out$sequence, cds)
})

test_that("sequences shorter than 20 bp are returned unchanged", {
  out <- remove_pirna_sites("ATGGCTAAATAA", "CAI1", 2)
  expect_identical(out$sequence, "ATGGCTAAATAA")
})

test_that("a threshold of 0 or 1 needs no change beyond exact matches", {
  out <- remove_pirna_sites(pirna_cds(), "CAI1", 1)
  expect_true(all(out$final_dists >= 1, na.rm = TRUE))
  expect_equal(protein_of(out$sequence), protein_of(pirna_cds()))
})

test_that("the returned distances describe the returned sequence", {
  out <- remove_pirna_sites(pirna_cds(), "CAI1", 3)
  expect_equal(out$final_dists, get_pirna_distances(out$sequence))
  expect_equal(out$initial_dists, get_pirna_distances(pirna_cds()))
})

test_that("other methods keep their random behaviour", {
  outs <- vapply(1:10, function(s) {
    set.seed(s); remove_pirna_sites(pirna_cds(), "Ubiq", 3)$sequence
  }, character(1))
  expect_gt(length(unique(outs)), 1)
})

test_that("the full CAI1 pipeline with piRNA and enzyme shields is repeatable", {
  set.seed(1)
  a <- optimize_transgene(pirna_cds(), opt_method = "CAI1", enzymes = "XhoI",
                          remove_pirnas = TRUE, target_min_hamming = 3)
  set.seed(2)
  b <- optimize_transgene(pirna_cds(), opt_method = "CAI1", enzymes = "XhoI",
                          remove_pirnas = TRUE, target_min_hamming = 3)
  expect_identical(a$final_seq, b$final_seq)
  expect_true(all(a$final_pirna_dists >= 3, na.rm = TRUE))
  expect_equal(protein_of(a$optimized_cds), protein_of(pirna_cds()))
})
