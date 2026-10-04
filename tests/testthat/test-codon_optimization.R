test_that("probabilistic_recode keeps the protein, length and alphabet", {
  cds <- fixture_cds()
  for (s in test_seeds) {
    set.seed(s)
    out <- probabilistic_recode(cds, "Ubiq")
    expect_equal(protein_of(out), protein_of(cds), info = paste("seed", s))
    expect_equal(nchar(out), nchar(cds), info = paste("seed", s))
    expect_match(out, "^[ACGT]+$", info = paste("seed", s))
  }
})

test_that("probabilistic_recode works with every usage column", {
  cds <- fixture_cds()
  cols <- c("Ubiq", "CAI1", "BringMansWeights", "Chance", "aCDS", "HighRCSU")
  for (col in cols) {
    set.seed(1)
    out <- probabilistic_recode(cds, col)
    expect_equal(protein_of(out), protein_of(cds), info = col)
  }
})

test_that("probabilistic_recode actually changes the DNA", {
  cds <- fixture_cds()
  outs <- vapply(test_seeds, function(s) {
    set.seed(s)
    probabilistic_recode(cds, "Ubiq")
  }, character(1))
  expect_true(any(outs != cds))
})

test_that("probabilistic_recode is reproducible for a fixed seed", {
  set.seed(42); a <- probabilistic_recode(fixture_cds(), "Ubiq")
  set.seed(42); b <- probabilistic_recode(fixture_cds(), "Ubiq")
  expect_identical(a, b)
})

test_that("probabilistic_recode leaves Met and Trp codons alone", {
  # ATG and TGG are the only codons for M and W, so they cannot change.
  for (s in test_seeds) {
    set.seed(s)
    out <- probabilistic_recode("ATGTGGATGTGGATGTGG", "Ubiq")
    expect_identical(out, "ATGTGGATGTGGATGTGG")
  }
})

test_that("probabilistic_recode keeps a stop codon as a stop codon", {
  # The stop codon itself can be swapped (TAA/TAG/TGA), but must remain a stop.
  cds <- fixture_cds()
  for (s in test_seeds) {
    set.seed(s)
    out <- probabilistic_recode(cds, "Ubiq")
    expect_true(is_stop_codon(last_codon(out)), info = paste("seed", s))
  }
})

test_that("probabilistic_recode uppercases input and strips junk characters", {
  set.seed(1)
  out <- probabilistic_recode("atg gct\naaa taa", "Ubiq")
  expect_match(out, "^[ACGT]+$")
  expect_equal(nchar(out), 12)
})

test_that("sample_new_codon returns a different synonymous codon", {
  for (s in test_seeds) {
    set.seed(s)
    new <- sample_new_codon("A", "GCT", "Ubiq")
    expect_false(new == "GCT")
    expect_equal(protein_of(new), "A")
  }
})

test_that("sample_new_codon returns the old codon for single-codon amino acids", {
  expect_equal(sample_new_codon("M", "ATG", "Ubiq"), "ATG")
  expect_equal(sample_new_codon("W", "tgg", "Ubiq"), "TGG")
})

test_that("sample_new_codon always returns uppercase", {
  set.seed(1)
  expect_match(sample_new_codon("K", "aaa", "Ubiq"), "^[ACGT]{3}$")
})

test_that("CAI1 recoding is fully deterministic and keeps the supplied stop codon", {
  # In the CAI1 column every codon except the best one has weight 0, and all
  # three stop codons have weight 0, so there is nothing to sample: each amino
  # acid becomes its best codon and the stop codon is left as supplied.
  cds <- fixture_cds()
  for (stop_codon in c("TAA", "TAG", "TGA")) {
    input <- paste0(substr(cds, 1, nchar(cds) - 3), stop_codon)
    outs <- vapply(test_seeds, function(s) {
      set.seed(s); probabilistic_recode(input, "CAI1")
    }, character(1))
    expect_length(unique(outs), 1)
    expect_equal(last_codon(outs[1]), stop_codon)
    expect_equal(protein_of(outs[1]), protein_of(input))
  }
})

test_that("CAI1 recoding gives every amino acid its best codon", {
  cf <- WormTransgeneBuilder::codon_freqs
  out <- probabilistic_recode(fixture_cds(), "CAI1")
  cods <- substring(out, seq(1, nchar(out), 3), seq(3, nchar(out), 3))
  body <- cods[-length(cods)]  # skip the stop codon
  expect_true(all(cf$CAI1[match(body, toupper(cf$Codon))] == 1))
  expect_equal(calc_cai(out, "BringMansWeights"), 1)
})
