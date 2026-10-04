# GLO (germline optimization, Dickinson et al.) is a greedy search: it tries a random
# synonymous codon at a random position and keeps the change only if the summed
# 12-mer scores of the affected windows go up. It is random, so tests check
# properties that must hold for any seed rather than exact sequences.

glo_seeds <- 1:5

test_that("the GLO fixture is a valid CDS", {
  expect_true(check_cds_validity(glo_cds())$valid)
})

test_that("GLO keeps the protein, the length and the alphabet", {
  for (s in glo_seeds) {
    set.seed(s)
    out <- optimize_glo(glo_cds())
    expect_equal(protein_of(out), protein_of(glo_cds()), info = paste("seed", s))
    expect_equal(nchar(out), nchar(glo_cds()), info = paste("seed", s))
    expect_match(out, "^[ACGT]+$", info = paste("seed", s))
  }
})

test_that("GLO never lowers the GLO score, and improves it", {
  before <- calc_glo_score(glo_cds())
  after <- vapply(glo_seeds, function(s) { set.seed(s); calc_glo_score(optimize_glo(glo_cds())) }, numeric(1))
  # Each accepted change raises the summed window scores, so the total cannot fall.
  expect_true(all(after >= before))
  expect_gt(median(after), before)
})

test_that("GLO keeps the last codon a stop codon", {
  for (s in glo_seeds) {
    set.seed(s)
    expect_true(is_stop_codon(last_codon(optimize_glo(glo_cds()))), info = paste("seed", s))
  }
})

test_that("GLO is reproducible for a fixed seed and varies between seeds", {
  set.seed(7); a <- optimize_glo(glo_cds())
  set.seed(7); b <- optimize_glo(glo_cds())
  expect_identical(a, b)
  outs <- vapply(glo_seeds, function(s) { set.seed(s); optimize_glo(glo_cds()) }, character(1))
  expect_gt(length(unique(outs)), 1)
})

test_that("GLO leaves Met and Trp codons alone", {
  # They have a single codon each, so no synonymous change exists.
  cds <- "ATGTGGATGTGGATGTGG"
  set.seed(1)
  expect_identical(optimize_glo(cds), cds)
})

test_that("sequences shorter than 4 codons are returned exactly as given", {
  expect_identical(optimize_glo("atggctaaa"), "atggctaaa")
})

test_that("lowercase input is accepted and comes back in uppercase", {
  set.seed(1)
  out <- optimize_glo(tolower(glo_cds()))
  expect_match(out, "^[ACGT]+$")
  expect_equal(protein_of(out), protein_of(glo_cds()))
})

test_that("calc_glo_score matches a direct lookup of each window in the score table", {
  # calc_glo_score finds each window's score by arithmetic on its bases. Check it against a
  # separate lookup that reads each window as a base-4 number. (More checks of the score
  # table are in test-glo_scores.R.)
  cds <- glo_cds()
  starts <- seq(1, nchar(cds) - 11, by = 3)
  windows <- substring(cds, starts, starts + 11)
  scores <- glo_score_table()[strtoi(chartr("ACGT", "0123", windows), base = 4L) + 1L]
  expect_equal(calc_glo_score(cds), sum(scores))
})
