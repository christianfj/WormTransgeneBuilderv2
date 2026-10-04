test_that("calc_gc_content returns the GC fraction", {
  expect_equal(calc_gc_content("ATGC"), 0.5)
  expect_equal(calc_gc_content("AAAA"), 0)
  expect_equal(calc_gc_content("GGCC"), 1)
})

test_that("calc_gc_content is case-insensitive and ignores junk characters", {
  expect_equal(calc_gc_content("atgc"), 0.5)
  expect_equal(calc_gc_content("AT GC\n"), 0.5)
})

test_that("calc_gc_content returns NA for empty or non-character input", {
  expect_true(is.na(calc_gc_content("")))
  expect_true(is.na(calc_gc_content(123)))
})

test_that("calc_cai is 1 for a sequence of only maximally-used codons", {
  # Met has a single codon, so its relative weight is 1 in every column.
  expect_equal(calc_cai("ATGATG", "Ubiq"), 1)
})

test_that("calc_cai equals the geometric mean of the codon weights", {
  cf <- WormTransgeneBuilder::codon_freqs
  w <- cf$Ubiq[match(c("ATG", "GCT", "GCT"), toupper(cf$Codon))]
  expect_equal(calc_cai("ATGGCTGCT", "Ubiq"), exp(mean(log(w))))
})

test_that("calc_cai ignores a dangling partial codon and handles short input", {
  expect_equal(calc_cai("ATGGCTGCTAA", "Ubiq"), calc_cai("ATGGCTGCT", "Ubiq"))
  expect_true(is.na(calc_cai("AT")))
})

test_that("calc_glo_score returns 0 for sequences under 4 codons", {
  expect_equal(calc_glo_score("ATGATG"), 0)
})

test_that("calc_glo_score returns a finite number for a normal CDS", {
  score <- calc_glo_score(fixture_cds())
  expect_true(is.numeric(score))
  expect_true(is.finite(score))
})
