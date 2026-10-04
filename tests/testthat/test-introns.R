# 2 introns need a reasonably long CDS. In this one every odd codon from the
# 3rd onwards is AAG (the preferred intron position), the others are GCT.
intron_cds <- function() paste0("ATG", strrep("GCTAAG", 100), "TAA")
intron_a <- "gtaagtttcag"

test_that("removing the lowercase introns gives back the original CDS", {
  for (type in c("early", "equidistant")) {
    res <- insert_introns(intron_cds(), c(intron_a, intron_a), type)
    expect_identical(gsub("[acgt]", "", res$seq), intron_cds(), info = type)
  }
})

test_that("the total length grows by exactly the intron lengths", {
  res <- insert_introns(intron_cds(), c(intron_a, intron_a))
  expect_equal(nchar(res$seq), nchar(intron_cds()) + 2 * nchar(intron_a))
})

test_that("introns are inserted in lowercase, the CDS stays uppercase", {
  res <- insert_introns(intron_cds(), c("GTAAGTTTCAG", "GTAAGTTTCAG"))
  expect_equal(sum(strsplit(res$seq, "")[[1]] %in% letters),
               2 * nchar(intron_a))
})

test_that("insertion points sit on codon boundaries, in increasing order", {
  for (type in c("early", "equidistant")) {
    res <- insert_introns(intron_cds(), c(intron_a, intron_a), type)
    expect_length(res$insertions, 2)
    expect_true(all(res$insertions %% 3 == 0), info = type)
    expect_false(is.unsorted(res$insertions), info = type)
  }
})

test_that("early spacing places introns near codons 50 and 100", {
  res <- insert_introns(intron_cds(), c(intron_a, intron_a), "early")
  expect_true(all(abs(res$insertions / 3 - c(50, 100)) <= 15))
})

test_that("equidistant spacing places introns near 1/3 and 2/3 of the CDS", {
  res <- insert_introns(intron_cds(), c(intron_a, intron_a), "equidistant")
  total_codons <- nchar(intron_cds()) / 3
  targets <- round(total_codons * c(1, 2) / 3)
  expect_true(all(abs(res$insertions / 3 - targets) <= 15))
})

test_that("introns prefer to follow an AAG codon", {
  res <- insert_introns(intron_cds(), c(intron_a, intron_a), "early")
  preceding <- substring(intron_cds(), res$insertions - 2, res$insertions)
  expect_equal(preceding, c("AAG", "AAG"))
})

test_that("an empty intron list is a no-op", {
  res <- insert_introns(intron_cds(), character(0))
  expect_identical(res$seq, intron_cds())
  expect_length(res$insertions, 0)
})

test_that("a single intron is inserted", {
  res <- insert_introns(intron_cds(), intron_a, "early")
  expect_length(res$insertions, 1)
  expect_identical(gsub("[acgt]", "", res$seq), intron_cds())
})
