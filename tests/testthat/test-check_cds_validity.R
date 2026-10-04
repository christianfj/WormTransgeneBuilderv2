test_that("a valid CDS passes with no warnings", {
  res <- check_cds_validity("ATGAAATAA")
  expect_true(res$valid)
  expect_length(res$warnings, 0)
  expect_true(check_cds_validity(fixture_cds())$valid)
})

test_that("lowercase and non-nucleotide characters are ignored", {
  expect_true(check_cds_validity("atg aaa\ntaa")$valid)
})

test_that("sequences shorter than 6 bp are rejected", {
  res <- check_cds_validity("AT")
  expect_false(res$valid)
  expect_match(res$warnings, "at least 6 bp")
})

test_that("a missing stop codon is reported", {
  res <- check_cds_validity("ATGAAA")
  expect_false(res$valid)
  expect_match(res$warnings, "valid stop codon")
})

test_that("a missing start codon is reported", {
  res <- check_cds_validity("CTGAAATAA")
  expect_false(res$valid)
  expect_match(paste(res$warnings, collapse = " "), "start codon")
})

test_that("a length that is not a multiple of 3 is reported", {
  res <- check_cds_validity("ATGAAAATAA")
  expect_false(res$valid)
  expect_match(paste(res$warnings, collapse = " "), "multiple of 3")
})

test_that("an internal stop codon is reported with its position", {
  res <- check_cds_validity("ATGTAAAAATAA")
  expect_false(res$valid)
  expect_match(res$warnings, "Internal stop codon found at position: 4")
})

test_that("all three stop codons are accepted at the end", {
  for (stop_codon in c("TAA", "TAG", "TGA")) {
    expect_true(check_cds_validity(paste0("ATGAAA", stop_codon))$valid)
  }
})
