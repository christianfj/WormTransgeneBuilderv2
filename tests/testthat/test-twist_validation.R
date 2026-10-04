# A 60 bp sequence with about 50% GC and no homopolymers or extreme windows.
balanced_seq <- function() strrep("ATGC", 15)

test_that("a balanced sequence passes with no warnings", {
  res <- validate_for_synthesis(balanced_seq())
  expect_true(res$pass)
  expect_length(res$warnings, 0)
})

test_that("empty or NULL input fails with an error message", {
  expect_false(validate_for_synthesis("")$pass)
  expect_match(validate_for_synthesis("")$warnings, "empty")
  expect_false(validate_for_synthesis(NULL)$pass)
})

test_that("a homopolymer of 14 bp fails but 13 bp does not", {
  gc_padding <- strrep("GCAT", 5)  # starts and ends without an "A" run
  bad  <- validate_for_synthesis(paste0(gc_padding, strrep("A", 14), gc_padding))
  ok   <- validate_for_synthesis(paste0(gc_padding, strrep("A", 13), gc_padding))
  expect_false(bad$pass)
  expect_match(paste(bad$warnings, collapse = " "), "Homopolymer of A")
  expect_true(ok$pass)
})

test_that("each of the four bases is checked for homopolymers", {
  for (base in c("A", "T", "G", "C")) {
    seq <- paste0(strrep("ATGC", 10), strrep(base, 14), strrep("ATGC", 10))
    res <- validate_for_synthesis(seq)
    expect_match(paste(res$warnings, collapse = " "),
                 paste("Homopolymer of", base), info = base)
  }
})

test_that("global GC below 25% or above 75% fails", {
  low  <- validate_for_synthesis(strrep("AT", 30))
  high <- validate_for_synthesis(strrep("GC", 30))
  expect_false(low$pass)
  expect_false(high$pass)
  expect_match(paste(low$warnings, collapse = " "), "Global GC")
  expect_match(paste(high$warnings, collapse = " "), "Global GC")
})

test_that("a high-GC 50 bp window fails even when global GC is acceptable", {
  # 50 bp of pure GC inside a long AT-rich sequence: global GC is ~50%.
  seq <- paste0(strrep("AT", 25), strrep("GC", 25), strrep("AT", 25))
  res <- validate_for_synthesis(seq)
  expect_false(res$pass)
  expect_match(paste(res$warnings, collapse = " "), "High local GC")
})

test_that("non-nucleotide characters are stripped before checking", {
  expect_true(validate_for_synthesis(paste0(balanced_seq(), " \n"))$pass)
})
