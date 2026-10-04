# Two different things must not be confused: a method may be *repeatable* (same
# input, same result) but the choice of synthetic introns must always be random.
# The server once fixed the global random seed for some methods, which made the
# "random" synthetic introns identical on every run.

long_cds <- function() paste0("ATG", strrep("GCTAAG", 30), "CTCGAG", strrep("GCTAAG", 30), "TAA")

# ---- with_local_seed() -----------------------------------------------------

test_that("with_local_seed gives the same draws for the same seed", {
  a <- with_local_seed(123, runif(5))
  b <- with_local_seed(123, runif(5))
  expect_identical(a, b)
  expect_false(identical(a, with_local_seed(124, runif(5))))
})

test_that("with_local_seed returns the value of the code it runs", {
  expect_equal(with_local_seed(1, 40 + 2), 42)
})

test_that("with_local_seed leaves the random stream exactly as it found it", {
  set.seed(99); expected <- runif(3)
  set.seed(99); with_local_seed(12345, runif(10)); actual <- runif(3)
  expect_identical(actual, expected)
})

test_that("with_local_seed restores the random stream even if the code fails", {
  set.seed(99); expected <- runif(3)
  set.seed(99)
  expect_error(with_local_seed(12345, stop("boom")), "boom")
  expect_identical(runif(3), expected)
})

test_that("with_local_seed with a NULL seed just runs the code without touching the seed", {
  set.seed(99); expected <- runif(3)
  set.seed(99); res <- with_local_seed(NULL, 7); actual <- runif(3)
  expect_equal(res, 7)
  expect_identical(actual, expected)
})

# ---- synthetic introns in the app ------------------------------------------

test_that("synthetic introns differ from run to run for every method", {
  for (method in c("7", "100", "1")) {  # High expression, No optimization, Ubiquitous
    picks <- vapply(1:5, function(i) {
      paste(introns_in(run_server_once(method, long_cds())), collapse = "|")
    }, character(1))
    expect_gt(length(unique(picks)), 1, label = paste("distinct intron choices, method", method))
  }
})

test_that("switching between methods does not make the intron choice repeat", {
  # A leaked seed from one method would show up as identical introns later on.
  picks <- vapply(rep(c("100", "7"), 3), function(m) {
    paste(introns_in(run_server_once(m, long_cds())), collapse = "|")
  }, character(1))
  expect_gt(length(unique(picks)), 3)
})

test_that("the app does not leave the global random seed fixed", {
  run_server_once("7", long_cds())
  a <- runif(4)
  run_server_once("7", long_cds())
  b <- runif(4)
  expect_false(identical(a, b))
  run_server_once("100", long_cds())
  c1 <- runif(4)
  run_server_once("100", long_cds())
  c2 <- runif(4)
  expect_false(identical(c1, c2))
})

test_that("two introns are chosen, both from the synthetic library, and are different", {
  lib <- safe_load_csv("Synthetic_Introns.csv")$Sequence
  picked <- introns_in(run_server_once("7", long_cds()))
  expect_length(picked, 2)
  expect_true(all(picked %in% lib))
  expect_false(picked[1] == picked[2])
})

# ---- repeatability that is meant to stay ------------------------------------

test_that("the High expression codon choices are repeatable even though the introns are not", {
  runs <- lapply(1:3, function(i) run_server_once("7", long_cds()))
  stripped <- vapply(runs, function(s) gsub("[acgt]", "", s), character(1))
  expect_length(unique(stripped), 1)
})

test_that("No optimization keeps its repeatable shields while the introns still vary", {
  # Random shield choices (Ubiq-weighted) use a fixed, temporary seed in this method.
  runs <- lapply(1:4, function(i) run_server_once("100", long_cds()))
  stripped <- vapply(runs, function(s) gsub("[acgt]", "", s), character(1))
  intron_sets <- vapply(runs, function(s) paste(introns_in(s), collapse = "|"), character(1))
  expect_length(unique(stripped), 1)
  expect_gt(length(unique(intron_sets)), 1)
})
