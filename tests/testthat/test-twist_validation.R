# Twist Bioscience's published gene synthesis guidelines, read on 2026-10-14 (docs/METHODS.md, section 8). Two levels: standard guidelines
# and high complexity.

# A non-repeating sequence with about 50% GC and no extreme windows (random, but the same every time).
balanced_seq <- function(n = 400) {
  withr::with_seed(11, paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = ""))
}
all_text <- function(res) paste(res$warnings, collapse = " ")

# The ccdB gene (the sequence supplied by the authors), 306 bases (101 amino acids).
ccdb_dna <- "atgcagtttaaggtttacacctataaaagagagagccgttatcgtctgtttgtggatgtacagagtgatattattgacacgcccgggcgacggatggtgatccccctggccagtgcacgtctgctgtcagataaagtctcccgtgaactttacccggtggtgcatatcggggatgaaagctggcgcatgatgaccaccgatatggccagtgtgccggtctccgttatcggggaagaagtggctgatctcagccaccgcgaaaatgacatcaaaaacgccattaacctgatgttctggggaatataa"

test_that("a balanced sequence passes with no warnings at either level", {
  res <- validate_for_synthesis(balanced_seq())
  expect_true(res$pass)
  expect_false(res$high_complexity)
  expect_length(res$warnings, 0)
  expect_length(res$standard_issues, 0)
  expect_length(res$high_complexity_issues, 0)
})

test_that("empty or NULL input fails with an error message", {
  expect_false(validate_for_synthesis("")$pass)
  expect_match(validate_for_synthesis("")$warnings, "empty")
  expect_false(validate_for_synthesis(NULL)$pass)
})

test_that("lower case bases (introns, UTRs, promoters) are checked, not dropped", {
  expect_match(all_text(validate_for_synthesis(tolower(paste0(balanced_seq(200), strrep("A", 14), balanced_seq(200))))), "Homopolymer of A")
})

test_that("non-nucleotide characters are stripped before checking", {
  expect_true(validate_for_synthesis(paste0(balanced_seq(), " \n"))$pass)
})

# ---- homopolymers: 14 or more is a standard issue, more than 30 is high complexity ------------------------

test_that("a homopolymer of 14 bp is a standard issue but 13 bp is fine", {
  pad <- balanced_seq(200)
  bad <- validate_for_synthesis(paste0(pad, strrep("A", 14), pad))
  ok <- validate_for_synthesis(paste0(pad, strrep("A", 13), pad))
  expect_false(bad$pass)
  expect_false(bad$high_complexity)
  expect_match(paste(bad$standard_issues, collapse = " "), "Homopolymer of A")
  expect_true(ok$pass)
})

test_that("a homopolymer of 30 bp is a standard issue and 31 bp is high complexity", {
  pad <- balanced_seq(200)
  thirty <- validate_for_synthesis(paste0(pad, strrep("T", 30), pad))
  thirty_one <- validate_for_synthesis(paste0(pad, strrep("T", 31), pad))
  expect_false(thirty$high_complexity)
  expect_match(paste(thirty$standard_issues, collapse = " "), "Homopolymer of T")
  expect_true(thirty_one$high_complexity)
  expect_match(paste(thirty_one$high_complexity_issues, collapse = " "), "Homopolymer of T of 31")
  expect_length(thirty_one$standard_issues, 0)
})

test_that("each of the four bases is checked for homopolymers", {
  for (base in c("A", "T", "G", "C")) {
    res <- validate_for_synthesis(paste0(balanced_seq(200), strrep(base, 14), balanced_seq(200)))
    expect_match(paste(res$standard_issues, collapse = " "), paste("Homopolymer of", base), info = base)
  }
})

# ---- GC content -------------------------------------------------------------------------------------------

test_that("global GC below 25% or above 75% is a standard issue", {
  # Complex enough (random) so that no repeat is involved; AT-only and GC-only sequences are also extreme in windows (high complexity).
  low <- withr::with_seed(2, paste(sample(c("A", "T", "A", "T", "A", "T", "A", "C", "G"), 400, replace = TRUE), collapse = ""))
  high <- withr::with_seed(3, paste(sample(c("G", "C", "G", "C", "G", "C", "G", "A", "T"), 400, replace = TRUE), collapse = ""))
  expect_match(paste(validate_for_synthesis(low)$standard_issues, collapse = " "), "Global GC")
  expect_match(paste(validate_for_synthesis(high)$standard_issues, collapse = " "), "Global GC")
  expect_equal(calc_gc_content(low) < 0.25, TRUE)
  expect_equal(calc_gc_content(high) > 0.75, TRUE)
})

test_that("a 50 bp window between 80% and 90% GC is a standard warning, above 90% or below 10% is high complexity", {
  at <- balanced_seq(300)
  make_window <- function(gc_count) {   # 50 bases with the stated number of G/C, mixed
    paste(withr::with_seed(5, sample(c(rep("G", gc_count %/% 2), rep("C", gc_count - gc_count %/% 2), rep("A", (50 - gc_count) %/% 2), rep("T", 50 - gc_count - (50 - gc_count) %/% 2)))), collapse = "")
  }
  at2 <- withr::with_seed(12, paste(sample(c("A", "C", "G", "T"), 300, replace = TRUE), collapse = ""))   # a different flank, so there is no repeat
  w85 <- validate_for_synthesis(paste0(at, make_window(43), at2))    # 86%
  w96 <- validate_for_synthesis(paste0(at, make_window(48), at2))    # 96%
  w4 <- validate_for_synthesis(paste0(at, make_window(2), at2))      # 4%
  expect_match(paste(w85$standard_issues, collapse = " "), "more than 80% GC")
  expect_false(w85$high_complexity)
  expect_true(w96$high_complexity)
  expect_match(paste(w96$high_complexity_issues, collapse = " "), "more than 90% GC")
  expect_true(w4$high_complexity)
  expect_match(paste(w4$high_complexity_issues, collapse = " "), "less than 10% GC")
})

# ---- repeats -----------------------------------------------------------------------------------------------

test_that("a direct repeat of 201 bp is high complexity but 200 bp is not", {
  a <- balanced_seq(300); unit201 <- substr(balanced_seq(500), 1, 201)
  bad <- validate_for_synthesis(paste0(unit201, a, unit201))
  ok <- validate_for_synthesis(paste0(substr(unit201, 1, 200), a, substr(unit201, 1, 200)))
  expect_true(bad$high_complexity)
  expect_match(paste(bad$high_complexity_issues, collapse = " "), "direct repeat")
  expect_false(ok$high_complexity)
})

# ---- length and ccdB -----------------------------------------------------------------------------------------

test_that("the length limits are 300 bp, 5,000 bp (gene fragments) and 7,000 bp (clonal genes)", {
  expect_match(all_text(validate_for_synthesis(balanced_seq(299))), "below the 300 bp minimum")
  expect_false(grepl("Length", all_text(validate_for_synthesis(balanced_seq(300)))))
  expect_match(all_text(validate_for_synthesis(balanced_seq(5001))), "5,000 bp maximum for gene fragments")
  expect_match(all_text(validate_for_synthesis(balanced_seq(7001))), "7,000 bp maximum")
  expect_false(grepl("Length", all_text(validate_for_synthesis(balanced_seq(5000)))))
})

test_that("ccdB is found as DNA, on the other strand, and when recoded with different codons", {
  expect_match(all_text(validate_for_synthesis(paste0(balanced_seq(200), ccdb_dna, balanced_seq(200)))), "ccdB")
  expect_match(all_text(validate_for_synthesis(paste0(balanced_seq(200), reverse_complement(toupper(ccdb_dna)), balanced_seq(200)))), "ccdB")
  withr::with_seed(4, recoded <- probabilistic_recode(toupper(ccdb_dna), "Chance"))
  expect_false(identical(recoded, toupper(ccdb_dna)))
  expect_match(all_text(validate_for_synthesis(paste0(balanced_seq(200), recoded, balanced_seq(200)))), "ccdB")
})

test_that("ordinary sequences are not reported as ccdB", {
  expect_false(grepl("ccdB", all_text(validate_for_synthesis(balanced_seq(3000)))))
  expect_false(grepl("ccdB", all_text(validate_for_synthesis(fixture_cds()))))
})
