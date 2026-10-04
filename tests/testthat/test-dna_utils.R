# Small DNA helpers in base R (R/utils_dna.R). They replace Biostrings and stringr, so the app
# no longer depends on Bioconductor. The reference values below are written out independently,
# from the standard IUPAC table, and the helpers were also compared with Biostrings on thousands
# of random inputs when they were written.

iupac_letters_all <- strsplit("ACGTRYSWKMBDHVN", "")[[1]]

# ---- reverse_complement -------------------------------------------------------------------

test_that("each IUPAC letter complements as the standard table says", {
  expected <- c(A = "T", C = "G", G = "C", T = "A", R = "Y", Y = "R", S = "S", W = "W",
                K = "M", M = "K", B = "V", D = "H", H = "D", V = "B", N = "N")
  for (l in iupac_letters_all) expect_identical(reverse_complement(l), unname(expected[l]), info = l)
})

test_that("a sequence is complemented and reversed", {
  expect_identical(reverse_complement("AAAC"), "GTTT")
  expect_identical(reverse_complement("GAGACC"), "GGTCTC")   # the reverse-strand form of the BsaI site
  expect_identical(reverse_complement("ATGGCT"), "AGCCAT")
})

test_that("palindromic sites are their own reverse complement", {
  for (s in c("GGATCC", "GAATTC", "AAGCTT", "CTCGAG", "GGYRCC", "CYCGRG")) {
    expect_identical(reverse_complement(s), s, info = s)
  }
})

test_that("reversing twice gives the original", {
  set.seed(3)
  x <- vapply(1:200, function(i) paste(sample(iupac_letters_all, sample(1:30, 1), TRUE), collapse = ""), "")
  expect_identical(reverse_complement(reverse_complement(x)), x)
})

test_that("reverse_complement works on several sequences at once and accepts lowercase", {
  expect_identical(reverse_complement(c("AAAC", "gagacc", "")), c("GTTT", "GGTCTC", ""))
  expect_length(reverse_complement(character(0)), 0)
})

test_that("letters that are not IUPAC DNA codes give a clear error", {
  expect_error(reverse_complement("ACGU"), "IUPAC")
  expect_error(reverse_complement("AC GT"), "IUPAC")
  expect_error(reverse_complement(NA_character_), "IUPAC")
})

# ---- count_base ------------------------------------------------------------------------------

test_that("count_base counts a letter, ignoring case", {
  expect_equal(count_base("ATGCGGC", "G"), 3)
  expect_equal(count_base("ATGCGGC", "C"), 2)
  expect_equal(count_base("atgcggc", "G"), 3)
  expect_equal(count_base("AAAA", "G"), 0)
  expect_equal(count_base("", "G"), 0)
})

# ---- approximate_match_starts -------------------------------------------------------------------

test_that("with no mismatches allowed it finds exact matches, overlapping ones included", {
  expect_identical(approximate_match_starts("ACA", "ACACACA", 0), c(1L, 3L, 5L))
  expect_identical(approximate_match_starts("GGATCC", "AAGGATCCAA", 0), 3L)
  expect_identical(approximate_match_starts("GGATCC", "AAAAAAAAAA", 0), integer(0))
})

test_that("it finds windows with up to the allowed number of mismatches", {
  subject <- "AAAAGGATCCAAAAGGTTCCAAAAGCTTCCAAAA"
  expect_identical(approximate_match_starts("GGATCC", subject, 0), 5L)
  expect_identical(approximate_match_starts("GGATCC", subject, 1), c(5L, 15L))       # GGTTCC has one mismatch
  expect_identical(approximate_match_starts("GGATCC", subject, 2), c(5L, 15L, 25L))  # GCTTCC has two
})

test_that("an N in the sequence is compatible with any base, and an N in the pattern too", {
  expect_identical(approximate_match_starts("GGATCC", "AAGGNTCCAA", 0), 3L)
  expect_identical(approximate_match_starts("GGNTCC", "AAGGATCCAA", 0), 3L)
  expect_identical(approximate_match_starts("GGATCC", "AAGGRTCCAA", 0), 3L)   # R = A or G includes A
  expect_identical(approximate_match_starts("GGATCC", "AAGGYTCCAA", 0), integer(0))  # Y = C or T excludes A
})

test_that("matches must lie entirely inside the sequence", {
  # A pattern hanging over an end is not a window of the sequence, however many mismatches are allowed.
  expect_identical(approximate_match_starts("GGATCC", "GATCC", 3), integer(0))
  expect_identical(approximate_match_starts("GGATCC", "AAAAGG", 2), integer(0))
  expect_true(all(approximate_match_starts("GGATCC", "ATATATATAT", 4) >= 1))
  x <- approximate_match_starts("GGATCC", "ATATATATAT", 4)
  expect_true(all(x + 5 <= 10))
})

test_that("a pattern longer than the sequence, or empty input, finds nothing", {
  expect_identical(approximate_match_starts("ACGTACGT", "ACG", 1), integer(0))
  expect_identical(approximate_match_starts("ACG", "", 1), integer(0))
  expect_identical(approximate_match_starts("", "ACG", 1), integer(0))
})

test_that("case does not matter, and characters that are not IUPAC letters match nothing", {
  expect_identical(approximate_match_starts("ggatcc", "aaGGATCCaa", 0), 3L)
  expect_identical(approximate_match_starts("GGATCC", "AAGGA-CCAA", 0), integer(0))
  expect_identical(approximate_match_starts("GGATCC", "AAGGA-CCAA", 1), 3L)
})

test_that("the result is a plain integer vector", {
  x <- approximate_match_starts("ACA", "ACACACA", 0)
  expect_type(x, "integer")
  expect_null(names(x))
})

test_that("a window of the sequence is reported when its mismatches are within the limit", {
  set.seed(9)
  for (i in 1:100) {
    pat <- paste(sample(c("A", "C", "G", "T"), 8, TRUE), collapse = "")
    sub <- paste(sample(c("A", "C", "G", "T"), 40, TRUE), collapse = "")
    k <- sample(0:3, 1)
    starts <- seq_len(nchar(sub) - 7)
    mism <- vapply(starts, function(s) sum(strsplit(pat, "")[[1]] != strsplit(substr(sub, s, s + 7), "")[[1]]), 0)
    expect_identical(approximate_match_starts(pat, sub, k), as.integer(starts[mism <= k]), info = paste(pat, sub, k))
  }
})

# ---- the viewer's use of it -----------------------------------------------------------------------------

test_that("the viewer highlights a piRNA-like target found with mismatches, as whole windows", {
  target <- "GGATCCGGATCCGGATCCAA"
  seq <- paste0("TTTTT", "GGATCCGGTTCCGGATCCAA", "TTTTT")    # one mismatch inside
  res <- add_matches(target, seq, "#f9e79f", "piRNA Target", "piRNA", max_mismatch = 2)
  expect_equal(res$starts, 6)
  expect_equal(res$ends, 25)
  expect_equal(res$tooltips, "piRNA Target")
})

test_that("the viewer does not highlight a target that only fits by hanging off the end", {
  target <- "GGATCCGGATCCGGATCCAA"
  seq <- paste0("TTTTTTTTTTTTTTTTTTTT", "GGATCCGGATCCGGAT")   # the target would overhang the end by 4
  res <- add_matches(target, seq, "#f9e79f", "piRNA Target", "piRNA", max_mismatch = 4)
  expect_length(res$starts, 0)
})

test_that("with no mismatches allowed the viewer still finds exact matches", {
  res <- add_matches("GGATCC", "AAGGATCCAA", "#000", "x", "y", max_mismatch = 0)
  expect_equal(res$starts, 3)
})
