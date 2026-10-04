# get_pirna_distances() and the random (non-CAI1) piRNA shield. The CAI1 version is
# tested in test-pirna_removal_cai.R.

# A CDS that contains a real piRNA target exactly (Hamming distance 0).
pirna_random_cds <- function() "ATGGATACCTCAACAAAAAACACATAA"

test_that("get_pirna_distances returns one distance per piRNA", {
  d <- get_pirna_distances(pirna_random_cds())
  expect_length(d, length(WormTransgeneBuilder::pirna_sequences))
})

test_that("a piRNA target inside the sequence has distance 0", {
  d <- get_pirna_distances(pirna_random_cds())
  target <- "GATACCTCAACAAAAAACAC"
  expect_equal(d[which(WormTransgeneBuilder::pirna_sequences == target)[1]], 0)
})

test_that("distances are whole numbers between 0 and 20", {
  d <- get_pirna_distances(pirna_random_cds())
  d <- d[!is.na(d)]
  expect_true(all(d == round(d)))
  expect_true(all(d >= 0 & d <= 20))
})

test_that("input shorter than 20 bp gives NA for every piRNA", {
  d <- get_pirna_distances("ATGGCTAAATAA")
  expect_length(d, length(WormTransgeneBuilder::pirna_sequences))
  expect_true(all(is.na(d)))
})

test_that("get_pirna_distances ignores case and junk characters", {
  expect_equal(get_pirna_distances(tolower(pirna_random_cds())),
               get_pirna_distances(pirna_random_cds()))
  expect_equal(get_pirna_distances(paste0(" ", substr(pirna_random_cds(), 1, 9), "\n",
                                          substr(pirna_random_cds(), 10, 27))),
               get_pirna_distances(pirna_random_cds()))
})

test_that("a longer sequence never has a larger minimum distance than one of its parts", {
  # Extra sequence can only add windows, so each piRNA's minimum can only fall.
  short <- get_pirna_distances(pirna_random_cds())
  long <- get_pirna_distances(paste0(pirna_random_cds(), "GCTAAAGGTCTTTCTGAA"))
  expect_true(all(long <= short, na.rm = TRUE))
})

test_that("the random shield lifts every piRNA above the threshold and keeps the protein", {
  for (target in c(2, 3)) {
    for (s in 1:5) {
      set.seed(s)
      out <- remove_pirna_sites(pirna_random_cds(), "Ubiq", target)
      info <- paste("target", target, "seed", s)
      expect_true(all(out$final_dists >= target, na.rm = TRUE), info = info)
      expect_equal(protein_of(out$sequence), protein_of(pirna_random_cds()), info = info)
      expect_equal(nchar(out$sequence), nchar(pirna_random_cds()), info = info)
      expect_match(out$sequence, "^[ACGT]+$", info = info)
    }
  }
})

test_that("the returned distances describe the returned sequence", {
  set.seed(1)
  out <- remove_pirna_sites(pirna_random_cds(), "Ubiq", 3)
  expect_equal(out$initial_dists, get_pirna_distances(pirna_random_cds()))
  expect_equal(out$final_dists, get_pirna_distances(out$sequence))
})

test_that("the random shield is reproducible for a fixed seed and varies between seeds", {
  set.seed(5); a <- remove_pirna_sites(pirna_random_cds(), "Ubiq", 3)$sequence
  set.seed(5); b <- remove_pirna_sites(pirna_random_cds(), "Ubiq", 3)$sequence
  expect_identical(a, b)
  outs <- vapply(1:6, function(s) { set.seed(s); remove_pirna_sites(pirna_random_cds(), "Ubiq", 3)$sequence }, character(1))
  expect_gt(length(unique(outs)), 1)
})

test_that("a sequence with no piRNA match is returned unchanged", {
  cds <- "ATGGCTAAAGGTCTTTCTGAACCTACTGTTAATCAAATTTATTGGTAA"
  expect_gte(min(get_pirna_distances(cds), na.rm = TRUE), 3)
  set.seed(1)
  out <- remove_pirna_sites(cds, "Ubiq", 3)
  expect_identical(out$sequence, cds)
})

test_that("a sequence shorter than 20 bp is returned unchanged", {
  out <- remove_pirna_sites("ATGGCTAAATAA", "Ubiq", 3)
  expect_identical(out$sequence, "ATGGCTAAATAA")
})

test_that("a threshold of 0 never changes anything", {
  # No distance can be below 0.
  out <- remove_pirna_sites(pirna_random_cds(), "Ubiq", 0)
  expect_identical(out$sequence, pirna_random_cds())
})
