# The piRNA distance engine (R/fct_pirna_engine.R) replaces a character-by-character
# comparison of every 20-base window with every piRNA by whole-number arithmetic. It must
# give exactly the same distances, so these tests keep the original stringdist method as a
# reference and compare the two on many kinds of input.

# ---- the original method, kept here as the reference ----------------------------------
reference_distance_matrix <- function(seq, piRNAs = WormTransgeneBuilder::pirna_sequences) {
  s <- toupper(gsub("[^ATGCNatgcn]", "", seq))
  starts <- seq_len(nchar(s) - 19)
  kmers <- substring(s, starts, starts + 19)
  stringdist::stringdistmatrix(kmers, piRNAs, method = "hamming", nthread = 2)
}
reference_get_distances <- function(seq, piRNAs = WormTransgeneBuilder::pirna_sequences) {
  if (nchar(toupper(gsub("[^ATGCNatgcn]", "", seq))) < 20) return(rep(NA, length(piRNAs)))
  res <- suppressWarnings(apply(reference_distance_matrix(seq, piRNAs), 2, min, na.rm = TRUE))
  res[is.infinite(res)] <- NA
  res
}

random_dna <- function(n, seed) {
  set.seed(seed)
  paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
}
# A sequence that contains real piRNA targets, some with a few bases changed.
planted_dna <- function(seed = 1) {
  p <- WormTransgeneBuilder::pirna_sequences
  set.seed(seed)
  chosen <- p[sample.int(length(p), 4)]
  bases <- c("A", "C", "G", "T")
  mutate <- function(x, n) {
    pos <- sample.int(20, n)
    ch <- strsplit(x, "")[[1]]
    ch[pos] <- vapply(ch[pos], function(b) sample(setdiff(bases, b), 1), "")
    paste(ch, collapse = "")
  }
  paste0(random_dna(37, seed), chosen[1], random_dna(23, seed + 1), mutate(chosen[2], 1),
         random_dna(31, seed + 2), mutate(chosen[3], 3), random_dna(29, seed + 3), chosen[4], random_dna(41, seed + 4))
}

# ---- building blocks -----------------------------------------------------------------------

test_that("the mismatch table counts the differing bases in a 10-base code", {
  tbl <- mismatch_count_table()
  expect_length(tbl, 4^10)
  expect_equal(tbl[1], 0)                 # code 0: no difference
  expect_equal(tbl[2], 1)                 # one base differs
  expect_equal(tbl[4 + 1], 1)             # a different single base
  expect_equal(tbl[4^10], 10)             # all ten bases differ (code 4^10 - 1)
  expect_equal(tbl[(1 + 4 * 3) + 1], 2)   # two bases differ
  expect_type(tbl, "integer")
})

test_that("the mismatch table agrees with counting the differences directly", {
  set.seed(1)
  codes <- sample(0:(4^10 - 1), 200)
  direct <- vapply(codes, function(x) sum(vapply(0:9, function(g) (x %/% 4^g) %% 4 != 0, TRUE)), numeric(1))
  expect_equal(mismatch_count_table()[codes + 1], direct)
})

test_that("20-mers are split into two 10-base numbers", {
  e <- encode_20mers(c("AAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTT", "AAAAAAAAAACCCCCCCCCC", "ACGTACGTACGTACGTACGT"))
  expect_equal(e$valid, rep(TRUE, 4))
  expect_equal(e$lo, c(0L, 4L^10 - 1L, 0L, strtoi("0123012301", base = 4L)))
  expect_equal(e$hi, c(0L, 4L^10 - 1L, strtoi("1111111111", base = 4L), strtoi("2301230123", base = 4L)))
})

test_that("anything that is not exactly 20 plain A/C/G/T bases is marked unusable", {
  e <- encode_20mers(c("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACG", "ACGTNCGTACGTACGTACGT", NA, "acgtacgtacgtacgtacgt", ""))
  expect_equal(e$valid, c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
  expect_true(all(is.na(e$lo[!e$valid])))
})

test_that("the piRNA database is encoded once and reused", {
  a <- pirna_reference(); b <- pirna_reference()
  expect_identical(a, b)
  expect_length(a$lo, length(WormTransgeneBuilder::pirna_sequences))
  expect_true(all(a$valid))
  expect_identical(a$source, WormTransgeneBuilder::pirna_sequences)
})

# ---- the distances ---------------------------------------------------------------------------

test_that("get_pirna_distances gives exactly the same result as the original method", {
  inputs <- list(
    random_20 = random_dna(20, 1), random_21 = random_dna(21, 2), random_100 = random_dna(100, 3),
    random_500 = random_dna(500, 4), planted_a = planted_dna(1), planted_b = planted_dna(2),
    fixture = "ATGGATACCTCAACAAAAAACACATAA", cds = glo_cds()
  )
  for (nm in names(inputs)) {
    expect_identical(get_pirna_distances(inputs[[nm]]), reference_get_distances(inputs[[nm]]), info = nm)
  }
})

test_that("results match for input with lowercase letters and stray characters", {
  s <- planted_dna(3)
  messy <- paste0(tolower(substr(s, 1, 60)), " \n", substr(s, 61, nchar(s)), "!!")
  expect_identical(get_pirna_distances(messy), reference_get_distances(messy))
})

test_that("results match when the sequence contains N", {
  s <- planted_dna(4)
  with_n <- paste0(substr(s, 1, 50), "NN", substr(s, 53, 120), "N", substr(s, 122, nchar(s)))
  expect_identical(get_pirna_distances(with_n), reference_get_distances(with_n))
  all_n <- strrep("N", 60)
  expect_identical(get_pirna_distances(all_n), reference_get_distances(all_n))
})

test_that("sequences shorter than 20 bases still give NA for every piRNA", {
  for (n in c(0, 1, 19)) {
    d <- get_pirna_distances(substr(strrep("ACGT", 10), 1, n))
    expect_identical(d, reference_get_distances(substr(strrep("ACGT", 10), 1, n)), info = n)
    expect_length(d, length(WormTransgeneBuilder::pirna_sequences))
  }
})

test_that("the distance result has the same type and shape as before", {
  d <- get_pirna_distances(planted_dna(1))
  expect_type(d, "double")
  expect_null(names(d))
  expect_null(dim(d))
})

# ---- the scan, with pairs and per-window minima ----------------------------------------------

scan_of <- function(seq, ...) {
  starts <- seq_len(nchar(seq) - 19)
  scan_kmers_against_pirnas(substring(seq, starts, starts + 19), ...)
}

test_that("the per-window minimum matches the original matrix", {
  s <- planted_dna(1)
  m <- reference_distance_matrix(s)
  expect_equal(scan_of(s)$row_min, apply(m, 1, min))
})

test_that("the (window, piRNA) pairs below a threshold match the original, in the same order", {
  s <- planted_dna(2)
  m <- reference_distance_matrix(s)
  for (thr in c(1, 2, 3, 5)) {
    hits <- which(m < thr, arr.ind = TRUE)
    p <- scan_of(s, threshold = thr)$pairs
    expect_identical(p$window, as.integer(hits[, 1]), info = thr)
    expect_identical(p$pirna, as.integer(hits[, 2]), info = thr)
  }
})

test_that("no pairs are returned when nothing is that close, and none when no threshold is asked for", {
  s <- random_dna(200, 9)
  expect_equal(nrow(scan_of(s, threshold = 0)$pairs), 0)
  expect_null(scan_of(s)$pairs)
})

test_that("a sequence with several matches gives pairs for each of them", {
  s <- planted_dna(1)
  p <- scan_of(s, threshold = 1)$pairs
  expect_gte(length(unique(p$pirna)), 2)
})

test_that("a custom reference with unusable entries behaves like the original method", {
  refs <- c(WormTransgeneBuilder::pirna_sequences[1:50], NA, "ACGT", "ACGTACGTACGTACGTACGTA", "ACGTNCGTACGTACGTACGT")
  s <- paste0(planted_dna(1), refs[3])
  ours <- scan_of(s, reference = make_pirna_reference(refs))
  res <- reference_get_distances(s, refs)
  expect_identical(ours$col_min, as.numeric(res))
  expect_true(all(is.na(ours$col_min[51:53])))   # missing, too short, too long: cannot be compared
  expect_false(is.na(ours$col_min[54]))          # 20 long with an N: compared, as the original does
  expect_equal(ours$row_min, apply(reference_distance_matrix(s, refs), 1, min, na.rm = TRUE))
})

test_that("an empty reference or no windows gives empty, well-formed results", {
  none <- scan_kmers_against_pirnas(character(0))
  expect_length(none$row_min, 0)
  expect_length(none$col_min, length(WormTransgeneBuilder::pirna_sequences))
  expect_true(all(is.na(none$col_min)))
  empty_ref <- scan_kmers_against_pirnas("ACGTACGTACGTACGTACGT", reference = make_pirna_reference(character(0)))
  expect_length(empty_ref$col_min, 0)
  expect_equal(empty_ref$row_min, Inf)
})

# ---- the shields give the same answers as before -------------------------------------------------
# Reference outputs produced by the previous implementation (before the engine existed).

pirna_fixture <- function() "ATGGATACCTCAACAAAAAACACATAA"

test_that("the random shield gives the same sequences as before for fixed seeds", {
  expected <- list(
    "2_1" = "ATGGACACTTCGACGAAGAATACATAA", "2_2" = "ATGGACACATCGACTAAGAATACATAA", "2_3" = "ATGGACACATCCACTAAGAATACATAA",
    "3_1" = "ATGGACACTTCGACGAAGAATACATAA", "3_2" = "ATGGACACATCGACTAAGAATACATAA", "3_3" = "ATGGACACATCCACTAAGAATACATAA")
  for (thr in c(2, 3)) for (s in 1:3) {
    set.seed(s)
    expect_identical(remove_pirna_sites(pirna_fixture(), "Ubiq", thr)$sequence, expected[[paste0(thr, "_", s)]],
                     info = paste("threshold", thr, "seed", s))
  }
})

test_that("the High-expression shield gives the same sequences as before", {
  best <- "ATGAAGGTCATCGTCATCTTCGACTAA"
  expected <- c("1" = "ATGAAGGTTATCGTCATCTTCGACTAA", "2" = "ATGAAGGTTATCGTTATCTTCGACTAA", "3" = "ATGAAGGTTATTGTTATCTTCGACTAA")
  for (thr in 1:3) expect_identical(remove_pirna_sites(best, "CAI1", thr)$sequence, expected[[as.character(thr)]], info = thr)
})

test_that("the shields give exactly the same results as before on the benchmark genes", {
  helpers <- testthat::test_path("..", "..", "dev", "benchmarks", "benchmark_helpers.R")
  skip_if(!file.exists(helpers), "benchmark helpers are not available")
  source(helpers, local = TRUE)
  g <- make_benchmark_cds(600); bg <- make_best_codon_benchmark_cds(600)

  set.seed(1); x1 <- remove_pirna_sites(g, "Ubiq", 3)
  expect_identical(rlang::hash(x1$sequence), "98327b40119c141c8e8845302e8b1314")
  expect_identical(rlang::hash(x1), "28b74990a261263e45070b60a71acbca")
  set.seed(2); x2 <- remove_pirna_sites(g, "Ubiq", 3)
  expect_identical(rlang::hash(x2$sequence), "3fe1aa7885a419f2359881c8816b8a77")
  expect_identical(rlang::hash(x2), "550ef93f774bf283e1d0fa4c4d9f9f46")
  xc <- remove_pirna_sites(bg, "CAI1", 3)
  expect_identical(rlang::hash(xc$sequence), "1c5cc88a0b9a3285f8477b1fde862a3c")
  expect_identical(rlang::hash(xc), "904c4ab56f91cbddd05b3b7a803bf0af")
  expect_identical(rlang::hash(get_pirna_distances(g)), "049a8815734d5e492c7977e47daeaac7")
})

test_that("the shield reports distances that match the original method", {
  set.seed(1)
  x <- remove_pirna_sites(pirna_fixture(), "Ubiq", 3)
  expect_identical(x$initial_dists, reference_get_distances(pirna_fixture()))
  expect_identical(x$final_dists, reference_get_distances(x$sequence))
})

# ---- several cores ---------------------------------------------------------------------------------

test_that("the result does not depend on the number of cores", {
  skip_on_os("windows")
  s <- paste0(planted_dna(1), planted_dna(2), planted_dna(3), random_dna(300, 5))   # about 900 windows
  starts <- seq_len(nchar(s) - 19)
  kmers <- substring(s, starts, starts + 19)
  one <- scan_kmers_against_pirnas(kmers, threshold = 4, cores = 1)
  for (k in c(2, 3, 7)) {
    many <- scan_kmers_against_pirnas(kmers, threshold = 4, cores = k)
    expect_identical(many$col_min, one$col_min, info = paste(k, "cores"))
    expect_identical(many$row_min, one$row_min, info = paste(k, "cores"))
    expect_identical(many$pairs, one$pairs, info = paste(k, "cores"))
  }
  expect_gt(nrow(one$pairs), 0)  # the comparison is only meaningful if there are pairs to compare
})

test_that("more cores than windows is harmless", {
  kmers <- substring(planted_dna(1), 1:5, 20:24)
  expect_identical(scan_kmers_against_pirnas(kmers, cores = 64)$col_min,
                   scan_kmers_against_pirnas(kmers, cores = 1)$col_min)
})

test_that("the multi-core scan still matches the original method", {
  skip_on_os("windows")
  s <- paste0(planted_dna(1), planted_dna(2), planted_dna(3), random_dna(300, 5))
  starts <- seq_len(nchar(s) - 19)
  many <- scan_kmers_against_pirnas(substring(s, starts, starts + 19), cores = 4)
  expect_identical(many$col_min, as.numeric(reference_get_distances(s)))
})

test_that("small scans use one core and large scans use several", {
  expect_equal(pirna_scan_cores(50), 1)
  expect_equal(pirna_scan_cores(299), 1)
  withr::local_options(wormtransgenebuilder.cores = 4)
  skip_on_os("windows")
  expect_equal(pirna_scan_cores(3000), 4)
  expect_equal(pirna_scan_cores(450), 3)          # never fewer than 150 windows per core
  withr::local_options(wormtransgenebuilder.cores = 1)
  expect_equal(pirna_scan_cores(3000), 1)
})

test_that("the number of CPU cores is looked up once, not on every scan", {
  # parallel::detectCores() asks the operating system and took about 7 ms per call; the
  # shield scans about 160 times per run, so it must not be asked each time.
  skip_on_os("windows")
  calls <- 0
  testthat::local_mocked_bindings(
    detectCores = function(...) { calls <<- calls + 1; 8L },
    .package = "parallel"
  )
  expect_equal(available_cores(refresh = TRUE), 7L)   # all cores but one
  for (i in 1:50) pirna_scan_cores(5000)
  for (i in 1:50) available_cores()
  expect_equal(calls, 1)
  withr::defer(available_cores(refresh = TRUE))       # do not leave the pretend value behind
})

test_that("small scans do not even ask how many cores there are", {
  skip_on_os("windows")
  asked <- 0
  testthat::local_mocked_bindings(available_cores = function(...) { asked <<- asked + 1; 8L })
  for (n in c(1, 22, 44, 299)) expect_equal(pirna_scan_cores(n), 1)
  expect_equal(asked, 0)
})

test_that("available_cores is at least 1, and one fewer than the machine has when it has several", {
  n <- available_cores(refresh = TRUE)
  expect_gte(n, 1)
  expect_equal(n, max(1L, parallel::detectCores() - 1L, na.rm = TRUE))
})

test_that("a scan on several cores uses fewer cores than there are when the sequence is short", {
  withr::local_options(wormtransgenebuilder.cores = 8)
  skip_on_os("windows")
  expect_lte(pirna_scan_cores(600), 4)
})

# ---- fewer passes over the database -------------------------------------------------------------

test_that("the random shield scans the piRNA database twice, where it used to do three full passes", {
  calls <- 0
  real_scan <- scan_kmers_against_pirnas
  testthat::local_mocked_bindings(
    scan_kmers_against_pirnas = function(...) { calls <<- calls + 1; real_scan(...) }
  )
  set.seed(1)
  remove_pirna_sites(pirna_fixture(), "Ubiq", 3)
  expect_equal(calls, 2)  # one scan of the input (distances and hits together), one of the result
})

test_that("get_pirna_distances scans once", {
  calls <- 0
  real_scan <- scan_kmers_against_pirnas
  testthat::local_mocked_bindings(
    scan_kmers_against_pirnas = function(...) { calls <<- calls + 1; real_scan(...) }
  )
  get_pirna_distances(pirna_fixture())
  expect_equal(calls, 1)
})
