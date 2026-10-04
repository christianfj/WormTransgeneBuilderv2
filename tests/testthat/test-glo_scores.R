# The GLO score of a 12-mer is found by arithmetic on its bases (R/fct_glo_scores.R), not by
# searching a dictionary whose names are the 12-mers. These tests check that the new score
# table is exactly the old dictionary's values, that a window's position is worked out
# correctly, and that calc_glo_score() and optimize_glo() give bit-for-bit the same answers
# as the versions that used the old dictionary.

# ---- the table --------------------------------------------------------------------------

test_that("the score table has one non-negative score for each of the 4^12 possible 12-mers", {
  tbl <- glo_score_table()
  expect_type(tbl, "double")
  expect_length(tbl, 4^12)
  expect_false(anyNA(tbl))
  expect_true(all(tbl >= 0))
  expect_null(names(tbl))
})

test_that("the score table is the same object every time, read once", {
  expect_identical(glo_score_table(), glo_score_table())
})

test_that("the score table holds exactly the values that were in the original dictionary", {
  # The original dictionary (glo_full, a named vector of all 4^12 12-mers) was removed once
  # the score table replaced it. This is the hash of its values, recorded when the table was
  # made (dev/make_glo_scores.R), and it was verified then that the table was identical to it,
  # entry for entry, and that its names were exactly the base-4 numbers 0, 1, 2, ...
  expect_identical(rlang::hash(glo_score_table()), "d6a73d13333dc7d215434c189192e18b")
})

# 195 scores read straight from the original dictionary by 12-mer name before it was removed:
# the first, second and last entries, the largest and smallest non-zero scores, and a random
# selection of others, about a third of them zero. They check independently that a 12-mer is
# matched to the right entry of the table.
known_glo_scores <- r"(AAAAAAAAAAAA=53.719999999999999
AAAAAAAAAAAC=0
AAAAAAAAAAAG=10.050000000000001
AAAAAAAAAAAT=15.24
AAAAAAAAAACA=0
AAAAAAAAAATT=16.57
AAAAAAAAACAA=74
AAAAAAAAATTT=23.780000000000001
AAAAAAAACAAA=0.27000000000000002
AAACAAGGACGG=0.02
AAACAGAGCCAA=51.960000000000001
AAACGTGATCAA=0.40000000000000002
AAAGTTTGTGAC=94.069999999999993
AAATTCGAGCAT=636.46000000000004
AACGCGTTTCCA=122.31
AACGCTAGATAG=0
AACTGGATGGCA=0.71999999999999997
AACTTGGCCCAC=185.81
AAGCACTCTACT=684.01999999999998
AAGCCACATCTT=8.7699999999999996
AAGCTTCTTGAT=71.780000000000001
AAGTCCGTTGAG=46995.82
AAGTGTTCTGAG=17.629999999999999
AAGTTGATACGG=0
AAGTTTGCTGAT=87.269999999999996
AATACCATTTCT=0
AATACGTTTTTT=0
AATAGTAGTATT=40.689999999999998
AATCAACTGATT=43.969999999999999
AATCGACACTTG=23.93
AATCTTCGACTG=58.420000000000002
AATGGTTTTTGG=21.98
ACAAATGACGTG=33.259999999999998
ACAAGTGAAGAC=58.420000000000002
ACACAAAAGAAG=24.510000000000002
ACACGGAATATA=0
ACACTAACAACG=0
ACAGACATAAAG=0
ACCCTTGAATAC=1.1799999999999999
ACCGCGCATGAG=0
ACCGGAACGATT=0.34000000000000002
ACGAGAATGCTG=38.740000000000002
ACGAGACGATGG=1.29
ACGAGATCCGTT=0.11
ACGAGCACCAGA=0
ACGTGGCTTACC=22.899999999999999
ACGTTCACGTCC=0
ACTATTATTCTA=3.3799999999999999
ACTGAAACTTTT=14.67
ACTGACAAAAAT=5.1399999999999997
AGAAAATATAAC=0.35999999999999999
AGAAATGTTATG=5.3700000000000001
AGAGACATAGTG=0
AGATGACACGGA=0
AGCCCAGATCTC=0.20999999999999999
AGCGGAACGTAC=6
AGGAATCCCACT=2.3900000000000001
AGGTGATTCAGG=0
AGTACTTCGCAC=0
AGTGCATTAGAT=2.1499999999999999
ATACAGTAAGAA=0
ATAGAGCTCCGG=0
ATCAAGCTTTTC=99.409999999999997
ATCCACGGACTA=21.41
ATCTTCCAGATT=20.789999999999999
ATGCTCGATCCT=0.14000000000000001
ATGGTTTTGAAA=72.650000000000006
ATTCAGTGCTTT=69.239999999999995
ATTCCTGAACTT=19.73
ATTCGAATCTTC=45.310000000000002
ATTCTTGCATAC=6.2999999999999998
ATTGATTTCTCA=71.640000000000001
ATTGCCTGGATC=0
ATTGTAGCAATC=4.0999999999999996
ATTTTTTTTTTT=3.2999999999999998
CAAAAAAAAAAA=99.209999999999994
CAAAGGAAGCAA=0
CAACTCCGTGCC=97.569999999999993
CAATTCCATCTT=3.9300000000000002
CACACATCTGGC=0
CACACTGTAATA=37.57
CACCACAATCAC=11.26
CAGGCCGCGCCA=0
CATCCGTCAATA=4.4000000000000004
CATCTGAAATCT=24.649999999999999
CATGAAAATCAG=12.640000000000001
CATGTGTGATCA=0
CATTCAAAACAG=0.31
CCAAAGCATGAA=17.859999999999999
CCATTGTTTGTT=5
CCGATGTACAGA=71.530000000000001
CCGGTCCTTGTG=0
CCGTATACACCA=14.44
CCTCACATCTAC=108.16
CCTTGGCCTGTA=9.8499999999999996
CGAACAGCTAAT=4.7199999999999998
CGAGGTCAGAAA=0
CGATACTACATC=476.52999999999997
CGCAAGTCACAA=18.109999999999999
CGCAGGAGTTAA=0
CGCCAGGATATG=0.65000000000000002
CGGACTGGTTAG=0
CGGCATTAGAGA=0
CGGGAATAATCA=0
CGTACCGCCGCG=11.41
CGTGGCATGTTC=76.379999999999995
CGTGGGGCAGCG=0
CTATCCATTCTA=51.270000000000003
CTATTTTTTACT=0
CTCAACCACGAT=19.890000000000001
CTCAGTTTCTCG=0
CTCCGCATGGCT=137
CTCCTTGCTCAA=0
CTGCTTCGCGTG=0
CTGGGCCTGGTC=0
CTGTTCGCGAGA=24.390000000000001
CTTAACGATTGT=180.24000000000001
CTTCAAGAGCTC=2102
CTTCTCTGCTAT=213.74000000000001
CTTGAACGATCG=6.6900000000000004
CTTTTAATGAAT=341.93000000000001
GAACCAGCAGAT=4.3099999999999996
GAACCCATCTAG=0
GAAGAAGGTCGA=2.3100000000000001
GAAGGCCCAGCT=15.289999999999999
GAATATGTAAAA=242.36000000000001
GACAATACTCTT=57.969999999999999
GACATCCTCAAT=22.199999999999999
GACATCGCATGT=0
GACGCGAACTTT=18.149999999999999
GACGCGGCAGGG=0
GAGAAGCGTATC=246.69999999999999
GAGCAATGTGCC=46.700000000000003
GAGGCGAGAGAT=0.75
GAGGGCTTGAAA=36.350000000000001
GATATCCACTAG=0
GATGAATTGAGG=0.72999999999999998
GATTAAATGGTG=0
GCCAAGGCAGGG=0
GCGCCGACCTGG=0
GCTAAGGTTCTC=11.24
GGACTAGGTGCA=15.460000000000001
GGAGTCACTATC=0
GGATGCGTTTGC=0
GGCGGTCGTCAT=17.699999999999999
GGCGTTGCTGTC=0.5
GGGCAGTCGTAG=0
GGTCCGGTATAG=0
GGTGTCGCGTTG=5.7400000000000002
GGTTACGATGAA=1.45
GTAAAAGCCATT=218.03999999999999
GTACTGCAAATG=281.43000000000001
GTAGCTTCTTTG=15.98
GTATTGAAAACC=0.39000000000000001
GTCACTGATGAT=91.980000000000004
GTGAGGGTATTG=0
GTGCCCAATTGG=0.23000000000000001
GTGTTTTCCTAC=0
GTTAAAAACTGC=57.289999999999999
GTTCCAATTAAT=1.72
GTTCTCCTCTGT=58.240000000000002
GTTTTAAACGCG=113.17
TACACCCTGCTT=0
TACATTCTCATT=0
TATAAGGAGAAC=79.269999999999996
TCAACGCCATCT=192.28999999999999
TCAGTTACCGGT=0.25
TCATTGAACGGA=151.09999999999999
TCCCCGTGGTTC=106.98
TCCCCTGGGGCA=0
TCCCGCACCCAA=41.68
TCGGTCCAACAT=13.32
TCGTTGGCAGTA=1.8999999999999999
TCTGTCAGTTTT=0.96999999999999997
TCTTGGCGGCAG=0
TCTTGTGGCTAT=0
TGACTAGAGAAC=0
TGACTCATACGC=0
TGACTTGGCTGC=0
TGCGCTAGGACC=0
TGCGTAGTTTTA=0
TGGCACTCGAAC=0
TGTGCATTCTGT=1.8999999999999999
TGTTAGGGTACT=0
TGTTCATACATG=0
TTAGCACCAGAT=80.430000000000007
TTAGCCATGAAA=3.3500000000000001
TTAGGATGGCAC=4.9500000000000002
TTCCGTGTACGC=70.430000000000007
TTGCGGACTGCG=0
TTGGCTCCTGGA=15.84
TTGGGTGGATCA=4.3700000000000001
TTTTTGATTCAA=218.91999999999999
TTTTTTTTTTTG=0
TTTTTTTTTTTT=57.030000000000001)"

known_glo_pairs <- function() {
  parts <- strsplit(strsplit(known_glo_scores, "\n")[[1]], "=", fixed = TRUE)
  data.frame(window = vapply(parts, `[`, "", 1), score = as.numeric(vapply(parts, `[`, "", 2)), stringsAsFactors = FALSE)
}

# An independent way to find a window's entry: read its bases as a base-4 number.
table_score_by_base4 <- function(window) {
  glo_score_table()[strtoi(chartr("ACGT", "0123", window), base = 4L) + 1L]
}

test_that("the score table gives the known score for each of 195 12-mers read from the original dictionary", {
  known <- known_glo_pairs()
  expect_equal(nrow(known), 195)
  expect_gt(sum(known$score != 0), 100)
  expect_identical(table_score_by_base4(known$window), known$score)
})

test_that("the window lookup gives those same scores, all in one long sequence", {
  known <- known_glo_pairs()
  chars <- strsplit(paste(known$window, collapse = ""), "")[[1]]
  ours <- glo_window_scores(glo_base_codes(chars), (seq_len(nrow(known)) - 1) * 12 + 1)
  expect_identical(ours, known$score)
})

test_that("the first, second and last 12-mers are entries 1, 2 and the last one", {
  known <- known_glo_pairs()
  expect_identical(glo_score_table()[1], known$score[known$window == "AAAAAAAAAAAA"])
  expect_identical(glo_score_table()[2], known$score[known$window == "AAAAAAAAAAAC"])
  expect_identical(glo_score_table()[4^12], known$score[known$window == "TTTTTTTTTTTT"])
})

# ---- finding a window -----------------------------------------------------------------------

test_that("bases are numbered A = 0, C = 1, G = 2, T = 3, and anything else has no number", {
  expect_identical(glo_base_codes(c("A", "C", "G", "T")), 0:3)
  expect_identical(glo_base_codes(c("N", "a", " ", "X")), rep(NA_integer_, 4))
})

test_that("the first and last 12-mers are the first and last entries", {
  tbl <- glo_score_table()
  first <- glo_window_scores(glo_base_codes(strsplit("AAAAAAAAAAAA", "")[[1]]), 1)
  second <- glo_window_scores(glo_base_codes(strsplit("AAAAAAAAAAAC", "")[[1]]), 1)
  last <- glo_window_scores(glo_base_codes(strsplit("TTTTTTTTTTTT", "")[[1]]), 1)
  expect_identical(first, tbl[1])
  expect_identical(second, tbl[2])
  expect_identical(last, tbl[4^12])
})

test_that("the vectorised window lookup agrees with the one-window-at-a-time base-4 lookup", {
  set.seed(11)
  windows <- vapply(1:3000, function(i) paste(sample(c("A", "C", "G", "T"), 12, replace = TRUE), collapse = ""), "")
  # All windows in one long sequence, one after another.
  chars <- strsplit(paste(windows, collapse = ""), "")[[1]]
  ours <- glo_window_scores(glo_base_codes(chars), (seq_along(windows) - 1) * 12 + 1)
  expect_identical(ours, table_score_by_base4(windows))
})

test_that("overlapping windows are scored one codon apart or one base apart correctly", {
  s <- strsplit(glo_cds(), "")[[1]]
  codes <- glo_base_codes(s)
  starts <- c(1, 2, 3, 4, 10)
  expected <- table_score_by_base4(vapply(starts, function(p) paste(s[p:(p + 11)], collapse = ""), ""))
  expect_identical(glo_window_scores(codes, starts), expected)
})

test_that("a window with a base that is not A/C/G/T gets the median score", {
  s <- strsplit("ACGTNCGTACGTACGTACGTAAAAAAAAAAAA", "")[[1]]
  scores <- glo_window_scores(glo_base_codes(s), c(1, 5, 9, 21))
  expect_equal(scores[1], WormTransgeneBuilder::glo_median_score)  # contains N
  expect_equal(scores[2], WormTransgeneBuilder::glo_median_score)  # starts on the N
  expect_false(is.na(scores[3]))
  expect_identical(scores[4], glo_score_table()[1])                # AAAAAAAAAAAA
})

test_that("no windows gives no scores", {
  expect_identical(glo_window_scores(0:3, integer(0)), numeric(0))
})

# ---- calc_glo_score and optimize_glo give the same answers as before ------------------------------
# Reference values produced by the previous implementation (binary search over the named
# dictionary), captured before it was replaced.

reference_glo <- r"(score|cds|0|1988.3799999999999
opt|cds|1|ATGGCTAAAGGCCTCTCGGAGCCAACTGTCAATCAAATCTACTGGCGGAAAGACGGAGAAGCATTCCCCAACTTCCTTGGACTTTGCTGGAAAGCCGAAGAGTGTTAA
opt|cds|2|ATGGCCAAGGGACTCTCCGAACCCACCGTTAATCAAATCTACTGGCGGAAAGACGGAGAAGCATTTCCAAATTTTCTGGGTCTCTGCTGGAAAGCCGAAGAGTGTTAA
opt|cds|3|ATGGCAAAGGGATTATCTGAACCGACCGTTAATCAAATCTACTGGAGAAAAGACGGAGAAGCATTTCCAAATTTTCTGGGCCTCTGTTGGAAGGCCGAAGAATGCTGA
score|lower|0|1988.3799999999999
opt|lower|1|ATGGCTAAAGGCCTCTCGGAGCCAACTGTCAATCAAATCTACTGGCGGAAAGACGGAGAAGCATTCCCCAACTTCCTTGGACTTTGCTGGAAAGCCGAAGAGTGTTAA
opt|lower|2|ATGGCCAAGGGACTCTCCGAACCCACCGTTAATCAAATCTACTGGCGGAAAGACGGAGAAGCATTTCCAAATTTTCTGGGTCTCTGCTGGAAAGCCGAAGAGTGTTAA
opt|lower|3|ATGGCAAAGGGATTATCTGAACCGACCGTTAATCAAATCTACTGGAGAAAAGACGGAGAAGCATTTCCAAATTTTCTGGGCCTCTGTTGGAAGGCCGAAGAATGCTGA
score|withN|0|1840.28
opt|withN|1|ATGGCTAAAGGCCTCTCGGAACCCACCGTTAATCAAATCTACTGGNNNAAAGACGGAGAAGCATTTCCAAATTTTTTGGGACTATGCTGGAAAGCCGAAGAGTGTTAA
opt|withN|2|ATGGCGAAGGGTTTGTCTGAACCCACCGTTAATCAAATCTACTGGNNNAAAGATGGTGAAGCTTTTCCAAATTTTCTAGGTCTCTGCTGGAAAGCCGAAGAGTGTTAA
opt|withN|3|ATGGCTAAAGGATTATCTGAACCGACCGTTAATCAAATCTACTGGNNNAAAGACGGAGAAGCATTTCCAAATTTTCTGGGCCTCTGTTGGAAGGCCGAAGAATGCTGA
score|odd_len|0|1988.3799999999999
opt|odd_len|1|ATGGCTAAAGGCCTCTCGGAGCCAACTGTCAATCAAATCTACTGGCGGAAAGACGGAGAAGCATTCCCCAACTTCCTTGGACTTTGCTGGAAAGCCGAAGAGTGTTAAAC
opt|odd_len|2|ATGGCCAAGGGACTCTCCGAACCCACCGTTAATCAAATCTACTGGCGGAAAGACGGAGAAGCATTTCCAAATTTTCTGGGTCTCTGCTGGAAAGCCGAAGAGTGTTAAAC
opt|odd_len|3|ATGGCAAAGGGATTATCTGAACCGACCGTTAATCAAATCTACTGGAGAAAAGACGGAGAAGCATTTCCAAATTTTCTGGGCCTCTGTTGGAAGGCCGAAGAATGCTGAAC
score|odd_len2|0|1988.3799999999999
opt|odd_len2|1|ATGGCTAAAGGCCTCTCGGAGCCAACTGTCAATCAAATCTACTGGCGGAAAGACGGAGAAGCATTCCCCAACTTCCTTGGACTTTGCTGGAAAGCCGAAGAGTGTTAAA
opt|odd_len2|2|ATGGCCAAGGGACTCTCCGAACCCACCGTTAATCAAATCTACTGGCGGAAAGACGGAGAAGCATTTCCAAATTTTCTGGGTCTCTGCTGGAAAGCCGAAGAGTGTTAAA
opt|odd_len2|3|ATGGCAAAGGGATTATCTGAACCGACCGTTAATCAAATCTACTGGAGAAAAGACGGAGAAGCATTTCCAAATTTTCTGGGCCTCTGTTGGAAGGCCGAAGAATGCTGAA
score|spaced|0|1988.3799999999999
opt|spaced|1|ATGGCTAAAGGCCTTTCTGAACCCACTGTA AATCAAATCTTCTCGCATAAAGATGGTGATCCTTCTCCAAGTTCTCTGGAACCCTGTTGGAGAGCCGAAGAATGCTTA
opt|spaced|2|ATGGCTAAAGGATTATCTGAACCGACTGTT AATCCAATCTTCTTGCGTAAAGATGGTGATCATTTTCCAAGTTTTCTGGAACTCTCTTGGAATCGAGAAGAATGCTGA
opt|spaced|3|ATGGCAAAGGGCCTCTCGGAGCCAACTGTC AATCAAATCTTCTCGCTTAAAGATGGTGATCATTTTCCAAGTTTTCGGGAACCCTGTTGGAGAGCCGAAGAATGCTTA
score|gene300|0|2089.4400000000001
opt|gene300|1|c3c99382ca3397e93aeda50a24856c53
opt|gene300|2|bdc89d690c9c61de0e667c73354cd3bd
opt|gene300|3|36923b632fbeb6e056e6c12e4d60b62d
score|gene3000|0|59509.669999999976
opt|gene3000|1|47083201c73ba061c46f71c3c72b73c3
opt|gene3000|2|ad8452252066d70d5e7d7b21cfdb60b2
opt|gene3000|3|4d444481189e3bf89d52d7735e53b4d7
score|short12|0|0
opt|short12|1|ATGGCTAAAGGT
opt|short12|2|ATGGCTAAAGGT
opt|short12|3|ATGGCTAAAGGT
score|short9|0|0
opt|short9|1|ATGGCTAAA
opt|short9|2|ATGGCTAAA
opt|short9|3|ATGGCTAAA
score|exact4|0|0
opt|exact4|1|ATGGCTAAAGGTTAA
opt|exact4|2|ATGGCTAAAGGTTAA
opt|exact4|3|ATGGCTAAAGGTTAA
score_of_opt|gene3000|1|325715.54999999993)"

glo_reference_inputs <- function() {
  cds <- glo_cds()
  inputs <- list(
    cds = cds, lower = tolower(cds),
    withN = paste0(substr(cds, 1, 45), "NNN", substr(cds, 49, nchar(cds))),
    odd_len = paste0(cds, "AC"), odd_len2 = paste0(cds, "A"),
    spaced = paste0(substr(cds, 1, 30), " ", substr(cds, 31, nchar(cds))),
    short12 = "ATGGCTAAAGGT", short9 = "ATGGCTAAA", exact4 = "ATGGCTAAAGGTTAA"
  )
  helpers <- testthat::test_path("..", "..", "dev", "benchmarks", "benchmark_helpers.R")
  if (file.exists(helpers)) {
    source(helpers, local = TRUE)
    inputs$gene300 <- make_benchmark_cds(300)
    inputs$gene3000 <- make_benchmark_cds(3000)
  }
  inputs
}

test_that("GLO scores and optimized sequences are exactly what the previous implementation gave", {
  inputs <- glo_reference_inputs()
  ref <- strsplit(strsplit(reference_glo, "\n")[[1]], "|", fixed = TRUE)
  checked <- 0
  for (r in ref) {
    kind <- r[1]; id <- r[2]; seed <- as.integer(r[3]); value <- r[4]
    if (!id %in% names(inputs)) next  # the large benchmark genes need dev/benchmarks
    x <- inputs[[id]]
    if (kind == "score") {
      expect_identical(sprintf("%.17g", calc_glo_score(x)), value, info = paste(kind, id))
    } else if (kind == "opt") {
      set.seed(seed)
      out <- optimize_glo(x)
      expect_identical(if (nchar(out) > 200) rlang::hash(out) else out, value, info = paste(kind, id, "seed", seed))
    } else if (kind == "score_of_opt") {
      set.seed(seed)
      expect_identical(sprintf("%.17g", calc_glo_score(optimize_glo(x))), value, info = paste(kind, id))
    }
    checked <- checked + 1
  }
  expect_gt(checked, 35)
})

test_that("calc_glo_score adds the window scores in order, as it always did", {
  cds <- glo_cds()
  starts <- seq(1, nchar(cds) - 11, by = 3)
  windows <- substring(cds, starts, starts + 11)
  total <- 0
  for (w in windows) total <- total + table_score_by_base4(w)
  expect_identical(calc_glo_score(cds), total)
})

test_that("calc_glo_score ignores case, spaces and other characters, and is 0 for fewer than 4 codons", {
  cds <- glo_cds()
  expect_identical(calc_glo_score(tolower(cds)), calc_glo_score(cds))
  expect_identical(calc_glo_score(paste0(substr(cds, 1, 30), " \n", substr(cds, 31, nchar(cds)))), calc_glo_score(cds))
  expect_identical(calc_glo_score("ATGGCTAAA"), 0)
  expect_identical(calc_glo_score(""), 0)
})
