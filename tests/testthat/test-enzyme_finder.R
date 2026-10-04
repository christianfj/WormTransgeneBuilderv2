# find_enzyme_site_ranges() searches with regular expressions built from each site. (It used to
# call Biostrings::matchPattern, which spent about 15 ms per pattern building R objects; the
# rule it followed was measured for all 15 x 15 letter pairs and reproduced. When these tests
# were written the two were compared on thousands of random cases and gave identical matches.)
#
# The reference below is written straight from the DEFINITION of the rule, independently of
# the regular expressions: a site letter matches a sequence letter when every base the sequence
# letter stands for is one of the bases the site letter stands for. So R (A or G) matches A, G and
# a sequence R, N (any base) matches everything including a sequence N, and a plain A does not
# match a sequence R or N.

iupac_letters <- strsplit("ACGTRYSWKMBDHVN", "")[[1]]

bases_of <- list(A = "A", C = "C", G = "G", T = "T", R = c("A", "G"), Y = c("C", "T"), S = c("C", "G"),
                 W = c("A", "T"), K = c("G", "T"), M = c("A", "C"), B = c("C", "G", "T"),
                 D = c("A", "G", "T"), H = c("A", "C", "T"), V = c("A", "C", "G"), N = c("A", "C", "G", "T"))

# compatible[site_letter, sequence_letter] follows the definition above.
compatible <- outer(iupac_letters, iupac_letters,
                    Vectorize(function(site, seq) all(bases_of[[seq]] %in% bases_of[[site]])))
dimnames(compatible) <- list(iupac_letters, iupac_letters)

reference_site_ranges <- function(dna_seq, sites) {
  subject <- strsplit(gsub("[^ACGTRYSWKMBDHVN]", "N", dna_seq), "")[[1]]
  starts <- integer(0)
  ends <- integer(0)
  for (pat in sites) {
    site <- strsplit(toupper(pat), "")[[1]]
    len <- length(site)
    if (len > length(subject)) next
    windows <- length(subject) - len + 1
    ok <- rep(TRUE, windows)
    for (j in seq_len(len)) ok <- ok & compatible[cbind(site[j], subject[j:(j + windows - 1)])]
    hit <- which(ok)
    starts <- c(starts, hit)
    ends <- c(ends, hit + len - 1)
  }
  data.frame(start = as.integer(starts), end = as.integer(ends))
}

random_string <- function(n, alphabet) paste(sample(alphabet, n, replace = TRUE), collapse = "")

# ---- the matching rule ------------------------------------------------------------------------

test_that("every site letter against every sequence letter matches as the definition says", {
  for (p in iupac_letters) for (s in iupac_letters) {
    expect_identical(find_enzyme_site_ranges(s, p), reference_site_ranges(s, p), info = paste("site", p, "sequence", s))
  }
})

test_that("the rule is: a site letter matches the letters whose bases it includes", {
  hit <- function(site, seq) nrow(find_enzyme_site_ranges(seq, site)) == 1
  expect_true(hit("R", "A"));  expect_true(hit("R", "G"));  expect_false(hit("R", "C"))
  expect_true(hit("Y", "C"));  expect_true(hit("Y", "T"));  expect_false(hit("Y", "A"))
  expect_true(hit("N", "A"));  expect_true(hit("N", "T"));  expect_true(hit("N", "N"))
  expect_true(hit("A", "A"));  expect_false(hit("A", "R"))   # an ambiguous base in the sequence is not a plain A
  expect_true(hit("R", "R"))                                 # ... but it is inside the set R
  expect_false(hit("A", "N"))                                # an unknown base does not match a specific base
  expect_true(hit("B", "Y"));  expect_false(hit("B", "A"))
})

test_that("BanI (GGYRCC) is found in the four sequences it recognises and no others", {
  for (s in c("GGCGCC", "GGTACC", "GGCACC", "GGTGCC")) expect_equal(nrow(find_enzyme_site_ranges(s, "GGYRCC")), 1, info = s)
  for (s in c("GGAGCC", "GGCCCC", "GGNGCC", "GGCNCC")) expect_equal(nrow(find_enzyme_site_ranges(s, "GGYRCC")), 0, info = s)
})

# ---- the same answer as the reference on many random cases ---------------------------------------------

test_that("random sequences and sites give exactly the same matches as the reference", {
  set.seed(2026)
  for (i in 1:600) {
    subject <- random_string(sample(0:70, 1), c(iupac_letters, iupac_letters[1:4], iupac_letters[1:4], "a", "-", "*"))
    sites <- vapply(seq_len(sample(1:4, 1)), function(j) random_string(sample(1:8, 1), iupac_letters), "")
    expect_identical(find_enzyme_site_ranges(subject, sites), reference_site_ranges(subject, sites),
                     info = paste("subject", subject, "sites", paste(sites, collapse = ",")))
  }
})

test_that("sequences of plain A/C/G/T with ambiguity sites match the reference", {
  set.seed(7)
  for (i in 1:300) {
    subject <- random_string(sample(20:120, 1), c("A", "C", "G", "T"))
    sites <- vapply(1:3, function(j) random_string(sample(4:7, 1), c("A", "C", "G", "T", "R", "Y", "N", "S", "W")), "")
    expect_identical(find_enzyme_site_ranges(subject, sites), reference_site_ranges(subject, sites))
  }
})

test_that("all the real enzyme sites and their reverse complements match the reference on a real-sized gene", {
  helpers <- testthat::test_path("..", "..", "dev", "benchmarks", "benchmark_helpers.R")
  skip_if(!file.exists(helpers), "benchmark helpers are not available")
  source(helpers, local = TRUE)
  sites <- unique(unlist(lapply(rownames(WormTransgeneBuilder::enzy), function(e) {
    s <- WormTransgeneBuilder::enzy[e, "Site"]
    c(s, reverse_complement(s))
  })))
  for (g in list(make_benchmark_cds(3000), make_best_codon_benchmark_cds(3000), glo_cds())) {
    expect_identical(find_enzyme_site_ranges(g, sites), reference_site_ranges(g, sites))
  }
})

# ---- details ---------------------------------------------------------------------------------------------

test_that("overlapping matches are all found, in order", {
  expect_identical(find_enzyme_site_ranges("AAAAAA", "AAA"), data.frame(start = 1:4, end = 3:6))
  expect_identical(find_enzyme_site_ranges("ACACACA", "ACA"), data.frame(start = c(1L, 3L, 5L), end = c(3L, 5L, 7L)))
})

test_that("matches are listed site by site, in the order the sites were given", {
  out <- find_enzyme_site_ranges("GGATCCGAATTC", c("GAATTC", "GGATCC"))
  expect_identical(out$start, c(7L, 1L))
  expect_identical(out$end, c(12L, 6L))
})

test_that("a site given twice is reported twice, as before", {
  expect_equal(nrow(find_enzyme_site_ranges("GGATCC", c("GGATCC", "GGATCC"))), 2)
})

test_that("start and end are integers, and no match gives an empty data frame", {
  hit <- find_enzyme_site_ranges("AAGGATCCAA", "GGATCC")
  expect_type(hit$start, "integer")
  expect_type(hit$end, "integer")
  none <- find_enzyme_site_ranges("AAAAAAAA", "GGATCC")
  expect_identical(none, data.frame(start = integer(0), end = integer(0)))
  expect_identical(find_enzyme_site_ranges("", "GGATCC"), data.frame(start = integer(0), end = integer(0)))
  expect_identical(find_enzyme_site_ranges("GGATCC", character(0)), data.frame(start = integer(0), end = integer(0)))
})

test_that("characters that are not IUPAC letters (including lowercase) are treated as unknown bases", {
  expect_equal(nrow(find_enzyme_site_ranges("ggatcc", "GGATCC")), 0)
  expect_equal(nrow(find_enzyme_site_ranges("GGA-CC", "GGATCC")), 0)
  expect_equal(nrow(find_enzyme_site_ranges("GGA*CC", "GGANCC")), 1)   # ... which only an N in the site matches
  expect_identical(find_enzyme_site_ranges("aGGATCCa", "GGATCC"), reference_site_ranges("aGGATCCa", "GGATCC"))
})

test_that("a site in lowercase is accepted, as it always was", {
  expect_identical(find_enzyme_site_ranges("AAGGATCCAA", "ggatcc"), reference_site_ranges("AAGGATCCAA", "ggatcc"))
})

test_that("a site with a character that is not an IUPAC letter gives a clear error", {
  expect_error(find_enzyme_site_ranges("AAGGATCCAA", "GGAXCC"), "IUPAC")
  expect_error(find_enzyme_site_ranges("AAGGATCCAA", "GG CC"), "IUPAC")
})

test_that("find_enzyme_site_starts is the start column", {
  expect_identical(find_enzyme_site_starts("AAGGATCCAAGGATCC", "GGATCC"), c(3L, 11L))
})

# ---- the whole search gives the same answers as before ---------------------------------------------------------
# Reference results produced by the previous implementation (which used Biostrings), captured before it
# was replaced: the hash of the result, and for the mixed-enzyme example the sequence itself.

test_that("restriction-site removal gives exactly the results it gave before, in both modes", {
  mix <- "ATGCTCGAGAAAGGTCTCGGCGCCGAATTCGGTACCAAGCTTGGATCCTCTAGACTGCAGTAA"
  all_enz <- original_21_enzymes  # the reference values below were captured with these 21
  expect_identical(remove_enzyme_sites(mix, all_enz, "CAI1"),
                   "ATGCTTGAGAAAGGACTCGGAGCCGAGTTCGGAACCAAGCTCGGATCTTCTCGTCTTCAGTAA")
  expected <- c("1" = "ATGCTTGAGAAAGGCCTCGGAGCCGAGTTCGGAACCAAACTTGGTTCCTCAAGATTGCAGTAA",
                "2" = "ATGTTAGAGAAAGGACTCGGAGCCGAGTTCGGAACCAAACTTGGTTCCTCCAGATTGCAGTAA",
                "3" = "ATGTTGGAGAAAGGACTCGGAGCCGAGTTCGGAACCAAACTTGGTTCCTCAAGACTTCAGTAA")
  for (s in 1:3) {
    set.seed(s)
    expect_identical(remove_enzyme_sites(mix, all_enz, "Ubiq"), expected[[as.character(s)]], info = paste("seed", s))
  }
})

test_that("restriction-site removal on the benchmark genes gives exactly the results it gave before", {
  helpers <- testthat::test_path("..", "..", "dev", "benchmarks", "benchmark_helpers.R")
  skip_if(!file.exists(helpers), "benchmark helpers are not available")
  source(helpers, local = TRUE)
  all_enz <- original_21_enzymes  # the reference values below were captured with these 21
  changed <- function(a, b) sum(strsplit(a, "")[[1]] != strsplit(b, "")[[1]])
  refs <- list(
    "300" = list(cai1_all = c("1adfb9bbadf85578cd428e64c3a61dc7", 5), xhoi = c("115aa5cccc41d5815fe77ce6f7116cf4", 1),
                 gene_cai1 = c("f0d6ad473922b9589cc935bbc3be63de", 2), rand1 = c("f0d6ad473922b9589cc935bbc3be63de", 2),
                 rand2 = c("f0d6ad473922b9589cc935bbc3be63de", 2), bsai = c("ddab16188f116c8c42bca9969f9bd1a6", 1)),
    "3000" = list(cai1_all = c("f778a1c819a1f2450be2505fc709259a", 27), xhoi = c("171bfce923f2f1e8ab214e9167b4e15a", 3),
                  gene_cai1 = c("39173d8be313b4f17abe686c4bd5026f", 25), rand1 = c("0e8a57f5f0b7d10386367cff0d3b81d4", 26),
                  rand2 = c("203f8f9e1b92747fb085a1f52019f2f5", 27), bsai = c("ad1ac3fa0947a99622fc9c9e22e2a43e", 1)))
  for (n in c(300, 3000)) {
    b <- make_best_codon_benchmark_cds(n); g <- make_benchmark_cds(n)
    ref <- refs[[as.character(n)]]
    check <- function(out, input, key) {
      expect_identical(rlang::hash(out), ref[[key]][1], info = paste(n, key))
      expect_equal(changed(out, input), as.numeric(ref[[key]][2]), info = paste(n, key))
    }
    check(remove_enzyme_sites(b, all_enz, "CAI1"), b, "cai1_all")
    check(remove_enzyme_sites(b, "XhoI", "CAI1"), b, "xhoi")
    check(remove_enzyme_sites(g, all_enz, "CAI1"), g, "gene_cai1")
    set.seed(1); check(remove_enzyme_sites(g, all_enz, "Ubiq"), g, "rand1")
    set.seed(2); check(remove_enzyme_sites(g, all_enz, "Ubiq"), g, "rand2")
    set.seed(1); check(remove_enzyme_sites(g, "BsaI", "Chance"), g, "bsai")
  }
})
