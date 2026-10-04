# add_matches() feeds the sequence viewer: it returns the start/end positions
# (1-based, on the string that is displayed) of each feature to highlight.

enzyme_matches <- function(site, seq, tooltip = "Enz site") {
  add_matches(site, seq, "#d7bde2", tooltip, "RE_site", is_enzyme = TRUE)
}

test_that("an exact enzyme site is highlighted with its start and end", {
  res <- enzyme_matches("GGTCTC", "AAGGTCTCTT")
  expect_equal(res$starts, 3)
  expect_equal(res$ends, 8)
  expect_equal(res$tooltips, "Enz site")
  expect_equal(res$feats, "RE_site")
})

test_that("a site on the reverse strand is highlighted and labelled (rev)", {
  res <- enzyme_matches("GGTCTC", "TTGAGACCAA")  # GAGACC = reverse complement
  expect_equal(res$starts, 3)
  expect_equal(res$ends, 8)
  expect_equal(res$tooltips, "Enz site (rev)")
})

test_that("a palindromic site is reported once, not once per strand", {
  res <- enzyme_matches("GAATTC", "AAGAATTCTT")
  expect_equal(res$starts, 3)
  expect_length(res$tooltips, 1)
})

test_that("no match returns empty results", {
  res <- enzyme_matches("GGTCTC", "AAAAAAAAAA")
  expect_length(res$starts, 0)
  expect_length(res$ends, 0)
})

test_that("ambiguity-code sites are highlighted (BanI = GGYRCC)", {
  for (site in c("GGCGCC", "GGTACC", "GGCACC", "GGTGCC")) {
    res <- enzyme_matches("GGYRCC", paste0("AA", site, "TT"), "BanI site")
    expect_equal(res$starts, 3, info = site)
    expect_equal(res$ends, 8, info = site)
  }
})

test_that("bases outside the ambiguity code are not highlighted", {
  expect_length(enzyme_matches("GGYRCC", "AAGGAGCCTT")$starts, 0)
  expect_length(enzyme_matches("GGYRCC", "AAGGCCCCTT")$starts, 0)
})

test_that("BsoBI (CYCGRG) is highlighted, including CTCGAG", {
  res <- enzyme_matches("CYCGRG", "AACTCGAGTT", "BsoBI site")
  expect_equal(res$starts, 3)
})

test_that("a non-palindromic ambiguous site is found on the reverse strand", {
  # GGTCTNN reverse-complements to NNAGACC.
  res <- enzyme_matches("GGTCTNN", "TTCCAGACCTT")
  expect_equal(res$starts, 3)
  expect_equal(res$tooltips, "Enz site (rev)")
})

test_that("overlapping sites are all highlighted", {
  res <- enzyme_matches("GGYRCC", "GGCGCCGGTACC")
  expect_equal(sort(res$starts), c(1, 7))
  res_ol <- enzyme_matches("AAA", "AAAAAA")
  expect_equal(sort(res_ol$starts[res_ol$tooltips == "Enz site"]), 1:4)
})

test_that("an N in the sequence does not create a false highlight", {
  expect_length(enzyme_matches("GGYRCC", "AAGGNGCCTT")$starts, 0)
})

test_that("stray letters or * in the displayed sequence do not error or shift positions", {
  # The app only strips characters outside A-Z and *, so these can reach here.
  res <- enzyme_matches("GGYRCC", "EF*GGCGCCQQ")
  expect_equal(res$starts, 4)
  expect_equal(res$ends, 9)
})

test_that("lowercase input is matched", {
  expect_equal(enzyme_matches("ggyrcc", "aaggcgcctt")$starts, 3)
})

test_that("non-enzyme features still use literal matching", {
  # is_enzyme = FALSE: the pattern is a fixed string, so Y is not an ambiguity code.
  res <- add_matches("GGYRCC", "AAGGCGCCTT", "#000", "x", "promoter")
  expect_length(res$starts, 0)
  res2 <- add_matches("GGCGCC", "AAGGCGCCTT", "#000", "x", "promoter")
  expect_equal(res2$starts, 3)
  expect_equal(res2$ends, 8)
})

test_that("non-enzyme features do not get a reverse-strand search", {
  res <- add_matches("GGTCTC", "TTGAGACCAA", "#000", "x", "promoter")
  expect_length(res$starts, 0)
})

test_that("every real enzyme in the table can be searched without error", {
  enz <- read.csv(system.file("extdata", "Enzymes.csv", package = "WormTransgeneBuilder"),
                  stringsAsFactors = FALSE)
  expect_gt(nrow(enz), 0)
  for (i in seq_len(nrow(enz))) {
    expect_no_error(enzyme_matches(enz$Site[i], "ATGGCTAAAGGTACCGAATTCTAA"))
  }
})
