# generate_sequence_viewer() returns HTML with a JavaScript "coverage" list. Each
# entry styles a block of bases; start is 0-based and end is exclusive. These
# tests read that list back out of the generated HTML.

viewer_html <- function(...) as.character(generate_sequence_viewer("v", ...))

coverage_entries <- function(html) {
  m <- regmatches(html, gregexpr("\\{start: [0-9]+, end: [0-9]+, color: [^}]*\\}", html))[[1]]
  data.frame(
    start = as.integer(sub(".*start: ([0-9]+),.*", "\\1", m)),
    end = as.integer(sub(".*end: ([0-9]+),.*", "\\1", m)),
    bgcolor = sub('.*bgcolor: "([^"]*)".*', "\\1", m),
    underscore = grepl("underscore: true", m, fixed = TRUE),
    tooltip = sub('.*tooltip: "([^"]*)".*', "\\1", m),
    stringsAsFactors = FALSE
  )
}
legend_entries <- function(html) {
  m <- regmatches(html, gregexpr('\\{name: "[^"]*", color: "[^"]*", underscore: (true|false)\\}', html))[[1]]
  m
}

seq20 <- strrep("ACGT", 5)

test_that("empty sequences give NULL", {
  expect_null(generate_sequence_viewer("v", ""))
  expect_null(generate_sequence_viewer("v", NULL))
})

test_that("without underline arguments nothing is underlined", {
  html <- viewer_html(seq20, 1, 3, "#a9dfbf", "Start")
  e <- coverage_entries(html)
  expect_equal(nrow(e), 1)
  expect_false(any(e$underscore))
  expect_false(grepl("Non-optimal", html, fixed = TRUE))
})

test_that("an underline on its own becomes one underlined block with a legend entry", {
  html <- viewer_html(seq20, underline_starts = 4, underline_ends = 6)
  e <- coverage_entries(html)
  expect_equal(nrow(e), 1)
  expect_equal(e$start, 3)
  expect_equal(e$end, 6)
  expect_true(e$underscore)
  expect_equal(e$tooltip, "Non-optimal codon")
  expect_true(any(grepl('name: "Non-optimal codon".*underscore: true', legend_entries(html))))
})

test_that("neighbouring underlined codons are merged into a single block", {
  html <- viewer_html(seq20, underline_starts = c(4, 7), underline_ends = c(6, 9))
  e <- coverage_entries(html)
  expect_equal(nrow(e), 1)
  expect_equal(c(e$start, e$end), c(3, 9))
})

test_that("separate underlined codons stay separate", {
  html <- viewer_html(seq20, underline_starts = c(4, 13), underline_ends = c(6, 15))
  e <- coverage_entries(html)
  expect_equal(nrow(e), 2)
  expect_true(all(e$underscore))
})

test_that("an underline over a coloured feature keeps that feature's colour", {
  # Feature covers bases 1-9; underline covers 4-6. The feature is split into
  # three blocks and only the middle one is underlined.
  html <- viewer_html(seq20, 1, 9, "#d2b4de", "Tag", underline_starts = 4, underline_ends = 6)
  e <- coverage_entries(html)
  expect_equal(nrow(e), 3)
  expect_equal(e$bgcolor, rep("#d2b4de", 3))
  expect_equal(e$underscore, c(FALSE, TRUE, FALSE))
  expect_equal(e$tooltip, rep("Tag", 3))
  expect_equal(c(e$start[2], e$end[2]), c(3, 6))
})

test_that("an underline that partly overlaps a feature covers both correctly", {
  html <- viewer_html(seq20, 1, 5, "#d2b4de", "Tag", underline_starts = 4, underline_ends = 8)
  e <- coverage_entries(html)
  expect_equal(e$underscore, c(FALSE, TRUE, TRUE))
  expect_equal(e$bgcolor, c("#d2b4de", "#d2b4de", "transparent"))
})

test_that("underlined bases outside the sequence are ignored", {
  html <- viewer_html(seq20, underline_starts = c(19, 30), underline_ends = c(25, 35))
  e <- coverage_entries(html)
  expect_equal(nrow(e), 1)
  expect_equal(c(e$start, e$end), c(18, 20))
})

test_that("features and underlines together keep the normal legend as well", {
  html <- viewer_html(seq20, 1, 3, "#a9dfbf", "Start", underline_starts = 7, underline_ends = 9)
  leg <- legend_entries(html)
  expect_length(leg, 2)
  expect_true(any(grepl('name: "Start".*underscore: false', leg)))
  expect_true(any(grepl('name: "Non-optimal codon".*underscore: true', leg)))
})

test_that("the underline label can be changed", {
  html <- viewer_html(seq20, underline_starts = 1, underline_ends = 3, underline_label = "Not best codon")
  expect_equal(coverage_entries(html)$tooltip, "Not best codon")
})
