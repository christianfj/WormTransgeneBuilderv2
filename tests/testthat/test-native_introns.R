# The native introns of the app (rpsa-1, rps-5, rps-20 in Introns.csv) against the gene models downloaded from the Alliance of Genome Resources
# (version 9.1.0, 2026-10-11; data-raw/intron_genes.txt). Written up in docs/DATA.md. Skipped when the file is not there.

test_that("the native introns of rpsa-1, rps-5 and rps-20 are the introns of the current gene models, in order", {
  path <- testthat::test_path("..", "..", "data-raw", "intron_genes.txt")
  skip_if(!file.exists(path), "the gene model file is not available (not running from the source tree)")
  x <- readLines(path)
  headers <- x[grepl("^>", x)]
  seqs <- x[!grepl("^>", x)]
  app <- utils::read.csv(system.file("extdata", "Introns.csv", package = "WormTransgeneBuilder"), stringsAsFactors = FALSE)
  expect_equal(sub("^>([^_]+)_.*", "\\1", headers), c("rpsa-1", "rps-5", "rps-20"))
  for (i in seq_along(seqs)) {
    gene <- sub("^>([^_]+)_.*", "\\1", headers[i])
    introns <- regmatches(seqs[i], gregexpr("[acgt]+", seqs[i]))[[1]]
    row <- app[app$Name == gene, ]
    expect_equal(nrow(row), 1, info = gene)
    expect_identical(introns, unlist(row[, c("Sequence_1", "Sequence_2", "Sequence_3")], use.names = FALSE), info = gene)
    expect_match(introns, "^gt.*ag$")
    cds <- toupper(gsub("[acgt]+", "", seqs[i]))
    expect_equal(nchar(cds) %% 3, 0, info = gene)
    expect_identical(substr(cds, 1, 3), "ATG")
    expect_true(substr(cds, nchar(cds) - 2, nchar(cds)) %in% c("TAA", "TAG", "TGA"), info = gene)
  }
})

test_that("the menu has the renamed gene rpsa-1 (formerly rps-0) and not the old name", {
  app <- utils::read.csv(system.file("extdata", "Introns.csv", package = "WormTransgeneBuilder"), stringsAsFactors = FALSE)
  expect_true("rpsa-1" %in% app$Name)
  expect_false("rps-0" %in% app$Name)
})
