# Removing the commonly used enzymes (common_enzymes in helper-fixtures.R; the set the laboratory
# avoids in a transgene) from a real-sized gene, in both modes. Removing every enzyme of the full
# 234-enzyme list would make no biological sense (many are 4-base cutters that occur every few
# hundred bases), so this is the realistic case to test.

common_sites <- function() {
  unique(unlist(lapply(common_enzymes, function(e) {
    s <- WormTransgeneBuilder::enzy[e, "Site"]
    c(s, reverse_complement(s))
  })))
}
sites_left <- function(dna) nrow(find_enzyme_site_ranges(dna, common_sites()))

bench_helpers <- function() {
  path <- testthat::test_path("..", "..", "dev", "benchmarks", "benchmark_helpers.R")
  skip_if(!file.exists(path), "benchmark helpers are not available")
  env <- new.env()
  source(path, local = env)
  env
}

test_that("every site of the common enzymes is removed from a typical gene, or a warning says which cannot be", {
  b <- bench_helpers()
  for (n in c(300, 3000)) {
    for (mode in c("Ubiq", "CAI1")) {
      gene <- if (mode == "CAI1") b$make_best_codon_benchmark_cds(n) else b$make_benchmark_cds(n)
      expect_gt(sites_left(gene), 0)               # the gene really contains sites to remove
      set.seed(1)
      warned <- FALSE
      out <- withCallingHandlers(remove_enzyme_sites(gene, common_enzymes, mode),
                                 warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") })
      info <- paste(n, "bp", mode)
      expect_equal(protein_of(out), protein_of(gene), info = info)
      expect_equal(nchar(out), nchar(gene), info = info)
      expect_lt(sites_left(out), sites_left(gene), label = info)
      if (!warned) expect_equal(sites_left(out), 0, info = info)   # no warning means none left
    }
  }
})

test_that("removing the common enzymes changes only a small part of the gene", {
  b <- bench_helpers()
  gene <- b$make_benchmark_cds(3000)
  set.seed(1)
  out <- suppressWarnings(remove_enzyme_sites(gene, common_enzymes, "Ubiq"))
  changed <- sum(strsplit(out, "")[[1]] != strsplit(gene, "")[[1]])
  expect_lt(changed / nchar(gene), 0.05)
})

test_that("High expression removes the common enzymes reproducibly", {
  b <- bench_helpers()
  gene <- b$make_best_codon_benchmark_cds(3000)
  set.seed(1); a <- suppressWarnings(remove_enzyme_sites(gene, common_enzymes, "CAI1"))
  set.seed(99); d <- suppressWarnings(remove_enzyme_sites(gene, common_enzymes, "CAI1"))
  expect_identical(a, d)
})

test_that("the common enzymes can be removed through the whole pipeline", {
  cds <- glo_cds()
  for (m in c("Ubiq", "CAI1", "GLO")) {
    set.seed(1)
    res <- suppressWarnings(optimize_transgene(cds, opt_method = m, enzymes = common_enzymes, remove_pirnas = FALSE))
    expect_equal(protein_of(res$optimized_cds), protein_of(cds), info = m)
    expect_equal(sites_left(res$optimized_cds), 0, info = m)
  }
})

test_that("a site that no synonymous change can break is reported, not hidden", {
  # HpyCH4III (ACNGT) in a Thr-Val pair: every Thr codon is ACN and every Val codon is GTN.
  cds <- "ATGGCTACTGTTAAATAA"
  expect_true(grepl("ACTGT", cds))
  expect_warning(out <- remove_enzyme_sites(cds, "HpyCH4III", "Ubiq"), "could not be removed")
  expect_true(grepl("AC[ACGT]GT", out))
  expect_warning(remove_enzyme_sites(cds, "HpyCH4III", "CAI1"), "could not be removed")
})

test_that("the app removes and highlights plain, ambiguity-code and run-of-any-base sites together", {
  # XhoI (CTCGAG), BanI (GGYRCC) and SfiI (GGCCNNNNNGGCC, a run of any five bases in the middle),
  # each present once in this CDS: ATG CTC GAG GGC GCC GGC CAA AAA GGC CGA TAA.
  cds <- "ATGCTCGAGGGCGCCGGCCAAAAAGGCCGATAA"
  enz <- c("XhoI", "BanI", "SfiI")
  sites <- unique(unlist(lapply(enz, function(e) { s <- WormTransgeneBuilder::enzy[e, "Site"]; c(s, reverse_complement(s)) })))
  expect_equal(nrow(find_enzyme_site_ranges(cds, sites)), 3)

  d <- run_server_downloads("1", cds, enzymes = enz)
  original <- parse_genbank(d$ori)
  expect_true(all(c("XhoI site", "BanI site", "SfiI site") %in% original$labels))   # all three are highlighted
  expect_equal(original$sequence, cds)

  optimized <- toupper(d$viewer_seq)
  expect_equal(protein_of(optimized), protein_of(cds))
  expect_equal(nrow(find_enzyme_site_ranges(optimized, sites)), 0)                   # all three are gone
})
