# The benchmark script (dev/benchmarks/) is only useful if its test genes are valid and
# reproducible, and if its compare script reads results correctly. These tests load the
# helper file from dev/, and skip themselves if it is not there (for example in an
# installed copy of the package, where dev/ is not included).

helpers_path <- testthat::test_path("..", "..", "dev", "benchmarks", "benchmark_helpers.R")
compare_path <- testthat::test_path("..", "..", "dev", "benchmarks", "compare.R")
skip_if_no_benchmarks <- function() skip_if(!file.exists(helpers_path), "benchmark helpers are not available")
if (file.exists(helpers_path)) source(helpers_path, local = TRUE)

bench_stops <- c("TAA", "TAG", "TGA")

test_that("benchmark genes are valid coding sequences of about the requested length", {
  skip_if_no_benchmarks()
  for (n in c(300, 3000)) {
    g <- make_benchmark_cds(n)
    expect_true(check_cds_validity(g)$valid, info = n)
    expect_lte(abs(nchar(g) - n), 6)
    expect_equal(nchar(g) %% 3, 0)
    cd <- codons_of(g)
    expect_false(any(cd[-length(cd)] %in% bench_stops))
  }
})

test_that("benchmark genes are reproducible, and change with the seed", {
  skip_if_no_benchmarks()
  expect_identical(make_benchmark_cds(600), make_benchmark_cds(600))
  expect_false(identical(make_benchmark_cds(600), make_benchmark_cds(600, seed = 1)))
  expect_identical(make_best_codon_benchmark_cds(600), make_best_codon_benchmark_cds(600))
})

test_that("making a benchmark gene does not disturb R's random numbers", {
  skip_if_no_benchmarks()
  set.seed(99); expected <- runif(3)
  set.seed(99); make_benchmark_cds(600); actual <- runif(3)
  expect_identical(actual, expected)
})

test_that("benchmark genes contain real piRNA targets and restriction sites", {
  skip_if_no_benchmarks()
  small <- make_benchmark_cds(300)
  expect_gte(sum(get_pirna_distances(small) == 0, na.rm = TRUE), 1)
  expect_true(grepl("GGTCTC", small))
  medium <- make_benchmark_cds(3000)
  expect_gte(sum(get_pirna_distances(medium) == 0, na.rm = TRUE), 3)
  expect_true(grepl("CTCGAG", medium))
})

test_that("the best-codon gene is unchanged by High-expression recoding and keeps its targets", {
  skip_if_no_benchmarks()
  g <- make_best_codon_benchmark_cds(3000)
  expect_true(check_cds_validity(g)$valid)
  expect_identical(probabilistic_recode(g, "CAI1"), g)
  expect_gte(sum(get_pirna_distances(g) == 0, na.rm = TRUE), 3)
  expect_true(grepl("CTCGAG", g))
})

test_that("piRNA blocks are 21 bases with no stop codon, and best-only blocks are all best codons", {
  skip_if_no_benchmarks()
  blocks <- pirna_blocks(limit = 50)
  expect_true(all(nchar(blocks) == 21))
  expect_false(any(unlist(lapply(blocks, codons_of)) %in% bench_stops))
  best <- pirna_blocks(best_only = TRUE, limit = 2)
  cf <- WormTransgeneBuilder::codon_freqs
  best_codons <- toupper(cf$Codon[cf$CAI1 == 1])
  expect_true(all(unlist(lapply(best, codons_of)) %in% best_codons))
})

test_that("fingerprints are stable, and different for different results", {
  skip_if_no_benchmarks()
  x <- list(a = 1:3, b = "ATG")
  expect_identical(fingerprint(x), fingerprint(list(a = 1:3, b = "ATG")))
  expect_false(identical(fingerprint(x), fingerprint(list(a = 1:3, b = "ATC"))))
  expect_type(fingerprint(x), "character")
})

test_that("every benchmark case is well formed and the cheap ones run", {
  skip_if_no_benchmarks()
  cases <- benchmark_cases()
  ids <- vapply(cases, `[[`, "", "id")
  expect_equal(anyDuplicated(ids), 0)
  for (cs in cases) {
    expect_true(all(cs$sizes %in% names(benchmark_sizes)), info = cs$id)
    expect_true(is.function(cs$fn), info = cs$id)
  }
  x <- list(cds = make_benchmark_cds(300), best = make_best_codon_benchmark_cds(300))
  cheap <- c("calc_gc_content", "calc_cai", "validate_for_synthesis", "insert_introns_early",
             "insert_introns_equidistant", "generate_sequence_viewer", "export_genbank")
  for (cs in cases[ids %in% cheap]) expect_no_error(cs$fn(x))
})

test_that("every benchmark case has a size that the standard suite runs", {
  skip_if_no_benchmarks()
  cases <- benchmark_cases()
  expect_true(all(vapply(cases, function(cs) any(c("small", "medium") %in% cs$sizes), logical(1))))
})

test_that("compare.R reports speed-ups, memory and changed output", {
  skip_if(!file.exists(compare_path), "compare script is not available")
  make_results <- function(label, times, hashes) {
    data.frame(label = label, suite = "quick", timestamp = "2026-01-01 00:00:00", git_commit = "abc1234",
               git_state = "clean", package_version = "2.0.0", r_version = "4.5.3", platform = "test",
               cores = 8, rnafold = TRUE, id = c("case_a", "case_b"), size = "small", size_bp = 300,
               iterations = 3, median_s = times, min_s = times, mem_alloc_mb = c(10, 20), n_gc = 0,
               output_hash = hashes, status = "ok", stringsAsFactors = FALSE)
  }
  before <- withr::local_tempfile(fileext = ".csv"); after <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(make_results("before", c(2, 1), c("h1", "h2")), before, row.names = FALSE)
  utils::write.csv(make_results("after", c(1, 1), c("h1", "CHANGED-HASH")), after, row.names = FALSE)
  out <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(compare_path), shQuote(before), shQuote(after)),
                 stdout = TRUE, stderr = TRUE)
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "case_a \\[small\\]")
  expect_match(txt, "2.00x")            # 2 s -> 1 s
  expect_match(txt, "same")             # case_a output unchanged
  expect_match(txt, "CHANGED")          # case_b output changed
  expect_match(txt, "10.0 -> 10.0 MB")
  expect_no_match(txt, "WARNING")
})

test_that("compare.R warns when the two runs used different machines", {
  skip_if(!file.exists(compare_path), "compare script is not available")
  base <- data.frame(label = "x", suite = "quick", timestamp = "t", git_commit = "a", git_state = "clean",
                     package_version = "2.0.0", r_version = "4.5.3", platform = "mac", cores = 8, rnafold = TRUE,
                     id = "case_a", size = "small", size_bp = 300, iterations = 3, median_s = 1, min_s = 1,
                     mem_alloc_mb = 1, n_gc = 0, output_hash = "h", status = "ok", stringsAsFactors = FALSE)
  other <- base; other$cores <- 4
  f1 <- withr::local_tempfile(fileext = ".csv"); f2 <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(base, f1, row.names = FALSE); utils::write.csv(other, f2, row.names = FALSE)
  out <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(compare_path), shQuote(f1), shQuote(f2)), stdout = TRUE, stderr = TRUE)
  expect_match(paste(out, collapse = "\n"), "WARNING: the two runs used a different")
})

test_that("the benchmark script records a meaningful, repeatable fingerprint of each result", {
  # Runs the real script on a few cheap cases, twice. A fingerprint that were the same for
  # every case (as happened when the output was not captured properly) would make every
  # future comparison say "same output" whatever had changed.
  script_path <- testthat::test_path("..", "..", "dev", "benchmarks", "benchmark.R")
  skip_if(!file.exists(script_path), "benchmark script is not available")
  skip_on_cran()
  run_once <- function() {
    outdir <- withr::local_tempdir(.local_envir = parent.frame())
    out <- system2(file.path(R.home("bin"), "Rscript"),
                   c(shQuote(script_path), "--suite", "quick", "--label", "t", "--iterations", "1",
                     "--only", shQuote("calc_gc_content|calc_glo_score|insert_introns_early"),
                     "--outdir", shQuote(outdir)), stdout = TRUE, stderr = TRUE)
    file <- file.path(outdir, "bench_t.csv")
    expect_true(file.exists(file), info = paste(tail(out, 5), collapse = "\n"))
    utils::read.csv(file, stringsAsFactors = FALSE)
  }
  first <- run_once()
  second <- run_once()

  expect_equal(first$id, c("calc_gc_content", "calc_glo_score", "insert_introns_early"))
  expect_true(all(first$status == "ok"))
  expect_false(anyNA(first$output_hash))
  expect_equal(length(unique(first$output_hash)), 3)          # different cases, different fingerprints
  expect_identical(first$output_hash, second$output_hash)      # same code, same fingerprints
  expect_true(all(first$median_s > 0))
  expect_true(all(c("git_commit", "r_version", "cores", "package_version") %in% names(first)))
})
