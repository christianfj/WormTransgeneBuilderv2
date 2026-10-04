#!/usr/bin/env Rscript
# Benchmark WormTransgeneBuilder: how long each step takes and how much memory it uses,
# on fixed, reproducible genes. Run it from anywhere:
#
#   Rscript dev/benchmarks/benchmark.R                       # standard suite (300 bp and 3 kb)
#   Rscript dev/benchmarks/benchmark.R --suite quick         # 300 bp only, a couple of minutes
#   Rscript dev/benchmarks/benchmark.R --suite full          # adds the 20 kb size (slow: many minutes)
#   Rscript dev/benchmarks/benchmark.R --label before-speedup
#   Rscript dev/benchmarks/benchmark.R --only "pirna|enzymes"   # only cases whose name matches
#   Rscript dev/benchmarks/benchmark.R --iterations 10          # repeat every case 10 times
#   Rscript dev/benchmarks/benchmark.R --outdir /some/folder    # write the results elsewhere
#
# Results go to dev/benchmarks/results/bench_<label>.csv. Compare two runs with
#   Rscript dev/benchmarks/compare.R results/bench_before.csv results/bench_after.csv
#
# Every row also stores a fingerprint (hash) of what the code returned, under a fixed
# random seed. If a later version is faster but the fingerprint changed, the output
# changed too: that is expected only when an algorithm was changed on purpose.
#
# Timings depend on the computer and on what else it is doing, so only compare runs made
# on the same machine, ideally with nothing else running.

args <- commandArgs(trailingOnly = TRUE)
arg <- function(name, default = NULL) {
  i <- match(paste0("--", name), args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- dirname(normalizePath(script_path))
root <- normalizePath(file.path(script_dir, "..", ".."))
setwd(root)

for (pkg in c("bench", "devtools", "rlang")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Please install the '", pkg, "' package first.", call. = FALSE)
}
suppressMessages(devtools::load_all(root, quiet = TRUE))
source(file.path(script_dir, "benchmark_helpers.R"))

suite <- arg("suite", "standard")
stopifnot(suite %in% c("quick", "standard", "full"))
label <- arg("label", format(Sys.Date(), "%Y%m%d"))
only <- arg("only", NULL)
iter_override <- arg("iterations", NULL)

sizes_in_suite <- switch(suite, quick = "small", standard = c("small", "medium"), full = c("small", "medium", "large"))
default_iterations <- c(small = if (suite == "quick") 3 else 5, medium = 3, large = 1)

run_git <- function(...) tryCatch(trimws(system2("git", c(...), stdout = TRUE, stderr = FALSE)), error = function(e) NA_character_, warning = function(w) NA_character_)
meta <- data.frame(
  label = label, suite = suite, timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  git_commit = run_git("rev-parse", "--short", "HEAD")[1],
  git_state = if (length(run_git("status", "--porcelain", "--untracked-files=no")) == 0) "clean" else "uncommitted-changes",
  package_version = as.character(utils::packageVersion("WormTransgeneBuilder")),
  r_version = paste(R.version$major, R.version$minor, sep = "."), platform = R.version$platform,
  cores = parallel::detectCores(), rnafold = Sys.which("RNAfold") != "",
  stringsAsFactors = FALSE
)
cat(sprintf("Benchmark '%s' (%s suite) at commit %s, %s, R %s, %d cores\n\n",
            label, suite, meta$git_commit, meta$git_state, meta$r_version, meta$cores))

# ---- start-up cost, in fresh R sessions ------------------------------------------------
# Two views. "dev": devtools::load_all() (what you use while developing). "installed": the
# package installed into a temporary library, which is how the app runs in Docker; there the
# GLO score table is only read the first time it is used (before the optimization of the GLO
# lookup this row timed loading the old 54 MB glo_full dictionary, which took about 9 s).
time_fresh_sessions <- function(code, repeats = 3, lib = NULL) {
  prefix <- if (is.null(lib)) "" else sprintf(".libPaths(c(%s, .libPaths()));", deparse(lib))
  t(vapply(seq_len(repeats), function(i) {
    out <- system2(file.path(R.home("bin"), "Rscript"), c("-e", shQuote(paste(prefix, code))), stdout = TRUE, stderr = FALSE)
    as.numeric(strsplit(trimws(out[length(out)]), " ")[[1]])
  }, numeric(2)))
}
startup_rows <- function(ids, res, repeats) {
  data.frame(id = ids, size = NA_character_, size_bp = NA_real_, iterations = repeats,
             median_s = apply(res, 2, median), min_s = apply(res, 2, min), mem_alloc_mb = NA_real_,
             n_gc = NA_real_, output_hash = NA_character_, status = "ok", stringsAsFactors = FALSE)
}
time_startup <- function(repeats = 3, installed = TRUE) {
  dev <- time_fresh_sessions(paste(
    't0 <- Sys.time(); suppressMessages(devtools::load_all(quiet = TRUE)); t1 <- Sys.time();',
    'n <- length(WormTransgeneBuilder:::glo_score_table()); t2 <- Sys.time();',
    'cat(as.numeric(difftime(t1, t0, units = "secs")), as.numeric(difftime(t2, t1, units = "secs")), "\n")'), repeats)
  out <- startup_rows("startup_dev_load_all", dev[, 1, drop = FALSE], repeats)
  if (installed) {
    lib <- tempfile("benchlib"); dir.create(lib)
    on.exit(unlink(lib, recursive = TRUE), add = TRUE)
    cat("  installing the package into a temporary library ...\n")
    ok <- tryCatch({
      suppressWarnings(utils::install.packages(root, repos = NULL, type = "source", lib = lib, quiet = TRUE,
                                               INSTALL_opts = c("--no-docs", "--no-multiarch")))
      requireNamespace("WormTransgeneBuilder", lib.loc = lib, quietly = TRUE) &&
        dir.exists(file.path(lib, "WormTransgeneBuilder"))
    }, error = function(e) FALSE)
    if (ok) {
      inst <- time_fresh_sessions(paste(
        't0 <- Sys.time(); suppressMessages(library(WormTransgeneBuilder)); t1 <- Sys.time();',
        'n <- length(WormTransgeneBuilder:::glo_score_table()); t2 <- Sys.time();',
        'cat(as.numeric(difftime(t1, t0, units = "secs")), as.numeric(difftime(t2, t1, units = "secs")), "\n")'), repeats, lib)
      out <- rbind(out, startup_rows(c("startup_installed_library_call", "startup_installed_first_use_of_GLO_dictionary"), inst, repeats))
    } else cat("  (could not install the package; skipping the installed-package start-up measurement)\n")
  }
  out
}

rows <- list()
if (is.null(only) || grepl(only, "startup")) {
  cat("Measuring start-up in fresh R sessions ...\n")
  rows[[length(rows) + 1]] <- time_startup(installed = suite != "quick")
  print(rows[[1]][, c("id", "median_s")], row.names = FALSE)
  cat("\n")
}

cases <- benchmark_cases()
if (!is.null(only)) cases <- Filter(function(cs) grepl(only, cs$id), cases)
genes <- lapply(benchmark_sizes[sizes_in_suite], function(n) list(cds = make_benchmark_cds(n), best = make_best_codon_benchmark_cds(n)))

jobs <- do.call(rbind, lapply(cases, function(cs) {
  sz <- intersect(sizes_in_suite, cs$sizes)
  if (length(sz) == 0) NULL else data.frame(id = cs$id, size = sz, stringsAsFactors = FALSE)
}))
if (is.null(jobs)) jobs <- data.frame(id = character(0), size = character(0))
started <- Sys.time()
for (j in seq_len(nrow(jobs))) {
  cs <- cases[[match(jobs$id[j], vapply(cases, `[[`, "", "id"))]]
  size <- jobs$size[j]
  iters <- if (is.null(iter_override)) default_iterations[[size]] else as.integer(iter_override)
  cat(sprintf("[%2d/%d] %-32s %-7s x%d ... ", j, nrow(jobs), cs$id, size, iters))
  x <- genes[[size]]
  res <- tryCatch({
    # bench::mark does not hand back what the code returned, so the output is kept here
    # (the seed is fixed, so every iteration returns the same thing).
    kept <- new.env()
    run_mark <- function(memory) suppressWarnings(bench::mark(
      { kept$output <- bench_with_seed(1L, cs$fn(x)); invisible(NULL) },
      iterations = iters, check = FALSE, memory = memory, filter_gc = FALSE
    ))
    # bench's memory profiler cannot follow the worker processes that the multi-core piRNA
    # scan starts and stops with "Memory profiling failed". In that case time it without
    # memory profiling and record the memory as unavailable.
    memory_ok <- TRUE
    m <- tryCatch(run_mark(TRUE), error = function(e) {
      if (!grepl("Memory profiling failed", conditionMessage(e))) stop(e)
      memory_ok <<- FALSE
      run_mark(FALSE)
    })
    data.frame(id = cs$id, size = size, size_bp = nchar(x$cds), iterations = iters,
               median_s = as.numeric(m$median), min_s = as.numeric(m$min),
               mem_alloc_mb = if (memory_ok) as.numeric(m$mem_alloc) / 1024^2 else NA_real_, n_gc = sum(m$n_gc),
               output_hash = fingerprint(kept$output), status = "ok", stringsAsFactors = FALSE)
  }, error = function(e) {
    data.frame(id = cs$id, size = size, size_bp = nchar(x$cds), iterations = iters, median_s = NA_real_, min_s = NA_real_,
               mem_alloc_mb = NA_real_, n_gc = NA_real_, output_hash = NA_character_,
               status = paste("error:", conditionMessage(e)), stringsAsFactors = FALSE)
  })
  cat(if (res$status == "ok") sprintf("median %s\n", format(bench::as_bench_time(res$median_s))) else paste0(res$status, "\n"))
  rows[[length(rows) + 1]] <- res
}

results <- cbind(meta, do.call(rbind, rows), row.names = NULL)
outdir <- arg("outdir", file.path(script_dir, "results"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
out <- file.path(outdir, paste0("bench_", label, ".csv"))
utils::write.csv(results, out, row.names = FALSE)
cat(sprintf("\nFinished in %.1f minutes. Results written to %s\n", as.numeric(difftime(Sys.time(), started, units = "mins")),
            sub(paste0(root, "/"), "", out, fixed = TRUE)))
if (any(results$status != "ok")) cat("Some cases failed:\n", paste(results$id[results$status != "ok"], collapse = ", "), "\n")
