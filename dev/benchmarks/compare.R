#!/usr/bin/env Rscript
# Compare two benchmark runs:
#   Rscript dev/benchmarks/compare.R results/bench_before.csv results/bench_after.csv
#
# For every case that appears in both files it shows the median time in each, the
# speed-up (before / after, so 2 means twice as fast), the change in memory, and whether
# the output fingerprint changed. A changed fingerprint means the code now returns
# something different under the same seed: fine if you changed an algorithm on purpose,
# a problem if you only meant to make it faster.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript compare.R <before.csv> <after.csv>", call. = FALSE)
a <- utils::read.csv(args[1], stringsAsFactors = FALSE)
b <- utils::read.csv(args[2], stringsAsFactors = FALSE)

cat(sprintf("before: %s  (%s, commit %s, %s, R %s, %d cores)\n", a$label[1], a$timestamp[1], a$git_commit[1], a$git_state[1], a$r_version[1], a$cores[1]))
cat(sprintf("after:  %s  (%s, commit %s, %s, R %s, %d cores)\n", b$label[1], b$timestamp[1], b$git_commit[1], b$git_state[1], b$r_version[1], b$cores[1]))
if (a$platform[1] != b$platform[1] || a$cores[1] != b$cores[1] || a$r_version[1] != b$r_version[1]) {
  cat("WARNING: the two runs used a different platform, R version or number of cores; timings are not comparable.\n")
}
cat("\n")

key <- function(d) paste(d$id, ifelse(is.na(d$size), "", d$size), sep = "|")
a$key <- key(a); b$key <- key(b)
m <- merge(a[, c("key", "id", "size", "median_s", "mem_alloc_mb", "output_hash", "status")],
           b[, c("key", "median_s", "mem_alloc_mb", "output_hash", "status")],
           by = "key", suffixes = c("_before", "_after"))
m <- m[order(match(m$id, unique(a$id)), match(m$size, c("small", "medium", "large"))), ]

fmt_time <- function(s) ifelse(is.na(s), "-", ifelse(s < 1e-3, sprintf("%.0f us", s * 1e6), ifelse(s < 1, sprintf("%.1f ms", s * 1e3), sprintf("%.2f s", s))))
speedup <- m$median_s_before / m$median_s_after
out <- data.frame(
  case = ifelse(is.na(m$size), m$id, paste0(m$id, " [", m$size, "]")),
  before = fmt_time(m$median_s_before), after = fmt_time(m$median_s_after),
  speedup = ifelse(is.na(speedup), "-", sprintf("%.2fx", speedup)),
  memory = ifelse(is.na(m$mem_alloc_mb_before) | is.na(m$mem_alloc_mb_after), "-",
                  sprintf("%.1f -> %.1f MB", m$mem_alloc_mb_before, m$mem_alloc_mb_after)),
  output = ifelse(is.na(m$output_hash_before) | is.na(m$output_hash_after), "-",
                  ifelse(m$output_hash_before == m$output_hash_after, "same", "CHANGED")),
  stringsAsFactors = FALSE)
print(out, row.names = FALSE, right = FALSE)

only_a <- setdiff(a$key, b$key); only_b <- setdiff(b$key, a$key)
if (length(only_a)) cat("\nOnly in 'before':", paste(only_a, collapse = ", "), "\n")
if (length(only_b)) cat("\nOnly in 'after':", paste(only_b, collapse = ", "), "\n")
if (any(out$output == "CHANGED")) cat("\nNOTE: output CHANGED for", sum(out$output == "CHANGED"), "case(s). Check that this is intended.\n")
