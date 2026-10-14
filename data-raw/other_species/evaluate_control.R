# Step 4 of the brief: how well does the calling method rediscover the known C. elegans piRNAs (the rebuilt database)?
# Usage: Rscript data-raw/other_species/evaluate_control.R <prefix of the control library>
args <- commandArgs(trailingOnly = TRUE)
x <- readRDS(paste0(args[1], ".candidates.rds"))
cand <- x$cand[!x$cand$miRNA, ]
db <- read.csv("data-raw/staging/piRNA_Database_staging.csv", na.strings = "", stringsAsFactors = FALSE)
known_t1 <- db$Full_21mer[db$Type == "type 1"]
known_all <- db$Full_21mer
known_t1_t <- known_t1[substr(known_t1, 1, 1) == "T"]
say <- function(...) cat(sprintf(...), "\n")
say("known type 1 piRNAs: %d (%d start with T); whole database: %d", length(known_t1), length(known_t1_t), length(known_all))
say("known type 1 with T-start seen at all in the library (>=1 read, exact match): %d", sum(known_t1_t %in% x$cand$Full_21mer))
say("candidates (T-start, 21 nt, exact, not miRNA): %d", nrow(cand))
rows <- list()
for (rule in c("reads >= 1", "reads >= 2", "RPM >= 1", "RPM >= 5")) {
  keep <- switch(rule, "reads >= 1" = cand$reads >= 1, "reads >= 2" = cand$reads >= 2, "RPM >= 1" = cand$RPM >= 1, "RPM >= 5" = cand$RPM >= 5)
  for (tier in c("all calls", "with Ruby motif (tier A)", "without motif (tier B)")) {
    sel <- keep & switch(tier, "all calls" = TRUE, "with Ruby motif (tier A)" = cand$motif, "without motif (tier B)" = !cand$motif)
    calls <- cand$Full_21mer[sel]
    rows[[length(rows) + 1]] <- data.frame(rule = rule, tier = tier, calls = length(calls),
      in_type1 = sum(calls %in% known_t1), in_database = sum(calls %in% known_all),
      precision_type1 = round(100 * mean(calls %in% known_t1), 1), precision_database = round(100 * mean(calls %in% known_all), 1),
      recall_type1_T = round(100 * sum(known_t1_t %in% calls) / length(known_t1_t), 1))
  }
}
print(do.call(rbind, rows), row.names = FALSE)
# for the reads >= 2 / tier A call, how do false calls look?
a <- cand[cand$reads >= 2 & cand$motif, ]
say("\ntier A, >= 2 reads: %d calls; in database %d; not in database %d. Of the not-in-database calls: median RPM %.2f, multi-locus %.1f%%",
    nrow(a), sum(a$Full_21mer %in% known_all), sum(!a$Full_21mer %in% known_all), median(a$RPM[!a$Full_21mer %in% known_all]), 100 * mean(a$n_loci[!a$Full_21mer %in% known_all] > 1))
