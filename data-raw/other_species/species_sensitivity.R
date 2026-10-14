# How strict should the tiers be? For remanei and brenneri: (1) how often does the Ruby-motif rule fire at random genome positions
# (C. elegans: 2.7%), and (2) how many tier A / B entries survive stricter thresholds. Usage: Rscript data-raw/other_species/species_sensitivity.R
source("data-raw/other_species/common.R")
S <- "/Volumes/CFJ_backup1S1TB/pirna_other_species/"
for (sp in c("remanei", "brenneri")) {
  g <- read_fasta(paste0(S, "genomes/", sp, ".fa"))
  g <- g[nchar(g) > 5000]
  set.seed(2)
  n <- 20000
  ch <- sample(names(g), n, TRUE, prob = nchar(g))
  st <- vapply(ch, function(c) sample(100:(nchar(g[[c]]) - 200), 1), 1L)
  sd <- sample(c("+", "-"), n, TRUE)
  up <- upstream_of(g, ch, st, st + 20L, sd)
  cat(sprintf("%s: random positions with the motif: %.1f%%  (genome GC %.1f%%)\n", sp, 100 * mean(has_ruby_core(up), na.rm = TRUE),
      100 * mean(unlist(strsplit(paste(substr(g, 1, 2e5), collapse = ""), "")) %in% c("G", "C"))))
  d <- read.csv(sprintf("data-raw/staging/other_species/piRNA_Database_%s_staging.csv.gz", sp), stringsAsFactors = FALSE)
  d$adult_or_best <- pmax(d$RPM_adult, d$RPM_embryo, d$RPM_male)
  cat("  entries by tier and rule (best library):\n")
  for (rule in c("current", "RPM >= 1", "RPM >= 5", ">= 5 reads in the adult library", "adult library only, >= 2 reads")) {
    keep <- switch(rule, "current" = TRUE, "RPM >= 1" = d$adult_or_best >= 1, "RPM >= 5" = d$adult_or_best >= 5,
                   ">= 5 reads in the adult library" = d$reads_adult >= 5, "adult library only, >= 2 reads" = d$reads_adult >= 2)
    cat(sprintf("   %-34s A %6d   B %6d   multi-locus %5d\n", rule, sum(keep & d$Tier == "A"), sum(keep & d$Tier == "B"), sum(keep & d$n_genomic_hits > 1)))
  }
  cat(sprintf("  tier A with >= 2 reads in at least two libraries: %d of %d\n", sum(d$Tier == "A" & ((d$reads_adult >= 2) + (d$reads_embryo >= 2) + (d$reads_male >= 2)) >= 2), sum(d$Tier == "A")))
}
