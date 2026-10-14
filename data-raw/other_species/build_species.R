# Steps 6-8 of data-raw/pirna_other_species_brief.md for C. remanei (PB4641) and C. brenneri (PB2801): combines the candidate
# tables of the three libraries of a species, assigns tiers and writes the staging database.
#   Tier A: Ruby upstream motif (GTTTC core at -50..-35) AND at least 1 RPM in some library  -> the only tier kept
#   Tier B (no motif, at least 1 RPM) was dropped by the owner's decision (2026-10-06): only 77% precision in the C. elegans control.
# The 1 RPM cut was chosen by the owner after the C. elegans control (evaluate_control.R): at >= 1 RPM the calls had 88% precision
# (motif 94%, no motif 77%) against 34% for >= 2 reads, at a recall of about 21% of the known type 1 piRNAs. Everything below
# 1 RPM is left out. All entries are low confidence: one library per condition,
# no replicates, no PRG-1 IP, and the assemblies are scaffold-level.
# Usage: Rscript data-raw/other_species/build_species.R remanei|brenneri
source("data-raw/other_species/common.R")
sp <- commandArgs(trailingOnly = TRUE)[1]
W <- "/Volumes/CFJ_backup1S1TB/pirna_other_species/collapsed/"
cfg <- list(
  remanei = list(species = "Caenorhabditis remanei", strain = "PB4641", code = "Cre", assembly = "PRJNA53967 (ASM14951v1; WormBase ParaSite WBPS19)",
                 libs = c(adult = "SRR585541", embryo = "SRR585542", male = "SRR585543")),
  brenneri = list(species = "Caenorhabditis brenneri", strain = "PB2801", code = "Cbr", assembly = "PRJNA20035 (C_brenneri-6.0.1b; WormBase ParaSite WBPS19)",
                  libs = c(adult = "SRR585544", embryo = "SRR585545", male = "SRR585546")))[[sp]]
cfg$code <- c(remanei = "Cre", brenneri = "Cbn")[[sp]]
parts <- lapply(names(cfg$libs), function(k) {
  x <- readRDS(paste0(W, cfg$libs[[k]], ".candidates.rds"))
  cat(sprintf("%s (%s): reads %d, mapped %d, candidates (not miRNA) %d\n", k, cfg$libs[[k]], x$total_reads, x$total_mapped, sum(!x$cand$miRNA)))
  d <- x$cand[!x$cand$miRNA, c("Full_21mer", "reads", "RPM", "n_loci", "motif", "Chrom", "Start", "End", "Strand")]
  names(d)[2:3] <- paste0(c("reads_", "RPM_"), k); d
})
all_seqs <- unique(unlist(lapply(parts, `[[`, "Full_21mer")))
m <- data.frame(Full_21mer = all_seqs, stringsAsFactors = FALSE)
for (i in seq_along(parts)) {
  p <- parts[[i]]; k <- names(cfg$libs)[i]; j <- match(m$Full_21mer, p$Full_21mer)
  m[[paste0("reads_", k)]] <- ifelse(is.na(j), 0L, p[[paste0("reads_", k)]][j])
  m[[paste0("RPM_", k)]] <- ifelse(is.na(j), 0, p[[paste0("RPM_", k)]][j])
  for (col in c("n_loci", "motif", "Chrom", "Start", "End", "Strand")) {
    if (is.null(m[[col]])) m[[col]] <- p[[col]][j] else m[[col]] <- ifelse(is.na(m[[col]]), p[[col]][j], m[[col]])
  }
}
m$max_reads <- pmax(m$reads_adult, m$reads_embryo, m$reads_male)
m$max_RPM <- pmax(m$RPM_adult, m$RPM_embryo, m$RPM_male)
m$Tier <- ifelse(m$max_RPM >= 1 & m$motif, "A", NA)
cat(sprintf("\nunique candidate sequences over the three libraries: %d\n", nrow(m)))
print(table(Tier = ifelse(is.na(m$Tier), "left out", m$Tier), useNA = "no"))
kept <- m[!is.na(m$Tier), ]
kept <- kept[order(kept$Tier, -kept$max_RPM), ]
kept$Chrom <- as.character(kept$Chrom)
sd <- ifelse(kept$Strand == "+", "plus", "minus")
out <- data.frame(
  Name = sprintf("%s_21U_%s_%d_%s", cfg$code, kept$Chrom, kept$Start, sd),
  Target_20mer = revcomp(substr(kept$Full_21mer, 2, 21)),
  Full_21mer = kept$Full_21mer,
  Abundance = sprintf("%s RPM (small-RNA-seq, best library)", signif(kept$max_RPM, 3)),
  Enrichment = "No Data",
  Species = cfg$species, Strain = cfg$strain, Assembly = cfg$assembly,
  Confidence = "low", Tier = kept$Tier, Motif_Ruby = kept$motif,
  Evidence = "Shi et al. 2013 (GEO GSE41461) reads, called here: 21 nt, 5' T, exact match, not a miRNA, >= 1 RPM",
  Chrom = kept$Chrom, Start = kept$Start, End = kept$End, Strand = kept$Strand, n_genomic_hits = kept$n_loci,
  RPM_adult = kept$RPM_adult, RPM_embryo = kept$RPM_embryo, RPM_male = kept$RPM_male,
  reads_adult = kept$reads_adult, reads_embryo = kept$reads_embryo, reads_male = kept$reads_male,
  stringsAsFactors = FALSE)
cat(sprintf("duplicated Name %d, Full_21mer %d, Target_20mer %d\n", sum(duplicated(out$Name)), sum(duplicated(out$Full_21mer)), sum(duplicated(out$Target_20mer))))
# several sequences can start at the same position on the same strand only if they are identical, so names are unique; if a
# Target_20mer is shared (two sequences differing only in the first base), keep the one with more reads and report it.
dup <- duplicated(out$Target_20mer) | duplicated(out$Target_20mer, fromLast = TRUE)
if (any(dup)) {
  cat(sprintf("%d entries share a Target_20mer with another entry (they differ only in their first base); keeping the most abundant of each\n", sum(dup)))
  out <- out[order(-(out$reads_adult + out$reads_embryo + out$reads_male)), ]; out <- out[!duplicated(out$Target_20mer), ]
  out <- out[order(out$Tier, -(out$RPM_adult + out$RPM_embryo + out$RPM_male)), ]
}
cat(sprintf("final: %d entries (tier A %d, tier B %d); multi-locus %d; with an adult-library read: %d\n", nrow(out), sum(out$Tier == "A"), sum(out$Tier == "B"), sum(out$n_genomic_hits > 1), sum(out$reads_adult > 0)))
cat(sprintf("compare with Shi et al. 2013 (observed piRNAs): remanei 23.7k, brenneri 29.8k\n"))
utils::write.csv(out, gzfile(sprintf("data-raw/staging/other_species/piRNA_Database_%s_staging.csv.gz", sp)), row.names = FALSE, na = "")   # gzip keeps the repository small
