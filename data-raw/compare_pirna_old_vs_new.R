# Phase 2/3 of the piRNA rebuild: compare the piRNAs extracted from WormBase (WS298) with the earlier,
# cleaned database. Usage, from the project root:
#   Rscript data-raw/compare_pirna_old_vs_new.R path/to/piRNA_sequence.csv   (old file: Name, Target_20mer, Full_21mer, no header)
old_file <- commandArgs(trailingOnly = TRUE)[1]
new <- utils::read.csv("data-raw/staging/pirna_extracted.csv", stringsAsFactors = FALSE)
old <- utils::read.csv(old_file, header = FALSE, col.names = c("Name", "Target_20mer", "Full_21mer"), stringsAsFactors = FALSE)
revcomp <- function(x) vapply(strsplit(x, ""), function(b) paste(rev(chartr("ACGT", "TGCA", b)), collapse = ""), "")
say <- function(...) cat(sprintf(...), "\n")

old$type2 <- grepl(" \\(type 2\\)$", old$Name)
old$base_name <- sub(" \\(type 2\\)$", "", old$Name)
say("OLD: %d rows (%d type 1, %d type 2); NEW: %d", nrow(old), sum(!old$type2), sum(old$type2), nrow(new))
say("old internal check Target == revcomp(Full[2..21]): %d of %d", sum(old$Target_20mer == revcomp(substr(old$Full_21mer, 2, 21))), nrow(old))
say("NEW internal check Target == revcomp(Full[2..21]): %d of %d", sum(new$Target_20mer == revcomp(substr(new$Full_21mer, 2, 21))), nrow(new))

cat("\n== Names ==\n")
t1 <- old[!old$type2, ]
say("old type 1 names unique: %s", !anyDuplicated(t1$Name))
say("old type 1 in new (by name): %d; old type 1 NOT in new: %d; new NOT in old type 1: %d",
    sum(t1$Name %in% new$Name), sum(!t1$Name %in% new$Name), sum(!new$Name %in% t1$Name))
if (any(!t1$Name %in% new$Name)) { cat("old type 1 not in new:\n"); print(t1[!t1$Name %in% new$Name, c("Name", "Full_21mer")]) }
if (any(!new$Name %in% t1$Name)) { cat("new not in old type 1:\n"); print(new[!new$Name %in% t1$Name, c("Name", "WBGene", "chrom", "start", "strand", "Full_21mer")]) }

cat("\n== Sequences of the names in both ==\n")
m <- merge(t1, new, by = "Name", suffixes = c("_old", "_new"))
say("same Full_21mer: %d of %d; same Target_20mer: %d of %d", sum(m$Full_21mer_old == m$Full_21mer_new), nrow(m),
    sum(m$Target_20mer_old == m$Target_20mer_new), nrow(m))
d <- m[m$Full_21mer_old != m$Full_21mer_new, ]
if (nrow(d)) { say("differing rows: %d", nrow(d)); print(utils::head(d[, c("Name", "Full_21mer_old", "Full_21mer_new", "chrom", "start", "strand")], 15))
  say("of these, differ by one extra/missing base shift?  old[2..21]==new[1..20]: %d; old[1..20]==new[2..21]: %d",
      sum(substr(d$Full_21mer_old, 2, 21) == substr(d$Full_21mer_new, 1, 20)), sum(substr(d$Full_21mer_old, 1, 20) == substr(d$Full_21mer_new, 2, 21)))
  say("mismatching bases per differing row (Hamming): %s", paste(names(table(mapply(function(a, b) sum(strsplit(a, "")[[1]] != strsplit(b, "")[[1]]), d$Full_21mer_old, d$Full_21mer_new))), table(mapply(function(a, b) sum(strsplit(a, "")[[1]] != strsplit(b, "")[[1]]), d$Full_21mer_old, d$Full_21mer_new)), sep = ":", collapse = "  ")) }

cat("\n== 5' base ==\n")
say("new: %d of %d start with T (%.2f%%); non-T: %d", sum(substr(new$Full_21mer, 1, 1) == "T"), nrow(new), 100 * mean(substr(new$Full_21mer, 1, 1) == "T"), sum(substr(new$Full_21mer, 1, 1) != "T"))
say("first-base table (new): %s", paste(names(table(substr(new$Full_21mer, 1, 1))), table(substr(new$Full_21mer, 1, 1)), sep = "=", collapse = " "))
old_nonT <- t1$Name[substr(t1$Full_21mer, 1, 1) != "T"]; new_nonT <- new$Name[substr(new$Full_21mer, 1, 1) != "T"]
say("old type 1 non-T: %d; new non-T: %d; in both: %d", length(old_nonT), length(new_nonT), length(intersect(old_nonT, new_nonT)))
say("old type 2 non-T: %d of %d", sum(substr(old$Full_21mer[old$type2], 1, 1) != "T"), sum(old$type2))

cat("\n== The old 'type 2' entries ==\n")
t2 <- old[old$type2, ]
say("type 2: %d rows; their base names are in new: %d; their Full_21mer equals some new Full_21mer: %d; equals the type-1 sibling's: %d",
    nrow(t2), sum(t2$base_name %in% new$Name), sum(t2$Full_21mer %in% new$Full_21mer),
    sum(t2$Full_21mer == t1$Full_21mer[match(t2$base_name, t1$Name)], na.rm = TRUE))
sib <- t1[match(t2$base_name, t1$Name), ]
say("type 2 Target_20mer equals its type-1 sibling's Target: %d", sum(t2$Target_20mer == sib$Target_20mer, na.rm = TRUE))
say("type 2 Full_21mer reverse-complement of sibling's Full_21mer: %d", sum(t2$Full_21mer == revcomp(sib$Full_21mer), na.rm = TRUE))
sh <- function(a, b, k) substr(a, 1 + k, 21) == substr(b, 1, 21 - k)
say("type 2 overlaps its sibling shifted by 1..5 bases (either direction): %s",
    paste(sapply(1:5, function(k) sum(sh(t2$Full_21mer, sib$Full_21mer, k) | sh(sib$Full_21mer, t2$Full_21mer, k), na.rm = TRUE)), collapse = " "))
print(utils::head(data.frame(name = t2$Name, type2 = t2$Full_21mer, sibling = sib$Full_21mer, chrom = new$chrom[match(t2$base_name, new$Name)]), 6))

cat("\n== Duplicates ==\n")
say("new: duplicated Full_21mer %d; duplicated Target_20mer %d; duplicated names %d", sum(duplicated(new$Full_21mer)), sum(duplicated(new$Target_20mer)), sum(duplicated(new$Name)))
say("old: duplicated Target_20mer %d", sum(duplicated(old$Target_20mer)))
say("old Target_20mer present among new Target_20mer: %d of %d (all old rows); type 1: %d of %d",
    sum(old$Target_20mer %in% new$Target_20mer), nrow(old), sum(t1$Target_20mer %in% new$Target_20mer), nrow(t1))
