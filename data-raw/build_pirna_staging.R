# Phase 4 and 6 of the piRNA rebuild: build the staging database. Usage, from the project root (after
# extract_pirna_from_gff3.R and validate_new21U_and_type2.R):
#   Rscript data-raw/build_pirna_staging.R path/to/piRNA_abundance.csv path/to/piRNA_sequence.csv
# (the first file is the Billi et al. 2013 male/female table: Name, Abundance, Enrichment, no header; the second the
# earlier database, used only for the old "type 2" entries.)
#
# What goes in (decisions of the owner, 2026-10-05 and 2026-10-06):
#  * Type 1: all 15,363 piRNAs annotated in WormBase WS298 (the authoritative list);
#  * Type 2: all 2,485 type 2 entries of the earlier database. This is the type II set of pirScan (Wu et al. 2018,
#    NAR 46:W43), "collected from PRG-1 IP data with at least 1 read per million reads" in Gu et al. 2012 and Tang et al.
#    2016 (Cell 151:1488, Cell 164:974); the app offers them as an option ("Type" column);
#  * New: the 466 new 21U loci that Claycomb lab (Seroussi et al. 2023) found enriched in PRG-1 immunoprecipitation, each
#    with an internal label (new21U_<chrom>_<leftmost base>_<plus|minus>), NOT a 21ur- gene name. 30 of them are
#    identical to type 2 entries and stay as type 2 entries (with their Claycomb location filled in);
#  * Abundance = PRG-1 IP abundance (RPM, average of the two replicates) from Claycomb; "not detected in PRG-1 IP" when a
#    type 1 piRNA is not in their table or has no reads in the IP; type 2 entries not in their table say so instead;
#  * Enrichment = the male/female enrichment of Billi et al. 2013, "No Data" where Billi has no entry or "n.d.".
args <- commandArgs(trailingOnly = TRUE)
source("data-raw/fasta_helpers.R")   # for revcomp()
say <- function(...) cat(sprintf(...), "\n")
read <- function(f, ...) utils::read.csv(f, stringsAsFactors = FALSE, na.strings = c("NA", ""), ...)

new <- read("data-raw/staging/pirna_extracted.csv")
s3 <- read("data-raw/inputs/claycomb_supp3_PRG1_piRNA_IP.csv")
s4 <- read("data-raw/inputs/claycomb_supp4_prg1_mutant_depleted_piRNA.csv")
c466 <- read("data-raw/staging/claycomb_466_with_sequences.csv")
hits <- read("data-raw/staging/type2_genome_hits.csv")
billi <- utils::read.csv(args[1], header = FALSE, col.names = c("Name", "Abundance_Billi", "Enrichment"), stringsAsFactors = FALSE)
old <- utils::read.csv(args[2], header = FALSE, col.names = c("Name", "Target_20mer", "Full_21mer"), stringsAsFactors = FALSE)

# ---- Claycomb PRG-1 IP (Supp 3) and prg-1 mutant depletion (Supp 4), by WBGene ---------------------------------
abundance_text <- function(rpm) ifelse(is.na(rpm) | rpm <= 0.01, "not detected in PRG-1 IP", sprintf("%s RPM (PRG-1 IP)", signif(rpm, 3)))
m3 <- match(new$WBGene, s3$WBGene)
new$PRG1_IP_avg_RPM <- s3$avg_IP[m3]
new$PRG1_IP_avg_enrichment <- s3$avg_enrichment[m3]
new$PRG1_IP_enriched <- as.logical(s3$is_target[m3])
# Supp 4 has two rows for 534 genes with different values; keep the most significant (smallest padj) and say how many rows there were.
s4 <- s4[order(is.na(s4$padj), s4$padj, -abs(s4$log2FoldChange)), ]
m4 <- match(new$WBGene, s4$WBGene)
new$PRG1_depleted_log2FC <- s4$log2FoldChange[m4]
new$PRG1_depleted_padj <- s4$padj[m4]
say("WS298 piRNAs: in Supp 3: %d; detected in the IP (avg RPM > 0.01): %d; PRG-1-enriched (is_target): %d; in Supp 4: %d",
    sum(!is.na(m3)), sum(new$PRG1_IP_avg_RPM > 0.01, na.rm = TRUE), sum(new$PRG1_IP_enriched, na.rm = TRUE), sum(!is.na(m4)))

# ---- Billi male/female enrichment, by name ---------------------------------------------------------------------
new$Enrichment <- billi$Enrichment[match(new$Name, billi$Name)]
say("Billi entry found for %d WS298 piRNAs; 'No Data' for %d. Billi names dropped (not in WS298): %d", sum(!is.na(new$Enrichment)),
    sum(is.na(new$Enrichment)), sum(!billi$Name %in% new$Name))
say("Billi lists \"n.d.\" (no enrichment determined) for %d of those; shown as \"No Data\" too", sum(new$Enrichment %in% "n.d."))
new$Enrichment[is.na(new$Enrichment) | new$Enrichment == "n.d."] <- "No Data"
new$Source <- "WormBase WS298"

# ---- type 2 entries (all 2,485) and the 466 new PRG-1 IP loci ---------------------------------------------------------
t2 <- old[grepl(" \\(type 2\\)$", old$Name), ]
stopifnot(all(t2$Target_20mer == revcomp(substr(t2$Full_21mer, 2, 21))))
n_hits <- table(hits$Name)[t2$Name]
single <- hits[hits$Name %in% names(n_hits)[n_hits == 1], ]
pos_of <- function(col) single[[col]][match(t2$Name, single$Name)]
supported <- c466[c466$identical_21mer_type2, ]
t2$Claycomb_location <- supported$location[match(t2$Full_21mer, supported$Full_21mer)]
t2$PRG1_IP_avg_RPM <- supported$avg_IP[match(t2$Full_21mer, supported$Full_21mer)]
t2$PRG1_IP_avg_enrichment <- supported$avg_enrichment[match(t2$Full_21mer, supported$Full_21mer)]
t2$Claycomb_class <- ifelse(is.na(t2$Claycomb_location), NA, "same sequence as an old type 2 entry")
type2 <- data.frame(Name = t2$Name, Type = "type 2", WBGene = NA_character_, chrom = pos_of("chrom"), start = as.integer(pos_of("start")),
                    end = as.integer(pos_of("end")), strand = pos_of("strand"), n_genomic_hits = as.integer(n_hits),
                    Full_21mer = t2$Full_21mer, Target_20mer = t2$Target_20mer, PRG1_IP_avg_RPM = t2$PRG1_IP_avg_RPM,
                    PRG1_IP_avg_enrichment = t2$PRG1_IP_avg_enrichment, PRG1_IP_enriched = ifelse(is.na(t2$Claycomb_location), NA, TRUE),
                    PRG1_depleted_log2FC = NA_real_, PRG1_depleted_padj = NA_real_, Enrichment = "No Data",
                    Source = "Mello lab type 2 (Gu et al. 2012); PRG-1 IP >= 1 RPM (Wu et al. 2018)",
                    Claycomb_location = t2$Claycomb_location, Claycomb_class = t2$Claycomb_class, stringsAsFactors = FALSE)
# the 30 with a Claycomb locus: use that locus, which is unambiguous
sup <- !is.na(type2$Claycomb_location)
pp <- do.call(rbind, strsplit(type2$Claycomb_location[sup], ":"))
type2$chrom[sup] <- pp[, 1]; type2$start[sup] <- as.integer(pp[, 2]); type2$end[sup] <- as.integer(pp[, 2]) + 20L; type2$strand[sup] <- pp[, 3]
say("type 2 entries: %d (%d with one genomic hit, %d with several, %d with a Claycomb PRG-1 IP locus)", nrow(type2), sum(n_hits == 1), sum(n_hits > 1), sum(sup))

c466$class <- ifelse(c466$identical_21mer_type2, "same sequence as an old type 2 entry",
              ifelse(c466$overlaps_annotated, "overlaps (shifted from) an annotated piRNA", "independent of annotated piRNAs and of type 2"))
nw <- c466[!c466$identical_21mer_type2, ]
pl <- do.call(rbind, strsplit(nw$location, ":"))
newloci <- data.frame(Name = sprintf("new21U_%s_%s_%s", pl[, 1], pl[, 2], ifelse(pl[, 3] == "+", "plus", "minus")), Type = "PRG-1 IP new",
                      WBGene = NA_character_, chrom = pl[, 1], start = as.integer(pl[, 2]), end = as.integer(pl[, 2]) + 20L, strand = pl[, 3],
                      n_genomic_hits = NA_integer_, Full_21mer = nw$Full_21mer, Target_20mer = revcomp(substr(nw$Full_21mer, 2, 21)),
                      PRG1_IP_avg_RPM = nw$avg_IP, PRG1_IP_avg_enrichment = nw$avg_enrichment, PRG1_IP_enriched = TRUE,
                      PRG1_depleted_log2FC = NA_real_, PRG1_depleted_padj = NA_real_, Enrichment = "No Data",
                      Source = "Claycomb lab, Seroussi et al. 2023, Supp 6: new 21U enriched in PRG-1 IP; internal label, no gene name",
                      Claycomb_location = nw$location, Claycomb_class = nw$class, stringsAsFactors = FALSE)
# Some sequences sit at several loci (multi-copy 21U genes on chrIV). The app screens sequences, so a sequence is one entry:
# the label of its first locus, all loci in Claycomb_location, no single coordinates, the largest IP value of its loci.
say("new loci: %d of the 466 (the other %d are type 2 entries); distinct sequences: %d", nrow(newloci), nrow(c466) - nrow(newloci),
    length(unique(newloci$Full_21mer)))
newloci <- newloci[order(newloci$chrom, newloci$start), ]
by_seq <- split(seq_len(nrow(newloci)), factor(newloci$Full_21mer, levels = unique(newloci$Full_21mer)))
collapse <- function(i) {
  x <- newloci[i[1], ]
  x$n_genomic_hits <- length(i)
  x$Claycomb_location <- paste(newloci$Claycomb_location[i], collapse = ";")
  if (length(i) > 1) {
    x$chrom <- NA_character_; x$start <- NA_integer_; x$end <- NA_integer_; x$strand <- NA_character_
    x$PRG1_IP_avg_RPM <- max(newloci$PRG1_IP_avg_RPM[i]); x$PRG1_IP_avg_enrichment <- max(newloci$PRG1_IP_avg_enrichment[i])
    x$Claycomb_class <- paste(unique(newloci$Claycomb_class[i]), collapse = "; ")
  }
  x
}
newloci <- do.call(rbind, lapply(by_seq, collapse))
rownames(newloci) <- NULL
say("new PRG-1 IP entries: %d", nrow(newloci))

new$Type <- "type 1"; new$n_genomic_hits <- NA_integer_; new$Claycomb_location <- NA_character_; new$Claycomb_class <- NA_character_
cols <- names(type2)
ord <- function(d) d[order(as.integer(sub("^21ur-([0-9]+).*", "\\1", d$Name))), ]
db <- rbind(ord(new)[, cols], ord(type2), newloci[order(newloci$Name), ])
abund <- abundance_text(db$PRG1_IP_avg_RPM)
abund[db$Type == "type 2" & is.na(db$PRG1_IP_avg_RPM)] <- "not in Claycomb PRG-1 IP data"
db$Abundance <- abund
db$WS_version <- "WS298"
db <- db[, c("Name", "Target_20mer", "Full_21mer", "Abundance", "Enrichment", "Type", "Source", "WS_version", "WBGene", "chrom", "start", "end",
             "strand", "n_genomic_hits", "PRG1_IP_avg_RPM", "PRG1_IP_avg_enrichment", "PRG1_IP_enriched", "PRG1_depleted_log2FC",
             "PRG1_depleted_padj", "Claycomb_location", "Claycomb_class")]
names(db)[names(db) == "chrom"] <- "Chrom"; names(db)[names(db) == "start"] <- "Start"; names(db)[names(db) == "end"] <- "End"; names(db)[names(db) == "strand"] <- "Strand"
utils::write.csv(db, "data-raw/staging/piRNA_Database_staging.csv", row.names = FALSE, na = "")
say("\nWrote %d rows: %d type 1, %d type 2, %d new PRG-1 IP loci", nrow(db), sum(db$Type == "type 1"), sum(db$Type == "type 2"), sum(db$Type == "PRG-1 IP new"))
say("duplicated Full_21mer %d, Target_20mer %d, Name %d; first base T: %.2f%%", sum(duplicated(db$Full_21mer)), sum(duplicated(db$Target_20mer)),
    sum(duplicated(db$Name)), 100 * mean(substr(db$Full_21mer, 1, 1) == "T"))
print(table(db$Type, ifelse(db$Abundance == "not detected in PRG-1 IP", "not detected in PRG-1 IP", ifelse(db$Abundance == "not in Claycomb PRG-1 IP data", "not in Claycomb data", "RPM value"))))
utils::write.csv(c466, "data-raw/staging/claycomb_466_classified.csv", row.names = FALSE, quote = FALSE)
cat("\nClaycomb Supp 6 loci by class:\n"); print(table(c466$class))

# ---- how the old type 2 entries relate to the annotation (all 2,485) ------------------------------------------------
ann <- new[, c("chrom", "start", "end")]
hits$overlaps_annotated <- vapply(seq_len(nrow(hits)), function(i) any(ann$chrom == hits$chrom[i] & ann$start <= hits$end[i] & ann$end >= hits$start[i]), NA)
say("\nold type 2: %d sequences, %d genomic hits; sequences with more than one hit: %d; sequences with a hit that overlaps an annotated WS298 piRNA: %d",
    length(unique(hits$Name)), nrow(hits), sum(table(hits$Name) > 1), length(unique(hits$Name[hits$overlaps_annotated])))
