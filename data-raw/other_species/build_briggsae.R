# C. briggsae (AF16) piRNA database from Pastore, Hertz & Tang 2022 (RNA Biol 19:1186-1197), Table S1: 25,883 piRNA loci
# defined by PRG-1 immunoprecipitation (GEO GSE203297). The 21-mer of every locus is checked against the CB4 genome
# (WormBase ParaSite WBPS19, PRJNA10731) before it is written. The table's `Start` is 0-based: a locus covers bases
# Start+1 .. End (1-based).
# Usage, from the project root:
#   Rscript data-raw/other_species/build_briggsae.R <briggsae genomic.fa (uncompressed)>
source("data-raw/other_species/common.R")
genome <- read_fasta(commandArgs(trailingOnly = TRUE)[1])
d <- read.csv("data-raw/inputs/pastore2022_TableS1_cbr_piRNA_loci.csv", stringsAsFactors = FALSE, check.names = FALSE)
names(d) <- c("Chrom", "start0", "End", "Name", "Full_21mer", "Strand", "IP_RPM")
stopifnot(all(d$Chrom %in% names(genome)), all(d$End - d$start0 == 21), all(nchar(d$Full_21mer) == 21))
d$Start <- d$start0 + 1L

# every sequence must equal the genome at its locus (read on the locus's strand)
fetch <- function(ch, s, e, sd) { x <- substr(genome[[ch]], s, e); if (sd == "-") revcomp(x) else x }
genomic <- mapply(fetch, d$Chrom, d$Start, d$End, d$Strand, USE.NAMES = FALSE)
stopifnot(identical(genomic, d$Full_21mer))
cat("all", nrow(d), "sequences equal the genome at their coordinates\n")

up <- upstream_of(genome, d$Chrom, d$Start, d$End, d$Strand)
out <- data.frame(
  Name = d$Name,
  Target_20mer = revcomp(substr(d$Full_21mer, 2, 21)),
  Full_21mer = d$Full_21mer,
  Abundance = sprintf("%s RPM (PRG-1 IP)", signif(d$IP_RPM, 3)),
  Enrichment = "No Data",
  Species = "Caenorhabditis briggsae", Strain = "AF16", Assembly = "CB4 (WormBase ParaSite WBPS19, PRJNA10731)",
  Confidence = "high",
  Evidence = "PRG-1 IP, Pastore et al. 2022 (GEO GSE203297), Table S1",
  Motif_Ruby = has_ruby_core(up),
  Chrom = d$Chrom, Start = d$Start, End = d$End, Strand = d$Strand, PRG1_IP_RPM = d$IP_RPM,
  stringsAsFactors = FALSE)
stopifnot(!anyNA(out[, c("Name", "Target_20mer", "Full_21mer")]))
cat("duplicated Name", sum(duplicated(out$Name)), "| Full_21mer", sum(duplicated(out$Full_21mer)), "| Target_20mer", sum(duplicated(out$Target_20mer)), "\n")
cat(sprintf("Ruby motif flag: %d of %d (%.1f%%); on chromosomes I-V,X: %.1f%%; on unplaced contigs: %.1f%%\n", sum(out$Motif_Ruby), nrow(out), 100 * mean(out$Motif_Ruby),
    100 * mean(out$Motif_Ruby[out$Chrom %in% c("I", "II", "III", "IV", "V", "X")]), 100 * mean(out$Motif_Ruby[!out$Chrom %in% c("I", "II", "III", "IV", "V", "X")])))
utils::write.csv(out, gzfile("data-raw/staging/other_species/piRNA_Database_briggsae_staging.csv.gz"), row.names = FALSE, na = "")
cat("wrote data-raw/staging/other_species/piRNA_Database_briggsae_staging.csv.gz\n")
