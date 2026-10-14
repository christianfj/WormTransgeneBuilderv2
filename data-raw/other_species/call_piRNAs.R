# Steps 4-6 of data-raw/pirna_other_species_brief.md for ONE library: turns the bowtie output of map_library.sh into a table of
# piRNA candidates (21 nt, 5' T, exact match to the genome, not a miRNA), with reads per million mapped reads (RPM; a read
# that maps to several loci is divided by the number of loci, as in Shi et al. 2013) and the Ruby upstream motif flag.
# Usage, from the project root:
#   Rscript data-raw/other_species/call_piRNAs.R <prefix> <genome.fa> <miRBase species codes, comma separated> <Shi S3 species code or "none">
# e.g.   ... /path/SRR585535 /path/c_elegans.fa cel Cel
source("data-raw/other_species/common.R")
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]; genome_file <- args[2]; mir_codes <- strsplit(args[3], ",")[[1]]; s3_code <- args[4]
W <- "/Volumes/CFJ_backup1S1TB/pirna_other_species"

collapsed <- read.delim(paste0(prefix, ".collapsed.tsv"), header = FALSE, col.names = c("seq", "reads"), stringsAsFactors = FALSE)
mapped <- read.delim(paste0(prefix, ".mapped_any.tsv"), header = FALSE, col.names = c("seq", "reads"), stringsAsFactors = FALSE)
total_reads <- sum(collapsed$reads)
total_mapped <- sum(mapped$reads)
cat(sprintf("reads kept after trimming: %d; mapped with 0 mismatches: %d (%.1f%%)\n", total_reads, total_mapped, 100 * total_mapped / total_reads))

# ---- miRNA exclusion: any 19-base stretch shared with a mature miRNA of the species -----------------------------------------
mature <- readLines(file.path(W, "tables/mirbase_mature.fa"))
hdr <- grep("^>", mature)
codes <- sub("^>([a-z]+)-.*", "\\1", mature[hdr])
mir <- toupper(gsub("U", "T", mature[hdr + 1]))[codes %in% mir_codes]
if (s3_code != "none") {   # new miRNAs identified by Shi et al. 2013 (Supplementary Table S3)
  s3 <- readLines(file.path(W, "tables/shi_S3_mature.txt"))
  mir <- c(mir, toupper(gsub("U", "T", s3[startsWith(s3, paste0(s3_code, "\t"))] |> sub(pattern = "^[^\t]+\t", replacement = ""))))
}
mir <- unique(mir[nchar(mir) >= 19])
kmers19 <- function(x) unlist(lapply(x, function(s) substring(s, 1:(nchar(s) - 18), 19:nchar(s))))
mir19 <- unique(c(kmers19(mir)))
cat(sprintf("miRNA sequences used for exclusion: %d (%d distinct 19-mers)\n", length(mir), length(mir19)))
is_mirna <- function(seqs) vapply(seqs, function(s) any(substring(s, 1:3, 19:21) %in% mir19), NA, USE.NAMES = FALSE)

# ---- alignments of the T-start 21-mers ---------------------------------------------------------------------------------------
aln <- read.delim(paste0(prefix, ".t21.aln"), header = FALSE, col.names = c("seq", "strand", "chrom", "offset"), stringsAsFactors = FALSE)
n_loci <- table(aln$seq)
cand <- data.frame(Full_21mer = names(n_loci), n_loci = as.integer(n_loci), stringsAsFactors = FALSE)
cand$reads <- collapsed$reads[match(cand$Full_21mer, collapsed$seq)]
cand$RPM <- cand$reads / cand$n_loci / total_mapped * 1e6
cand$miRNA <- is_mirna(cand$Full_21mer)
cat(sprintf("T-start 21-mers mapping exactly: %d unique sequences (%d reads); %d map to more than one locus; flagged as miRNA: %d\n",
    nrow(cand), sum(cand$reads), sum(cand$n_loci > 1), sum(cand$miRNA)))

# ---- Ruby motif for each locus (a sequence counts as having the motif if any of its loci has it) --------------------------------
genome <- read_fasta(genome_file)
aln$start <- aln$offset + 1L; aln$end <- aln$offset + 21L
up <- upstream_of(genome, aln$chrom, aln$start, aln$end, aln$strand)
aln$motif <- has_ruby_core(up)
cand$motif <- as.logical(tapply(aln$motif, factor(aln$seq, levels = cand$Full_21mer), any))
first <- !duplicated(aln$seq)
f <- aln[first, ]
cand$Chrom <- f$chrom[match(cand$Full_21mer, f$seq)]; cand$Start <- f$start[match(cand$Full_21mer, f$seq)]
cand$End <- f$end[match(cand$Full_21mer, f$seq)]; cand$Strand <- f$strand[match(cand$Full_21mer, f$seq)]
cand <- cand[order(-cand$RPM), ]
saveRDS(list(cand = cand, total_reads = total_reads, total_mapped = total_mapped), paste0(prefix, ".candidates.rds"))
utils::write.csv(cand[cand$reads >= 2, ], paste0(prefix, ".candidates_ge2reads.csv"), row.names = FALSE)
cat(sprintf("candidates (not miRNA): %d; with >= 2 reads: %d; >= 1 RPM: %d; with Ruby motif: %d (%.1f%%)\n", sum(!cand$miRNA), sum(!cand$miRNA & cand$reads >= 2),
    sum(!cand$miRNA & cand$RPM >= 1), sum(!cand$miRNA & cand$motif), 100 * mean(cand$motif[!cand$miRNA])))
