# Calibrates the Ruby-motif flag on the C. elegans database: what fraction of type 1, type 2 and random 21-mers have a
# GTTTC core at positions -50..-35? Usage: Rscript data-raw/other_species/calibrate_motif.R <C. elegans genomic.fa>
source("data-raw/other_species/common.R")
g <- read_fasta(commandArgs(trailingOnly = TRUE)[1])
db <- read.csv("data-raw/staging/piRNA_Database_staging.csv", na.strings = "", stringsAsFactors = FALSE)
one <- db[!is.na(db$Chrom) & !is.na(db$Start), ]
up <- upstream_of(g, one$Chrom, one$Start, one$End, one$Strand)
one$motif <- has_ruby_core(up)
cat("with a position:\n"); print(tapply(one$motif, one$Type, function(x) sprintf("%d of %d (%.1f%%)", sum(x), length(x), 100 * mean(x))))
t1 <- one[one$Type == "type 1", ]
cat(sprintf("type 1, with T first: %.1f%%; non-T first: %.1f%%\n", 100 * mean(t1$motif[substr(t1$Full_21mer, 1, 1) == "T"]), 100 * mean(t1$motif[substr(t1$Full_21mer, 1, 1) != "T"])))
set.seed(1)
rs <- one[sample(nrow(one), 5000), ]
shift <- sample(1000:50000, 5000, TRUE) * sample(c(-1, 1), 5000, TRUE)
st <- pmax(100, pmin(rs$Start + shift, nchar(g[rs$Chrom]) - 200))
ru <- upstream_of(g, rs$Chrom, st, st + 20, rs$Strand)
cat(sprintf("random positions (same chromosomes and strands, 1-50 kb away): %.1f%% have the motif\n", 100 * mean(has_ruby_core(ru), na.rm = TRUE)))
