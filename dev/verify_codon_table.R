# Checks where the columns of inst/extdata/Codon_frequencies.csv come from, by recomputing them from public data.
# Run from the project root (needs internet for the first run; about 40 MB are downloaded into `work_dir`):
#   Rscript dev/verify_codon_table.R [work_dir]
#
# What it recomputes, and what the result was when it was written (2026-10-07):
#  * Chance            1 / (number of synonymous codons): exact.
#  * Ubiq              codon usage (pooled, relative frequency of each codon among the codons of its amino acid) of the 500 genes with the
#                      highest mean expression across the five tissues of Serizay et al. 2020 (Genome Res 30:1752), among the genes the paper
#                      classified as ubiquitous ("Ubiq."), in the longest isoform of each gene (Supplemental Table S2, sheet RNA_metrics,
#                      column Annotation; expression = mean of the five tissue TPM columns): matches to within 1e-10.
#  * CAI1, BringMansWeights   come from the web tool "C. elegans Codon Adapter" (https://worm.mpi-cbg.de/codons/cgi-bin/optimize.py, Redemann et al.
#                      2011) and cannot be recomputed here; they can only be checked against the tool's output.
#
# Five more columns were removed from the table on 2026-10-07 because the app did not use them (they are in the git history; see
# docs/DATA.md): Germ, Neur, Soma (the same rule as Ubiq for the germline-specific, neuron-specific and soma-restricted classes: matched to
# within 0.0007), aCDS (codon usage of all annotated coding sequences, every isoform: within 0.0006) and HighRCSU (relative synonymous codon
# usage; not reproducible exactly, closest to the 200 most highly expressed genes). Run this script on the commit before the removal to
# reproduce those numbers.
cf <- utils::read.csv("inst/extdata/Codon_frequencies.csv", stringsAsFactors = FALSE)
work <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(work)) work <- file.path(tempdir(), "codon_check")
dir.create(work, showWarnings = FALSE, recursive = TRUE)
options(timeout = 3600)   # the supplement is 29 MB from a slow server
non_stop <- cf$Amino != "X"
say <- function(...) cat(sprintf(...), "\n")
rel_freq <- function(v) v / ave(v, cf$Amino, FUN = sum)

# ---- downloads ---------------------------------------------------------------------------------------------------
cds_gz <- file.path(work, "caenorhabditis_elegans.CDS_transcripts.fa.gz")
if (!file.exists(cds_gz)) utils::download.file(
  "https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/current/species/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS19.CDS_transcripts.fa.gz",
  cds_gz, mode = "wb", quiet = TRUE)
suppl <- file.path(work, "serizay_supplement.zip")
if (!file.exists(suppl)) utils::download.file("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7706728/supplementaryFiles", suppl, mode = "wb", quiet = TRUE)
xlsx <- file.path(work, "Table_S2.xlsx")
if (!file.exists(xlsx)) { utils::unzip(suppl, "supp_gr.265934.120_Supplemental_Table_S2.xlsx", exdir = work); file.rename(file.path(work, "supp_gr.265934.120_Supplemental_Table_S2.xlsx"), xlsx) }

# ---- a minimal xlsx reader (base R only): returns the cells of one sheet as a data frame of character columns ----------------------
read_xlsx_sheet <- function(path, sheet) {
  dir <- file.path(dirname(path), paste0("xlsx_", tools::file_path_sans_ext(basename(path)))); utils::unzip(path, exdir = dir)
  strings <- paste(readLines(file.path(dir, "xl", "sharedStrings.xml"), warn = FALSE), collapse = "")
  si <- regmatches(strings, gregexpr("<si>.*?</si>", strings))[[1]]
  strings <- vapply(si, function(x) gsub("&amp;", "&", paste(gsub("<[^>]*>", "", regmatches(x, gregexpr("<t[^>]*>[^<]*</t>", x))[[1]]), collapse = "")), "", USE.NAMES = FALSE)
  xml <- paste(readLines(file.path(dir, "xl", "worksheets", sprintf("sheet%d.xml", sheet)), warn = FALSE), collapse = "")
  rows <- regmatches(xml, gregexpr("<row [^>]*>.*?</row>", xml))[[1]]
  col_index <- function(ref) vapply(strsplit(gsub("[0-9]", "", ref), ""), function(ch) sum((utf8ToInt(paste(ch, collapse = "")) - 64) * 26^(rev(seq_along(ch)) - 1)), 0)
  out <- lapply(rows, function(r) {
    cells <- regmatches(r, gregexpr("<c [^>]*?(/>|>.*?</c>)", r))[[1]]
    ref <- sub('.* r="([A-Z]+)[0-9]+".*', "\\1", cells); type <- ifelse(grepl(' t="s"', cells), "s", "n")
    val <- ifelse(grepl("<v>", cells), sub(".*<v>([^<]*)</v>.*", "\\1", cells), NA_character_)
    val[type == "s" & !is.na(val)] <- strings[as.integer(val[type == "s" & !is.na(val)]) + 1]
    v <- rep(NA_character_, max(col_index(ref))); v[col_index(ref)] <- val; v
  })
  width <- max(lengths(out)); do.call(rbind, lapply(out, function(v) c(v, rep(NA_character_, width - length(v)))))
}
s2 <- read_xlsx_sheet(xlsx, 2)
s2 <- s2[!startsWith(ifelse(is.na(s2[, 1]), "", s2[, 1]), "##"), ]
colnames(s2) <- s2[1, ]; s2 <- as.data.frame(s2[-1, ], stringsAsFactors = FALSE)
tissues <- c("Germline", "Neurons", "Muscle", "Hypod.", "Intest.")
s2$mean_tpm <- rowMeans(sapply(s2[, tissues], as.numeric))
say("Serizay Table S2: %d genes; classes: %s", nrow(s2), paste(names(sort(table(s2$Annotation), decreasing = TRUE))[1:5], collapse = ", "))

# ---- coding sequences ---------------------------------------------------------------------------------------------------------
lines <- readLines(gzfile(cds_gz)); head_at <- grep("^>", lines); last <- c(head_at[-1] - 1, length(lines))
seq <- mapply(function(a, b) toupper(paste(lines[(a + 1):b], collapse = "")), head_at, last)
gene <- sub(".*gene=([^ ]+).*", "\\1", lines[head_at])
ok <- nchar(seq) %% 3 == 0 & grepl("^ATG", seq) & !grepl("[^ACGT]", seq)
codons_of <- function(x) substring(x, seq(1, nchar(x) - 5, 3), seq(3, nchar(x) - 3, 3))            # without the stop codon
count_codons <- function(s) { tab <- table(factor(tolower(unlist(lapply(s, codons_of))), levels = cf$Codon)); as.numeric(tab[cf$Codon]) }
longest <- function(genes) { i <- which(ok & gene %in% genes); i <- i[order(gene[i], -nchar(seq[i]))]; i[!duplicated(gene[i])] }

report <- function(label, column, computed) {
  d <- abs(computed - cf[[column]])[non_stop]
  say("%-60s %-9s mean |difference| %.5f   max %.5f   correlation %.4f", label, column, mean(d), max(d), cor(computed[non_stop], cf[[column]][non_stop]))
}
report("equal weight for each synonymous codon", "Chance", 1 / ave(rep(1, nrow(cf)), cf$Amino, FUN = length))
for (cl in list(c("Ubiq.", "Ubiq"))) {
  g <- s2$WBID[s2$Annotation == cl[1]]
  g <- g[order(-s2$mean_tpm[match(g, s2$WBID)])][1:500]
  i <- longest(g)
  report(sprintf("top 500 %s genes by mean TPM (%d found)", cl[1], length(i)), cl[2], rel_freq(count_codons(seq[i])))
}
