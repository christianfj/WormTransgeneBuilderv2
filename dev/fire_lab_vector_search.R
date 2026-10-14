# Parses the Addgene GenBank file of the whole Fire lab vector kit (data-raw/fire_lab_vectors/, 286 records) into a table of
# vectors (name, length, sequence) and searches them for the Fire lab introns of the app. Run from the project root:
#   Rscript dev/fire_lab_vector_search.R
source("dev/fire_lab_parse.R")
vec <- read_fire_lab_vectors()
tab <- data.frame(id = vapply(vec, `[[`, "", "id"), name = names(vec), len = nchar(vapply(vec, `[[`, "", "seq")), stringsAsFactors = FALSE)
seqs <- lapply(vec, `[[`, "seq")
rc <- reverse_complement_dna
count <- function(pat, s) { n <- 0; for (p in unique(c(pat, rc(pat)))) { m <- gregexpr(p, s, fixed = TRUE)[[1]]; if (m[1] > 0) n <- n + length(m) }; n }

introns <- utils::read.csv("inst/extdata/Introns.csv", stringsAsFactors = FALSE)
app <- toupper(unlist(introns[introns$Name == "Fire_lab_introns", c("Sequence_1", "Sequence_2", "Sequence_3")]))
delta_variant <- "GTAAGTTTAAACTTGGACTTACTAACTAACGGATTATATTTAAATTTTCAG"   # the third gfp intron in the pPD95.75 file the authors supplied (ApE, 2018)
targets <- c(alpha = app[[1]], beta = app[[2]], app3 = app[[3]], delta_variant = delta_variant)

cat("Records read:", nrow(tab), "; total length", sum(tab$len), "bp\n\n")
hits <- sapply(targets, function(t) vapply(seqs, count, 0, pat = t))
colnames(hits) <- names(targets)
rownames(hits) <- tab$name
print(table(Records_with_intron = apply(hits > 0, 2, sum)))
cat("\nNumber of vectors containing each intron exactly (either strand):\n"); print(colSums(hits > 0))
for (n in colnames(hits)) { cat("\n", n, ": ", sep = ""); h <- rownames(hits)[hits[, n] > 0]; cat(length(h), "vectors:", paste(head(h, 60), collapse = ", "), if (length(h) > 60) "...", "\n") }

# The third gfp intron (the one after ...GAATGGAATCAAAGTT) in every vector that has the Fire lab gfp coding region
anchor <- "AAAAGAATGGAATCAAAGTT"
get <- function(s) { for (z in c(s, rc(s))) { p <- regexpr(anchor, z, fixed = TRUE)[1]; if (p > 0) return(substr(z, p + nchar(anchor), p + nchar(anchor) + 50)) }; NA_character_ }
third <- vapply(seqs, get, "")
cat("\nThe 51 bases after the exon that precedes the third gfp intron, in", sum(!is.na(third)), "vectors:\n")
tt <- sort(table(third[!is.na(third)]), decreasing = TRUE)
for (i in seq_along(tt)) cat(sprintf("  %3d vectors: %s\n      e.g. %s\n", tt[i], names(tt)[i], paste(head(names(third)[which(third == names(tt)[i])], 8), collapse = ", ")))
