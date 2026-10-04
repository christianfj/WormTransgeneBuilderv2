# Converts the owner's tab-separated promoter list (Name, Gene, Sequence; group headings are rows
# with an empty Sequence) into inst/extdata/Promoters.csv. Run from the project root:
#   Rscript dev/make_promoters.R path/to/promoters_expanded.csv
# Small tidying done here: heading rows keep all their text in Name, and a row whose Name also holds the
# gene in brackets ("pAFD (WBGene00001535, gcy-8)") is split into Name and Gene.
src <- commandArgs(trailingOnly = TRUE)[1]
p <- read.delim(src, colClasses = "character", na.strings = character(0), quote = "", check.names = FALSE)
stopifnot(identical(names(p), c("Name", "Gene", "Sequence")))
p[] <- lapply(p, function(x) trimws(gsub("\\s+", " ", x)))
# A row missing its Gene tab (pAFD) has the sequence in the Gene column: move it across.
shifted <- p$Sequence == "" & grepl("^[ACGTacgt]{50,}$", p$Gene)
p$Sequence[shifted] <- p$Gene[shifted]
p$Gene[shifted] <- ""
heading <- p$Sequence == ""
p$Name[heading] <- trimws(paste(p$Name[heading], p$Gene[heading]))
p$Gene[heading] <- ""
split_gene <- !heading & p$Gene == "" & grepl(" \\(.*\\)$", p$Name)
p$Gene[split_gene] <- sub("^.*? (\\(.*\\))$", "\\1", p$Name[split_gene])
p$Name[split_gene] <- sub(" \\(.*\\)$", "", p$Name[split_gene])
stopifnot(!anyDuplicated(p$Name), all(grepl("^[ACGTacgt]+$", p$Sequence[!heading])), heading[1])
write.csv(p, "inst/extdata/Promoters.csv", row.names = FALSE)
cat(sum(heading), "groups,", sum(!heading), "promoters\n")
