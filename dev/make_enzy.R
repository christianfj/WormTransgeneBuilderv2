# Build data/enzy.rda, the table of restriction enzymes the removal code looks up by name,
# from inst/extdata/Enzymes.csv (the list the app also shows in its enzyme menu).
#
# Run from the project root:   Rscript dev/make_enzy.R
#
# Enzymes.csv has one row per enzyme with four columns:
#   Enzyme  the enzyme's name (unique)
#   Site    its recognition sequence, written with IUPAC letters (for example GGYRCC for BanI,
#           or GACNNNNNGTC where the middle bases can be anything)
#   GG      TRUE for the Golden Gate enzymes (BsaI, BsmBI, BbsI, SapI, ...)
#   Common  TRUE for the enzymes most often avoided in a transgene; they are listed first in the
#           app's enzyme menu (a Golden Gate enzyme is listed under Golden Gate instead)
# The list is New England Biolabs' full enzyme catalogue (234 enzymes). To add or change an
# enzyme, edit that file and run this script; the tests fail if the two disagree.
#
# Methylation (Dam/Dcm sensitivity) is deliberately not recorded: a synthesized transgene is not
# methylated, so every site is treated as a site whatever the later cloning host does.

enzymes <- read.csv(file.path("inst", "extdata", "Enzymes.csv"), stringsAsFactors = FALSE)

stopifnot(
  identical(names(enzymes), c("Enzyme", "Site", "GG", "Common")),
  !anyNA(enzymes),
  !anyDuplicated(enzymes$Enzyme),
  all(grepl("^[A-Za-z0-9_]+$", enzymes$Enzyme)),
  all(grepl("^[ACGTRYSWKMBDHVN]+$", enzymes$Site)),
  is.logical(enzymes$GG),
  is.logical(enzymes$Common)
)

# One row per enzyme, found by name: enzy["BsaI", "Site"]
enzy <- enzymes
rownames(enzy) <- enzy$Enzyme

# The same format as the other small datasets (bzip2, format version 3).
save(enzy, file = file.path("data", "enzy.rda"), compress = "bzip2", version = 3)
cat(sprintf("Wrote data/enzy.rda: %d enzymes, columns %s.\n", nrow(enzy), paste(names(enzy), collapse = ", ")))
