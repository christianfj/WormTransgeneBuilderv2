# Build inst/extdata/glo_scores.rds, the compact GLO score table the app actually reads.
#
# Run from the project root:   Rscript dev/make_glo_scores.R
#
# Background: data/glo_full.rda is a named numeric vector with one entry for every possible
# 12-mer of A/C/G/T: 4^12 = 16,777,216 entries, in alphabetical order (A < C < G < T). Because
# it is complete and in that order, the name of entry k is just k - 1 written in base 4
# (AAAAAAAAAAAA = 0, AAAAAAAAAAAC = 1, ..., TTTTTTTTTTTT = 4^12 - 1). The names therefore carry
# no information, but as text they take 1.15 GB of the 1.28 GB the object needs in memory and
# most of the 9 to 10 seconds it takes to load. The score table keeps only the values, in the
# same order, and a 12-mer's position is worked out from its bases (see R/fct_glo_scores.R).
#
# The values are stored exactly (as doubles), compressed with gzip: 12 MB on disk, 128 MB in
# memory, about 0.2 s to load.
#
# The script checks that the names really are the base-4 numbers 0, 1, 2, ..., that no value
# is missing, and that reading the new file back gives exactly the original values.

# data/glo_full.rda has since been removed from the repository (nothing used it once the score
# table existed). Its last version is kept in git history, commit bdee5b0, and is read from there.
glo_rda <- tempfile(fileext = ".rda")
system2("git", c("show", "bdee5b0:data/glo_full.rda"), stdout = glo_rda)
env <- new.env()
load(glo_rda, envir = env)
glo_full <- env$glo_full

names_glo <- names(glo_full)
stopifnot(length(names_glo) == 4^12,
          all(grepl("^[ACGT]{12}$", names_glo)),
          identical(strtoi(chartr("ACGT", "0123", names_glo), base = 4L), 0:(4^12 - 1)))

scores <- unname(glo_full)
stopifnot(is.double(scores), !anyNA(scores), length(scores) == 4^12)

out <- file.path("inst", "extdata", "glo_scores.rds")
saveRDS(scores, out, compress = "gzip")

check <- readRDS(out)
stopifnot(identical(check, scores))
cat(sprintf("Wrote %s: %d scores, %.1f MB on disk; identical to the values in glo_full.\n",
            out, length(check), file.size(out) / 1024^2))
cat("Hash of the score vector:", rlang::hash(check), "\n")
