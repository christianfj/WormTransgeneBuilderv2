# Promotes the rebuilt piRNA databases from data-raw/staging/ to the files the app reads. Run from the project root:
#   Rscript dev/make_pirna_data.R
#
#   inst/extdata/piRNA_Database.csv                      C. elegans (WormBase WS298 type 1, type 2, new PRG-1 IP entries)
#   inst/extdata/piRNA_Database_<species>.csv.gz         C. briggsae, C. remanei, C. brenneri
#   data/pirna_sequences.rda                             the Target_20mer column of the C. elegans file, in the same order
#
# The databases are described in docs/PIRNA_DATABASE.md and docs/PIRNA_OTHER_SPECIES.md; the staging files are built by the scripts
# in data-raw/ and data-raw/other_species/. The order of the C. elegans rows must match pirna_sequences (the app finds a piRNA's
# row from its position), which is checked here and by the tests.
cel <- utils::read.csv("data-raw/staging/piRNA_Database_staging.csv", na.strings = "", stringsAsFactors = FALSE)
stopifnot(identical(names(cel)[1:5], c("Name", "Target_20mer", "Full_21mer", "Abundance", "Enrichment")),
          !anyNA(cel$Target_20mer), !anyDuplicated(cel$Target_20mer), !anyDuplicated(cel$Name), all(grepl("^[ACGT]{20}$", cel$Target_20mer)))
file.copy("data-raw/staging/piRNA_Database_staging.csv", "inst/extdata/piRNA_Database.csv", overwrite = TRUE)
pirna_sequences <- cel$Target_20mer
save(pirna_sequences, file = "data/pirna_sequences.rda", compress = "bzip2", version = 3)
for (sp in c("briggsae", "remanei", "brenneri")) {
  src <- sprintf("data-raw/staging/other_species/piRNA_Database_%s_staging.csv.gz", sp)
  d <- utils::read.csv(src, na.strings = "", stringsAsFactors = FALSE)
  stopifnot(!anyNA(d$Target_20mer), !anyDuplicated(d$Target_20mer), all(grepl("^[ACGT]{20}$", d$Target_20mer)))
  file.copy(src, sprintf("inst/extdata/piRNA_Database_%s.csv.gz", sp), overwrite = TRUE)
  cat(sp, nrow(d), "entries\n")
}
cat("C. elegans", nrow(cel), "entries\n")
