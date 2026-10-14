# Checks on the staging piRNA databases for other Caenorhabditis species (data-raw/staging/other_species/, see
# data-raw/pirna_other_species_brief.md). Skipped when a staging file is not in this copy of the project.

staging_species <- function(file) {
  path <- testthat::test_path("..", "..", "data-raw", "staging", "other_species", file)
  skip_if(!file.exists(path), paste(file, "is not in this copy of the project"))
  utils::read.csv(path, stringsAsFactors = FALSE, na.strings = "")
}

expect_clean_database <- function(d) {
  expect_false(anyNA(d[, c("Name", "Target_20mer", "Full_21mer", "Abundance", "Species", "Strain", "Assembly", "Confidence")]))
  expect_true(all(grepl("^[ACGT]{21}$", d$Full_21mer)))
  expect_true(all(grepl("^[ACGT]{20}$", d$Target_20mer)))
  expect_identical(d$Target_20mer, reverse_complement(substring(d$Full_21mer, 2, 21)))
  expect_equal(anyDuplicated(d$Name), 0)
  expect_equal(anyDuplicated(d$Full_21mer), 0)
  expect_equal(anyDuplicated(d$Target_20mer), 0)
  expect_true(all(d$End - d$Start == 20))
  expect_equal(anyDuplicated(d[, c("Chrom", "Start", "Strand")]), 0)
  expect_true(all(d$Strand %in% c("+", "-")))
}

test_that("the C. briggsae staging database holds the 25,883 PRG-1 IP loci of Pastore et al. 2022", {
  d <- staging_species("piRNA_Database_briggsae_staging.csv.gz")
  expect_equal(nrow(d), 25883)
  expect_clean_database(d)
  expect_true(all(substr(d$Full_21mer, 1, 1) == "T"))
  expect_true(all(d$Species == "Caenorhabditis briggsae"))
  expect_true(all(d$Strain == "AF16"))
  expect_true(all(d$Confidence == "high"))
  expect_match(d$Abundance, "^[0-9.e+-]+ RPM \\(PRG-1 IP\\)$")
  expect_true(all(d$Enrichment == "No Data"))
  expect_type(d$Motif_Ruby, "logical")
})

test_that("the briggsae loci are on the expected chromosomes and the first locus is as in the paper's table", {
  d <- staging_species("piRNA_Database_briggsae_staging.csv.gz")
  expect_gt(sum(d$Chrom == "IV"), 17000)
  first <- d[d$Name == "Cbr-21UR-1", ]
  expect_identical(first$Full_21mer, "TGACAATTAGTTGAAGTGTGC")
  expect_identical(c(first$Chrom, first$Strand), c("I", "-"))
  expect_equal(c(first$Start, first$End), c(10643792, 10643812))   # the table's 0-based start 10643791 + 1
})

# ---- C. remanei and C. brenneri: low-confidence sets called from Shi et al. 2013 reads --------------------------------------

test_that("the remanei and brenneri databases are clean, low confidence, and tiered A or B", {
  for (sp in c("remanei", "brenneri")) {
    d <- staging_species(sprintf("piRNA_Database_%s_staging.csv.gz", sp))
    expect_clean_database(d)
    expect_true(all(substr(d$Full_21mer, 1, 1) == "T"), info = sp)
    expect_true(all(d$Confidence == "low"), info = sp)
    expect_true(all(d$Tier == "A"), info = sp)
    expect_true(all(d$Enrichment == "No Data"), info = sp)
    expect_gt(nrow(d), 10000)
  }
})

test_that("every remanei and brenneri entry has the Ruby motif and at least 1 RPM in some library (tier A only; tier B was dropped)", {
  for (sp in c("remanei", "brenneri")) {
    d <- staging_species(sprintf("piRNA_Database_%s_staging.csv.gz", sp))
    expect_true(all(d$Motif_Ruby), info = sp)
    expect_true(all(pmax(d$RPM_adult, d$RPM_embryo, d$RPM_male) >= 1), info = sp)
  }
})

test_that("each species names its strain and assembly, and names never use a formal 21ur- gene name", {
  r <- staging_species("piRNA_Database_remanei_staging.csv.gz")
  b <- staging_species("piRNA_Database_brenneri_staging.csv.gz")
  expect_true(all(r$Strain == "PB4641")); expect_true(all(grepl("PRJNA53967", r$Assembly)))
  expect_true(all(b$Strain == "PB2801")); expect_true(all(grepl("PRJNA20035", b$Assembly)))
  expect_match(r$Name, "^Cre_21U_"); expect_match(b$Name, "^Cbn_21U_")
})
