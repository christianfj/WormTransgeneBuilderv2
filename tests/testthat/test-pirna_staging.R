# Checks on the rebuilt piRNA database while it is still in data-raw/staging/ (see data-raw/pirna_rebuild_brief.md).
# When the staging file is promoted to inst/extdata/piRNA_Database.csv these checks move to that file.
# Skipped when the staging file is not there (for example in an installed package).

staging_file <- testthat::test_path("..", "..", "data-raw", "staging", "piRNA_Database_staging.csv")
staging <- function() {
  skip_if(!file.exists(staging_file), "the staging piRNA database is not in this copy")
  utils::read.csv(staging_file, stringsAsFactors = FALSE, na.strings = "")
}

test_that("the staging database has the columns the app uses plus the provenance columns", {
  d <- staging()
  expect_identical(names(d)[1:5], c("Name", "Target_20mer", "Full_21mer", "Abundance", "Enrichment"))
  expect_true(all(c("Type", "Source", "WS_version", "WBGene", "Chrom", "Start", "End", "Strand", "n_genomic_hits", "PRG1_IP_avg_RPM",
                    "PRG1_IP_avg_enrichment", "PRG1_IP_enriched", "PRG1_depleted_log2FC", "PRG1_depleted_padj",
                    "Claycomb_location", "Claycomb_class") %in% names(d)))
  expect_true(all(d$Type %in% c("type 1", "type 2", "PRG-1 IP new")))
})

test_that("it holds 15,363 type 1 (WS298), 2,485 type 2 and 359 new PRG-1 IP entries", {
  d <- staging()
  expect_equal(sum(d$Type == "type 1"), 15363)
  expect_equal(sum(d$Type == "type 2"), 2485)
  expect_equal(sum(d$Type == "PRG-1 IP new"), 359)
  expect_equal(nrow(d), 18207)
  expect_true(all(d$WS_version == "WS298"))
  expect_true(all(d$Source[d$Type == "type 1"] == "WormBase WS298"))
  expect_match(d$Source[d$Type == "type 2"], "Wu et al. 2018")      # the justification for the type 2 set
})

test_that("type 2 entries are named '21ur-N (type 2)'; new PRG-1 IP entries get internal labels, never a 21ur- gene name", {
  d <- staging()
  expect_match(d$Name[d$Type == "type 1"], "^21ur-[0-9]+$")
  expect_match(d$Name[d$Type == "type 2"], "^21ur-[0-9]+ \\(type 2\\)$")
  expect_match(d$Name[d$Type == "PRG-1 IP new"], "^new21U_(I|II|III|IV|V|X)_[0-9]+_(plus|minus)$")
  expect_false(any(grepl("^21ur-", d$Name[d$Type == "PRG-1 IP new"])))
})

test_that("all 466 Claycomb new-21U loci are in the database (359 entries, 30 of them as type 2 entries)", {
  d <- staging()
  locs <- unlist(strsplit(d$Claycomb_location[!is.na(d$Claycomb_location)], ";"))
  expect_equal(length(unique(locs)), 466)
  expect_equal(sum(!is.na(d$Claycomb_location) & d$Type == "type 2"), 30)
  expect_equal(sum(!is.na(d$Claycomb_location) & d$Type == "PRG-1 IP new"), 359)
  # a sequence that sits at several loci is one entry that lists them all
  multi <- d[d$Type == "PRG-1 IP new" & d$n_genomic_hits > 1, ]
  expect_equal(nrow(multi), 58)
  expect_equal(sum(multi$n_genomic_hits), 135)
  expect_equal(lengths(strsplit(multi$Claycomb_location, ";")), multi$n_genomic_hits)
})

test_that("every entry is a clean 21-mer with a matching 20-base target", {
  d <- staging()
  expect_false(anyNA(d[, c("Name", "Target_20mer", "Full_21mer", "Abundance", "Enrichment")]))
  expect_true(all(grepl("^[ACGT]{21}$", d$Full_21mer)))
  expect_true(all(grepl("^[ACGT]{20}$", d$Target_20mer)))
  expect_identical(d$Target_20mer, reverse_complement(substring(d$Full_21mer, 2, 21)))
})

test_that("names and sequences are unique", {
  d <- staging()
  expect_equal(anyDuplicated(d$Name), 0)
  expect_equal(anyDuplicated(d$Full_21mer), 0)
  expect_equal(anyDuplicated(d$Target_20mer), 0)
})

test_that("579 type 1 piRNAs do not start with T; all type 2 and new entries do", {
  d <- staging()
  t1 <- d[d$Type == "type 1", ]
  expect_equal(sum(substr(t1$Full_21mer, 1, 1) != "T"), 579)
  expect_gt(mean(substr(t1$Full_21mer, 1, 1) == "T"), 0.96)
  expect_true(all(substr(d$Full_21mer[d$Type != "type 1"], 1, 1) == "T"))
})

test_that("Abundance is the PRG-1 IP RPM from the Claycomb tables, 'not detected in PRG-1 IP', or (type 2 only) 'not in Claycomb PRG-1 IP data'", {
  d <- staging()
  detected <- !is.na(d$PRG1_IP_avg_RPM) & d$PRG1_IP_avg_RPM > 0.01
  expect_equal(sum(detected), 13587 + 30 + 359)
  expect_match(d$Abundance[detected], "^[0-9.e+-]+ RPM \\(PRG-1 IP\\)$")
  expect_equal(as.numeric(sub(" RPM.*", "", d$Abundance[detected])), signif(d$PRG1_IP_avg_RPM[detected], 3))
  expect_identical(d$Abundance[d$Type == "type 1" & !detected], rep("not detected in PRG-1 IP", 1776))
  expect_identical(d$Abundance[d$Type == "type 2" & !detected], rep("not in Claycomb PRG-1 IP data", 2455))
  expect_false(any(d$Type == "PRG-1 IP new" & !detected))
})

test_that("Enrichment is Billi et al.'s male/female call, or 'No Data' (751 not in Billi, 1,203 n.d., all type 2 and new entries)", {
  d <- staging()
  expect_equal(sum(d$Enrichment == "No Data"), 751 + 1203 + 2485 + 359)
  called <- d$Enrichment[d$Enrichment != "No Data"]
  expect_match(called, "^(male|female|not enriched) \\(score: -?[0-9]+\\)$")
  expect_true(all(d$Enrichment[d$Type != "type 1"] == "No Data"))
})

test_that("the anchor piRNA 21ur-5581 (WBGene00165728) is correct", {
  d <- staging()
  a <- d[d$Name == "21ur-5581", ]
  expect_equal(nrow(a), 1)
  expect_identical(a$WBGene, "WBGene00165728")
  expect_identical(c(a$Chrom, a$Strand), c("IV", "-"))
  expect_equal(c(a$Start, a$End), c(13843030, 13843050))
  expect_equal(a$PRG1_IP_avg_RPM, (1.356 + 2.174) / 2, tolerance = 0.01)
})

test_that("21ur-2868, which WormBase no longer annotates, is not in the database", {
  expect_false("21ur-2868" %in% staging()$Name)
})

test_that("every type 1 entry has a WBGene and its genomic interval is 21 bases", {
  d <- staging()
  expect_true(all(is.na(d$WBGene[d$Type != "type 1"])))
  ws <- d[d$Type == "type 1", ]
  expect_match(ws$WBGene, "^WBGene[0-9]{8}$")
  expect_true(all(ws$End - ws$Start == 20))
  expect_true(all(ws$Strand %in% c("+", "-")))
})
