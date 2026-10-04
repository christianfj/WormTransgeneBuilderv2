# The enzyme dropdown is filled from inst/extdata/Enzymes.csv, but the removal
# step looks enzymes up by name in the internal `enzy` table. If the two ever
# disagree, an enzyme can be chosen in the app and then silently ignored.

read_enzyme_csv <- function() {
  read.csv(system.file("extdata", "Enzymes.csv", package = "WormTransgeneBuilder"),
           stringsAsFactors = FALSE)
}

test_that("every enzyme offered in the app exists in the removal table", {
  csv <- read_enzyme_csv()
  expect_equal(setdiff(csv$Enzyme, rownames(WormTransgeneBuilder::enzy)), character(0))
})

test_that("the removal table and the CSV list the same enzymes", {
  csv <- read_enzyme_csv()
  expect_setequal(rownames(WormTransgeneBuilder::enzy), csv$Enzyme)
})

test_that("the CSV and the removal table agree on every recognition site", {
  csv <- read_enzyme_csv()
  expect_equal(toupper(WormTransgeneBuilder::enzy[csv$Enzyme, "Site"]),
               toupper(csv$Site))
})

test_that("every recognition site is a valid IUPAC DNA string", {
  expect_true(all(grepl("^[ACGTRYSWKMBDHVN]+$", toupper(read_enzyme_csv()$Site))))
})

test_that("the default enzyme (BsaI) is in the list", {
  expect_true("BsaI" %in% read_enzyme_csv()$Enzyme)
})

test_that("every enzyme in the table can actually be removed by name", {
  # A site placed in frame would be ideal, but here we only check the lookup does not error, keeps
  # the length and keeps the protein. A few sites genuinely cannot be broken by a synonymous change
  # (HpyCH4III, ACNGT, in the fixture's Thr-Val pair ACT GTT: every Thr codon is ACN and every Val
  # codon is GTN), and the function then warns, which is correct, so warnings are not errors here.
  csv <- read_enzyme_csv()
  for (enz in csv$Enzyme) {
    set.seed(1)
    out <- suppressWarnings(remove_enzyme_sites(fixture_cds(), enz, "Ubiq"))
    expect_equal(nchar(out), nchar(fixture_cds()), info = enz)
    expect_equal(protein_of(out), protein_of(fixture_cds()), info = enz)
  }
})

# ---- the table is New England Biolabs' full list, without methylation -----------------------------

test_that("the enzyme table has exactly the columns Enzyme, Site, GG and Common", {
  # Dam/Dcm methylation is deliberately not recorded: a synthesized transgene is not methylated, so
  # a site is a site whatever the later cloning host does.
  expect_identical(names(read_enzyme_csv()), c("Enzyme", "Site", "GG", "Common"))
  expect_identical(names(WormTransgeneBuilder::enzy), c("Enzyme", "Site", "GG", "Common"))
  expect_false("DamDcm" %in% c(names(read_enzyme_csv()), names(WormTransgeneBuilder::enzy)))
})

test_that("enzyme names are unique, and there are hundreds of enzymes", {
  csv <- read_enzyme_csv()
  expect_equal(anyDuplicated(csv$Enzyme), 0)
  expect_gt(nrow(csv), 200)
  expect_match(csv$Enzyme, "^[A-Za-z0-9_]+$")
})

test_that("the Golden Gate flag is logical and marks the expected enzymes", {
  csv <- read_enzyme_csv()
  expect_type(csv$GG, "logical")
  expect_false(anyNA(csv$GG))
  expect_true(all(c("BsaI", "BsmBI", "BbsI", "SapI", "BspQI", "Esp3I", "PaqCI") %in% csv$Enzyme[csv$GG]))
  expect_false("EcoRI" %in% csv$Enzyme[csv$GG])
})

test_that("every site is at least 4 bases long", {
  expect_true(all(nchar(read_enzyme_csv()$Site) >= 4))
})

test_that("the 21 enzymes of the earlier list are still there with the same sites", {
  expected <- c(ApaI = "GGGCCC", ApaLI = "GTGCAC", BamHI = "GGATCC", BanI = "GGYRCC", BsaI = "GGTCTC",
                BsmBI = "CGTCTC", BsoBI = "CYCGRG", DraI = "TTTAAA", EcoRI = "GAATTC", EcoRV = "GATATC",
                HindIII = "AAGCTT", KpnI = "GGTACC", NcoI = "CCATGG", NdeI = "CATATG", NheI = "GCTAGC",
                PstI = "CTGCAG", PvuII = "CAGCTG", SacII = "CCGCGG", SapI = "GCTCTTC", XbaI = "TCTAGA",
                XhoI = "CTCGAG")
  expect_setequal(names(expected), original_21_enzymes)
  expect_identical(WormTransgeneBuilder::enzy[names(expected), "Site"], unname(expected))
})

test_that("the commonly used enzymes are all in the table", {
  expect_true(all(common_enzymes %in% rownames(WormTransgeneBuilder::enzy)))
  expect_length(common_enzymes, 24)
})

test_that("enzymes with an N in the site (a run of any bases) are present and searched", {
  # For example AhdI (GACNNNNNGTC) and SfiI (GGCCNNNNNGGCC); the finder handles N in a site.
  n_sites <- WormTransgeneBuilder::enzy$Site[grepl("N", WormTransgeneBuilder::enzy$Site)]
  expect_gt(length(n_sites), 20)
  expect_equal(nrow(find_enzyme_site_ranges("AAGACAAAAAGTCAA", "GACNNNNNGTC")), 1)
})
