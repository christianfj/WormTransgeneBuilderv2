# export_genbank() writes the GenBank (.gb) file that users open in ApE or Benchling.
# These tests read the file back with parse_genbank() (in helper-fixtures.R) and check
# that what comes out is what went in: the sequence, the feature coordinates and the
# comments. The last group runs the real download buttons.

gb_seq <- function() "ATGGCTAAAGGTTAA"

# ---- the header --------------------------------------------------------------

test_that("the file starts with LOCUS, ends with '//' and has the standard sections", {
  gb <- export_genbank("test", gb_seq())
  expect_true(startsWith(gb[1], "LOCUS"))
  expect_equal(gb[length(gb)], "//")
  expect_true(all(c("FEATURES             Location/Qualifiers", "ORIGIN") %in% gb))
  expect_true(any(startsWith(gb, "DEFINITION")))
  expect_true(any(startsWith(gb, "ORGANISM")))
})

test_that("the LOCUS line has the name, the length, the type and a date", {
  locus <- export_genbank("my gene", gb_seq())[1]
  expect_match(locus, "my gene")
  expect_match(locus, "15bp ds-DNA")
  expect_match(locus, "linear")
  expect_match(locus, "[0-9]{2}-[A-Z][a-z]{2}-[0-9]{4}$")
})

test_that("genbank_date writes a two-digit day, a three-letter month and the year", {
  expect_equal(genbank_date(as.Date("2026-09-29")), "29-Sep-2026")
  expect_equal(genbank_date(as.Date("2026-09-09")), "09-Sep-2026")  # single-digit day
  expect_equal(genbank_date(as.Date("2026-01-01")), "01-Jan-2026")
  expect_equal(genbank_date(as.Date("2026-12-31")), "31-Dec-2026")
  expect_equal(genbank_date(as.Date("2028-02-29")), "29-Feb-2028")
})

test_that("genbank_date is correct for every day of a month", {
  for (d in 1:9) {
    expect_match(genbank_date(as.Date(sprintf("2026-03-%02d", d))), sprintf("^0%d-Mar-2026$", d))
  }
})

test_that("genbank_date uses English month names whatever the computer's language", {
  old <- Sys.getlocale("LC_TIME")
  on.exit(Sys.setlocale("LC_TIME", old), add = TRUE)
  suppressWarnings(Sys.setlocale("LC_TIME", "de_DE.UTF-8"))
  expect_equal(genbank_date(as.Date("2026-03-15")), "15-Mar-2026")
  expect_equal(genbank_date(as.Date("2026-05-15")), "15-May-2026")
  expect_equal(genbank_date(as.Date("2026-10-15")), "15-Oct-2026")
})

test_that("the LOCUS date is today's date", {
  expect_true(endsWith(export_genbank("t", gb_seq())[1], genbank_date(Sys.Date())))
})

# ---- comments ----------------------------------------------------------------

test_that("scores appear as comments only when they are supplied", {
  none <- export_genbank("t", gb_seq())
  expect_false(any(grepl("Codon Adaptation Index|GLO Score|COMMENT     GC", none)))
  gb <- export_genbank("t", gb_seq(), cai_score = 0.96, gc_score = "40%", glo_score = "1,234")
  expect_true("COMMENT     Codon Adaptation Index 0.96" %in% gb)
  expect_true("COMMENT     GC 40%" %in% gb)
  expect_true("COMMENT     GLO Score 1,234" %in% gb)
})

test_that("the input sequence and extra comments are written, one per line", {
  gb <- export_genbank("t", gb_seq(), input_seq = "ATGGCT",
                       extracomments = c("Optimization routine: CAI1", "Remove RE sites: Yes"))
  expect_true("COMMENT     Input sequence: ATGGCT" %in% gb)
  expect_true(any(grepl("^COMMENT +Optimization routine: CAI1$", gb)))
  expect_true(any(grepl("^COMMENT +Remove RE sites: Yes$", gb)))
})

# ---- features ------------------------------------------------------------------

test_that("each feature gets a location, a label, colours and the ApE arrow format", {
  gb <- export_genbank("t", gb_seq(), c(1, 13), c(3, 15), c("#a9dfbf", "#e6b0aa"),
                       c("Start", "Stop"), c("misc_feature", "misc_feature"))
  p <- parse_genbank(gb)
  expect_equal(p$starts, c(1, 13))
  expect_equal(p$ends, c(3, 15))
  expect_equal(p$labels, c("Start", "Stop"))
  expect_equal(vapply(p$features, `[[`, "", "color"), c("#a9dfbf", "#e6b0aa"))
  expect_equal(sum(grepl("/locus_tag=", gb)), 2)
  expect_equal(sum(grepl("/ApEinfo_revcolor=", gb)), 2)
  expect_equal(sum(grepl("/ApEinfo_graphicformat=", gb)), 2)
})

test_that("piRNA features are written on the complement strand", {
  gb <- export_genbank("t", strrep("ACGT", 10), c(1, 5), c(3, 24), c("#111", "#222"),
                       c("a", "b"), c("misc_feature", "piRNA"))
  p <- parse_genbank(gb)
  expect_equal(vapply(p$features, `[[`, TRUE, "complement"), c(FALSE, TRUE))
  expect_true(any(grepl("piRNA +complement\\(5\\.\\.24\\)", gb)))
})

test_that("a file with no features still has an empty FEATURES section", {
  gb <- export_genbank("t", gb_seq())
  feat <- which(startsWith(gb, "FEATURES"))
  expect_equal(gb[feat + 1], "ORIGIN")
})

test_that("features and their qualifiers keep the input order", {
  labels <- c("first", "second", "third")
  gb <- export_genbank("t", strrep("A", 30), c(1, 10, 20), c(5, 15, 25), c("#1", "#2", "#3"),
                       labels, rep("misc_feature", 3))
  expect_equal(parse_genbank(gb)$labels, labels)
})

test_that("mismatched feature vectors give a clear error", {
  expect_error(export_genbank("t", gb_seq(), 1, 3, "#a9dfbf", "Start"), "same length")
  expect_error(export_genbank("t", gb_seq(), c(1, 2), 3, "#a9dfbf", "Start", "misc_feature"), "same length")
})

# ---- the sequence ---------------------------------------------------------------

test_that("the sequence comes back exactly, whatever its length", {
  for (n in c(1, 9, 10, 11, 59, 60, 61, 119, 120, 121, 130, 600, 1234)) {
    s <- paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
    p <- parse_genbank(export_genbank("t", s))
    expect_identical(p$sequence, s, info = paste(n, "bp"))
    expect_match(p$locus, paste0(n, "bp ds-DNA"), info = paste(n, "bp"))
  }
})

test_that("ORIGIN lines hold 60 bases in blocks of 10, numbered by their first base", {
  s <- strrep("ACGTACGTAC", 13)  # 130 bp
  gb <- export_genbank("t", s)
  origin <- gb[(which(gb == "ORIGIN") + 1):(length(gb) - 1)]
  expect_length(origin, 3)
  expect_equal(as.integer(sub("^ *([0-9]+) .*", "\\1", origin)), c(1, 61, 121))
  expect_equal(nchar(gsub("[^A-Z]", "", origin)), c(60, 60, 10))
  expect_match(origin[1], "^        1 ([ACGT]{10} ){5}[ACGT]{10}$")
})

test_that("lowercase introns stay lowercase and the case is preserved", {
  s <- "ATGGCTgtaagtttcagAAATAA"
  expect_identical(parse_genbank(export_genbank("t", s))$sequence, s)
})

test_that("NULL, non-character and empty sequences give an empty result", {
  expect_length(export_genbank("t", NULL), 0)
  expect_length(export_genbank("t", 123), 0)
  expect_length(export_genbank("t", ""), 0)
})

# ---- the file the download buttons produce ---------------------------------------

test_that("the optimized file matches what the viewer shows", {
  cds <- paste0("ATG", strrep("GCTAAG", 30), "CTCGAG", strrep("GCTAAG", 30), "TAA")
  d <- run_server_downloads("7", cds, list(checkIntron = TRUE, intropt = "Synthetic"))
  p <- parse_genbank(d$opt)

  expect_identical(p$sequence, d$viewer_seq)
  expect_match(p$locus, paste0(nchar(p$sequence), "bp ds-DNA"))
  expect_equal(d$opt_name, "Optimized_My gene.gb")
  expect_match(p$locus, "My gene")
})

test_that("every feature in the optimized file lies inside the sequence", {
  cds <- paste0("ATG", strrep("GCTAAG", 30), "CTCGAG", strrep("GCTAAG", 30), "TAA")
  d <- run_server_downloads("7", cds, list(checkIntron = TRUE, intropt = "Synthetic"))
  p <- parse_genbank(d$opt)
  expect_gt(length(p$starts), 0)
  expect_true(all(p$starts >= 1))
  expect_true(all(p$ends <= nchar(p$sequence)))
  expect_true(all(p$starts <= p$ends))
})

test_that("Start, Stop and intron features point at the right bases in the file", {
  cds <- paste0("ATG", strrep("GCTAAG", 30), "CTCGAG", strrep("GCTAAG", 30), "TAA")
  d <- run_server_downloads("7", cds, list(checkIntron = TRUE, intropt = "Synthetic"))
  p <- parse_genbank(d$opt)
  region <- function(i) substr(p$sequence, p$starts[i], p$ends[i])

  expect_equal(region(which(p$labels == "Start")), "ATG")
  expect_true(is_stop_codon(region(which(p$labels == "Stop"))))

  introns <- which(p$labels == "Intron")
  expect_length(introns, 2)
  lib <- safe_load_csv("Synthetic_Introns.csv")$Sequence
  for (i in introns) {
    expect_match(region(i), "^[acgt]+$")
    expect_true(region(i) %in% lib)
  }
})

test_that("the original file holds the input sequence and marks the restriction site", {
  cds <- paste0("ATG", strrep("GCTAAG", 30), "CTCGAG", strrep("GCTAAG", 30), "TAA")
  d <- run_server_downloads("7", cds)
  p <- parse_genbank(d$ori)
  expect_identical(p$sequence, cds)
  expect_equal(d$ori_name, "Original_My gene.gb")
  site <- which(p$labels == "XhoI site")
  expect_length(site, 1)
  expect_equal(substr(p$sequence, p$starts[site], p$ends[site]), "CTCGAG")
  expect_equal(p$keys[site], "RE_site")
  expect_true("COMMENT     Original unmodified sequence" %in% d$ori |
                any(grepl("Original unmodified sequence", d$ori)))
})

test_that("the comments in the file record the settings and the scores", {
  cds <- paste0("ATG", strrep("GCTAAG", 30), "CTCGAG", strrep("GCTAAG", 30), "TAA")
  d <- run_server_downloads("7", cds, list(checkIntron = TRUE, intropt = "Synthetic"))
  p <- parse_genbank(d$opt)
  expect_true(any(grepl("Optimization routine: CAI1", p$comments)))
  expect_true(any(grepl("Restriction sites removed: XhoI", p$comments)))
  expect_true(any(grepl("Introns added: 2 Random Synthetic Introns", p$comments)))
  expect_true(any(grepl("^Input sequence: ATGGCTAAG", p$comments)))
  expect_true(any(grepl("^Codon Adaptation Index [0-9.]+$", p$comments)))
  expect_true(any(grepl("^GC [0-9]+%$", p$comments)))
  expect_true(any(grepl("^GLO Score [0-9,]+$", p$comments)))
})

test_that("tags, promoter and UTR features cover the right bases in the file", {
  cds <- paste0("ATG", strrep("GCTAAG", 30), "TAA")
  tags <- safe_load_csv("Tags.csv")
  prom <- safe_load_csv("Promoters.csv")
  utr5 <- safe_load_csv("5pUTRs.csv")
  d <- run_server_downloads("100", cds, list(
    checkTags = TRUE, n_tags = tags$Name[1], c_tags = tags$Name[4],
    checkPromoters = TRUE, oppromo = "pPVQ",
    checkUTRs = TRUE, p5UTR = "Fire lab synthetic spliced", p3UTR = "None",
    checkfouras = TRUE), enzymes = NULL)
  p <- parse_genbank(d$opt)
  region <- function(label) {
    i <- which(p$labels == label)
    expect_length(i, 1)
    substr(p$sequence, p$starts[i], p$ends[i])
  }
  expect_equal(toupper(region("SV40 NLS")), toupper(trimws(tags$Sequence[1])))
  expect_equal(toupper(region(gsub("<[^>]+>", "", tags$Name[4]))), toupper(trimws(tags$Sequence[4])))
  expect_equal(region("pPVQ"), tolower(prom$Sequence[prom$Name == "pPVQ"]))
  expect_equal(region("Fire lab synthetic spliced"), tolower(utr5$Sequence[utr5$Name == "Fire lab synthetic spliced"]))
  expect_true(grepl("aaaaATG", p$sequence, fixed = TRUE))
})
