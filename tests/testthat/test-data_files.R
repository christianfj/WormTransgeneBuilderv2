# "Data contract" tests. The app reads its lookup tables from two places: CSV
# files in inst/extdata (introns, tags, UTRs, promoters, ...) and pre-built
# .rda files in data/ (codon frequencies, piRNAs, enzymes). These tests check
# that the two stay in step and that the CSVs contain what the code expects.

extdata_csv <- function(name) {
  read.csv(system.file("extdata", name, package = "WormTransgeneBuilder"),
           stringsAsFactors = FALSE)
}

usage_columns <- c("Ubiq", "CAI1", "BringMansWeights", "Chance", "aCDS", "HighRCSU")

# ---- 1. Codon table --------------------------------------------------------

test_that("codon_freqs.rda matches Codon_frequencies.csv", {
  csv <- extdata_csv("Codon_frequencies.csv")
  rda <- as.data.frame(WormTransgeneBuilder::codon_freqs, row.names = NULL)
  expect_equal(rda, csv[, names(rda)], ignore_attr = TRUE)
})

test_that("the codon table has 64 unique codons, each with an amino acid", {
  cf <- WormTransgeneBuilder::codon_freqs
  expect_equal(nrow(cf), 64)
  expect_equal(anyDuplicated(toupper(cf$Codon)), 0)
  expect_false(anyNA(cf$Amino))
  expect_match(toupper(cf$Codon), "^[ACGT]{3}$")
})

test_that("stop codons are stored under amino acid 'X', not '*'", {
  # Recorded because it is easy to trip over: code that looks for "*" in this
  # table will find nothing.
  cf <- WormTransgeneBuilder::codon_freqs
  stops <- toupper(cf$Codon)[cf$Amino == "X"]
  expect_setequal(stops, c("TAA", "TAG", "TGA"))
  expect_false("*" %in% cf$Amino)
})

test_that("every usage column has non-negative, non-missing weights", {
  cf <- WormTransgeneBuilder::codon_freqs
  for (col in usage_columns) {
    expect_false(anyNA(cf[[col]]), info = col)
    expect_true(all(cf[[col]] >= 0), info = col)
  }
})

test_that("the usage columns offered in the app give every amino acid a positive weight", {
  # The app offers Ubiq, Chance and CAI1 (GLO uses Ubiq for its shields).
  # BringMansWeights is not offered in the app and gives the stop codons a
  # weight of 0, so it is deliberately not checked here.
  cf <- WormTransgeneBuilder::codon_freqs
  for (col in c("Ubiq", "Chance", "CAI1")) {
    per_aa <- tapply(cf[[col]], cf$Amino, function(x) any(x > 0))
    expect_true(all(per_aa), info = col)
  }
})

# ---- 2. piRNA database -----------------------------------------------------

test_that("pirna_sequences.rda matches Target_20mer in piRNA_Database.csv, in order", {
  # The app finds a piRNA's row in the CSV by matching its sequence against the
  # .rda, and distances are indexed by .rda position, so the order must agree.
  csv <- extdata_csv("piRNA_Database.csv")
  expect_identical(WormTransgeneBuilder::pirna_sequences, csv$Target_20mer)
})

test_that("every piRNA has a sequence, and it is a 20-mer of A/C/G/T", {
  p <- WormTransgeneBuilder::pirna_sequences
  expect_false(anyNA(p))
  expect_match(p, "^[ACGT]{20}$")
})

# 301 entries once had no sequence (NA in both Target_20mer and Full_21mer). They were
# removed on purpose: they are not in another list of piRNAs and not on WormBase, so
# they are not validated. A Hamming distance to a missing sequence is itself missing,
# so they could never have been flagged or removed anyway. This keeps them out.
test_that("no piRNA in the database is missing its name or a sequence", {
  csv <- extdata_csv("piRNA_Database.csv")
  expect_false(anyNA(csv[, c("Name", "Target_20mer", "Full_21mer")]))
  expect_true(all(nzchar(csv$Target_20mer)))
})

test_that("Full_21mer is the reverse complement of Target_20mer with one extra base in front", {
  # A strong check on every row: it would notice a shifted, swapped or damaged entry.
  csv <- extdata_csv("piRNA_Database.csv")
  rc <- reverse_complement(csv$Target_20mer)
  expect_true(all(nchar(csv$Full_21mer) == 21))
  expect_true(all(endsWith(csv$Full_21mer, rc)))
})

test_that("piRNA names and target sequences are each unique", {
  csv <- extdata_csv("piRNA_Database.csv")
  expect_equal(anyDuplicated(csv$Name), 0)
  expect_equal(anyDuplicated(WormTransgeneBuilder::pirna_sequences), 0)
})

test_that("the piRNA table has the columns the app reads", {
  csv <- extdata_csv("piRNA_Database.csv")
  expect_true(all(c("Name", "Target_20mer", "Full_21mer", "Abundance", "Enrichment") %in% names(csv)))
})

test_that("every piRNA can be flagged: get_pirna_distances gives a real distance for all of them", {
  d <- get_pirna_distances(fixture_cds())
  expect_length(d, length(WormTransgeneBuilder::pirna_sequences))
  expect_false(anyNA(d))
})

# ---- 3. Introns (natural and synthetic) ------------------------------------

test_that("natural introns are lowercase ACGT, begin with gt and end with ag", {
  introns <- extdata_csv("Introns.csv")
  seqs <- unlist(introns[, grep("^Sequence_", names(introns))])
  seqs <- seqs[nzchar(seqs)]
  expect_gt(length(seqs), 0)
  expect_match(seqs, "^[acgt]+$")
  expect_match(seqs, "^gt")
  expect_match(seqs, "ag$")
})

test_that("intron set names are unique", {
  expect_equal(anyDuplicated(extdata_csv("Introns.csv")$Name), 0)
})

test_that("synthetic introns are read by the app with a 'Sequence' column", {
  synth <- safe_load_csv("Synthetic_Introns.csv")
  expect_false(is.null(synth))
  expect_true("Sequence" %in% colnames(synth))
  expect_gte(nrow(synth), 3)
})

test_that("synthetic introns are unique lowercase ACGT that begin with gt and end with ag", {
  s <- safe_load_csv("Synthetic_Introns.csv")$Sequence
  expect_match(s, "^[acgt]+$")
  expect_match(s, "^gt")
  expect_match(s, "ag$")
  expect_equal(anyDuplicated(s), 0)
})

# ---- 4. Promoters, UTRs and tags -------------------------------------------

test_that("promoters and UTRs have unique names and valid DNA", {
  for (f in c("Promoters.csv", "5pUTRs.csv", "3pUTRs.csv")) {
    d <- extdata_csv(f)
    if (f == "Promoters.csv") {
      expect_identical(names(d), c("Name", "Gene", "Sequence"))
    } else {
      expect_true(all(c("Name", "UI_HTML", "Sequence") %in% names(d)), info = f)
    }
    expect_equal(anyDuplicated(d$Name), 0, info = f)
    # Rows of Promoters.csv without a sequence are group headings, not promoters.
    real <- d$Sequence[d$Name != "None" & nzchar(d$Sequence)]
    expect_gt(length(real), 0)
    expect_match(real, "^[ACGTacgt]+$", info = f)
  }
})

test_that("the UTR lists offer a 'None' choice", {
  expect_true("None" %in% extdata_csv("5pUTRs.csv")$Name)
  expect_true("None" %in% extdata_csv("3pUTRs.csv")$Name)
})

test_that("no tag contains an in-frame stop codon", {
  # A tag is inserted into a coding sequence, so a stop codon inside it would end
  # the protein early. (The GGSGG linker once ended in TAA; it must not again.)
  tags <- extdata_csv("Tags.csv")
  s <- toupper(trimws(tags$Sequence))
  has_stop <- vapply(s, function(x) {
    any(substring(x, seq(1, nchar(x), 3), seq(3, nchar(x), 3)) %in% c("TAA", "TAG", "TGA"))
  }, logical(1))
  expect_equal(gsub("<[^>]+>", "", tags$Name[has_stop]), character(0))
})

test_that("the GGSGG linker encodes exactly the five residues GGSGG", {
  tags <- extdata_csv("Tags.csv")
  linker <- toupper(trimws(tags$Sequence[grepl("GGSGG", tags$Name)]))
  expect_length(linker, 1)
  expect_equal(protein_of(linker), "GGSGG")
})

test_that("tags are unique, valid DNA, and a multiple of 3 long so they stay in frame", {
  tags <- extdata_csv("Tags.csv")
  expect_true(all(c("Name", "Sequence") %in% names(tags)))
  expect_equal(anyDuplicated(tags$Name), 0)
  s <- trimws(tags$Sequence)  # the file has a space after each comma
  expect_match(s, "^[ACGTacgt]+$")
  expect_true(all(nchar(s) %% 3 == 0))
})

# ---- 5. File encoding ------------------------------------------------------

test_that("no CSV in inst/extdata starts with a byte-order mark (BOM)", {
  # A BOM can turn the first column name into "X.U.FEFF.Sequence" on some R
  # versions, which would silently break the column lookup.
  files <- list.files(system.file("extdata", package = "WormTransgeneBuilder"),
                      pattern = "\\.csv$", full.names = TRUE)
  expect_gt(length(files), 0)
  for (f in files) {
    first3 <- readBin(f, "raw", 3)
    expect_false(identical(first3, as.raw(c(0xef, 0xbb, 0xbf))), info = basename(f))
  }
})

# ---- 6. The datasets that ship with the package --------------------------------

test_that("the package ships exactly the datasets the code uses", {
  # intron_seqs and aa_to_codon_f were removed because nothing used them; the app reads
  # introns from Introns.csv. A dataset added later should be added here on purpose.
  lazy <- asNamespace("WormTransgeneBuilder")[[".__NAMESPACE__."]][["lazydata"]]
  expect_setequal(ls(lazy), c("codon_freqs", "enzy", "glo_median_score", "pirna_sequences"))
})
