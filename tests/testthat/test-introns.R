# 2 introns need a reasonably long CDS. In this one every odd codon from the
# 3rd onwards is AAG (the preferred intron position), the others are GCT.
intron_cds <- function() paste0("ATG", strrep("GCTAAG", 100), "TAA")
intron_a <- "gtaagtttcag"

test_that("removing the lowercase introns gives back the original CDS", {
  for (type in c("early", "equidistant")) {
    res <- insert_introns(intron_cds(), c(intron_a, intron_a), type)
    expect_identical(gsub("[acgt]", "", res$seq), intron_cds(), info = type)
  }
})

test_that("the total length grows by exactly the intron lengths", {
  res <- insert_introns(intron_cds(), c(intron_a, intron_a))
  expect_equal(nchar(res$seq), nchar(intron_cds()) + 2 * nchar(intron_a))
})

test_that("introns are inserted in lowercase, the CDS stays uppercase", {
  res <- insert_introns(intron_cds(), c("GTAAGTTTCAG", "GTAAGTTTCAG"))
  expect_equal(sum(strsplit(res$seq, "")[[1]] %in% letters),
               2 * nchar(intron_a))
})

test_that("insertion points sit on codon boundaries, in increasing order", {
  for (type in c("early", "equidistant")) {
    res <- insert_introns(intron_cds(), c(intron_a, intron_a), type)
    expect_length(res$insertions, 2)
    expect_true(all(res$insertions %% 3 == 0), info = type)
    expect_false(is.unsorted(res$insertions), info = type)
  }
})

test_that("early spacing places introns near codons 50 and 100", {
  res <- insert_introns(intron_cds(), c(intron_a, intron_a), "early")
  expect_true(all(abs(res$insertions / 3 - c(50, 100)) <= 15))
})

test_that("equidistant spacing places introns near 1/3 and 2/3 of the CDS", {
  res <- insert_introns(intron_cds(), c(intron_a, intron_a), "equidistant")
  total_codons <- nchar(intron_cds()) / 3
  targets <- round(total_codons * c(1, 2) / 3)
  expect_true(all(abs(res$insertions / 3 - targets) <= 15))
})

test_that("introns prefer to follow an AAG codon", {
  res <- insert_introns(intron_cds(), c(intron_a, intron_a), "early")
  preceding <- substring(intron_cds(), res$insertions - 2, res$insertions)
  expect_equal(preceding, c("AAG", "AAG"))
})

test_that("an empty intron list is a no-op", {
  res <- insert_introns(intron_cds(), character(0))
  expect_identical(res$seq, intron_cds())
  expect_length(res$insertions, 0)
})

test_that("a single intron is inserted", {
  res <- insert_introns(intron_cds(), intron_a, "early")
  expect_length(res$insertions, 1)
  expect_identical(gsub("[acgt]", "", res$seq), intron_cds())
})

# ---- introns that need not be in the reading frame (in_frame = FALSE) ---------------------------------
# Here AAG never forms a codon: the codons are GCA AGC GCA AGC ..., and AAG straddles the boundary of
# two codons (it ends at base 8, 14, 20, ... which is not a multiple of 3).
straddle_cds <- function() paste0("ATG", strrep("GCAAGC", 100), "TAA")

test_that("with in_frame = FALSE the intron follows an AAG even when AAG is not a codon", {
  for (type in c("early", "equidistant")) {
    res <- insert_introns(straddle_cds(), c(intron_a, intron_a), type, in_frame = FALSE)
    expect_length(res$insertions, 2)
    expect_true(all(res$insertions %% 3 == 2), info = type)
    expect_identical(substring(straddle_cds(), res$insertions - 2, res$insertions), c("AAG", "AAG"), info = type)
  }
})

test_that("with in_frame = TRUE the same sequence gets introns between codons", {
  res <- insert_introns(straddle_cds(), c(intron_a, intron_a), "early", in_frame = TRUE)
  expect_true(all(res$insertions %% 3 == 0))
})

test_that("in_frame = TRUE is the default", {
  expect_identical(insert_introns(straddle_cds(), c(intron_a, intron_a)),
                   insert_introns(straddle_cds(), c(intron_a, intron_a), in_frame = TRUE))
})

test_that("with in_frame = FALSE the CDS is unchanged by the introns and they stay near their targets", {
  for (type in c("early", "equidistant")) {
    res <- insert_introns(intron_cds(), c(intron_a, intron_a), type, in_frame = FALSE)
    expect_identical(gsub("[acgt]", "", res$seq), intron_cds(), info = type)
    expect_equal(nchar(res$seq), nchar(intron_cds()) + 2 * nchar(intron_a))
    expect_false(is.unsorted(res$insertions))
    framed <- insert_introns(intron_cds(), c(intron_a, intron_a), type, in_frame = TRUE)
    expect_true(all(abs(res$insertions - framed$insertions) <= 3 * 15 * 2), info = type)
  }
})

test_that("with in_frame = FALSE no intron is placed inside the start or stop codon region", {
  res <- insert_introns(straddle_cds(), rep(intron_a, 6), "early", in_frame = FALSE)
  expect_true(all(res$insertions >= 15))
  expect_true(all(res$insertions <= nchar(straddle_cds()) - 15))
  expect_true(all(diff(res$insertions) >= 30))
})

test_that("the 'Force in reading frame' tick box in the app decides whether introns follow codons", {
  cds <- straddle_cds()
  position_of_intron <- function(framed) {
    html <- NULL
    shiny::testServer(app_server, {
      session$setInputs(
        intypeinput = 1, nameinput = "t", seqDNA = cds, seqPROT = "", selectCAI = "100",
        checkEnzySites = FALSE, Oenzymes = NULL, checkPirna = FALSE, selectPiMM = 3, checkTags = FALSE,
        checkIntron = TRUE, intropt = "rps-0", num_synth_introns = 1, checkUTRs = FALSE,
        checkPromoters = FALSE, checkboxRibo = FALSE, checkfouras = FALSE, checkTwisty = FALSE,
        checkintframe = framed, intdistop = "1")
      session$setInputs(actionSeq = 1)
      html <<- as.character(output$opt_sequence_viewer$html)
    })
    seq <- gsub("\\n", "", sub(".*new Sequence\\('([^']*)'\\).*", "\\1", gsub("\\n", "", html)))
    first_lower <- regexpr("[acgt]", seq)[[1]]
    first_lower - 1                      # bases before the first intron
  }
  expect_equal(position_of_intron(TRUE) %% 3, 0)
  expect_equal(position_of_intron(FALSE) %% 3, 2)
})
