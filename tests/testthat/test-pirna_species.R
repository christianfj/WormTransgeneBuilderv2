# The piRNA databases of the supported species are independent files; the app screens against the species the user selects
# (default: C. elegans with its type 2 piRNAs). See docs/PIRNA_DATABASE.md and docs/PIRNA_OTHER_SPECIES.md.

all_species <- c("elegans", "briggsae", "remanei", "brenneri")

# ---- the files that ship with the package ---------------------------------------------------------------------

test_that("every species database in inst/extdata is clean: 21-mers starting with T and their 20-base targets", {
  for (id in all_species) {
    d <- read_pirna_species(id)
    expect_false(anyNA(d[, c("Name", "Target_20mer", "Full_21mer")]), info = id)
    expect_true(all(grepl("^[ACGT]{21}$", d$Full_21mer)), info = id)
    expect_identical(d$Target_20mer, reverse_complement(substring(d$Full_21mer, 2, 21)), info = id)
    expect_equal(anyDuplicated(d$Name), 0, info = id)
    expect_equal(anyDuplicated(d$Target_20mer), 0, info = id)
    expect_true(all(c("Species", "Type", "Name", "Target_20mer", "Full_21mer", "Abundance", "Enrichment") %in% names(d)), info = id)
  }
})

test_that("the databases have the expected sizes", {
  n <- vapply(all_species, function(id) nrow(read_pirna_species(id)), 1L)
  expect_identical(unname(n), c(18207L, 25883L, 11801L, 17991L))
  expect_equal(sum(read_pirna_species("elegans")$Type == "type 2"), 2485)
})

test_that("pirna_sequences is the Target_20mer column of the C. elegans database, in the same order", {
  expect_identical(WormTransgeneBuilder::pirna_sequences, read_pirna_species("elegans")$Target_20mer)
})

test_that("the databases of different species are independent: practically no sequence is shared", {
  targets <- lapply(all_species, function(id) read_pirna_species(id)$Target_20mer)
  shared <- sum(duplicated(unlist(targets)))
  expect_lt(shared, 10)   # three coincidences among 74,000 targets at the time of writing
  expect_equal(nrow(select_pirna_database(all_species)), sum(lengths(targets)))   # nothing is merged or dropped
})

# ---- selecting ------------------------------------------------------------------------------------------------

test_that("the default selection is C. elegans with its type 2 piRNAs", {
  d <- select_pirna_database()
  expect_equal(nrow(d), 18207)
  expect_true(all(d$Species == "C. elegans"))
  expect_setequal(unique(d$Type), c("type 1", "type 2", "PRG-1 IP new"))
  expect_identical(d$Target_20mer, WormTransgeneBuilder::pirna_sequences)
})

test_that("the type 2 option removes exactly the type 2 entries of C. elegans and nothing else", {
  expect_equal(nrow(select_pirna_database("elegans", include_type2 = FALSE)), 18207 - 2485)
  expect_false("type 2" %in% select_pirna_database("elegans", include_type2 = FALSE)$Type)
  both <- select_pirna_database(c("elegans", "briggsae"), include_type2 = FALSE)
  expect_equal(nrow(both), 18207 - 2485 + 25883)   # the other species have no type 2 entries
})

test_that("selecting several species concatenates their databases, C. elegans first and then in a fixed order", {
  d <- select_pirna_database(c("brenneri", "briggsae", "elegans"))   # the order of the request does not matter
  expect_identical(unique(d$Species), c("C. elegans", "C. briggsae", "C. brenneri"))
  expect_equal(nrow(d), 18207 + 25883 + 17991)
  expect_identical(d$Target_20mer[seq_len(18207)], WormTransgeneBuilder::pirna_sequences)
})

test_that("each species is labelled with its confidence", {
  d <- select_pirna_database(all_species)
  expect_match(d$Type[d$Species == "C. briggsae"], "high confidence")
  expect_match(d$Type[d$Species %in% c("C. remanei", "C. brenneri")], "low confidence")
})

test_that("selecting no species, or an unknown one, is an error", {
  expect_error(select_pirna_database(character(0)), "at least one species")
  expect_error(select_pirna_database("elegans_x"), "at least one species")
})

# ---- the shield against the selected species ------------------------------------------------------------------

# A coding sequence that contains one target exactly, in frame (20 bases plus a padding base: 7 codons without a stop).
cds_with_target <- function(target) {
  stops <- c("TAA", "TAG", "TGA")
  for (pad in c("A", "C", "G", "T")) {
    block <- paste0(target, pad)
    codons <- substring(block, seq(1, 19, 3), seq(3, 21, 3))
    if (!any(codons %in% stops)) return(paste0("ATG", block, "GCTAAGTAA"))
  }
  NA_character_
}

test_that("a piRNA of a selected species is shielded, and the distances have one entry per selected piRNA", {
  d <- select_pirna_database(c("elegans", "briggsae"))
  b <- d[d$Species == "C. briggsae", ]
  target <- NA; cds <- NA
  for (t in b$Target_20mer[1:200]) { cds <- cds_with_target(t); if (!is.na(cds)) { target <- t; break } }
  skip_if(is.na(cds), "no usable target")
  idx <- which(d$Target_20mer == target)[1]
  before <- get_pirna_distances(cds, d$Target_20mer)
  expect_length(before, nrow(d))
  expect_equal(before[idx], 0)
  set.seed(1)
  res <- optimize_transgene(cds, opt_method = "Ubiq", enzymes = NULL, remove_pirnas = TRUE, target_min_hamming = 4,
                            pirna_targets = d$Target_20mer)
  expect_length(res$initial_pirna_dists, nrow(d))
  expect_equal(res$initial_pirna_dists[idx], 0)
  expect_identical(protein_of(res$optimized_cds), protein_of(cds))
  expect_gte(res$final_pirna_dists[idx], 4)
})

test_that("a briggsae piRNA is not looked for when only C. elegans is selected", {
  b <- read_pirna_species("briggsae")
  cel <- WormTransgeneBuilder::pirna_sequences
  target <- NA; cds <- NA
  for (t in b$Target_20mer[1:500]) { cds <- cds_with_target(t); if (!is.na(cds)) { target <- t; break } }
  skip_if(is.na(cds), "no usable target")
  # With C. elegans only, the design is left alone unless the same window happens to be close to a C. elegans piRNA.
  d_cel <- get_pirna_distances(cds, cel)
  d_cbr <- get_pirna_distances(cds, b$Target_20mer)
  expect_equal(min(d_cbr), 0)
  expect_gt(min(d_cel), 0)   # no exact C. elegans match: the briggsae piRNA is a different sequence
})

test_that("distances against a subset equal the corresponding entries of the distances against the full set", {
  d <- select_pirna_database(c("elegans", "remanei"))
  cds <- paste0("ATG", strrep("GCTAAGGATCCT", 12), "TAA")
  full <- get_pirna_distances(cds, d$Target_20mer)
  keep <- which(d$Type != "type 2")
  expect_identical(get_pirna_distances(cds, d$Target_20mer[keep]), full[keep])
})

# ---- the app --------------------------------------------------------------------------------------------------

run_pirna_app <- function(extra = list(), seq = "ATGGATACCTCAACAAAAAACACATAA") {
  # "21ur-1 (type 2)" has the target GATACCTCAACAAAAAACAC, which this sequence contains exactly.
  table_html <- NULL
  out <- NULL
  shiny::testServer(app_server, {
    inputs <- c(list(
      intypeinput = 1, nameinput = "t", seqDNA = seq, seqPROT = "", selectCAI = "100",
      checkEnzySites = FALSE, Oenzymes = NULL, checkPirna = TRUE, selectPiMM = 4, checkTags = FALSE,
      checkIntron = FALSE, intropt = "rpsa-1", num_synth_introns = 1, checkUTRs = FALSE, checkPromoters = FALSE,
      checkboxRibo = FALSE, checkfouras = FALSE, checkTwisty = FALSE, checkintframe = TRUE, intdistop = "1"), extra)
    inputs <- inputs[!duplicated(names(inputs), fromLast = TRUE)]
    do.call(session$setInputs, inputs)
    session$setInputs(actionSeq = 1)
    table_html <<- tryCatch(as.character(output$pirnaTable), error = function(e) NA_character_)
    out <<- tryCatch(readLines(output$downloadOptGenbank), error = function(e) NA_character_)
  })
  list(table = table_html, genbank = out)
}

test_that("by default the app avoids C. elegans piRNAs including type 2 and says so in the run log", {
  r <- run_pirna_app(list(pirnaSpecies = "elegans", pirnaType2 = TRUE))
  expect_match(r$table, "21ur-1 \\(type 2\\)", fixed = FALSE)
  expect_match(r$table, "<th[^>]*> Type </th>")
  expect_false(grepl("<th[^>]*> Species </th>", r$table))      # a single species: no Species column
  expect_true(any(grepl("piRNA species avoided: C. elegans", r$genbank, fixed = TRUE)))
  expect_true(any(grepl("Include C. elegans type 2 piRNAs: Yes", r$genbank, fixed = TRUE)))
})

test_that("with the type 2 box unticked the type 2 piRNA is not flagged", {
  r <- run_pirna_app(list(pirnaSpecies = "elegans", pirnaType2 = FALSE))
  expect_false(grepl("21ur-1 (type 2)", r$table, fixed = TRUE))
  expect_true(any(grepl("Include C. elegans type 2 piRNAs: No", r$genbank, fixed = TRUE)))
})

test_that("with several species selected the table has a Species column and the log lists them", {
  r <- run_pirna_app(list(pirnaSpecies = c("elegans", "briggsae", "remanei", "brenneri"), pirnaType2 = TRUE))
  expect_match(r$table, "<th[^>]*> Species </th>")
  expect_match(r$table, "C. elegans")
  expect_true(any(grepl("piRNA species avoided: C. elegans, C. briggsae, C. remanei, C. brenneri", r$genbank, fixed = TRUE)))
})

test_that("a species other than C. elegans is not screened unless it is selected", {
  b <- read_pirna_species("briggsae")
  cds <- NA
  for (t in b$Target_20mer[1:500]) { cds <- cds_with_target(t); if (!is.na(cds)) { target_name <- b$Name[b$Target_20mer == t]; break } }
  skip_if(is.na(cds), "no usable target")
  only_cel <- run_pirna_app(list(pirnaSpecies = "elegans"), seq = cds)
  with_cbr <- run_pirna_app(list(pirnaSpecies = c("elegans", "briggsae")), seq = cds)
  expect_false(grepl(target_name, only_cel$table, fixed = TRUE))
  expect_true(grepl(target_name, with_cbr$table, fixed = TRUE))
})

test_that("ticking 'Minimize piRNA homology' but no species gives an error message and no result", {
  r <- run_pirna_app(list(pirnaSpecies = NULL))
  expect_true(is.na(r$table))
})

test_that("the species box is in the interface: four species, only C. elegans ticked, plus the type 2 box ticked", {
  html <- as.character(app_ui(NULL))
  expect_equal(length(gregexpr('name="pirnaSpecies"', html, fixed = TRUE)[[1]]), 4)
  for (v in c("elegans", "briggsae", "remanei", "brenneri")) expect_true(grepl(sprintf('value="%s"', v), html, fixed = TRUE), info = v)
  expect_equal(length(gregexpr('name="pirnaSpecies"[^>]*checked', html)[[1]]), 1)
  expect_true(grepl('id="pirnaType2"[^>]*checked', html))
})
