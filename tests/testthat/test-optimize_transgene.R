# End-to-end tests of optimize_transgene(): the whole assembly line from input
# checking through tags, recoding, shields, introns, start/UTRs/promoter to the
# final sequence and its coordinates. Random methods are checked over several
# seeds against properties that must always hold, not exact sequences.

BSAI <- "GGTCTC|GAGACC"
methods_all <- c("None", "Ubiq", "Chance", "CAI1", "GLO")
introns_2 <- c("gtaagtttcag", "gtaagtttcag")

run_pipeline <- function(cds = glo_cds(), method = "None", seed = 1, ...) {
  args <- list(dna_seq = cds, opt_method = method, enzymes = NULL, remove_pirnas = FALSE)
  extra <- list(...)
  args[names(extra)] <- extra
  set.seed(seed)
  do.call(optimize_transgene, args)
}

test_that("the fixture is a valid CDS containing a BsaI site", {
  expect_true(check_cds_validity(glo_cds())$valid)
  expect_true(grepl(BSAI, glo_cds()))
})

# ---- every method keeps the biology intact ----------------------------------

test_that("every method keeps the protein, length, alphabet, start and a stop codon", {
  for (m in methods_all) {
    for (s in 1:2) {
      res <- run_pipeline(method = m, seed = s)
      info <- paste(m, "seed", s)
      expect_equal(protein_of(res$optimized_cds), protein_of(glo_cds()), info = info)
      expect_equal(nchar(res$optimized_cds), nchar(glo_cds()), info = info)
      expect_match(res$optimized_cds, "^[ACGT]+$", info = info)
      expect_equal(substr(res$optimized_cds, 1, 3), "ATG", info = info)
      expect_true(is_stop_codon(last_codon(res$optimized_cds)), info = info)
      expect_equal(nchar(res$optimized_cds) %% 3, 0, info = info)
    }
  }
})

test_that("with no optimization and no shields the CDS is returned unchanged", {
  res <- run_pipeline(method = "None")
  expect_identical(res$optimized_cds, glo_cds())
  expect_identical(res$final_seq, glo_cds())
})

test_that("every method removes the requested restriction site", {
  for (m in methods_all) {
    for (s in 1:2) {
      res <- run_pipeline(method = m, seed = s, enzymes = "BsaI")
      expect_false(grepl(BSAI, res$optimized_cds), info = paste(m, "seed", s))
      expect_equal(protein_of(res$optimized_cds), protein_of(glo_cds()), info = paste(m, "seed", s))
    }
  }
})

test_that("several enzymes can be removed in one run", {
  cds <- "ATGCTCGAGAAAGGTCTCGGCGCCTAA"  # XhoI, BsaI and BanI sites
  res <- run_pipeline(cds, method = "Ubiq", enzymes = c("XhoI", "BsaI", "BanI"))
  expect_false(grepl("CTCGAG|GGTCTC|GAGACC|GG[CT][AG]CC", res$optimized_cds))
  expect_equal(protein_of(res$optimized_cds), protein_of(cds))
})

test_that("the piRNA shield works inside the pipeline for every method", {
  cds <- "ATGGATACCTCAACAAAAAACACATAA"  # contains a real piRNA target
  for (m in methods_all) {
    res <- run_pipeline(cds, method = m, remove_pirnas = TRUE, target_min_hamming = 3)
    expect_equal(protein_of(res$optimized_cds), protein_of(cds), info = m)
    expect_length(res$initial_pirna_dists, length(WormTransgeneBuilder::pirna_sequences))
    expect_length(res$final_pirna_dists, length(WormTransgeneBuilder::pirna_sequences))
    expect_equal(min(res$initial_pirna_dists, na.rm = TRUE), 0, info = m)
    expect_true(all(res$final_pirna_dists >= 3, na.rm = TRUE), info = m)
  }
})

test_that("without the piRNA shield the final distances are still reported", {
  res <- run_pipeline(method = "Ubiq", remove_pirnas = FALSE)
  expect_length(res$final_pirna_dists, length(WormTransgeneBuilder::pirna_sequences))
})

test_that("GLO does not lower the GLO score when no shield changes the sequence", {
  for (s in 1:3) {
    res <- run_pipeline(method = "GLO", seed = s)
    expect_gte(calc_glo_score(res$optimized_cds), calc_glo_score(glo_cds()))
  }
})

# ---- repeatability ----------------------------------------------------------

test_that("High expression gives the same answer whatever the seed", {
  a <- run_pipeline(method = "CAI1", seed = 1, enzymes = "BsaI")
  b <- run_pipeline(method = "CAI1", seed = 2, enzymes = "BsaI")
  expect_identical(a$final_seq, b$final_seq)
})

test_that("random methods are reproducible for a fixed seed and differ between seeds", {
  a <- run_pipeline(method = "Ubiq", seed = 4)
  b <- run_pipeline(method = "Ubiq", seed = 4)
  expect_identical(a$final_seq, b$final_seq)
  outs <- vapply(1:5, function(s) run_pipeline(method = "Ubiq", seed = s)$final_seq, character(1))
  expect_gt(length(unique(outs)), 1)
})

# ---- assembly and coordinates ------------------------------------------------

test_that("promoter, UTRs and the consensus start are assembled in order", {
  res <- run_pipeline(method = "None", add_consensus_start = TRUE,
                      promoter_seq = "ACGTacgt", utr5_seq = "tttt", utr3_seq = "CCCC")
  expect_true(startsWith(res$final_seq, "acgtacgttttt"))          # promoter + 5' UTR
  expect_true(startsWith(res$final_seq, "acgtacgtttttaaaaATG"))  # ... then the consensus start
  expect_true(endsWith(res$final_seq, "TAAcccc"))
})

test_that("the CDS coordinates point at the start and stop codons", {
  res <- run_pipeline(method = "None", add_consensus_start = TRUE, introns = introns_2,
                      promoter_seq = "acgtacgt", utr5_seq = "tttt", utr3_seq = "cccc")
  expect_equal(res$cds_start, nchar("acgtacgt") + nchar("tttt") + nchar("aaaa") + 1)
  expect_equal(toupper(substr(res$final_seq, res$cds_start, res$cds_start + 2)), "ATG")
  expect_true(is_stop_codon(toupper(substr(res$final_seq, res$cds_end - 2, res$cds_end))))
  expect_equal(nchar(res$final_seq), res$cds_end + nchar("cccc"))
})

test_that("the consensus start replaces the ATG rather than adding a second one", {
  res <- run_pipeline(method = "None", add_consensus_start = TRUE)
  expect_true(startsWith(res$final_seq, "aaaaATG"))
  expect_equal(nchar(res$final_seq), nchar(glo_cds()) + 4)
  expect_equal(res$cds_start, 5)
})

test_that("tag coordinates point at the tags, even with introns inserted", {
  tag_n <- "CCAAAGAAGAAGCGTAAG"
  tag_c <- "GACTACAAGGAC"
  res <- run_pipeline(method = "None", introns = introns_2, tag_n_seqs = tag_n, tag_c_seqs = tag_c,
                      add_consensus_start = TRUE, promoter_seq = "acgtacgt", utr5_seq = "tttt")
  expect_equal(toupper(substr(res$final_seq, res$tag_n_starts, res$tag_n_ends)), tag_n)
  expect_equal(toupper(substr(res$final_seq, res$tag_c_starts, res$tag_c_ends)), tag_c)
})

test_that("two tags on one end are placed one after the other", {
  tags <- c("CCAAAGAAGAAGCGTAAG", "GACTACAAGGAC")
  res <- run_pipeline(method = "None", tag_n_seqs = tags)
  expect_equal(toupper(substring(res$final_seq, res$tag_n_starts, res$tag_n_ends)), tags)
  expect_equal(res$tag_n_starts[2], res$tag_n_ends[1] + 1)
})

test_that("tags stay in frame and keep their amino acids under every method", {
  tag_n <- "CCAAAGAAGAAGCGTAAG"
  for (m in methods_all) {
    res <- run_pipeline(method = m, tag_n_seqs = tag_n)
    expect_equal(nchar(res$optimized_cds) %% 3, 0, info = m)
    expect_equal(protein_of(substr(res$final_seq, res$tag_n_starts, res$tag_n_ends)),
                 protein_of(tag_n), info = m)
  }
})

test_that("every tag in the app can be added at either end without changing the protein", {
  tags <- safe_load_csv("Tags.csv")
  for (i in seq_len(nrow(tags))) {
    tag <- trimws(tags$Sequence[i])
    name <- gsub("<[^>]+>", "", tags$Name[i])
    n_end <- run_pipeline(method = "None", tag_n_seqs = tag)
    c_end <- run_pipeline(method = "None", tag_c_seqs = tag)
    for (res in list(n_end, c_end)) {
      prot <- protein_of(res$optimized_cds)
      expect_equal(nchar(res$optimized_cds) %% 3, 0, info = name)
      # Exactly one stop codon, and it is the last residue.
      expect_equal(lengths(regmatches(prot, gregexpr("\\*", prot))), 1, info = name)
      expect_true(endsWith(prot, "*"), info = name)
      # The protein is the original protein plus the tag's residues.
      expect_equal(nchar(prot), nchar(protein_of(glo_cds())) + nchar(tag) / 3, info = name)
    }
  }
})

test_that("introns are inserted in lowercase and can be removed to recover the CDS", {
  for (spacing in c("early", "equidistant")) {
    res <- run_pipeline(method = "None", introns = introns_2, spacing_type = spacing)
    expect_identical(toupper(gsub("[acgt]", "", res$final_seq)), toupper(glo_cds()), info = spacing)
    expect_equal(nchar(res$final_seq), nchar(glo_cds()) + 2 * nchar(introns_2[1]), info = spacing)
    expect_equal(res$intron_lengths, c(11, 11), info = spacing)
    expect_length(res$intron_insertions, 2)
  }
})

test_that("intron input is cleaned: case, non-DNA characters, empty and NA entries", {
  res <- run_pipeline(method = "None", introns = c("GTAAGtttcag", NA, "", "gt-aa xg"))
  expect_equal(res$intron_lengths, c(11, 5))  # "gt-aa xg" cleans to "gtaag"
  expect_identical(gsub("[acgt]", "", res$final_seq), glo_cds())
})

test_that("an empty intron list adds nothing", {
  res <- run_pipeline(method = "None", introns = character(0))
  expect_identical(res$final_seq, glo_cds())
  expect_length(res$intron_insertions, 0)
})

# ---- input handling ------------------------------------------------------------

test_that("original_seq is returned exactly as supplied", {
  messy <- "atg gct\naaa taa"
  res <- run_pipeline(messy, method = "None")
  expect_identical(res$original_seq, messy)
  expect_identical(res$optimized_cds, "ATGGCTAAATAA")
})

test_that("invalid DNA stops with an explanatory error", {
  expect_error(optimize_transgene("GCTAAATAA"), "Invalid CDS.*start codon")
  expect_error(optimize_transgene("ATGAAA"), "Invalid CDS.*stop codon")
  expect_error(optimize_transgene("ATGTAAAAATAA"), "Invalid CDS.*Internal stop")
  expect_error(optimize_transgene("ATGAAAATAA"), "Invalid CDS.*multiple of 3")
})

test_that("protein input is back-translated, given a stop codon and recoded", {
  prot <- "MAKGLSEPTVNQIYW"
  res <- run_pipeline(prot, method = "Ubiq", is_protein = TRUE)
  expect_equal(protein_of(res$optimized_cds), paste0(prot, "*"))
  expect_true(is_stop_codon(last_codon(res$optimized_cds)))
  expect_null(res$initial_pirna_dists)
  expect_match(res$log[1], "Protein")
})

test_that("protein input already ending in a stop does not get a second one", {
  res <- run_pipeline("MAKGL*", method = "None", is_protein = TRUE)
  expect_equal(protein_of(res$optimized_cds), "MAKGL*")
})

test_that("unknown residues in a protein become NNN", {
  res <- run_pipeline("MAXK", method = "None", is_protein = TRUE, enzymes = NULL)
  expect_true(grepl("NNN", res$optimized_cds))
})

test_that("the result has all the documented parts", {
  res <- run_pipeline(method = "None")
  expect_true(all(c("original_seq", "optimized_cds", "final_seq", "log", "initial_pirna_dists",
                    "final_pirna_dists", "cds_start", "cds_end", "tag_n_starts", "tag_n_ends",
                    "tag_c_starts", "tag_c_ends", "intron_insertions", "intron_lengths") %in% names(res)))
  expect_match(res$log[1], "DNA")
})

# ---- ribosome-binding-site optimization ----------------------------------------

test_that("RBS optimization adds the consensus start and keeps the protein", {
  skip_if(Sys.which("RNAfold") == "", "RNAfold (ViennaRNA) is not on the PATH")
  for (m in c("Ubiq", "CAI1", "GLO", "None")) {
    res <- run_pipeline(method = m, rbs_opt = TRUE)
    expect_true(startsWith(res$final_seq, "aaaaATG"), info = m)
    expect_equal(res$cds_start, 5, info = m)
    expect_equal(protein_of(res$optimized_cds), protein_of(glo_cds()), info = m)
    expect_equal(nchar(res$optimized_cds), nchar(glo_cds()), info = m)
  }
})

test_that("RBS optimization only recodes the first 39 bp when the method is None", {
  skip_if(Sys.which("RNAfold") == "", "RNAfold (ViennaRNA) is not on the PATH")
  res <- run_pipeline(method = "None", rbs_opt = TRUE)
  expect_identical(substr(res$optimized_cds, 40, nchar(res$optimized_cds)),
                   substr(glo_cds(), 40, nchar(glo_cds())))
})

test_that("RBS optimization is skipped for a CDS shorter than 39 bp, but the start is still added", {
  res <- run_pipeline("ATGGCTAAAGGTTAA", method = "None", rbs_opt = TRUE)
  expect_true(startsWith(res$final_seq, "aaaaATG"))
  expect_identical(res$optimized_cds, "ATGGCTAAAGGTTAA")
})

# ---- everything at once ----------------------------------------------------------

test_that("all options together produce a consistent construct", {
  skip_if(Sys.which("RNAfold") == "", "RNAfold (ViennaRNA) is not on the PATH")
  tag_n <- "CCAAAGAAGAAGCGTAAG"
  res <- run_pipeline(method = "Ubiq", seed = 3, enzymes = "BsaI", remove_pirnas = TRUE,
                      target_min_hamming = 2, introns = introns_2, rbs_opt = TRUE,
                      promoter_seq = "acgtacgt", utr5_seq = "tttt", utr3_seq = "cccc",
                      tag_n_seqs = tag_n)
  expect_false(grepl(BSAI, res$optimized_cds))
  expect_equal(protein_of(res$optimized_cds),
               protein_of(paste0(substr(glo_cds(), 1, 3), tag_n, substr(glo_cds(), 4, nchar(glo_cds())))))
  expect_equal(toupper(substr(res$final_seq, res$cds_start, res$cds_start + 2)), "ATG")
  expect_equal(protein_of(substr(res$final_seq, res$tag_n_starts, res$tag_n_ends)), protein_of(tag_n))
  expect_true(startsWith(res$final_seq, "acgtacgtttttaaaaATG"))
  expect_true(endsWith(res$final_seq, "cccc"))
})

# ---- "High expression" with the ribosome-binding-site step -----------------------------------------------
# The start is chosen among versions drawn with the graded weights (never the input codons), and the
# rest of the gene is recoded to the best codons.

test_that("High expression + RBS always recodes the start, keeps the protein, and the rest is all best codons", {
  skip_if(Sys.which("RNAfold") == "", "RNAfold (ViennaRNA) is not on the PATH")
  worst <- paste0("ATG", "AAA", "AAT", "GCA", "CTA", "GAA", "ACA", "TCA", "CCA", "GTA", "GGA", "ATA", "CAA")
  cds <- paste0(worst, strrep("GCAAAACTAGAA", 10), "TAA")
  for (s in 1:4) {
    set.seed(s)
    res <- optimize_transgene(cds, opt_method = "CAI1", enzymes = NULL, remove_pirnas = FALSE, rbs_opt = TRUE)
    out <- res$optimized_cds
    expect_identical(protein_of(out), protein_of(cds), info = s)
    expect_false(identical(substr(out, 1, 39), worst), info = s)
    expect_identical(substr(out, 40, nchar(out)), probabilistic_recode(substr(cds, 40, nchar(cds)), "CAI1"), info = s)
  }
})

test_that("High expression + RBS without RNAfold still recodes the start (deterministically)", {
  testthat::local_mocked_bindings(find_rnafold = function(...) "")
  worst <- paste0("ATG", "AAA", "AAT", "GCA", "CTA", "GAA", "ACA", "TCA", "CCA", "GTA", "GGA", "ATA", "CAA")
  cds <- paste0(worst, strrep("GCAAAACTAGAA", 10), "TAA")
  res <- suppressWarnings(optimize_transgene(cds, opt_method = "CAI1", enzymes = NULL, remove_pirnas = FALSE, rbs_opt = TRUE))
  expect_identical(res$optimized_cds, probabilistic_recode(cds, "CAI1"))
})

test_that("in the app, High expression + RBS is repeatable (a fixed seed is used)", {
  skip_if(Sys.which("RNAfold") == "", "RNAfold (ViennaRNA) is not on the PATH")
  cds <- paste0("ATG", strrep("GCAAAACTAGAA", 12), "TAA")
  run <- function() {
    html <- NULL
    shiny::testServer(app_server, {
      session$setInputs(
        intypeinput = 1, nameinput = "t", seqDNA = cds, seqPROT = "", selectCAI = "7",
        checkEnzySites = FALSE, Oenzymes = NULL, checkPirna = FALSE, selectPiMM = 3, checkTags = FALSE,
        checkIntron = FALSE, intropt = "rps-0", num_synth_introns = 1, checkUTRs = FALSE,
        checkPromoters = FALSE, checkboxRibo = TRUE, checkfouras = FALSE, checkTwisty = FALSE,
        checkintframe = TRUE, intdistop = "1")
      session$setInputs(actionSeq = 1)
      html <<- as.character(output$opt_sequence_viewer$html)
    })
    gsub("\n", "", sub(".*new Sequence\\('([^']*)'\\).*", "\\1", gsub("\n", "", html)))
  }
  expect_identical(run(), run())
})
