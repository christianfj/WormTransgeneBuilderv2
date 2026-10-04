# Shared helpers for the test suite. testthat loads every helper-*.R file
# automatically before running the tests.

# A short, valid CDS (starts with ATG, ends with a stop, no internal stops).
# Built from a list of codons so it is easy to read and edit.
fixture_cds <- function() {
  codons <- c("ATG", "GCT", "AAA", "GGT", "CTT", "TCT", "GAA", "CCT",
              "ACT", "GTT", "AAT", "CAA", "ATT", "TAT", "TGG", "TAA")
  paste(codons, collapse = "")
}

# The standard genetic code, written out here on purpose so that the tests do not rely on the
# app's own codon table (if that table assigned a codon to the wrong amino acid, a test that used
# it to translate would never notice). Codons are in the usual T, C, A, G order.
standard_genetic_code <- local({
  bases <- c("T", "C", "A", "G")
  codons <- paste0(rep(bases, each = 16), rep(rep(bases, each = 4), 4), rep(bases, 16))
  amino_acids <- strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "")[[1]]
  setNames(amino_acids, codons)
})

# Translate DNA to a protein string, used to check that recoding never changes
# the amino acids. Stop codons appear as "*", a codon with an unknown base as "X",
# and a last incomplete codon is ignored.
protein_of <- function(dna) {
  dna <- toupper(dna)
  n <- nchar(dna) %/% 3
  if (n == 0) return("")
  codons <- substring(dna, seq(1, by = 3, length.out = n), seq(3, by = 3, length.out = n))
  amino_acids <- unname(standard_genetic_code[codons])
  amino_acids[is.na(amino_acids)] <- "X"
  paste(amino_acids, collapse = "")
}

# Seeds used to repeat randomised tests. The functions are random, so instead
# of checking one exact output we check that a property holds for every seed.
test_seeds <- 1:10

is_stop_codon <- function(codon) codon %in% c("TAA", "TAG", "TGA")

last_codon <- function(dna) substr(dna, nchar(dna) - 2, nchar(dna))

# Run the real Shiny server once with simulated inputs and return the optimized
# sequence exactly as the viewer receives it (introns are lowercase). Used to test
# behaviour that depends on the server code, such as which random choices it makes.
run_server_once <- function(method, seq, enzymes = "XhoI", introns = "Synthetic",
                            n_introns = 2, pirna = FALSE) {
  html <- NULL
  shiny::testServer(app_server, {
    session$setInputs(
      intypeinput = 1, nameinput = "t", seqDNA = seq, seqPROT = "", selectCAI = method,
      checkEnzySites = !is.null(enzymes), Oenzymes = enzymes, checkPirna = pirna, selectPiMM = 3,
      checkTags = FALSE, checkIntron = !is.null(introns), intropt = introns,
      num_synth_introns = n_introns, checkUTRs = FALSE, checkPromoters = FALSE,
      checkboxRibo = FALSE, checkfouras = FALSE, checkTwisty = FALSE,
      checkintframe = TRUE, intdistop = "early"
    )
    session$setInputs(actionSeq = 1)
    html <<- as.character(output$opt_sequence_viewer$html)
  })
  gsub("\n", "", sub(".*new Sequence\\('([^']*)'\\).*", "\\1", gsub("\n", "", html)))
}

# The lowercase runs of a final sequence, i.e. the inserted introns.
introns_in <- function(final_seq) regmatches(final_seq, gregexpr("[acgt]+", final_seq))[[1]]

# A 36-codon CDS (ATG ... TAA) used by the GLO, RBS and pipeline tests. It contains a
# BsaI site (GGTCTC, in frame at codons 3-4) and has no internal stop codon.
glo_cds <- function() {
  paste0("ATGGCTAAAGGTCTCTCTGAACCTACTGTTAATCAAATTTATTGGCGTAAAGATGGTGAAGCTTTT",
         "CCAAATTTTCTGGGACTCTGTTGGAAAGCCGAAGAATGTTAA")
}

# Run the real Shiny server once and return what the two GenBank download buttons
# produce, together with the sequence the viewer shows for the same run.
run_server_downloads <- function(method, seq, extra = list(), enzymes = "XhoI") {
  out <- NULL
  shiny::testServer(app_server, {
    inputs <- c(list(
      intypeinput = 1, nameinput = "My gene", seqDNA = seq, seqPROT = "", selectCAI = method,
      checkEnzySites = !is.null(enzymes), Oenzymes = enzymes, checkPirna = FALSE, selectPiMM = 3,
      checkTags = FALSE, checkIntron = FALSE, intropt = "Synthetic", num_synth_introns = 2,
      checkUTRs = FALSE, checkPromoters = FALSE, checkboxRibo = FALSE, checkfouras = FALSE,
      checkTwisty = FALSE, checkintframe = TRUE, intdistop = "early"
    ), extra)
    inputs <- inputs[!duplicated(names(inputs), fromLast = TRUE)]
    do.call(session$setInputs, inputs)
    session$setInputs(actionSeq = 1)
    html <- gsub("\n", "", as.character(output$opt_sequence_viewer$html))
    opt_path <- output$downloadOptGenbank
    ori_path <- output$downloadOriGenbank
    out <<- list(viewer_seq = sub(".*new Sequence\\('([^']*)'\\).*", "\\1", html),
                 opt = readLines(opt_path), ori = readLines(ori_path),
                 opt_name = basename(opt_path), ori_name = basename(ori_path))
  })
  out
}

# Read a GenBank file (as produced by export_genbank) back into its parts.
parse_genbank <- function(lines) {
  origin <- which(lines == "ORIGIN")
  end <- which(lines == "//")
  seq_lines <- lines[(origin + 1):(end - 1)]
  sequence <- paste(gsub("[^A-Za-z]", "", seq_lines), collapse = "")

  feat_start <- which(startsWith(lines, "FEATURES"))
  feat_lines <- lines[(feat_start + 1):(origin - 1)]
  key_idx <- which(grepl("^     [^ ]", feat_lines))
  features <- lapply(key_idx, function(i) {
    m <- regmatches(feat_lines[i], regexec("^     ([^ ]+) +(complement\\()?([0-9]+)\\.\\.([0-9]+)\\)?$", feat_lines[i]))[[1]]
    q <- feat_lines[(i + 1):min(length(feat_lines), i + 5)]
    list(key = m[2], complement = nzchar(m[3]), start = as.integer(m[4]), end = as.integer(m[5]),
         label = sub('.*/ApEinfo_label="([^"]*)".*', "\\1", grep("ApEinfo_label", q, value = TRUE)),
         color = sub('.*/ApEinfo_fwdcolor="([^"]*)".*', "\\1", grep("ApEinfo_fwdcolor", q, value = TRUE)))
  })
  list(
    locus = lines[startsWith(lines, "LOCUS")],
    comments = sub("^COMMENT +", "", lines[startsWith(lines, "COMMENT")]),
    sequence = sequence,
    features = features,
    keys = vapply(features, `[[`, "", "key"),
    starts = vapply(features, `[[`, 0L, "start"),
    ends = vapply(features, `[[`, 0L, "end"),
    labels = vapply(features, `[[`, "", "label")
  )
}

# The 21 enzymes the enzyme table held before it became New England Biolabs' full list (234
# enzymes). Reference results captured back then are tied to these 21, so the tests that use
# them name them explicitly instead of taking "everything in the table".
original_21_enzymes <- c("ApaI", "ApaLI", "BamHI", "BanI", "BsaI", "BsmBI", "BsoBI", "DraI", "EcoRI",
                         "EcoRV", "HindIII", "KpnI", "NcoI", "NdeI", "NheI", "PstI", "PvuII", "SacII",
                         "SapI", "XbaI", "XhoI")

# The enzymes commonly avoided in a transgene (the set the laboratory uses in ApE). Removing every
# enzyme of the full list makes no biological sense, so tests use this set. BanIII is left out: it
# is not a New England Biolabs enzyme, and its site (ATCGAT) is the same as that of ClaI and BspDI.
common_enzymes <- c("ApaI", "ApaLI", "AvaI", "AvaII", "BamHI", "BanI", "BanII", "BclI", "BsaI", "BsmBI",
                    "DraI", "EcoRI", "EcoRV", "HindIII", "KpnI", "NciI", "NdeI", "PmlI", "PstI", "PvuII",
                    "SacII", "SmaI", "XbaI", "XhoI")
