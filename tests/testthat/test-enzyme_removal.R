# A CDS with a BsaI site (GGTCTC) in frame at codons 3-4.
bsai_cds <- function() "ATGGCTGGTCTCAAAGGTCTTTCTGAACCTACTGTTAATCAAATTTATTGGTAA"

has_bsai <- function(dna) grepl("GGTCTC|GAGACC", toupper(dna))

test_that("the BsaI fixture really contains a site", {
  expect_true(has_bsai(bsai_cds()))
})

test_that("BsaI sites are removed and the protein is unchanged", {
  for (s in test_seeds) {
    set.seed(s)
    out <- remove_enzyme_sites(bsai_cds(), "BsaI", "Ubiq")
    expect_false(has_bsai(out), info = paste("seed", s))
    expect_equal(protein_of(out), protein_of(bsai_cds()), info = paste("seed", s))
    expect_equal(nchar(out), nchar(bsai_cds()), info = paste("seed", s))
  }
})

test_that("a site spanning a codon boundary is removed", {
  # AGG TCT CAA contains GGTCTC across codons 1-3.
  cds <- "ATGAGGTCTCAATAA"
  expect_true(has_bsai(cds))
  for (s in test_seeds) {
    set.seed(s)
    out <- remove_enzyme_sites(cds, "BsaI", "Ubiq")
    expect_false(has_bsai(out), info = paste("seed", s))
    expect_equal(protein_of(out), protein_of(cds), info = paste("seed", s))
  }
})

test_that("the reverse-complement strand (GAGACC) is also removed", {
  cds <- "ATGGAGACCAAATAA"
  expect_true(has_bsai(cds))
  for (s in test_seeds) {
    set.seed(s)
    out <- remove_enzyme_sites(cds, "BsaI", "Ubiq")
    expect_false(has_bsai(out), info = paste("seed", s))
    expect_equal(protein_of(out), protein_of(cds), info = paste("seed", s))
  }
})

test_that("several enzymes can be removed at once", {
  # BsaI (GGTCTC) and BsmBI (CGTCTC), one of each.
  cds <- "ATGGCTGGTCTCAAACGTCTCGAATAA"
  set.seed(1)
  out <- remove_enzyme_sites(cds, c("BsaI", "BsmBI"), "Ubiq")
  expect_false(grepl("GGTCTC|GAGACC|CGTCTC|GAGACG", out))
  expect_equal(protein_of(out), protein_of(cds))
})

test_that("a sequence without the site is returned unchanged", {
  cds <- fixture_cds()
  expect_false(has_bsai(cds))
  expect_identical(remove_enzyme_sites(cds, "BsaI", "Ubiq"), cds)
})

test_that("an empty enzyme list leaves the sequence unchanged", {
  expect_identical(remove_enzyme_sites(bsai_cds(), character(0), "Ubiq"),
                   bsai_cds())
})

test_that("output is uppercase and junk characters are stripped", {
  set.seed(1)
  out <- remove_enzyme_sites("atg gct\naaa taa", "BsaI", "Ubiq")
  expect_identical(out, "ATGGCTAAATAA")
})

# ---- IUPAC ambiguity codes -------------------------------------------------

test_that("find_enzyme_site_starts expands IUPAC codes in the site", {
  # BanI = GGYRCC: Y = C/T, R = A/G. Matches GGCGCC, GGTACC, GGCACC, GGTGCC.
  for (site in c("GGCGCC", "GGTACC", "GGCACC", "GGTGCC")) {
    expect_equal(find_enzyme_site_starts(paste0("AA", site, "TT"), "GGYRCC"), 3L,
                 info = site)
  }
})

test_that("find_enzyme_site_starts does not match bases outside the code", {
  # Y is C or T only, R is A or G only.
  expect_length(find_enzyme_site_starts("AAGGAGCCTT", "GGYRCC"), 0)  # A at the Y position
  expect_length(find_enzyme_site_starts("AAGGCCCCTT", "GGYRCC"), 0)  # C at the R position
})

test_that("an N in the searched sequence does not match a specific base", {
  expect_length(find_enzyme_site_starts("AAGGNGCCTT", "GGYRCC"), 0)
  expect_length(find_enzyme_site_starts("AAGGCNCCTT", "GGYRCC"), 0)
})

test_that("an N in the site matches any base", {
  expect_equal(find_enzyme_site_starts("TTGACAAAGTCTT", "GACNNNGTC"), 3L)
})

test_that("find_enzyme_site_starts finds overlapping and repeated sites", {
  expect_equal(sort(find_enzyme_site_starts("GGCGCCGGTACC", "GGYRCC")), c(1L, 7L))
  expect_equal(sort(find_enzyme_site_starts("AAAAAA", "AAA")), 1:4)
})

test_that("the reverse complement of an ambiguous, non-palindromic site is found", {
  # GACNNNGTC is palindromic; use a non-palindromic ambiguous site, GGTCTNN
  # whose reverse complement is NNAGACC.
  site <- "GGTCTNN"
  rc <- reverse_complement(site)
  expect_equal(rc, "NNAGACC")
  expect_equal(find_enzyme_site_starts("TTCCAGACCTT", c(site, rc)), 3L)
})

test_that("BanI sites are now removed and the protein is unchanged", {
  # Codons GGC GCC = Gly Ala, which contain the BanI site GGCGCC.
  cds <- "ATGGGCGCCAAATAA"
  for (s in test_seeds) {
    set.seed(s)
    out <- remove_enzyme_sites(cds, "BanI", "Ubiq")
    expect_false(grepl("GG[CT][AG]CC", out), info = paste("seed", s))
    expect_equal(protein_of(out), protein_of(cds), info = paste("seed", s))
  }
})

test_that("BsoBI (CYCGRG) sites are removed, including CTCGAG", {
  cds <- "ATGCTCGAGAAATAA"  # Leu Glu: CTCGAG
  expect_true(grepl("C[CT]CG[AG]G", cds))
  for (s in test_seeds) {
    set.seed(s)
    out <- remove_enzyme_sites(cds, "BsoBI", "Ubiq")
    expect_false(grepl("C[CT]CG[AG]G", out), info = paste("seed", s))
    expect_equal(protein_of(out), protein_of(cds), info = paste("seed", s))
  }
})

test_that("ambiguous and exact enzymes can be removed together", {
  cds <- "ATGGGCGCCAAAGGTCTCGAATAA"  # BanI site then BsaI site
  set.seed(1)
  out <- remove_enzyme_sites(cds, c("BanI", "BsaI"), "Ubiq")
  expect_false(grepl("GG[CT][AG]CC|GGTCTC|GAGACC", out))
  expect_equal(protein_of(out), protein_of(cds))
})

# ---- Warnings when a site cannot be removed or an enzyme is unknown ---------------

test_that("an unknown enzyme name gives a warning and is otherwise ignored", {
  expect_warning(out <- remove_enzyme_sites(bsai_cds(), "NotAnEnzyme", "Ubiq"), "NotAnEnzyme")
  expect_identical(out, bsai_cds())
})

test_that("known enzymes are still removed when an unknown name is in the list", {
  set.seed(1)
  expect_warning(out <- remove_enzyme_sites(bsai_cds(), c("BsaI", "Nope"), "Ubiq"), "Nope")
  expect_false(has_bsai(out))
  expect_equal(protein_of(out), protein_of(bsai_cds()))
})

test_that("the warning names each unknown enzyme once", {
  w <- tryCatch(remove_enzyme_sites(bsai_cds(), c("Nope", "BsaI", "Nope", "Zilch"), "Ubiq"),
                warning = function(w) conditionMessage(w))
  expect_match(w, "Nope, Zilch")
  expect_equal(lengths(regmatches(w, gregexpr("Nope", w))), 1)
})

test_that("unknown enzymes also give a warning in High expression (CAI1) mode", {
  expect_warning(out <- remove_enzyme_sites(bsai_cds(), "NotAnEnzyme", "CAI1"), "NotAnEnzyme")
  expect_identical(out, bsai_cds())
})

# ---- Sites whose first codon cannot be changed ---------------------------------------
# A BsaI site (GGTCTC) that starts on the second base of the Trp codon TGG: the codons are
# TGG | TCT | CAA. The method used to change only the first codon a site touches, and Trp
# has no synonyms, so it could not remove this site. It now falls back to the other codons
# the site covers, but only when the first codon cannot change.
trp_site_cds <- function() "ATGGCTTGGTCTCAATAA"

test_that("a site starting inside a Trp codon is now removed by changing a neighbour", {
  expect_true(has_bsai(trp_site_cds()))
  for (method in c("Ubiq", "Chance")) {
    for (s in test_seeds) {
      set.seed(s)
      out <- remove_enzyme_sites(trp_site_cds(), "BsaI", method)
      info <- paste(method, "seed", s)
      expect_false(has_bsai(out), info = info)
      expect_equal(protein_of(out), protein_of(trp_site_cds()), info = info)
      expect_equal(nchar(out), nchar(trp_site_cds()), info = info)
      expect_equal(substr(out, 7, 9), "TGG", info = info)  # the Trp codon cannot change
    }
  }
})

test_that("that fix gives no warning, and only codons the site covers are changed", {
  for (s in test_seeds) {
    set.seed(s)
    expect_no_warning(out <- remove_enzyme_sites(trp_site_cds(), "BsaI", "Ubiq"))
    changed <- which(strsplit(out, "")[[1]] != strsplit(trp_site_cds(), "")[[1]])
    expect_true(all(changed >= 7 & changed <= 15), info = paste("seed", s))  # codons 3-5
  }
})

test_that("the fallback still chooses at random among the changes that work", {
  outs <- vapply(1:30, function(s) { set.seed(s); remove_enzyme_sites(trp_site_cds(), "BsaI", "Chance") }, character(1))
  expect_gt(length(unique(outs)), 1)
})

test_that("the fallback works in the full pipeline for methods that do not recode the site", {
  for (m in c("None", "GLO", "Ubiq", "Chance")) {
    for (s in 1:3) {
      set.seed(s)
      expect_no_warning(res <- optimize_transgene(trp_site_cds(), opt_method = m, enzymes = "BsaI", remove_pirnas = FALSE))
      expect_false(has_bsai(res$optimized_cds), info = paste(m, "seed", s))
      expect_equal(protein_of(res$optimized_cds), protein_of(trp_site_cds()), info = paste(m, "seed", s))
    }
  }
})

# ---- Reference outputs: nothing that already worked may change ------------------------
# These strings were produced by the previous version of the algorithm (before the
# fallback existed). Whenever the first codon of a site can be changed, the new version
# must draw exactly the same random numbers and give exactly the same sequence.
reference_cases <- list(
  bsai_ubiq = list(cds = "ATGGCTAAAGGTCTCTCTGAACCTACTGTTAATCAAATTTATTGGCGTAAAGATGGTGAAGCTTTTCCAAATTTTCTGGGACTCTGTTGGAAAGCCGAAGAATGTTAA", enz = "BsaI", col = "Ubiq"),
  bsai_chance = list(cds = "ATGGCTAAAGGTCTCTCTGAACCTACTGTTAATCAAATTTATTGGCGTAAAGATGGTGAAGCTTTTCCAAATTTTCTGGGACTCTGTTGGAAAGCCGAAGAATGTTAA", enz = "BsaI", col = "Chance"),
  multi_ubiq = list(cds = "ATGCTCGAGAAAGGTCTCGGCGCCTAA", enz = c("XhoI", "BsaI", "BanI"), col = "Ubiq"),
  multi_chance = list(cds = "ATGCTCGAGAAAGGTCTCGGCGCCTAA", enz = c("XhoI", "BsaI", "BanI"), col = "Chance"),
  revcomp = list(cds = "ATGGAGACCAAAGAGACCAAATAA", enz = "BsaI", col = "Ubiq")
)
reference_outputs <- rbind(
  c("bsai_ubiq", 1, "ATGGCTAAAGGACTCTCTGAACCTACTGTTAATCAAATTTATTGGCGTAAAGATGGTGAAGCTTTTCCAAATTTTCTGGGACTCTGTTGGAAAGCCGAAGAATGTTAA"),
  c("bsai_ubiq", 2, "ATGGCTAAAGGACTCTCTGAACCTACTGTTAATCAAATTTATTGGCGTAAAGATGGTGAAGCTTTTCCAAATTTTCTGGGACTCTGTTGGAAAGCCGAAGAATGTTAA"),
  c("bsai_ubiq", 3, "ATGGCTAAAGGACTCTCTGAACCTACTGTTAATCAAATTTATTGGCGTAAAGATGGTGAAGCTTTTCCAAATTTTCTGGGACTCTGTTGGAAAGCCGAAGAATGTTAA"),
  c("bsai_chance", 1, "ATGGCTAAAGGCCTCTCTGAACCTACTGTTAATCAAATTTATTGGCGTAAAGATGGTGAAGCTTTTCCAAATTTTCTGGGACTCTGTTGGAAAGCCGAAGAATGTTAA"),
  c("bsai_chance", 2, "ATGGCTAAAGGCCTCTCTGAACCTACTGTTAATCAAATTTATTGGCGTAAAGATGGTGAAGCTTTTCCAAATTTTCTGGGACTCTGTTGGAAAGCCGAAGAATGTTAA"),
  c("bsai_chance", 3, "ATGGCTAAAGGCCTCTCTGAACCTACTGTTAATCAAATTTATTGGCGTAAAGATGGTGAAGCTTTTCCAAATTTTCTGGGACTCTGTTGGAAAGCCGAAGAATGTTAA"),
  c("multi_ubiq", 1, "ATGCTTGAGAAAGGACTCGGAGCCTAA"),
  c("multi_ubiq", 2, "ATGCTTGAGAAAGGACTCGGAGCCTAA"),
  c("multi_ubiq", 3, "ATGCTTGAGAAAGGACTCGGAGCCTAA"),
  c("multi_chance", 1, "ATGCTTGAGAAAGGGCTCGGAGCCTAA"),
  c("multi_chance", 2, "ATGCTGGAGAAAGGACTCGGAGCCTAA"),
  c("multi_chance", 3, "ATGCTGGAGAAAGGACTCGGGGCCTAA"),
  c("revcomp", 1, "ATGGAAACCAAAGAAACCAAATAA"),
  c("revcomp", 2, "ATGGAAACCAAAGAAACCAAATAA"),
  c("revcomp", 3, "ATGGAAACCAAAGAAACCAAATAA")
)

test_that("outputs are identical to the previous algorithm whenever the first codon can be changed", {
  for (i in seq_len(nrow(reference_outputs))) {
    id <- reference_outputs[i, 1]; seed <- as.integer(reference_outputs[i, 2])
    cs <- reference_cases[[id]]
    set.seed(seed)
    out <- remove_enzyme_sites(cs$cds, cs$enz, cs$col)
    expect_identical(out, reference_outputs[i, 3], info = paste(id, "seed", seed))
  }
})

# ---- A site that truly cannot be removed -----------------------------------------------
# The enzyme table has no site made only of Met and Trp codons, so these tests give the
# internal function a made-up one: TGGATG (Trp Met) and its reverse complement CATCCA.
unfixable_sites <- c("TGGATG", "CATCCA")
unfixable_cds <- function() "ATGTGGATGTGGTAA"

test_that("a site made only of Met and Trp codons gives one clear warning", {
  expect_warning(out <- remove_enzyme_sites_weighted(unfixable_cds(), unfixable_sites, "Ubiq"),
                 "1 restriction site\\(s\\) could not be removed.*position\\(s\\) 4")
  expect_identical(out, unfixable_cds())
})

test_that("the warning is given once, not once per attempt", {
  n <- 0
  withCallingHandlers(
    remove_enzyme_sites_weighted(unfixable_cds(), unfixable_sites, "Ubiq"),
    warning = function(w) { n <<- n + 1; invokeRestart("muffleWarning") }
  )
  expect_equal(n, 1)
})

test_that("only the sites that remain are reported when others were removed", {
  # A removable BsaI site (codons 4-5) plus the unremovable TGGATG.
  cds <- "ATGGCTAAAGGTCTCTCTGAATGGATGTGGTAA"
  sites <- c("GGTCTC", "GAGACC", unfixable_sites)
  set.seed(1)
  expect_warning(out <- remove_enzyme_sites_weighted(cds, sites, "Ubiq"),
                 "1 restriction site\\(s\\) could not be removed")
  expect_false(has_bsai(out))
  expect_true(grepl("TGGATG", out))
})

test_that("no warning is given when every site is removed", {
  for (method in c("Ubiq", "Chance", "CAI1")) {
    set.seed(1)
    expect_no_warning(remove_enzyme_sites(bsai_cds(), "BsaI", method))
  }
})

test_that("no warning is given for an empty enzyme list or a sequence without sites", {
  expect_no_warning(remove_enzyme_sites(bsai_cds(), character(0), "Ubiq"))
  expect_no_warning(remove_enzyme_sites(fixture_cds(), "BsaI", "Ubiq"))
})

test_that("removal that cannot progress stops early instead of repeating 100 times", {
  calls <- 0
  real_ranges <- find_enzyme_site_ranges
  testthat::local_mocked_bindings(
    find_enzyme_site_ranges = function(dna_seq, sites) { calls <<- calls + 1; real_ranges(dna_seq, sites) }
  )
  suppressWarnings(remove_enzyme_sites_weighted(unfixable_cds(), unfixable_sites, "Ubiq"))
  expect_lt(calls, 10)
})
