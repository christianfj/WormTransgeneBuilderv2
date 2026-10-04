# "High expression" (CAI1) mode: the sequence is recoded to the best codon for
# every amino acid, so when a restriction site is removed the change should cost
# as little CAI as possible. The rule is deterministic, so no seeds are needed
# for the results themselves, and repeated runs must give identical output.
#
# CAI is scored with the graded BringMansWeights column (best codon = 1).

cai_w <- function(codon) {
  cf <- WormTransgeneBuilder::codon_freqs
  cf$BringMansWeights[match(toupper(codon), toupper(cf$Codon))]
}

codons_of <- function(dna) substring(dna, seq(1, nchar(dna), 3), seq(3, nchar(dna), 3))

# Brute-force reference: try every single synonymous codon change, keep those
# that leave no site, and return the best CAI-loss score among them.
best_possible_score <- function(dna, site_regex) {
  cf <- WormTransgeneBuilder::codon_freqs
  cods <- codons_of(dna)
  best <- -Inf
  for (i in seq_along(cods)) {
    aa <- cf$Amino[match(cods[i], toupper(cf$Codon))]
    if (is.na(aa)) next
    for (alt in setdiff(toupper(cf$Codon[cf$Amino == aa]), cods[i])) {
      trial <- cods; trial[i] <- alt
      if (!grepl(site_regex, paste(trial, collapse = ""))) {
        best <- max(best, log(max(cai_w(alt), 1e-6)) - log(max(cai_w(cods[i]), 1e-6)))
      }
    }
  }
  best
}

realised_score <- function(before, after) {
  b <- codons_of(before); a <- codons_of(after)
  changed <- which(a != b)
  sum(log(pmax(cai_w(a[changed]), 1e-6)) - log(pmax(cai_w(b[changed]), 1e-6)))
}

xhoi <- "CTCGAG"

test_that("a site made only of best codons is removed instead of raising an error", {
  cds <- "ATGCTCGAGAAGTAA"  # CTC GAG: both are the best codons; CTCGAG is XhoI
  out <- remove_enzyme_sites(cds, "XhoI", "CAI1")
  expect_false(grepl(xhoi, out))
  expect_equal(protein_of(out), protein_of(cds))
  expect_equal(nchar(out), nchar(cds))
})

test_that("the CAI1 site removal is deterministic", {
  cds <- "ATGCTCGAGAAGTAA"
  set.seed(1); a <- remove_enzyme_sites(cds, "XhoI", "CAI1")
  set.seed(2); b <- remove_enzyme_sites(cds, "XhoI", "CAI1")
  set.seed(99); d <- remove_enzyme_sites(cds, "XhoI", "CAI1")
  expect_identical(a, b)
  expect_identical(a, d)
})

test_that("the change chosen has the smallest possible CAI loss", {
  cds <- "ATGCTCGAGAAGTAA"
  out <- remove_enzyme_sites(cds, "XhoI", "CAI1")
  expect_equal(realised_score(cds, out), best_possible_score(cds, xhoi))
  # Exactly one codon should have changed.
  expect_equal(sum(codons_of(out) != codons_of(cds)), 1)
})

test_that("CAI drops below 1 after a site had to be changed, but only slightly", {
  cds <- "ATGCTCGAGAAGTAA"
  out <- remove_enzyme_sites(cds, "XhoI", "CAI1")
  expect_equal(calc_cai(cds, "BringMansWeights"), 1)
  expect_lt(calc_cai(out, "BringMansWeights"), 1)
  expect_gt(calc_cai(out, "BringMansWeights"), 0.9)
})

test_that("the cheapest codon anywhere in the site is changed, not just the first one", {
  # In CTC GAG the Glu change (GAG -> GAA) and the Leu change (CTC -> CTT) cost
  # different amounts. Whichever is cheaper must be the one chosen.
  cds <- "ATGCTCGAGAAGTAA"
  out <- remove_enzyme_sites(cds, "XhoI", "CAI1")
  changed <- which(codons_of(out) != codons_of(cds))
  cost <- function(i, alt) log(cai_w(alt)) - log(cai_w(codons_of(cds)[i]))
  expect_gte(cost(changed, codons_of(out)[changed]),
             min(cost(2, "CTT"), cost(3, "GAA")))
  expect_equal(realised_score(cds, out), best_possible_score(cds, xhoi))
})

test_that("a site starting inside a Trp codon is still removed by changing a neighbour", {
  # TGG TCT C.. contains BsaI GGTCTC starting on the G of TGG. Trp has only one
  # codon, so the fix has to come from another codon of the site.
  cds <- "ATGGCTTGGTCTCAATAA"
  expect_true(grepl("GGTCTC", cds))
  out <- remove_enzyme_sites(cds, "BsaI", "CAI1")
  expect_false(grepl("GGTCTC|GAGACC", out))
  expect_equal(protein_of(out), protein_of(cds))
  expect_equal(realised_score(cds, out), best_possible_score(cds, "GGTCTC|GAGACC"))
})

test_that("ambiguity-code enzymes are handled in CAI1 mode", {
  cds <- "ATGGGCGCCAAATAA"  # GGC GCC = BanI GGYRCC
  out <- remove_enzyme_sites(cds, "BanI", "CAI1")
  expect_false(grepl("GG[CT][AG]CC", out))
  expect_equal(protein_of(out), protein_of(cds))
})

test_that("several sites and several enzymes are all removed", {
  cds <- "ATGCTCGAGAAGGGTCTCGGCGCCTAA"  # XhoI + BsaI + BanI
  out <- remove_enzyme_sites(cds, c("XhoI", "BsaI", "BanI"), "CAI1")
  expect_false(grepl("CTCGAG|GGTCTC|GAGACC|GG[CT][AG]CC", out))
  expect_equal(protein_of(out), protein_of(cds))
})

test_that("a sequence with no site is returned unchanged", {
  cds <- "ATGGCTAAAGGTTAA"
  expect_identical(remove_enzyme_sites(cds, "XhoI", "CAI1"), cds)
})

test_that("a site that cannot be removed gives a warning and leaves the sequence intact", {
  # A custom site made only of Met/Trp codons cannot be broken synonymously.
  cds <- "ATGTGGATGTGGTAA"
  expect_warning(
    out <- remove_enzyme_sites_max_cai(cds, "TGGATG", "BringMansWeights"),
    "could not be removed"
  )
  expect_identical(out, cds)
})

test_that("the full pipeline with CAI1 and enzyme removal works and is repeatable", {
  cds <- "ATGCTCGAGAAGTAA"
  set.seed(1)
  a <- optimize_transgene(cds, opt_method = "CAI1", enzymes = "XhoI", remove_pirnas = FALSE)
  set.seed(2)
  b <- optimize_transgene(cds, opt_method = "CAI1", enzymes = "XhoI", remove_pirnas = FALSE)
  expect_false(grepl(xhoi, a$final_seq))
  expect_identical(a$final_seq, b$final_seq)
  expect_equal(protein_of(a$optimized_cds), protein_of(cds))
  expect_lt(calc_cai(a$optimized_cds, "BringMansWeights"), 1)
})

test_that("other methods keep their random, weighted behaviour", {
  cds <- "ATGCTCGAGAAGTAA"
  outs <- vapply(1:10, function(s) {
    set.seed(s); remove_enzyme_sites(cds, "XhoI", "Ubiq")
  }, character(1))
  expect_false(any(grepl(xhoi, outs)))
  expect_gt(length(unique(outs)), 1)  # random: not always the same answer
})
