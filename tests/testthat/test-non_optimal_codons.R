# For the "High expression" (CAI1) method the viewer underlines every codon that
# is not the best codon for its amino acid (CAI1 weight not 1). These tests check
# that the right codons are found and that their positions are mapped correctly
# onto the final sequence, which also contains tags, introns, UTRs and a promoter.

cai1_weight <- function(codon) {
  cf <- WormTransgeneBuilder::codon_freqs
  cf$CAI1[match(toupper(codon), toupper(cf$Codon))]
}
codons_in <- function(dna) substring(dna, seq(1, nchar(dna) - 2, 3), seq(3, nchar(dna), 3))

test_that("only codons that are not the best codon are reported", {
  # ATG (Met, best) | CTT (Leu, not best) | GAG (Glu, best) | AAA (Lys, not best) | TAA (stop)
  r <- non_optimal_codon_ranges("ATGCTTGAGAAATAA")
  expect_equal(r$start, c(4, 10))
  expect_equal(r$end, c(6, 12))
})

test_that("a sequence made only of best codons has nothing to underline", {
  r <- non_optimal_codon_ranges("ATGCTCGAGAAGTAA")
  expect_equal(nrow(r), 0)
})

test_that("stop codons are never reported, whichever stop is used", {
  for (stop_codon in c("TAA", "TAG", "TGA")) {
    expect_equal(nrow(non_optimal_codon_ranges(paste0("ATGCTC", stop_codon))), 0, info = stop_codon)
  }
})

test_that("positions are shifted by the start of the CDS in the final sequence", {
  r <- non_optimal_codon_ranges("ATGCTTGAGAAATAA", cds_start = 11)
  expect_equal(r$start, c(14, 20))
  expect_equal(r$end, c(16, 22))
})

test_that("positions are shifted by every intron inserted before them", {
  # Two introns (10 and 20 bp) inserted after bases 3 and 9 of the CDS.
  r <- non_optimal_codon_ranges("ATGCTTGAGAAATAA", cds_start = 1,
                                intron_insertions = c(3, 9), intron_lengths = c(10, 20))
  expect_equal(r$start, c(4 + 10, 10 + 30))
  expect_equal(r$end, c(6 + 10, 12 + 30))
})

test_that("codons containing N, and an incomplete last codon, are ignored", {
  expect_equal(nrow(non_optimal_codon_ranges("ATGNNNCTC")), 0)
  r <- non_optimal_codon_ranges("ATGCTTGA")
  expect_equal(r$start, 4)
})

test_that("empty and very short input give an empty result", {
  expect_equal(nrow(non_optimal_codon_ranges("")), 0)
  expect_equal(nrow(non_optimal_codon_ranges("AT")), 0)
})

test_that("the pipeline returns the intron positions needed for the mapping", {
  cds <- paste0("ATG", strrep("GCTAAG", 60), "TAA")
  res <- optimize_transgene(cds, opt_method = "CAI1", enzymes = NULL, remove_pirnas = FALSE,
                            introns = c("gtaagtttcag", "gtaagtttcag"))
  expect_length(res$intron_insertions, 2)
  expect_equal(res$intron_lengths, c(11, 11))
  res0 <- optimize_transgene(cds, opt_method = "CAI1", enzymes = NULL, remove_pirnas = FALSE)
  expect_length(res0$intron_insertions, 0)
  expect_length(res0$intron_lengths, 0)
})

test_that("underlined ranges land on the right codons in the final sequence", {
  # 120 codons with XhoI sites (CTCGAG) before, between and after the introns.
  # In CAI1 mode every codon is the best codon except those changed to remove a site.
  filler <- function(n) strrep("GCTAAG", n / 2)
  cds <- paste0("ATG", filler(20), "CTCGAG", filler(80), "CTCGAG", filler(80), "CTCGAG", filler(20), "TAA")
  expect_equal(nchar(cds) %% 3, 0)

  res <- optimize_transgene(
    cds, opt_method = "CAI1", enzymes = "XhoI", remove_pirnas = FALSE,
    introns = c("gtaagtttcag", "gtaagtttcag", "gtaagtttcag"),
    tag_n_seqs = "CCAAAGAAGAAGCGTAAG", tag_c_seqs = "GACTACAAGGAC",
    add_consensus_start = TRUE, utr5_seq = "ttttttt", promoter_seq = "acgtacgtacgt"
  )
  ranges <- non_optimal_codon_ranges(res$optimized_cds, res$cds_start,
                                     res$intron_insertions, res$intron_lengths)

  # One changed codon per site.
  expected <- codons_in(res$optimized_cds)
  expected <- expected[!is.na(cai1_weight(expected)) & cai1_weight(expected) != 1]
  expect_equal(nrow(ranges), 3)
  expect_equal(length(expected), 3)

  # Each range, read from the FINAL sequence, spells exactly that codon.
  found <- toupper(substring(res$final_seq, ranges$start, ranges$end))
  expect_equal(found, expected)
  expect_true(all(cai1_weight(found) == 0))
  expect_true(all(ranges$end - ranges$start == 2))
  # Ranges are in increasing order and do not overlap.
  expect_false(is.unsorted(ranges$start))
})
