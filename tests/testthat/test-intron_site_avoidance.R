# Introns must not create a recognition site of an enzyme the user wants removed, neither inside the intron
# nor where it joins the exon (insert_introns(avoid_sites =, pool =), and optimize_transgene(intron_pool =)).

# As in test-introns.R: every odd codon from the 3rd on is AAG (the preferred intron position).
intron_cds <- function() paste0("ATG", strrep("GCTAAG", 100), "TAA")

site_in <- function(seq, sites) any(vapply(sites, function(s) nrow(find_enzyme_site_ranges(toupper(seq), s)) > 0, TRUE))

bsai <- c("GGTCTC", "GAGACC")
bad_intron <- "gtaagtttaaggtctcaaattttcag"      # contains GGTCTC
good_a <- "gtaagtttaacctgcaaattttcag"
good_b <- "gtaagtttcatgcaacaattttcag"
good_c <- "gtgagtttcgcaaccaaattttcag"

test_that("without avoid_sites nothing changes and no random numbers are used", {
  set.seed(1); a <- insert_introns(intron_cds(), c(bad_intron, good_a))
  set.seed(1); b <- insert_introns(intron_cds(), c(bad_intron, good_a), pool = c(good_b, good_c))
  expect_identical(a$seq, b$seq)
  expect_identical(a$introns, c(bad_intron, good_a))
  expect_length(a$unresolved, 0)
})

test_that("an intron with a site is swapped for a spare one that has none", {
  set.seed(1)
  res <- insert_introns(intron_cds(), c(bad_intron, good_a), avoid_sites = bsai, pool = c(bad_intron, good_a, good_b, good_c))
  expect_false(site_in(res$seq, bsai))
  expect_length(res$unresolved, 0)
  expect_true(res$introns[1] %in% c(good_b, good_c))      # not the bad one, and not the one already used at position 2
  expect_identical(res$introns[2], good_a)
  expect_identical(gsub("[acgt]", "", res$seq), intron_cds())
})

test_that("a site formed only where the intron joins the exon is found and avoided", {
  # The exon before the intron ends in AAG, the intron starts gtaagttt: AAGGTAAGTTT spans the junction.
  site <- "AAGGTAAGTTT"
  set.seed(1)
  res <- insert_introns(intron_cds(), c(good_a, good_a), avoid_sites = site, pool = c(good_c))
  expect_true(site_in(insert_introns(intron_cds(), c(good_a, good_a))$seq, site))   # the junction does make the site
  expect_length(res$unresolved, 1)                                                   # only one spare intron for two positions
  expect_true(good_c %in% res$introns)
})

test_that("an intron that cannot be replaced is reported, not hidden", {
  res <- insert_introns(intron_cds(), c(bad_intron, good_a), avoid_sites = bsai)
  expect_identical(res$unresolved, 1L)
  expect_identical(res$introns, c(bad_intron, good_a))
})

test_that("the pipeline swaps introns and warns when it cannot", {
  cds <- glo_cds()
  set.seed(2)
  out <- optimize_transgene(cds, opt_method = "Ubiq", enzymes = "BsaI", remove_pirnas = FALSE,
                            introns = c(bad_intron, good_a), intron_pool = c(good_b, good_c), force_reading_frame = TRUE)
  expect_false(site_in(out$final_seq, bsai))
  expect_false(bad_intron %in% out$introns_used)
  expect_equal(out$intron_lengths, nchar(out$introns_used))
  set.seed(2)
  expect_warning(optimize_transgene(cds, opt_method = "Ubiq", enzymes = "BsaI", remove_pirnas = FALSE,
                                    introns = c(bad_intron, good_a)), "creates a recognition site")
})
