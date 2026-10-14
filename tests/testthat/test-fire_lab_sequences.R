# The Fire lab introns and the 5' synthetic spliced UTR of the app against pPD95.75 as deposited at Addgene (plasmid 1494; record "pPD95_75" in
# data-raw/fire_lab_vectors/, parsed by dev/fire_lab_parse.R) and against the authors' own pPD95.75 file (data-raw/fire_lab/pPD95.75_excerpt.txt).
# Written up in docs/FIRE_LAB_SEQUENCES.md. Skipped when the files are not there (not running from the source tree).

fire_csv <- function(name) utils::read.csv(system.file("extdata", name, package = "WormTransgeneBuilder"), stringsAsFactors = FALSE)

fire_root <- function() {
  root <- normalizePath(testthat::test_path("..", ".."))
  skip_if(!file.exists(file.path(root, "data-raw", "fire_lab_vectors", "fire-lab-c-elegans-vector-kit-genbank-sequences.txt")),
          "the Fire lab vector file is not available (not running from the source tree)")
  root
}

addgene_ppd95_75 <- function() {
  root <- fire_root()
  source(file.path(root, "dev", "fire_lab_parse.R"), local = TRUE)
  list(seq = read_fire_lab_vectors(file.path(root, "data-raw", "fire_lab_vectors", "fire-lab-c-elegans-vector-kit-genbank-sequences.txt"))[["pPD95_75"]]$seq,
       find = find_in_circular, root = root)
}

fire_introns <- function() {
  x <- fire_csv("Introns.csv")
  toupper(unlist(x[x$Name == "Fire_lab_introns", c("Sequence_1", "Sequence_2", "Sequence_3")], use.names = FALSE))
}

test_that("the Addgene record of pPD95.75 is read correctly (4,487 bp, gfp is a reading frame of ATG...TAG)", {
  a <- addgene_ppd95_75()
  expect_equal(nchar(a$seq), 4487)
  start <- a$find("ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTC", a$seq)
  expect_identical(start$strand, "+")
  expect_identical(substr(a$seq, start$start + 867, start$start + 869), "TAG")   # the gfp coding region is 870 bases
})

test_that("all three Fire lab introns of the app are in pPD95.75 (Addgene), at the three gfp intron positions", {
  a <- addgene_ppd95_75()
  start <- a$find("ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTC", a$seq)$start
  before <- vapply(fire_introns(), function(x) { h <- a$find(x, a$seq); expect_identical(h$strand, "+"); h$start - start }, 0)
  expect_equal(unname(before), c(171, 376, 591))        # bases of gfp before each intron: phases 0, 1, 0
})

test_that("the 5' UTR 'Fire lab synthetic spliced' is in pPD95.75 (Addgene), with intron A in the middle", {
  a <- addgene_ppd95_75()
  utr <- toupper(fire_csv("5pUTRs.csv")$Sequence[2])
  h <- a$find(utr, a$seq)
  expect_identical(h$strand, "+")
  expect_equal(nchar(utr), 92)
  expect_match(substr(utr, 26, 67), "^GT.*TTTTCAG$")                   # the 42-base intron A
  expect_equal(round(100 * mean(strsplit(utr, "")[[1]] %in% c("G", "C"))), 47)   # the menu says 47% GC
})

test_that("the authors' own pPD95.75 file differs from Addgene's only in the third gfp intron (9 bases)", {
  a <- addgene_ppd95_75()
  lines <- grep("^ +[0-9]+ ", readLines(file.path(a$root, "data-raw", "fire_lab", "pPD95.75_excerpt.txt")), value = TRUE)
  excerpt <- toupper(paste(gsub("[^A-Za-z]", "", sub("^ +[0-9]+ ", "", lines)), collapse = ""))
  start <- a$find(substr(excerpt, 1, 40), a$seq)$start
  same <- substr(paste0(a$seq, a$seq), start, start + nchar(excerpt) - 1)
  differing <- which(strsplit(excerpt, "")[[1]] != strsplit(same, "")[[1]])
  expect_length(differing, 9)
  third <- fire_introns()[3]
  intron_start <- regexpr(third, same, fixed = TRUE)[1]
  expect_true(all(differing >= intron_start & differing < intron_start + nchar(third)))
})

test_that("all three Fire lab introns have the Fire lab frame (12 + 7 + 11 + 7 + 14 = 51 bases)", {
  expect_match(tolower(fire_introns()), "^gtaagtttaaac[acgt]{7}tactaactaac[acgt]{7}atttaaattttcag$")
})

test_that("the three introns of a gfp optimized with the Codon Adapter (Redemann et al. 2011) are the three Fire lab introns of the app, in order", {
  path <- testthat::test_path("..", "..", "data-raw", "fire_lab", "redemann2011_gfp_output.fa")
  skip_if(!file.exists(path), "the Redemann output is not available (not running from the source tree)")
  s <- readLines(path)[2]
  runs <- regmatches(s, gregexpr("[acgt]+", s))[[1]]
  expect_identical(runs[1], "aaaa")                 # the consensus start
  expect_identical(toupper(runs[-1]), fire_introns())
  exon <- toupper(gsub("[acgt]+", "", s))
  expect_equal(nchar(exon) %% 3, 0)
  expect_identical(substr(exon, nchar(exon) - 2, nchar(exon)), "TAA")
})
