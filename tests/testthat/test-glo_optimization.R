# GLO (germline optimization, D. J. Dickinson; Fielmich et al. 2018) works with whole four-codon words: it scores every coding
# word, assembles a first sequence from the best words, and improves it word by word (see ?optimize_glo). The search is random, so
# most tests check properties that must hold for any seed. Two tests compare the result with the ORIGINAL Perl implementation
# (https://github.com/dannyhmg/germline), which was run unchanged on the same proteins: the scores it reached (per codon) are
# recorded below, and the app must reach them. See docs/GLO_VERIFICATION.md.

glo_seeds <- 1:5

test_that("the GLO fixture is a valid CDS", {
  expect_true(check_cds_validity(glo_cds())$valid)
})

test_that("GLO keeps the protein, the length and the alphabet", {
  for (s in glo_seeds) {
    set.seed(s)
    out <- optimize_glo(glo_cds())
    expect_equal(protein_of(out), protein_of(glo_cds()), info = paste("seed", s))
    expect_equal(nchar(out), nchar(glo_cds()), info = paste("seed", s))
    expect_match(out, "^[ACGT]+$", info = paste("seed", s))
  }
})

test_that("GLO never lowers the GLO score, and improves it a lot", {
  before <- calc_glo_score(glo_cds())
  after <- vapply(glo_seeds, function(s) { set.seed(s); calc_glo_score(optimize_glo(glo_cds())) }, numeric(1))
  expect_true(all(after >= before))
  expect_gt(min(after), 5 * before)
})

test_that("GLO keeps the last codon a stop codon", {
  for (s in glo_seeds) {
    set.seed(s)
    expect_true(is_stop_codon(last_codon(optimize_glo(glo_cds()))), info = paste("seed", s))
  }
})

test_that("GLO is reproducible for a fixed seed", {
  set.seed(7); a <- optimize_glo(glo_cds())
  set.seed(7); b <- optimize_glo(glo_cds())
  expect_identical(a, b)
})

test_that("GLO leaves Met and Trp codons alone", {
  # They have a single codon each, so no synonymous change exists.
  cds <- "ATGTGGATGTGGATGTGG"
  set.seed(1)
  expect_identical(optimize_glo(cds), cds)
})

test_that("sequences shorter than 4 codons are returned exactly as given", {
  expect_identical(optimize_glo("atggctaaa"), "atggctaaa")
})

test_that("lowercase input is accepted and comes back in uppercase", {
  set.seed(1)
  out <- optimize_glo(tolower(glo_cds()))
  expect_match(out, "^[ACGT]+$")
  expect_equal(protein_of(out), protein_of(glo_cds()))
})

test_that("calc_glo_score matches a direct lookup of each window in the score table", {
  # calc_glo_score finds each window's score by arithmetic on its bases. Check it against a
  # separate lookup that reads each window as a base-4 number. (More checks of the score
  # table are in test-glo_scores.R.)
  cds <- glo_cds()
  starts <- seq(1, nchar(cds) - 11, by = 3)
  windows <- substring(cds, starts, starts + 11)
  scores <- glo_score_table()[strtoi(chartr("ACGT", "0123", windows), base = 4L) + 1L]
  expect_equal(calc_glo_score(cds), sum(scores))
})

# ---- the same scores as the original implementation -------------------------------------------------------------------
# Scores per codon (the original's unit) reached by the unchanged original Perl code (OptimizerTools::optimize) for these proteins,
# the same for every random seed tried: 509.813429, 703.833824 and 668.897467.

original_glo_proteins <- c(
  small = "MAKGLSEPTVNQIYWRKDGEAFPNFLGLCWKAEEC",
  gfp = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFCYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
  random = "YSKPRCYGMRSAFCGTHRSDSFNQTKPLGQRRCDSRNCTYCIINFTNIPQSKYMRDPHEFHLTVSASNCYCSAQQRRMPCDWVYEHLVQIKDASHVLPPKWYAAGIGLSVYYCNIIIFGNPEHKNNLANNPSPYYSTAAKSTDLFYFRFL")
original_glo_scores <- c(small = 509.813429, gfp = 703.833824, random = 668.897467)

test_that("the heuristic method reaches the score the original implementation reaches, on three proteins and several seeds", {
  for (id in names(original_glo_proteins)) {
    dna <- back_translate(original_glo_proteins[[id]])
    for (s in 1:3) {
      set.seed(s)
      out <- optimize_glo(dna, method = "heuristic")
      expect_identical(protein_of(out), original_glo_proteins[[id]], info = id)
      expect_equal(calc_glo_score(out) / (nchar(out) / 3), original_glo_scores[[id]], tolerance = 1e-6, info = paste(id, "seed", s))
    }
  }
})

test_that("the search needs the four-codon words: a change of one codon at a time cannot reach the score", {
  # A first-guess-free, single-codon hill climb (what this function used to do) stalls far below the original's score.
  dna <- back_translate(original_glo_proteins[["gfp"]])
  expect_gt(calc_glo_score(optimize_glo(dna)) / 238, 1.5 * 391)   # the old search averaged about 391 per codon on this protein
})

test_that("stop codons and codons that are not plain A/C/G/T are left exactly as they were", {
  dna <- back_translate(original_glo_proteins[["small"]])
  withstop <- paste0(dna, "TGA")
  set.seed(1); out <- optimize_glo(withstop)
  expect_identical(substr(out, nchar(out) - 2, nchar(out)), "TGA")
  withn <- paste0(substr(dna, 1, 30), "NNN", substr(dna, 34, nchar(dna)))
  set.seed(1); out <- optimize_glo(withn)
  expect_identical(substr(out, 31, 33), "NNN")
  expect_identical(protein_of(gsub("NNN", "GCT", out)), protein_of(gsub("NNN", "GCT", withn)))
})

test_that("bases after the last complete codon are kept, and the search is exact on a second look", {
  dna <- paste0(back_translate(original_glo_proteins[["small"]]), "AC")
  set.seed(1); out <- optimize_glo(dna)
  expect_identical(substr(out, nchar(out) - 1, nchar(out)), "AC")
  # a sequence that is already a local optimum comes back unchanged
  set.seed(2); again <- optimize_glo(out)
  expect_identical(again, out)
})

test_that("a long sequence is optimized in reasonable time with every window scored exactly", {
  skip_on_cran()
  helpers <- testthat::test_path("..", "..", "dev", "benchmarks", "benchmark_helpers.R")
  skip_if(!file.exists(helpers), "benchmark helpers are not available")
  source(helpers, local = TRUE)
  g <- make_benchmark_cds(3000)
  set.seed(1)
  t <- system.time(out <- optimize_glo(g))[["elapsed"]]
  expect_lt(t, 30)
  expect_identical(protein_of(out), protein_of(g))
  expect_gt(calc_glo_score(out), 5 * calc_glo_score(g))
})

# ---- the exact method -----------------------------------------------------------------------------------------------------

test_that("the exact method finds the best of ALL codings (brute force on small proteins)", {
  syn <- glo_synonym_indices(); aas <- setdiff(names(syn)[lengths(syn) > 1 & lengths(syn) <= 4], "X")   # stop codons are never recoded
  set.seed(5)
  for (rep in 1:10) {
    p <- sample(aas, 7, TRUE)
    dna <- paste(vapply(p, function(x) glo_codon_names()[syn[[x]][1] + 1], ""), collapse = "")
    all <- expand.grid(lapply(p, function(x) syn[[x]]))
    best <- max(apply(all, 1, function(cd) sum(glo_window_s100(cd, 0:3))))
    out <- optimize_glo(dna)
    got <- sum(glo_window_s100(glo_codon_index(substring(out, seq(1, 19, 3), seq(3, 21, 3))), 0:3))
    expect_equal(got, best, info = paste(p, collapse = ""))
    expect_identical(protein_of(out), protein_of(dna))
  }
})

test_that("the exact optimum is at least as good as the heuristic, and as the original's scores", {
  for (id in names(original_glo_proteins)) {
    dna <- back_translate(original_glo_proteins[[id]])
    exact <- calc_glo_score(optimize_glo(dna)) / (nchar(dna) / 3)
    expect_gte(exact, original_glo_scores[[id]] - 1e-5)    # (the recorded scores are rounded to 6 decimals)
    set.seed(1)
    expect_gte(exact, calc_glo_score(optimize_glo(dna, method = "heuristic")) / (nchar(dna) / 3) - 1e-6)
  }
  expect_equal(calc_glo_score(optimize_glo(back_translate(original_glo_proteins[["gfp"]]))) / 238, 717.640, tolerance = 1e-5)   # recorded 2026-10-07
})

test_that("the exact score does not depend on the seed or on the tie-breaking jitter; diversity 0 is deterministic", {
  dna <- back_translate(original_glo_proteins[["random"]])
  scores <- vapply(c(0, 0.3, 3), function(d) { set.seed(11); calc_glo_score(optimize_glo(dna, diversity = d)) }, 0)
  expect_equal(max(scores) - min(scores), 0)
  expect_identical(optimize_glo(dna, diversity = 0), optimize_glo(dna, diversity = 0))
  set.seed(3); a <- optimize_glo(dna); set.seed(3); b <- optimize_glo(dna)
  expect_identical(a, b)
})

test_that("among all codings with the best GLO score the exact method returns the one with the highest CAI (brute force)", {
  syn <- glo_synonym_indices(); aas <- setdiff(names(syn)[lengths(syn) > 1 & lengths(syn) <= 4], "X")
  cf <- WormTransgeneBuilder::codon_freqs; u <- setNames(cf$Ubiq, toupper(cf$Codon))
  set.seed(8); checked <- 0
  for (rep in 1:40) {
    p <- sample(aas, 8, TRUE)
    all <- expand.grid(lapply(p, function(x) syn[[x]]))
    sc <- apply(all, 1, function(cd) sum(glo_window_s100(cd, 0:4)))
    top <- which(sc == max(sc))
    if (length(top) < 2) next                                                   # only cases with a tie say anything
    logw <- apply(all[top, , drop = FALSE], 1, function(cd) sum(log(u[glo_codon_names()[cd + 1]])))
    dna <- paste(vapply(p, function(x) glo_codon_names()[syn[[x]][1] + 1], ""), collapse = "")
    out <- optimize_glo(dna)
    got <- sum(log(u[substring(out, seq(1, 22, 3), seq(3, 24, 3))]))
    expect_equal(got, max(logw), info = paste(p, collapse = ""))
    checked <- checked + 1
  }
  expect_gt(checked, 5)
})

test_that("among equally scoring codings the exact method prefers the common codons", {
  # Positions that no high-scoring word depends on are free to take any codon; they should take the more frequent ones.
  dna <- back_translate(original_glo_proteins[["random"]])
  cf <- WormTransgeneBuilder::codon_freqs; u <- setNames(cf$Ubiq, toupper(cf$Codon))
  cod <- function(x) substring(x, seq(1, nchar(x) - 2, 3), seq(3, nchar(x), 3))
  mean_usage <- function(x) mean(u[cod(x)])
  set.seed(1)
  expect_gt(mean_usage(optimize_glo(dna)), mean_usage(optimize_glo(dna, method = "heuristic")) - 0.02)
})

test_that("the exact method keeps stop codons, N codons and leftover bases, and leaves short sequences alone", {
  dna <- back_translate(original_glo_proteins[["small"]])
  expect_identical(substr(optimize_glo(paste0(dna, "TGA")), nchar(dna) + 1, nchar(dna) + 3), "TGA")
  withn <- paste0(substr(dna, 1, 30), "NNN", substr(dna, 34, nchar(dna)))
  expect_identical(substr(optimize_glo(withn), 31, 33), "NNN")
  expect_identical(substr(optimize_glo(paste0(dna, "AC")), nchar(dna) + 1, nchar(dna) + 2), "AC")
  expect_identical(optimize_glo("ATGGCTAAA"), "ATGGCTAAA")
  expect_identical(optimize_glo("ATGGCTAAAGGT"), optimize_glo("ATGGCTAAAGGT", method = "heuristic"))
})

test_that("the exact method handles a 20 kb gene in about a second", {
  skip_on_cran()
  helpers <- testthat::test_path("..", "..", "dev", "benchmarks", "benchmark_helpers.R")
  skip_if(!file.exists(helpers), "benchmark helpers are not available")
  source(helpers, local = TRUE)
  g <- make_benchmark_cds(20000)
  t <- system.time(out <- optimize_glo(g))[["elapsed"]]
  expect_lt(t, 10)
  expect_identical(protein_of(out), protein_of(g))
})

# ---- avoiding restriction sites inside the exact optimization -----------------------------------------------------------------

test_that("with avoid_sites the exact method finds the best coding among those WITHOUT the sites (brute force)", {
  syn <- glo_synonym_indices()
  sites <- enzyme_search_sites(c("EcoRI", "BsaI"))
  proteins <- c("AEFKEFA", "GEFEFSA", "MGEFGEFK", "KEFRSEF")                      # E followed by F can give GAATTC (EcoRI)
  for (pr in proteins) {
    p <- strsplit(pr, "")[[1]]
    dna <- paste(vapply(p, function(x) glo_codon_names()[syn[[x]][1] + 1], ""), collapse = "")
    all <- expand.grid(lapply(p, function(x) syn[[x]]))
    strings <- apply(all, 1, function(cd) paste(glo_codon_names()[cd + 1], collapse = ""))
    has_site <- vapply(strings, function(s) nrow(find_enzyme_site_ranges(s, sites)) > 0, NA)
    totals <- apply(all, 1, function(cd) sum(glo_window_s100(cd, 0:(length(p) - 4))))
    expect_true(any(has_site) && any(!has_site), info = pr)                      # the test is only meaningful if sites can occur and can be avoided
    best_free <- max(totals[!has_site])
    out <- optimize_glo(dna, avoid_sites = sites)
    expect_equal(nrow(find_enzyme_site_ranges(out, sites)), 0, info = pr)
    got <- sum(glo_window_s100(glo_codon_index(substring(out, seq(1, 3 * length(p) - 2, 3), seq(3, 3 * length(p), 3))), 0:(length(p) - 4)))
    expect_equal(got, best_free, info = pr)
    expect_identical(protein_of(out), pr)
  }
})

test_that("avoid_sites makes a gene free of the sites at a small cost, and does nothing when there are none", {
  helpers <- testthat::test_path("..", "..", "dev", "benchmarks", "benchmark_helpers.R")
  skip_if(!file.exists(helpers), "benchmark helpers are not available")
  source(helpers, local = TRUE)
  g <- make_benchmark_cds(3000)
  sites <- enzyme_search_sites(c("BsaI", "BsmBI", "BbsI", "SapI", "EcoRI", "XhoI"))
  set.seed(1); plain <- optimize_glo(g); set.seed(1); free <- optimize_glo(g, avoid_sites = sites)
  expect_gt(nrow(find_enzyme_site_ranges(plain, sites)), 0)
  expect_equal(nrow(find_enzyme_site_ranges(free, sites)), 0)
  expect_identical(protein_of(free), protein_of(g))
  expect_lt(calc_glo_score(free), calc_glo_score(plain))                       # avoiding sites can only cost score
  expect_gt(calc_glo_score(free), 0.95 * calc_glo_score(plain))
  set.seed(1)
  expect_identical(optimize_glo(g, avoid_sites = "GGTCTCNNNNNNNNNNNNNNNNNNNNNNNN"), plain)   # a site that cannot occur changes nothing
})

# ---- GLO-aware shields ------------------------------------------------------------------------------------------------------

test_that("glo_change_gain equals the change of the whole-sequence GLO score", {
  g <- back_translate(original_glo_proteins[["random"]])
  set.seed(4); seq_vec <- strsplit(optimize_glo(g), "")[[1]]
  cod <- function(v) paste(v, collapse = "")
  syn <- glo_synonym_indices(); a <- glo_codon_amino()
  for (ci in c(1, 2, 4, 20, 75, 149, 150)) {
    old <- cod(seq_vec[((ci - 1) * 3 + 1):(ci * 3)])
    alt <- setdiff(names(a)[a == a[old]], old)[1]
    if (is.na(alt)) next
    trial <- seq_vec; trial[((ci - 1) * 3 + 1):(ci * 3)] <- strsplit(alt, "")[[1]]
    expect_equal(round(glo_change_gain(seq_vec, ci, alt)),
                 round(100 * (calc_glo_score(cod(trial)) - calc_glo_score(cod(seq_vec)))), info = paste("codon", ci))
  }
})

test_that("the GLO shields remove the sites and the piRNA matches and cost less GLO score than the Ubiquitous-weighted repairs", {
  helpers <- testthat::test_path("..", "..", "dev", "benchmarks", "benchmark_helpers.R")
  skip_if(!file.exists(helpers), "benchmark helpers are not available")
  source(helpers, local = TRUE)
  g <- make_benchmark_cds(3000); set.seed(1); ex <- optimize_glo(g)
  enz <- c("BsaI", "BsmBI", "BbsI", "SapI", "EcoRI", "XhoI"); sites <- enzyme_search_sites(enz)
  old <- benchmark_pirna_targets()
  glo_p <- remove_pirna_sites(ex, "GLO", 4, pirnas = old)$sequence
  glo_e <- remove_enzyme_sites(glo_p, enz, "GLO")
  set.seed(2)
  ubi_p <- remove_pirna_sites(ex, "Ubiq", 4, pirnas = old)$sequence
  ubi_e <- remove_enzyme_sites(ubi_p, enz, "Ubiq")
  expect_equal(nrow(find_enzyme_site_ranges(glo_e, sites)), 0)
  expect_true(all(get_pirna_distances(glo_p, old) >= 4, na.rm = TRUE))
  expect_identical(protein_of(glo_e), protein_of(g))
  expect_gt(calc_glo_score(glo_e), calc_glo_score(ubi_e))
  expect_gt(calc_glo_score(glo_p), calc_glo_score(ubi_p))
})

test_that("the whole pipeline with GLO, enzymes and the piRNA shield gives a clean sequence with a high GLO score", {
  helpers <- testthat::test_path("..", "..", "dev", "benchmarks", "benchmark_helpers.R")
  skip_if(!file.exists(helpers), "benchmark helpers are not available")
  source(helpers, local = TRUE)
  g <- make_benchmark_cds(3000)
  enz <- c("BsaI", "BsmBI", "BbsI", "SapI", "EcoRI", "XhoI")
  set.seed(1)
  r <- optimize_transgene(g, opt_method = "GLO", enzymes = enz, remove_pirnas = TRUE, target_min_hamming = 4)
  expect_equal(nrow(find_enzyme_site_ranges(r$optimized_cds, enzyme_search_sites(enz))), 0)
  expect_true(all(r$final_pirna_dists >= 4, na.rm = TRUE))
  expect_identical(protein_of(r$optimized_cds), protein_of(g))
  expect_gt(calc_glo_score(r$optimized_cds) / (nchar(g) / 3), 500)                # the exact optimum is about 530 per codon on this gene
})
