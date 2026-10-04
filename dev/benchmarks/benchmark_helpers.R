# Helpers for the benchmark script (benchmark.R). They live in their own file so the
# test suite can check that the benchmark inputs are valid genes.
#
# Everything here is deterministic: the same call always gives the same sequence, so
# timings from different days or different versions of the code are comparable.

benchmark_sizes <- c(small = 300, medium = 3000, large = 20000)  # bp; 20 kb is the app's input limit

# Run code with a fixed random seed and leave R's random state exactly as it was.
bench_with_seed <- function(seed, code) WormTransgeneBuilder:::with_local_seed(seed, code)

# A stable fingerprint of a result, used to tell whether an optimization changed the
# output. rlang::hash() ignores R's serialization header, so it does not depend on the
# R version.
fingerprint <- function(x) rlang::hash(x)

codons_of <- function(dna) substring(dna, seq(1, nchar(dna) - 2, 3), seq(3, nchar(dna), 3))

# piRNA 20-mers that can sit in frame inside a coding sequence: 20 bases plus one padding
# base (7 codons) with no stop codon. `best_only = TRUE` keeps only those made entirely of
# best (CAI1) codons, which survive "High expression" recoding intact.
pirna_blocks <- function(best_only = FALSE, limit = Inf) {
  cf <- WormTransgeneBuilder::codon_freqs
  stops <- c("TAA", "TAG", "TGA")
  best <- toupper(cf$Codon[cf$CAI1 == 1 & cf$Amino != "X"])
  p <- WormTransgeneBuilder::pirna_sequences
  p <- unique(p[!is.na(p)])
  blocks <- character(0)
  for (x in p) {
    for (pad in c("A", "C", "G", "T")) {
      cd <- codons_of(paste0(x, pad))
      if (any(cd %in% stops)) next
      if (best_only && !all(cd %in% best)) next
      blocks <- c(blocks, paste0(x, pad))
      break
    }
    if (length(blocks) >= limit) break
  }
  blocks
}

# Where to plant things in a gene body of `n_body` codons: `n_pirna` piRNA targets and
# `n_sites` restriction sites, spread evenly and interleaved so that none can overlap
# (neighbouring slots are far more than one 7-codon target apart).
plant_positions <- function(n_body, n_pirna, n_sites) {
  total <- n_pirna + n_sites
  slots <- round(seq(10, n_body - 20, length.out = total + 2))[-c(1, total + 2)]
  is_pirna <- rep(FALSE, total)
  is_pirna[unique(round(seq(1, total, length.out = n_pirna)))] <- TRUE
  list(pirna = slots[is_pirna], sites = slots[!is_pirna])
}

# A reproducible gene of about `n_bp` bases: ATG, random sense codons, TAA. Codons are
# drawn with the Ubiq usage weights. Real piRNA targets (one per ~900 bp, at least one)
# and restriction sites (BsaI and XhoI, one per ~1500 bp, at least one each) are planted
# in frame so that the shields have work to do.
make_benchmark_cds <- function(n_bp, seed = 2026) {
  cf <- WormTransgeneBuilder::codon_freqs
  sense <- cf[cf$Amino != "X", ]
  n_codons <- max(40, n_bp %/% 3)
  codons <- bench_with_seed(seed, {
    body <- sample(toupper(sense$Codon), n_codons - 2, replace = TRUE, prob = sense$Ubiq)
    blocks <- pirna_blocks(limit = 400)
    n_pirna <- max(1, n_codons %/% 300)
    picked <- blocks[sample.int(length(blocks), n_pirna)]
    list(body = body, picked = picked)
  })
  body <- codons$body
  n_sites <- max(1, n_codons %/% 500)
  pos <- plant_positions(length(body), length(codons$picked), n_sites)
  for (i in seq_along(pos$pirna)) {
    cd <- codons_of(codons$picked[i])
    body[pos$pirna[i]:(pos$pirna[i] + length(cd) - 1)] <- cd
  }
  site_codons <- rep(list(c("GGT", "CTC"), c("CTC", "GAG")), length.out = n_sites)  # BsaI, XhoI
  for (i in seq_along(pos$sites)) body[pos$sites[i]:(pos$sites[i] + 1)] <- site_codons[[i]]
  paste0("ATG", paste(body, collapse = ""), "TAA")
}

# The same idea for the High-expression shields: every codon is the best codon, so the
# gene is unchanged by CAI1 recoding, and the planted piRNA targets and XhoI sites are
# made only of best codons so they survive it.
make_best_codon_benchmark_cds <- function(n_bp, seed = 2026) {
  cf <- WormTransgeneBuilder::codon_freqs
  best <- toupper(cf$Codon[cf$CAI1 == 1 & cf$Amino != "X"])
  n_codons <- max(40, n_bp %/% 3)
  blocks <- pirna_blocks(best_only = TRUE)
  body <- bench_with_seed(seed, sample(best, n_codons - 2, replace = TRUE))
  n_pirna <- max(1, n_codons %/% 300)
  n_sites <- max(1, n_codons %/% 500)
  pos <- plant_positions(length(body), n_pirna, n_sites)
  for (i in seq_along(pos$pirna)) {
    cd <- codons_of(blocks[(i - 1) %% length(blocks) + 1])
    body[pos$pirna[i]:(pos$pirna[i] + length(cd) - 1)] <- cd
  }
  for (i in seq_along(pos$sites)) body[pos$sites[i]:(pos$sites[i] + 1)] <- c("CTC", "GAG")
  paste0("ATG", paste(body, collapse = ""), "TAA")
}

# ---- what gets timed ------------------------------------------------------------------

# Each case has an id, the sizes it runs at, and a function of `x`, a list holding the
# gene for that size: x$cds (random codons) and x$best (all best codons).
benchmark_cases <- function() {
  pkg <- asNamespace("WormTransgeneBuilder")
  load_csv <- get("safe_load_csv", pkg)
  # The 21 enzymes the table held before it became NEB's full list of 234; the "all21" cases keep
  # using exactly these so their results stay comparable with earlier runs. "Common" is the set the
  # laboratory avoids in a transgene (BanIII is not an NEB enzyme and is left out); removing all 234
  # enzymes makes no biological sense, so that is not benchmarked.
  all_enzymes <- c("ApaI", "ApaLI", "BamHI", "BanI", "BsaI", "BsmBI", "BsoBI", "DraI", "EcoRI", "EcoRV",
                   "HindIII", "KpnI", "NcoI", "NdeI", "NheI", "PstI", "PvuII", "SacII", "SapI", "XbaI", "XhoI")
  common_enzymes <- c("ApaI", "ApaLI", "AvaI", "AvaII", "BamHI", "BanI", "BanII", "BclI", "BsaI", "BsmBI",
                      "DraI", "EcoRI", "EcoRV", "HindIII", "KpnI", "NciI", "NdeI", "PmlI", "PstI", "PvuII",
                      "SacII", "SmaI", "XbaI", "XhoI")
  introns3 <- load_csv("Synthetic_Introns.csv")$Sequence[1:3]
  tags <- trimws(load_csv("Tags.csv")$Sequence)
  promoter <- Filter(nzchar, load_csv("Promoters.csv")$Sequence)[1]  # first real promoter (rows without a sequence are group headings)
  utr5 <- load_csv("5pUTRs.csv")$Sequence[2]
  utr3 <- load_csv("3pUTRs.csv")$Sequence[2]
  std <- c("small", "medium")
  all3 <- c("small", "medium", "large")

  pipe <- function(x, method, ...) {
    args <- list(dna_seq = x$cds, opt_method = method, enzymes = "BsaI", remove_pirnas = TRUE, target_min_hamming = 3)
    extra <- list(...)
    args[names(extra)] <- extra
    do.call(optimize_transgene, args)
  }
  features <- function(n) {
    st <- round(seq(1, n - 30, length.out = 10))
    list(starts = st, ends = st + 20, colors = rep("#a9dfbf", 10),
         tooltips = paste("Feature", 1:10), feats = rep("misc_feature", 10))
  }

  list(
    list(id = "calc_gc_content", sizes = all3, fn = function(x) calc_gc_content(x$cds)),
    list(id = "calc_cai", sizes = all3, fn = function(x) calc_cai(x$cds, "Ubiq")),
    list(id = "calc_glo_score", sizes = all3, fn = function(x) calc_glo_score(x$cds)),
    list(id = "validate_for_synthesis", sizes = all3, fn = function(x) validate_for_synthesis(x$cds)),

    list(id = "probabilistic_recode_ubiq", sizes = all3, fn = function(x) probabilistic_recode(x$cds, "Ubiq")),
    list(id = "optimize_glo", sizes = all3, fn = function(x) optimize_glo(x$cds)),

    list(id = "get_pirna_distances", sizes = all3, fn = function(x) get_pirna_distances(x$cds)),
    list(id = "pirna_shield_random", sizes = std, fn = function(x) remove_pirna_sites(x$cds, "Ubiq", 3)),
    list(id = "pirna_shield_cai1", sizes = std, fn = function(x) remove_pirna_sites(x$best, "CAI1", 3)),

    list(id = "enzymes_random_BsaI", sizes = all3, fn = function(x) remove_enzyme_sites(x$cds, "BsaI", "Ubiq")),
    list(id = "enzymes_random_all21", sizes = std, fn = function(x) remove_enzyme_sites(x$cds, all_enzymes, "Ubiq")),
    list(id = "enzymes_cai1_XhoI", sizes = std, fn = function(x) remove_enzyme_sites(x$best, "XhoI", "CAI1")),
    list(id = "enzymes_cai1_all21", sizes = std, fn = function(x) remove_enzyme_sites(x$best, all_enzymes, "CAI1")),
    list(id = "enzymes_random_common24", sizes = std, fn = function(x) remove_enzyme_sites(x$cds, common_enzymes, "Ubiq")),
    list(id = "enzymes_cai1_common24", sizes = std, fn = function(x) remove_enzyme_sites(x$best, common_enzymes, "CAI1")),

    list(id = "insert_introns_early", sizes = all3, fn = function(x) insert_introns(x$cds, tolower(introns3), "early")),
    list(id = "insert_introns_equidistant", sizes = all3, fn = function(x) insert_introns(x$cds, tolower(introns3), "equidistant")),

    list(id = "optimize_rbs_100_candidates", sizes = "small", fn = function(x) optimize_rbs(substr(x$cds, 1, 39), "Ubiq")),

    list(id = "generate_sequence_viewer", sizes = all3, fn = function(x) {
      f <- features(nchar(x$cds))
      as.character(generate_sequence_viewer("v", x$cds, f$starts, f$ends, f$colors, f$tooltips))
    }),
    list(id = "export_genbank", sizes = all3, fn = function(x) {
      f <- features(nchar(x$cds))
      export_genbank("bench", x$cds, f$starts, f$ends, f$colors, f$tooltips, f$feats, 0.9, "50%", "1,000", x$cds, "Benchmark")
    }),

    list(id = "pipeline_None", sizes = std, fn = function(x) pipe(x, "None")),
    list(id = "pipeline_Ubiq", sizes = std, fn = function(x) pipe(x, "Ubiq")),
    list(id = "pipeline_Chance", sizes = std, fn = function(x) pipe(x, "Chance")),
    list(id = "pipeline_CAI1", sizes = std, fn = function(x) pipe(x, "CAI1")),
    list(id = "pipeline_GLO", sizes = std, fn = function(x) pipe(x, "GLO")),
    list(id = "pipeline_Ubiq_all_options", sizes = std, fn = function(x) {
      pipe(x, "Ubiq", enzymes = c("BsaI", "XhoI"), introns = tolower(introns3), tag_n_seqs = tags[1],
           tag_c_seqs = tags[4], promoter_seq = promoter, utr5_seq = utr5, utr3_seq = utr3,
           add_consensus_start = TRUE, rbs_opt = Sys.which("RNAfold") != "")
    })
  )
}
