# 1. Load the raw laboratory CSVs
# (Codon_frequencies.csv has exactly the columns the app uses: Codon, Amino, Ubiq, Chance, CAI1, BringMansWeights. Five more columns,
# Germ, Neur, Soma, aCDS and HighRCSU, were removed on 2026-10-07 because nothing used them; they are in the git history, and
# docs/DATA.md says how they were made.)
codon_freqs <- read.csv("inst/extdata/Codon_frequencies.csv", stringsAsFactors = FALSE)

# 2. CRUCIAL: Set row names to the codons themselves
# This allows the optimizer to do fast lookups like: codon_freqs["ATG", "Amino"]
rownames(codon_freqs) <- toupper(as.character(codon_freqs$Codon))

# 3. Save to the data/ folder for the package
# (Two more datasets used to be built here, `intron_seqs` from Introns.csv and
# `aa_to_codon_f`, a list of codons per amino acid. Nothing used them, and the app reads
# introns straight from inst/extdata/Introns.csv, so they were removed. They can be
# recreated from git history or from these two sources if they are ever needed.)
usethis::use_data(codon_freqs, overwrite = TRUE)

## --- piRNA Data Processing ---
# Load the raw piRNA list
pirna_raw <- read.csv("data-raw/piRNA_sequence.csv", header = FALSE, stringsAsFactors = FALSE)

# The legacy code extracts the 20-mer sequences (column 2) to use for the string distance math
pirna_sequences <- as.character(pirna_raw[, 2])

# Drop entries without a sequence. 301 piRNAs had NA here; they are in no other piRNA list and
# not on WormBase, so they are not validated (and a missing sequence can never be matched).
# piRNA_Database.csv in inst/extdata has the same 301 rows removed, so the two stay in step.
pirna_sequences <- pirna_sequences[!is.na(pirna_sequences)]

# Save to the package data
usethis::use_data(pirna_sequences, overwrite = TRUE)

## --- Enzyme Data Processing ---
# data/enzy.rda is no longer built here. It is built from inst/extdata/Enzymes.csv (New England
# Biolabs' full list) by dev/make_enzy.R, which also checks the file.
