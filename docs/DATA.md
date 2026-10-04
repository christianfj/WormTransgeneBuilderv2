# Data tables

Everything the app uses lives in `inst/extdata/` (human-readable source tables) and `data/` (the same
information in R's binary format, built by the scripts in `dev/`). The **Source and licence** column must be
confirmed by the authors before the first release: the entries marked *to be completed* are not recorded
anywhere in the code, so they are not guessed here.

| Table | What it contains | Built into | Source and licence |
|---|---|---|---|
| `Codon_frequencies.csv` | One row per codon (64): amino acid and the weight columns listed below | `data/codon_freqs.rda` (via `inst/extdata/codon_data.R`) | `Ubiq`: labelled "Ubiquitous expression (Serizay)" in the app; `CAI1`: "High expression (Redemann et al., 2011)". Full references and licence *to be completed*. |
| `piRNA_Database.csv` | 17,849 *C. elegans* piRNAs: name, 20-base target sequence (reverse complement of the piRNA without its first base), 21-base piRNA, abundance, enrichment | `data/pirna_sequences.rda` (target sequences, same order) | Source publication and licence *to be completed*. 301 entries without a target sequence were removed on 2026-09-30 (they are in no other piRNA list and not on WormBase). |
| `glo_scores.rds` | A score for each of the 4^12 possible 12-mers, in alphabetical order | read directly (12 MB) | The GLO scores of Dickinson et al. (PubMed 30109984). The raw score file is not in the repository; `dev/make_glo_scores.R` documents how the compact table was made from the earlier object, and the tests keep 195 values read from the original as a check. Permission to redistribute *to be completed*. |
| `Enzymes.csv` | 234 restriction enzymes: name, recognition site (IUPAC codes allowed), `GG` (used for Golden Gate cloning), `Common` (commonly avoided in transgenes) | `data/enzy.rda` (`dev/make_enzy.R`) | New England Biolabs' enzyme list (the full catalogue supplied by the authors, 2026). Methylation sensitivity is not recorded. |
| `Promoters.csv` | 38 validated promoters (300 bp) in two groups (single neuron, tissue-specific) with their WormBase gene | read directly | Authors' own validated promoters. Rows with an empty sequence are group headings. `dev/make_promoters.R` converts the authors' tab-separated file. |
| `5pUTRs.csv`, `3pUTRs.csv` | One 5' UTR (Fire-lab synthetic spliced) and three 3' UTRs (*rps-1*, *rps-4*, *tbb-2*) | read directly | *To be completed.* |
| `Introns.csv` | Three native introns each for *rps-0*, *rps-5*, *rps-20* and the Fire-lab introns | read directly | *To be completed.* |
| `Synthetic_Introns.csv` | 1,000 synthetic introns | read directly | *To be completed* (how generated). |
| `Tags.csv` | SV40 NLS, *egl-13* NLS (Lyssenko et al., 2007), flexible linkers (one cited as Waldo et al., 1999) | read directly | As cited in the table; *to be completed*. |

## Codon table columns

`Codon`, `Amino` (one letter; `X` = stop). The weights are per amino acid:

* `Ubiq`: relative frequency of each synonymous codon (sums to 1 for each amino acid); used for the Ubiquitous method.
* `Chance`: equal weight for every synonymous codon; used for Random recoding.
* `CAI1`: 1 for the single preferred codon of each amino acid and 0 for the others; used for High expression.
* `BringMansWeights`: graded weights (largest = 1) used to compute the Codon Adaptation Index of the High-expression method and to rank changes by their CAI loss.
* `Germ`, `Neur`, `Soma`, `aCDS`, `HighRCSU`: present in the source file; `Germ`, `Neur`, `Soma` are not stored in the package data; `aCDS` and `HighRCSU` are stored but not used by any method in the app.

## Rebuilding

* `Rscript dev/make_enzy.R` rebuilds `data/enzy.rda` from `Enzymes.csv` (and checks its format).
* `Rscript dev/make_promoters.R file.tsv` rebuilds `Promoters.csv` from a tab-separated list.
* `data/glo_median_score.rda` is the median of all GLO scores; it is used as the score of a 12-mer that contains an `N`.
* `data/codon_freqs.rda`, `pirna_sequences.rda` and `glo_median_score.rda` were built by `inst/extdata/codon_data.R` and its companions from source files that are not all kept in the repository.

The automated tests (`tests/testthat/test-data_files.R`) check that the tables and the built data files agree.
