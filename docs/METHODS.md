# Methods

What each step of WormTransgeneBuilder does. The code for each step is named, so a statement here can be
checked against the program. Sources of the data tables are in [DATA.md](DATA.md).

*Beta release: the manuscript is in preparation. Full literature references will be added when it is
available; sources are named below as they appear in the program.*

## 1. Input and order of the steps

The input is a coding sequence (DNA, which must start with ATG, end with a stop codon, be a multiple of 3
long and have no internal stop) or a protein sequence, up to 20,000 characters. A protein is first
back-translated with one fixed codon per amino acid (`back_translate()`) and a stop codon is added if it
is missing; the next steps then recode it.

`optimize_transgene()` (R/fct_optimize_transgene.R) applies the steps in this order:

1. add N- and C-terminal tags (N-terminal tags directly after the ATG, C-terminal tags directly before the stop codon);
2. record the distance to piRNAs of the starting sequence;
3. recode codons (section 2), optionally optimizing the ribosome binding site (section 5);
4. remove piRNA homology (section 3);
5. remove restriction sites (section 4);
6. insert introns (section 6);
7. add the consensus start, 5' UTR, promoter and 3' UTR.

Every recoding step replaces a codon by a **synonymous** codon, so the encoded protein is the same as the
input. This is checked by the automated tests.

## 2. Codon optimization

Codon tables are in `inst/extdata/Codon_frequencies.csv` (columns in [DATA.md](DATA.md)). In the app:

| Menu entry | Internal name | What is done |
|---|---|---|
| No codon optimization | `None` | codons are left as given |
| Ubiquitous expression (Serizay) | `Ubiq` | each codon is replaced by a synonymous codon drawn at random with probability proportional to its frequency in the `Ubiq` column |
| Random recoding | `Chance` | each codon is replaced by a synonymous codon drawn with equal probability |
| High expression (Redemann et al., 2011) | `CAI1` | each codon is replaced by the single preferred codon of its amino acid (weight 1 in the `CAI1` column); no random choice, so the result is always the same |
| Germline optimized (GLO) | `GLO` | see below |

(`probabilistic_recode()`)

**Germline optimization (GLO).** The GLO score of a sequence is the sum, over all 12-base windows that start
at a codon boundary and advance one codon at a time, of an empirical score for that 12-mer (Dickinson et
al.; PubMed 30109984). The table has a score for each of the 4^12 = 16,777,216 possible 12-mers
(`inst/extdata/glo_scores.rds`; a window's score is found from its bases read as a base-4 number,
`R/fct_glo_scores.R`). `optimize_glo()` is a greedy local search: codons are visited in random order; a
random synonymous alternative is tried; it is kept if it raises the summed score of the (up to four)
12-mer windows that overlap that codon. This is repeated for up to five passes or until a pass changes
nothing. Because choices are random, two runs give different (similarly scoring) sequences.

**Codon Adaptation Index (CAI).** Reported with the result: the geometric mean of the per-codon weights
(`calc_cai()`); weights are `BringMansWeights` for the High-expression method and `Ubiq` otherwise.

## 3. piRNA homology removal

The reference is a set of 17,849 *C. elegans* piRNA target sequences of 20 bases
(`inst/extdata/piRNA_Database.csv`). For the designed sequence, every 20-base window is compared with every
piRNA and the **Hamming distance** (number of mismatching positions) is computed exactly
(`get_pirna_distances()`; the comparison uses 2-bit-encoded bases, XOR and a lookup table, which gives the
same numbers as comparing letters but is much faster, `R/fct_pirna_engine.R`). Windows are on the given
strand; the user chooses the minimum number of mismatches to be considered safe (0 to 4; default 3).

Windows closer than the threshold are repaired by synonymous recoding (`remove_pirna_sites()`):

* **Ubiquitous, Random, GLO and no optimization:** every codon under a flagged window is re-drawn once at
  random from its other synonymous codons (weights from the `Ubiq` column; equal weights for Random), and
  the distances of the final sequence are recomputed. This is a single pass, so a window can still be
  below the threshold afterwards; the distances of the final sequence are returned with the result.
* **High expression:** a deterministic search. For the first flagged window, every synonymous change to any
  codon under it is tried; only changes that increase that window's distance and make no other window
  worse are allowed; of those, the change that loses the least CAI (`log(w_new) - log(w_old)` with
  `BringMansWeights`) is applied. Ties are broken by larger gain, earlier position, then alphabetical codon.
  This repeats until no window is under the threshold or no change can help (a warning is then shown).

Only piRNAs in the database are considered; the 301 entries of the original list with no target sequence were removed
(see [DATA.md](DATA.md)).

## 4. Restriction-site removal

The enzyme table (234 enzymes from New England Biolabs, `inst/extdata/Enzymes.csv`) gives each recognition
site, which may contain IUPAC ambiguity codes (for example BanI, `GGYRCC`, or `N` for any base). Both strands
are searched (the site and its reverse complement). A site letter matches the sequence letters whose bases
it includes (`find_enzyme_site_ranges()`).

Sites are removed by synonymous changes (`remove_enzyme_sites()`).
*Random methods:* the codon at the start of each site is re-drawn using the method's weights; if it cannot
change (Met or Trp), the other codons the site covers are tried and a change that lowers the number of
sites is chosen at random (weighted). The search is repeated, for up to 100 passes, until no site is
left or a pass changes nothing. *High expression:* the change with the smallest CAI loss that removes a
site is applied, one site at a time (deterministic). If a site cannot be removed by any synonymous change
(for example `ACNGT` in a Thr-Val pair, where every Thr codon is `ACN` and every Val codon `GTN`), a
warning lists it. Methylation (Dam/Dcm) is not considered: synthesized DNA is not methylated.

## 5. Ribosome binding site (optional)

For the first 39 bases of the coding sequence, 100 alternative synonymous versions are generated by random
recoding with the `Ubiq` weights (or the `CAI1` column for High expression) and added to the input version
(`optimize_rbs()`). For High expression every draw gives the same sequence, so there are only two candidates. Each candidate is prefixed with `AAAA` and folded with
**ViennaRNA** (`RNAfold --noPS`); the candidate with the **highest minimum free energy** (least stable
structure) is kept, to leave the start of the mRNA accessible. The start is written as `aaaaATG`.
If `RNAfold` is not installed this step is skipped with a warning. The version of ViennaRNA used is shown
in the About tab and in the log (reference: Lorenz et al., ViennaRNA Package 2.0, Algorithms Mol Biol 2011).

## 6. Introns

Native introns (rps-0, rps-5, rps-20 and the Fire-lab introns) or 1 to 10 synthetic introns (picked at
random from a list of 1,000 sequences in `Synthetic_Introns.csv`) are inserted (`insert_introns()`).
Placement is "early" (first intron near codon 50, then every 50 codons) or "equidistant" (evenly across the
coding sequence). Each intron goes at the best nearby codon boundary (an intron is always inserted between two codons): within
±15 codons of the target, a preceding codon `AAG` scores 3, `CAG` scores 2, `xAG` scores 1; the highest
score wins, and the nearest wins ties. Introns are written in lower case. A warning is given when the
average exon would be shorter than 150 bp. (The menu has a "Force in reading frame" tick box; at present it
changes nothing, see [LIMITATIONS.md](LIMITATIONS.md).)

## 7. Tags, promoters, UTRs

Tags, promoters and UTRs are taken as they are from the data tables and are not recoded. Promoters are
added upstream of the 5' UTR; the 3' UTR follows the stop codon. In the GenBank file and the viewer each
element is labelled with its **name** only (for example `pADF`); the gene identifier shown in the menu is
not part of the annotation.

## 8. Synthesis checks (optional)

`validate_for_synthesis()` checks the final sequence against rules that follow Twist Bioscience's
guidelines: no homopolymer of 14 or more identical bases; global GC content between 25 % and 75 %; no
50-base window with more than 80 % GC. Passing these checks does not guarantee that a vendor will accept
the sequence; check the vendor's current guidelines.

## 9. Reproducibility of a design

Steps that use random choices (all methods except High expression) give different sequences each time.
For "No codon optimization" the seed is fixed (12345), so the same input gives the same output. The
GenBank file records the settings used in its comment lines. The versions of the software are shown in
the About tab. See [REPRODUCIBILITY.md](REPRODUCIBILITY.md).
