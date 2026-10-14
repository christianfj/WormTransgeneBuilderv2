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
at a codon boundary and advance one codon at a time, of an empirical score for that 12-mer (GLO, GermLine Optimized, developed by D. J. Dickinson and first described and used in Fielmich et al. 2018, eLife 7:e38198, doi:10.7554/eLife.38198; the published design assembles a coding sequence from 12-nucleotide words found in native germline-expressed genes). The table has a score for each of the 4^12 = 16,777,216 possible 12-mers
(`inst/extdata/glo_scores.rds`; a window's score is found from its bases read as a base-4 number,
`R/fct_glo_scores.R`). The sequence designed for a protein is the one with the **highest possible GLO score**, found exactly: the score is
a sum of terms that each involve four consecutive codons, so dynamic programming (as in the Viterbi algorithm) goes along the sequence, keeps for
every choice of the last three codons the best score reached so far, and traces back. No search is needed, so it takes a fraction of a second even for 20 kb. (The
original Perl implementation of D. J. Dickinson, https://github.com/dannyhmg/germline, reaches 0.1 to 2% less than this optimum with a randomized search over whole four-codon
words; that algorithm is kept as `optimize_glo(method = "heuristic")`, and the app reproduces the scores of the unchanged original code exactly,
[GLO_VERIFICATION.md](GLO_VERIFICATION.md).) Many codings share the best score because 82% of 12-mers score zero; among them the one with the **highest Codon Adaptation Index** (with the Ubiquitous codon table,
as the app reports it) is chosen. The CAI is the geometric mean of the codon weights, so this is a maximization of the sum of the log weights, and it can never override
the GLO score (GLO is the first criterion, CAI the second). In practice the best GLO coding is almost unique: with a random jitter on the tie-break two runs differ in about 0.2% of codons. The amino-acid
sequence, stop codons and any codon that is not plain A/C/G/T are never changed. Note that the original reports the score per codon, whereas the app reports the sum over the sequence.

**GLO together with the other steps.** When restriction sites are to be removed, the optimization takes them into account: every window of four codons that contains a site is ruled out, so the best
sequence *without* the sites is found directly (any site of up to 10 bases lies inside one such window; this covers 219 of the 234 enzymes and every commonly used one; the
longer ones, mostly with runs of N, are handled afterwards). The piRNA shield and whatever restriction sites remain are then removed one change at a time, each time choosing the
synonymous change that costs the least GLO score (only the four windows that contain the changed codon are affected), as the "High expression" method does with the Codon Adaptation Index.

**Codon Adaptation Index (CAI).** Reported with the result: the geometric mean of the per-codon weights
(`calc_cai()`); weights are `BringMansWeights` for the High-expression method and `Ubiq` otherwise.

## 3. piRNA homology removal

The *C. elegans* database (`inst/extdata/piRNA_Database.csv`) holds 18,207 entries: 15,363 type 1 piRNAs annotated in WormBase WS298, 2,485 type 2 piRNAs and 359 further sequences found in PRG-1 immunoprecipitation (see [PIRNA_DATABASE.md](PIRNA_DATABASE.md)). For each, the 21-base
piRNA (`Full_21mer`, 5' to 3'; 96.8 % begin with U/T) is stored together with its **target sequence**
(`Target_20mer`): the reverse complement of the piRNA without its first base. This target is the 20-base
stretch of the *sense* transcript to which the piRNA, which pairs with its target in the antisense
orientation, would bind. The tests check that this holds for every entry. Because the target is already
written in the sense orientation, windows of the designed (sense) sequence are compared with it directly. For the designed sequence, every 20-base window is compared with every
piRNA and the **Hamming distance** (number of mismatching positions) is computed exactly
(`get_pirna_distances()`; the comparison uses 2-bit-encoded bases, XOR and a lookup table, which gives the
same numbers as comparing letters but is much faster, `R/fct_pirna_engine.R`). Windows are on the given
strand; the user chooses the minimum number of mismatches to be considered safe (0 to 4; **default 4**, see below).

**Which species.** Each species has its own, independent database ([PIRNA_OTHER_SPECIES.md](PIRNA_OTHER_SPECIES.md)): *C. elegans* (default; the type 2 piRNAs can be switched off, and are on by default), *C. briggsae* (PRG-1 IP, high confidence), and *C. remanei* and *C. brenneri* (called from small-RNA sequencing, low confidence). The user ticks every species the transgene will be used in, for example *C. elegans* and *C. briggsae*; the windows are then compared with the union of the selected databases, and a window is flagged if it is within the threshold of any piRNA of any selected species. Nothing is combined unless the user selects it. The 4-mismatch default comes from experiments in *C. elegans* and is applied unchanged to the other species. Codon optimization, GLO scoring and the other steps always use *C. elegans* data: there is no codon-usage or GLO data for the other species, so only the piRNA screen is species-specific.

**Why mismatches are counted equally (and how this differs from pirScan).** The app scores a window by its plain Hamming
distance to a piRNA target and does not weight the seed (nucleotides 2 to 7) or wobble pairs differently. This is
deliberate. Priyadarshini et al. (2022, Nature Methods 19:187-194) tested synthetic guide piRNAs expressed from transgenes
that silence the endogenous gene *him-5* (male frequency as read-out). Silencing was relatively insensitive to one, two or
three mismatches in each of the six guide piRNAs and was abolished by four or five mismatches, "regardless of where
mismatches were located" (their Fig. 2a and Extended Data Fig. 5). In their regressions the number of mismatches explained
the silencing about as well as, or slightly better than, the position- and wobble-aware score of pirScan (R2 = 0.34 and 0.38
for the number of mismatches, against 0.30 and 0.35 for the pirScan-based score, for all and for the four best-expressed guide
piRNAs). Their web tool for designing piRNA-interference transgenes (www.wormbuilder.org/piRNAi) therefore uses the
Hamming distance, with four mismatches as its stringent setting. For the same reason the threshold here is a count of
mismatches: **a window with at least 4 mismatches to every piRNA is not expected to be targeted**, while windows with up to
3 mismatches may be.

**Scope of that evidence.** It comes from one target gene and from guide piRNAs expressed in extrachromosomal arrays, with the
leading U left unchanged, so it is evidence about the tolerance of piRNA-mediated silencing and is not a complete model of
piRNA targeting. Where users need the finer details for individual piRNAs, the position of each mismatch, seed pairing and
G:U pairs, they should use **pirScan** (Wu et al. 2018; http://cosbi4.ee.ncku.edu.tw/pirScan/), which is built for that. The
two tools can disagree for a given site: pirScan may predict a site that has 3 or fewer mismatches but a disrupted seed to
be weak, or a site with more mismatches but a perfect seed to be active. This app is a general-purpose design tool for
transgenes (codon optimization, shields, parts and synthesis checks together) and screens against a larger, more
conservative set of piRNAs (type 1, optional type 2 and the newly identified PRG-1-bound sequences), so its flags are not
meant to reproduce pirScan's predictions. **The default threshold is 4 mismatches**, chosen to match our experimental data: in Priyadarshini et al. four
mismatches were enough to abolish silencing while up to three were tolerated, so a design in which no 20-base window is
within 3 mismatches of any piRNA is not expected to be targeted. Billi et al. (2013) and Bagijn et al. (2012) likewise
describe 21U RNA targets as having up to three mismatches. Lower settings (3, 2, 1, 0) are available for users who want
fewer sequence changes, at the price of leaving windows that a piRNA could in principle target. A stricter setting
means more codons are changed in the shield, which can lower the codon-optimization score a little.

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

The enzyme table (234 enzymes sold by New England Biolabs according to REBASE (see `docs/DATA.md`), `inst/extdata/Enzymes.csv`) gives each recognition
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
recoding (`optimize_rbs()`). For Ubiquitous, Random and GLO the draws use the `Ubiq` weights and the input
version is a further candidate. For **High expression**, which has a single best codon per amino acid (so
all draws would be identical), the draws use the graded `BringMansWeights`: frequent codons stay likely, the
candidates differ, and the input codons are not a candidate, so the start is always recoded. These start
codons need not all be the best codon, so the viewer underlines them as non-optimal. If `RNAfold` is
missing, the High-expression start falls back to the single best-codon version. Each candidate is prefixed with `AAAA` and folded with
**ViennaRNA** (`RNAfold --noPS`); the candidate with the **highest minimum free energy** (least stable
structure) is kept, to leave the start of the mRNA accessible. The start is written as `aaaaATG`.
If `RNAfold` is not installed this step is skipped with a warning. The version of ViennaRNA used is shown
in the About tab and in the log (reference: Lorenz et al., ViennaRNA Package 2.0, Algorithms Mol Biol 2011).

## 6. Introns

Native introns (rpsa-1, rps-5, rps-20 and the Fire-lab introns) or 1 to 10 synthetic introns (picked at
random from a list of 1,000 sequences in `Synthetic_Introns.csv`, made by `dev/make_synthetic_introns.R`) are inserted (`insert_introns()`). The Fire-lab
introns and the 5' UTR with intron A are described in [FIRE_LAB_SEQUENCES.md](FIRE_LAB_SEQUENCES.md).
Placement is "early" (first intron near codon 41, then every 41 codons: one median exon of 123 bp) or "equidistant" (evenly across the
coding sequence). Each intron goes at the best place within 15 codons (45 bases) either side of its target. A place scores 3
if the three bases before the intron are `AAG`, 2 if they are `CAG` and 1 if they end in `AG`; the highest
score wins, and the place nearest the target wins ties. Introns are written in lower case. A warning is
given when the average exon would be shorter than 110 bp (a little under the median exon of 123 bp). The spacing and the warning follow the median, not the average, because the average is skewed by a few very large exons. The reference for exon and intron sizes is Spieth J, Lawson D, Davis P, Williams G, Howe K (2014), Overview of gene structure in *C. elegans*, WormBook, doi:10.1895/wormbook.1.65.2 (PMID 25368915; open access, CC BY): in WormBase release WS237 the average size of unique coding exons is 200.7 bp (median 123 bp, skewed by a few very large exons), the most common intron size is 47 bp (median 65 bp), and short introns are preferentially found in highly expressed genes. The text under the number of synthetic introns in the app quotes the median of 123 bp and the average of 200 bp.

**Restriction sites.** When enzymes are selected for removal, the inserted introns must not create a recognition site, either inside the intron or where it joins
the exon (the stretch checked is the intron with up to one base less than the longest site on each side). A synthetic intron that would create one is replaced by a random
spare one from the list of 1,000 that does not. The fixed intron sets (*rpsa-1*, *rps-5*, *rps-20*, Fire lab) cannot be changed: if one of them creates a site the app
gives a warning and leaves it.

The tick box **Force in reading frame** (ticked by default) decides where an intron may go. Ticked: only
between two codons (the intron follows a complete codon, phase 0), so the three bases before the intron are
a codon. Not ticked: after any base, so an intron can fall inside a codon (phase 1 or 2) and the
`AAG`/`CAG` motif can straddle a codon boundary. The first intron is at least 15 bases from the start and
every intron at least 15 bases from the end; introns are at least 30 bases apart.

## 7. Tags, promoters, UTRs

Tags, promoters and UTRs are taken as they are from the data tables and are not recoded. Promoters are
added upstream of the 5' UTR; the 3' UTR follows the stop codon. In the GenBank file and the viewer each
element is labelled with its **name** only (for example `pADF`); the gene identifier shown in the menu is
not part of the annotation.

## 8. Synthesis checks (optional)

`validate_for_synthesis()` checks the final sequence (including promoter, introns and UTRs) against Twist Bioscience's published gene synthesis guidelines, **as read on 2026-10-14**
from Twist's web pages (these change; check them before ordering):

* *Are there any sequence limitations/design guidelines for genes which I should follow?* (https://www.twistbioscience.com/faq/gene-synthesis/are-there-any-sequence-limitationsdesign-guidelines-genes-which-i-should-follow)
* *Why don't complex sequences qualify for the Express service?* (https://www.twistbioscience.com/faq/gene-synthesis/why-dont-complex-sequences-qualify-express-service)
* *Complex Gene Synthesis Services* (https://www.twistbioscience.com/products/genes/complex-genes)

Twist scores each design as Low, Medium or High complexity and decides itself. The app reports two levels separately.

**Standard guidelines** (Twist's rules for ordinary genes; a message means a problem may arise):

| Check | Rule |
|---|---|
| Homopolymer | 14 or more identical bases (Twist: avoid homopolymers of 14 bp or more) |
| Global GC content | outside 25 to 75 % (the authors' decision, 2026-10-14; the pages read disagreed on the upper limit) |
| Local GC | a 50-base window with more than 80 % GC (a softer warning than Twist's own limit of 90 %) |
| Length | below 300 bp; above 5,000 bp (the maximum for gene fragments); above 7,000 bp (the maximum for clonal genes) |
| ccdB | a coding region for ccdB (Twist: do not include ccdB). The sequence is the full ccdB coding sequence (101 amino acids), supplied by the authors; any stretch of 30 of its amino acids found in a translation of the sequence (all six frames) is reported, so a recoded copy is found too |

**High complexity** (not eligible for the Express service; the design would be made as a Complex Gene, which takes longer):

| Check | Rule |
|---|---|
| Homopolymer | any run of more than 30 identical bases |
| Local GC | a 50-base window below 10 % or above 90 % GC |
| Direct repeats | a stretch of more than 200 bases that occurs twice (non-overlapping, same strand) |

**Not checked:** short tandem repeats longer than 100 bases, hairpin-forming repeats longer than 100 bases, and Twist's machine-learning complexity score (which also uses the overall GC, the longest repeat, the length and the repeat density). Passing the checks does not
guarantee that Twist will accept or quickly make the sequence.

*Changed on 2026-10-14.* Earlier versions checked only a homopolymer of 14 or more, global GC and local GC above 80 % (as a failure), and, because of an error, ignored all lower-case bases (introns, UTRs and promoters were not checked).

## 9. Reproducibility of a design

Steps that use random choices (all methods except High expression) give different sequences each time.
For "No codon optimization", and for "High expression" with the ribosome-binding-site step, the seed is
fixed (12345), so the same input gives the same output. The
GenBank file records the settings used in its comment lines. The versions of the software are shown in
the About tab. See [REPRODUCIBILITY.md](REPRODUCIBILITY.md).
