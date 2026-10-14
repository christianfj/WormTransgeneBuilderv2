# GLO: comparison with the original implementation

*Written 2026-10-07; the problem described under "The optimization" was fixed the same day (see "After the fix" at the end). The GLO (GermLine Optimized) method of D. J. Dickinson was first described and used in Fielmich et al. 2018 (eLife 7:e38198,
doi:10.7554/eLife.38198). Its code is at https://github.com/dannyhmg/germline (Perl, MIT licence; frozen copy at
https://github.com/elifesciences-publications/germline). We compared it with this app. Scripts: `dev/verify_glo_against_original.R`,
`dev/glo_original_perl/`.*

## What was checked

| Check | Result |
|---|---|
| **The score table.** Every 12-mer's score in `inst/extdata/glo_scores.rds` against the original `data/sequence_lib_scores.db` (737 MB, dumped with `db_dump`) | **Identical**, all 16,777,216 values. 82% of the 12-mers score 0; the maximum is 46,995.8. |
| **The scoring function.** The original (`Seqscore.pm`) sums, for every fourth codon, seven overlapping 12-mer words with weights 0.25, 0.5, 0.75, 1, 0.75, 0.5, 0.25, and divides by the number of codons. | **Equivalent** to the app's sum over all 12-mer windows (every window ends up with total weight 1): the two agree to 1e-15 on random sequences. **The original reports the score per codon; the app reports the sum**, so the app's number is the original's multiplied by the number of codons. |
| **The input.** The original takes an amino-acid sequence and designs the DNA from scratch; the app recodes the DNA's own codons | Same task: the amino acids are fixed in both. |
| **The optimization.** The original (`OptimizerTools.pm`) and the app's `optimize_glo()` on the same proteins | **Different algorithms with very different results (below).** |

## The optimization

**Original.** For every group of four consecutive amino acids it scores *all* synonymous four-codon coding 12-mers (up to 6^4 = 1,296 of them). It
assembles a first sequence by placing the highest-scoring words first and filling the gaps with the best compatible words. It then improves the
sequence by choosing a random position, scoring every coding 12-mer there in the context of the neighbouring codons, taking the best one, and
keeping the change if the score of the whole sequence rises; this stops after about 4.6 times the length in successive failures, and it is repeated
10 times, keeping the best (`maxiter = -4.6052 / log(1 - 1/(L - 3))`).

**App.** It visits the codons in random order and tries *one* random synonymous alternative for a single codon, keeping it if the (up to four)
12-mer windows containing that codon score higher; up to five passes.

Because 82% of the possible 12-mers score zero, a change of a single codon rarely creates a high-scoring word: the best words need four codons to
change together. The app's search therefore stalls on plateaus the original search leaps over.

**Result** (score per codon, higher is better; the original algorithm ported to R, 10 rounds, and checked against the unchanged Perl code, see below):

| Protein | Original: first guess | Original: final | App, 5 runs (range, mean) | App's mean is lower by |
|---|---|---|---|---|
| GFP (S65C), 238 residues | 693.3 | **703.8** | 351 to 417 (391) | 44% |
| first 60 residues of GFP | 999.7 | **1000.0** | 298 to 379 (350) | 65% |
| random protein, 120 residues | 405.6 | **409.9** | 184 to 214 (193) | 53% |

The port was verified by running the **unchanged original Perl optimizer** (with small stand-ins for its database, random-number and progress modules)
on the GFP protein: its score, 703.834, equals the port's 703.8.

## What follows (before the fix)

* The app's "Germline optimized (GLO)" sequences scored about half of what the published method gives. They were improved over the input and
  kept the protein, but they were not the sequences the published algorithm would design.
* The score the app displays is the sum over the sequence; the original reports it per codon.

## After the fix (2026-10-07)

`optimize_glo()` was rewritten as a port of the original algorithm (see `R/fct_glo_optimization.R` and [METHODS.md](METHODS.md)). It keeps the
four-codon words and the first guess, and replaces the original's "stop after about 4.6 times the length in successive failures" by "revisit a position
only when a neighbouring codon has changed", with exact whole-number scores and, for the first guess, a heap. Results (score per codon, unchanged
original Perl code on the same proteins in brackets):

| Protein | App, any seed tried | Original |
|---|---|---|
| 35 residues | 509.813 | 509.813 |
| GFP (S65C), 238 residues | 703.834 | 703.834 |
| random protein, 150 residues | 668.897 | 668.897 |

The tests keep these three numbers (`tests/testthat/test-glo_optimization.R`) and fail if the refinement or the context scoring is broken (checked
by mutation). Time for the GLO step on the benchmark genes: 0.1 s (300 bp), 1.1 s (3 kb), 8.5 s (20 kb); the earlier, weaker search took 0.1 s, 0.3 s
and 0.8 s. The whole design with GLO and the piRNA shield (threshold 4) takes 1.0, 3.3 and 14 s.

Two things remain: the piRNA and restriction-site shields run *after* GLO and recode codons at random with the `Ubiq` weights, which lowers the GLO score
a little (about 539 to 528 per codon on a 20 kb gene); and the displayed score is the sum rather than the per-codon value of the original.

## The two repositories of the original author, and the number of rounds

Dan Dickinson has two public repositories: `dickinson-lab/germline` ("new", April 2018, one commit; the same files as the frozen copy
`elifesciences-publications/germline` we used: `OptimizerTools.pm`, `Seqscore.pm` and `OptimizerMain.pm` are byte-identical) and `dickinson-lab/germline_old`
(2014 to 2018, 358 commits). The optimizer in `germline_old` is the same algorithm in every version we read (from the first full version of
7 November 2014 to the last of October 2015): the same scoring, the same four-codon words, 10 rounds and the same stopping rule; the files differ from the "new"
repository only in two `use` lines and the name of the progress task. **We found no parameter for the number of iterations and no faster variant in either
repository**; a version with such a parameter (default 4, as remembered in the lab) is not in these repositories, and the original web server
(http://104.131.81.59/) no longer answers. The app's `optimize_glo()` has a `rounds` argument (10 as in the published code).

Whatever that parameter is, the number of rounds hardly matters (`dev/glo_rounds_sensitivity.R`; score per codon, mean of four seeds):

| | GFP (238 aa) | 999-aa benchmark protein |
|---|---|---|
| first guess only | 693.35 | 521.04 |
| original code, 1 round | 703.834 | 529.037 |
| original code, 4 rounds | 703.834 | 529.066 |
| original code, 10 rounds | 703.834 | 529.066 |
| original code, 4 rounds, stopping count x 0.25 | 703.608 | 529.075 |
| app, 4 rounds | 703.834 | 529.066 |
| app, 10 rounds | 703.834 | 529.094 |

The first guess reaches 98.5% of the final score; one round of refinement reaches essentially all of it. A faster variant that runs fewer rounds or stops sooner would
therefore give practically the same sequences' scores as the published one, and so does the app. If the lab's faster version differs in some other way, we would need its code to compare.

## The exact optimum and GLO-aware shields (implemented 2026-10-07)

The GLO score is a sum of terms that each depend on four consecutive codons only, so the best coding can be found exactly by dynamic programming (as in the
Viterbi algorithm): go along the sequence and, for every choice of the last three codons (at most 6^3 = 216), keep the best score reached so far, then trace back.
`optimize_glo()` now does this by default (`method = "exact"`; the port of the original search is `method = "heuristic"`). It was checked against brute force on small proteins
(the best of all codings every time).

| Sequence | Exact optimum (per codon) | Time | Port of the original (per codon) | Time | Original's shortfall |
|---|---|---|---|---|---|
| GFP, 238 aa | 717.640 | 0.3 s | 703.834 | 0.4 s | 1.9% |
| benchmark 300 bp | 717.495 | 0.0 s | 716.897 | 0.1 s | 0.08% |
| benchmark 3 kb | 530.400 | 0.1 s | 528.565 | 1.2 s | 0.35% |
| benchmark 20 kb | 541.592 | 0.5 s | 539.112 | 8.7 s | 0.46% |

Because 82% of 12-mers score zero, many codings can share the best score. A second criterion that can never override the score breaks the ties: the Codon Adaptation Index
(Ubiquitous table; the sum of the log weights), so the result is the highest-CAI coding among the best-scoring ones (checked by brute force); an optional random jitter (`diversity`, default 0) lets runs
differ. In practice the best GLO coding is almost unique: with jitter two runs differ in about 0.2% of codons. On the 3 kb benchmark gene the GLO optimum has a CAI of 0.433
(input 0.422, all-best-codon recoding 0.475): the GLO score fixes most codons, so the CAI left to choose is small.

**Combining with the shields.** (a) *Restriction sites* are avoided inside the optimization: windows of four codons that contain a site are ruled out (with a very large penalty, so a window is only
used if no alternative exists). Any site of up to 10 bases lies inside one such window, which covers 219 of the 234 enzymes and every commonly used one; the 15 longer sites (mostly with runs of N) are left to the
restriction-site step. (b) *piRNA matches* (20-base windows, up to 8 codons) cannot be handled this way, so after the optimization the shield changes one codon at a time, each time
the synonymous change that costs the least GLO score (`glo_change_gain()`), as the High-expression method does with the CAI; the same rule removes whatever restriction sites remain.
Results on the 3 kb and 20 kb benchmark genes (which contain planted piRNA targets and enzyme sites; six enzymes, piRNA threshold 4), GLO score per codon:

| | 3 kb | 20 kb |
|---|---|---|
| exact optimum, no shields | 530.40 | 541.59 |
| shields as before (Ubiquitous-weighted random changes) | 506.38 | 507.01 |
| piRNA shield GLO-aware, then enzymes repaired afterwards | 518.84 | 519.90 |
| **site avoidance inside the optimization + GLO-aware piRNA shield (the app)** | **520.94** | **524.18** |
| exact optimum with the sites avoided (before the piRNA shield) | 521.80 | 525.24 |

So the shields now cost 1.8% (3 kb) and 3.2% (20 kb) of the GLO score on these deliberately difficult genes, where before they cost about 5% and 6%; the piRNA shield alone costs about 0.2%.
The result contains none of the selected sites and no piRNA window below the threshold; the protein is unchanged (tests).
Time for the whole design with GLO, six enzymes and the piRNA shield: 3.4 s (3 kb) and 13 s (20 kb).

**Not done:** an exact treatment of the piRNA constraint (windows of 20 bases) and of the restriction sites longer than 10 bases, and the ribosome-binding-site step, which still recodes the first 13 codons with
the Ubiquitous weights before GLO designs the rest.

