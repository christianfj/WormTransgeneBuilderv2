# The Fire lab introns and the 5' synthetic spliced UTR

Three elements of the app come from the Fire lab *C. elegans* vector kits (Andrew Fire and colleagues, Carnegie Institution of
Washington; the kits are now distributed through Addgene, https://www.addgene.org/kits/firelab/):

* the three artificial introns of the intron menu entry **"3 synthetic introns Fire lab kit"** (`Fire_lab_introns` in `Introns.csv`);
* the **5' UTR "Fire lab synthetic spliced"** (`5pUTRs.csv`).

The kit documents used here are *Fire Lab Vector Kit, June 1995 (updated to August 1996)*, *Fire Lab 1997 Vector Supplement* and
*Fire Lab 1999 Vector Supplement* (cited below as the 1995, 1997 and 1999 documents). The sequences were compared with **pPD95.75**
(the S65C gfp vector with the unc-54 3' end that carries this 5' UTR and these introns; Addgene plasmid 1494, 4,487 bp) as it is in the file of the whole
kit from Addgene (`data-raw/fire_lab_vectors/fire-lab-c-elegans-vector-kit-genbank-sequences.txt`, 286 vectors) by `dev/verify_fire_lab_sequences.R`;
`dev/fire_lab_vector_search.R` searches all 286 vectors. The authors also supplied their own pPD95.75 file (ApE, 14 November 2018, 4,494 bp; positions
841 to 2160 are kept in `data-raw/fire_lab/pPD95.75_excerpt.txt`). A test (`tests/testthat/test-fire_lab_sequences.R`) repeats the comparisons.

**Result in one paragraph.** The 5' UTR and all three Fire lab introns of the app are in the Addgene pPD95.75, base for base. The third intron is also the third gfp
intron of 155 other vectors of the kit (156 of the 167 vectors that have the Fire lab gfp with its introns). The authors' own pPD95.75 file differs from Addgene's
in that intron only (9 bases); it has the variant that two other vectors of the kit carry (L3920 and pJKL464). The app follows the Addgene sequence.

## Why these introns exist

The Fire lab found that introns raise the expression of transgenes in *C. elegans*: with the *unc-54* promoter, lacZ reporters with 1 to 12 introns
gave 2 to 3 orders of magnitude more expression than reporters with one intron, and gfp reporters behaved in the same way (1995 document, section II;
the same effect with a quantitative CPRG assay, "several hundred fold" between 1 and 12 introns). The effect did not depend on the
reporter, the 3' UTR or the tissue, and it did not require the intron to be in the coding region (an intron in the 5' UTR or the 3' UTR also helped).
Their explanation is open: introns may help export of the transcript (unspliced transcripts of intron-poor constructs stayed in the nucleus in
in situ experiments), RNA processing, or stability. For spacing, their limited experience favoured introns spread evenly over the transcript.
These findings are the reason the app can insert introns, and why the Fire lab gfp coding region carries three of them.

## The three gfp introns (å, ß, ∂)

In the kit documents the introns of the gfp coding region are called **å, ß and ∂** (the 1995 document, "IVS" in the vector tables), each **51 bases** (1999 document).
All of the Fire lab synthetic introns of this family have the same frame, which is the degenerate "intron transposon" oligonucleotide of the 1995 document
(section VI): the splice sites and a few fixed bases, with two stretches of 7 variable bases, so that every intron differs from the others (the aim was to
lower the chance of recombination between introns in the same construct; flanking cryptic splice sites were mutated where possible):

```
gtaagtttaaac NNNNNNN tactaactaac NNNNNNN atttaaattttcag        (12 + 7 + 11 + 7 + 14 = 51 bases; the N's were 34% GC)
5' splice site                                      3' splice site (ttttcag)
```

Their positions in pPD95.75 (Addgene numbering; the gfp coding region starts at 134 and is 870 bases), and where they sit in the gfp reading frame:

| Intron | pPD95.75 (Addgene) | Sequence (the same as in the app) | GC | Phase in gfp |
|---|---|---|---|---|
| å | 305 to 355 | `gtaagtttaaacatatatatactaactaaccctgattatttaaattttcag` | 21.6% | 0 (between codons 57 and 58) |
| ß | 510 to 560 | `gtaagtttaaacagttcggtactaactaaccatacatatttaaattttcag` | 27.5% | 1 (after the first base of a codon) |
| ∂ | 725 to 775 | `gtaagtttaaacatgattttactaactaactaatctgatttaaattttcag` | 21.6% | 0 |

**What the app has.** `Introns.csv`, entry `Fire_lab_introns`: `Sequence_1` is å, `Sequence_2` is ß and `Sequence_3` is ∂, all identical to the Addgene pPD95.75. The menu label
"avg. length 51 bp, 24% GC" is right (23.5%). The three are in 188, 188 and 186 of the 286 vectors of the kit; among the 167 vectors that have the Fire lab gfp coding region with introns, 156 have this
∂ (e.g. pPD93.97, pPD95.02, pPD95.67 to pPD95.81, pPD96.02, pPD96.12 and the vectors with 4xNLS, let-858, myo-2 and so on).

**Redemann et al. (2011).** The *C. elegans* Codon Adapter inserts the Fire lab introns into the optimized sequence. For a gfp optimized with it (NLS-gfp-NLS, CAI 0.978, three introns, Rbs-optimized; the sequence is in `data-raw/fire_lab/redemann2011_gfp_output.fa`, obtained 2026-10-11) the three inserted introns are identical, in order, to `Sequence_1`, `Sequence_2` and `Sequence_3` of the app. So the introns of the paper are the standard kit introns, and the app and the paper agree. In that output the introns follow `...TTCAAG`, `...GACTTCAAG` and `...CTCTCCAAG`, each after a complete codon (phase 0 in the coding sequence), and the gfp is the S65C form of the Fire lab (`TFCYGVQCF`).

**The ∂ variant.** The authors' own pPD95.75 file has a different ∂ intron, `gtaagtttaaacttggacttactaactaacggattatatttaaattttcag`, which differs from the app's at 9 of 51 positions
(positions 13, 16 to 18, 31, 32 and 35 to 37 of the intron). Outside this intron the two files agree over the 1,320 bases compared. In the Addgene kit this variant is in only two vectors, **L3920** and **pJKL464**
(the CeTwist promoter vector of the 1999 kit). The 1997 document says (version 1.01, March 1997) that the sequence of ivs∂ as first posted had several errors, which were corrected and
the corrected sequences posted; the Addgene sequences date from 2005 and the documents do not say which of the two variants is the corrected one or which is the plasmid as it was built.
Both fit the Fire lab frame. If the authors' file came from sequencing their own copy of pPD95.75, it would be worth knowing whether the difference is real; the app uses the sequence that
nearly all the kit's gfp vectors have. The 7-base difference in length between the two files (4,494 and 4,487 bp) is not in the stretch compared and was not investigated.

All three contain the DraI site `TTTAAA` (the frame `gtaagtttaaac`/`atttaaatttt` has it twice); none contains a site of another commonly used or Golden Gate enzyme.
Fire lab introns are inserted as they are, so a DraI site selected for removal stays in them (the app warns about that; see "Restriction sites" in `docs/METHODS.md`, section 6).

## The 5' synthetic spliced UTR

`5pUTRs.csv`, entry "Fire lab synthetic spliced": **92 bases, 46.7% GC** (the menu says 47%). It is **positions 37 to 128 of the Addgene pPD95.75** (the same stretch as positions 1002 to 1093 in the authors' file), identical base for base, and consists of

| Part | Bases | Sequence |
|---|---|---|
| 5' exon fragment | 25 | `CCGGGATTGGCCAAAGGACCCAAAG` |
| **intron A** | 42 | `gtatgtttcgaatgatactaacataacatagaacattttcag` |
| 3' exon fragment | 25 | `GAGGACCCTTGGAGGGTACCGGTAG` |

Intron A is the original synthetic intron of the first Fire lab lacZ vectors (Fire, Harrison and Dixon 1990, Gene 93:189-198, reference 3 of the 1995 document), placed upstream of the
reporter coding region; the 1999 document gives its length as **42 bases** (the other introns are 51). Its 5' splice site is `AAG|gtatgt` and its 3' splice site `ttttcag|GAGG`.
The 1999 document warns that intron A can leave only a short exon at the 5' end of the mRNA in some fusions, and omits it in a subset of its constructs. In pPD95.75 the sequence that follows the UTR is `AAAAAATGAGT...`, i.e. `aaaaaa` and then the start of gfp.

The 1995 document describes the same arrangement for the ectopic expression vectors: a first intron in the 5' UTR and a last intron in the 3' UTR already gave a substantial stimulation of expression, so adding the UTR with intron A is
a way to give an otherwise intron-less construct its first intron.

**Restriction sites in the UTR.** The UTR is added as it is (it is not recoded), so its sites stay when enzymes are selected for removal. It contains sites of **AvaII, BanI, KpnI** (`GGTACC`) and **NciI** among the commonly used enzymes
(and also AgeI, `ACCGGT`, which is next to KpnI in the same stretch). In the vectors these sites sit at the junction with the reporter (the 1995 document: KpnI and AgeI are in the same frame in all vectors, so reporter coding regions can be swapped as AgeI-ApaI segments).

## Reference data the sequences were checked against

The comparisons above use files that are **not** in the repository (the first is third-party depositor data, 3 MB). To repeat them, put the file in the place shown and run
`Rscript dev/verify_fire_lab_sequences.R`; the tests skip themselves when it is missing.

| What | Details |
|---|---|
| Kit file (local only, `data-raw/fire_lab_vectors/`, listed in `.gitignore`) | `fire-lab-c-elegans-vector-kit-genbank-sequences.txt`: the GenBank-format sequences of the whole Fire Lab *C. elegans* Vector Kit as deposited at Addgene (https://www.addgene.org/kits/firelab/), 286 records (`ADDGENE_nnnn`, deposit date 01/03/2005), 3,053,316 bytes, SHA-256 `682ebe51fe808f178e8bd76aabb74407a484a06a8a6bd82eb6906952785ae376`. Obtained by the authors on 2026-10-08. |
| The vector used for the main comparison | pPD95.75, Addgene plasmid 1494 (record `ADDGENE_1494`, 4,487 bp). Individual sequences can also be downloaded from the plasmid's page at Addgene. |
| Kit documentation | *Fire Lab Vector Kit, June 1995 (updated to August 1996)*; *Fire Lab 1997 Vector Supplement, February 1997* (history, corrections and other info: version 1.01 of March 1997, the ivs∂ correction); *Fire Lab 1999 Vector Supplement, May 1999* (version 1.1, July 2000). Supplied by the authors as PDFs (`vec95-fire-lab-vector-kit-june-1995-documentation.pdf`, `vec97-fire-lab-1997-vector-supplement-documentation.pdf`, `vec99-fire-lab-1999-vector-supplement-documentation.pdf`); they are archived by the Fire lab and linked from the Addgene kit page. They are not copied into the repository. |
| Excerpt of the authors' own pPD95.75 file (in the repository) | `data-raw/fire_lab/pPD95.75_excerpt.txt`, positions 841 to 2160 of an ApE file dated 14 November 2018. **The authors later found a problem with that file** (its third gfp intron differs from the Addgene and kit sequence); the excerpt is kept only to document the difference that was investigated, and a test pins it. It should not be used as a reference. |
| Redemann et al. 2011 tool output (in the repository) | `data-raw/fire_lab/redemann2011_gfp_output.fa`: one gfp sequence optimized with the Codon Adapter (https://worm.mpi-cbg.de/codons/cgi-bin/optimize.py), supplied by the authors on 2026-10-11, used to confirm the introns. |
| Intron design chapter | Blumenthal T and Steward K (1997), RNA processing and gene structure, chapter 6 of *C. elegans II* (NCBI Bookshelf NBK20075), section "Cis-splicing in worms"; used for the synthetic intron pool (`dev/make_synthetic_introns.R`). |
