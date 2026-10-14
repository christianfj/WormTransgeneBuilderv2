# Limitations and known issues

This is a beta release. Please check designs before ordering synthesis.

**What the tool does not do**
* It does not predict expression; the codon schemes are published heuristics, not guarantees.
* piRNA screening compares 20-base windows of the given strand with the target sequences in the database. The mismatch threshold (0 to 4) is the user's choice; the default is 4, which follows our experiments with synthetic guide piRNAs (four mismatches abolished silencing of *him-5*, up to three did not; Priyadarshini et al. 2022) and is not a validated cut-off for every locus or every piRNA. Only the sense strand of the transcript is compared with the stored targets (the reverse complements of the piRNAs without their first base), because piRNAs pair with the sense transcript; matches on the opposite strand are irrelevant and are not considered. The viewer marks only these sense-strand matches for piRNAs (restriction sites are marked on both strands). Mismatches are counted equally at every position, with no special weight for the seed (nucleotides 2 to 7) or for G:U pairs. This is a deliberate simplification based on our own experiments with synthetic guide piRNAs (Priyadarshini et al. 2022: up to three mismatches tolerated, four or five abolish silencing, wherever they are), and it differs from pirScan (Wu et al. 2018; Zhang et al. 2018), which scores position and wobble pairing. The two tools can therefore flag different sites. Users who want the details for individual piRNAs and mismatch positions should use pirScan; see [METHODS.md](METHODS.md), section 3. The evidence behind the simplification comes from one target gene in extrachromosomal arrays.
* Methylation (Dam/Dcm) is ignored.
* The Twist check follows a short list of rules read from Twist's pages on 2026-10-14 (homopolymers, global and local GC, direct repeats, length, ccdB; not tandem repeats or hairpins, and not Twist's complexity score); the vendor's own checks decide whether a sequence can be made.
* Input is limited to 20,000 characters (20 kb).

**Cases where a request cannot be met**
* Some restriction sites cannot be removed without changing the protein (for example `ACNGT` in a Thr-Val pair). The app warns and lists them.
* Selecting very many enzymes at once is slow, leaves many unremovable sites and limits how far the sequence can be optimized. There is deliberately no "select all" button.
* The piRNA shield for the random methods makes a single pass, so a window can remain under the threshold; the High-expression shield stops with a warning when no synonymous change can help.
* With High expression, the ribosome-binding-site step recodes the first 13 codons using graded codon weights (not only the single best codon), so those codons are not all best codons.

**Known issues (not changed in this release)**
* With the GLO method the piRNA shield and the removal of restriction sites longer than 10 bases (and of any site that cannot be avoided) change one codon at a time, choosing the change that costs the least GLO score; this is not the exact optimum under those constraints (the loss is about 0.2% for the piRNA shield). The ribosome-binding-site step recodes the first 13 codons with the Ubiquitous weights, independently of GLO.
* Random methods give a different sequence each time, apart from "No codon optimization" (fixed seed) and High expression (no random choices).

**Maintenance.** The app is expected to be maintained until about mid-2027. The tagged releases are archived
(see the README) and the images can be run unchanged afterwards.
