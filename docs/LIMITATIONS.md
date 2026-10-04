# Limitations and known issues

This is a beta release. Please check designs before ordering synthesis.

**What the tool does not do**
* It does not predict expression; the codon schemes are published heuristics, not guarantees.
* piRNA screening compares 20-base windows of the given strand with the target sequences in the database. The mismatch threshold (0 to 4) is the user's choice; the default of 3 is a design setting, not a validated cut-off for every locus. Only the given (sense) strand is compared; reverse-strand matches are not considered. The viewer marks only forward matches for piRNAs (restriction sites are marked on both strands).
* Methylation (Dam/Dcm) is ignored.
* The Twist check follows a short list of rules (homopolymers, global and local GC); the vendor's own checks decide whether a sequence can be made.
* Input is limited to 20,000 characters (20 kb).

**Cases where a request cannot be met**
* Some restriction sites cannot be removed without changing the protein (for example `ACNGT` in a Thr-Val pair). The app warns and lists them.
* Selecting very many enzymes at once is slow, leaves many unremovable sites and limits how far the sequence can be optimized. There is deliberately no "select all" button.
* The piRNA shield for the random methods makes a single pass, so a window can remain under the threshold; the High-expression shield stops with a warning when no synonymous change can help.
* For the High-expression method the ribosome-binding-site step has only two candidates to choose from (the input start and its High-expression recoding), because that method has no alternative codons.

**Known issues found while documenting (not changed in this release)**
* The tick box "Force in reading frame" for introns currently has no effect: introns are always inserted between codons.
* Random methods give a different sequence each time, apart from "No codon optimization" (fixed seed) and High expression (no random choices).

**Maintenance.** The app is expected to be maintained until about mid-2027. The tagged releases are archived
(see the README) and the images can be run unchanged afterwards.
