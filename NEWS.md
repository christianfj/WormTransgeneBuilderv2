# WormTransgeneBuilder (development version)
* The piRNA database no longer contains 301 entries that had no target sequence (17,849 piRNAs remain). They were not present in another list of piRNAs or on WormBase, so they are not validated, and with no sequence they could never be flagged or removed by the piRNA shield.
* "High expression" recoding is now deterministic, and its restriction-site and piRNA shields choose the change with the smallest loss of Codon Adaptation Index. Non-optimal codons are underlined in the sequence viewer.
* Restriction-site removal understands IUPAC ambiguity codes (for example BanI), the enzyme list is read from `Enzymes.csv`, and a site that cannot be removed now gives a warning. The random methods can fall back to a neighbouring codon.
* Synthetic introns are picked at random again for every method.
* The GGSGG linker tag no longer ends in a stop codon.
* piRNA screening is several times faster (about 6 times for a 3 kb gene and 11 times for a 20 kb gene), with identical results.
* Germline optimization (GLO) is faster and uses far less memory, with identical results: its 12-mer scores are found by arithmetic instead of searching a dictionary of names. The first GLO request after the app starts no longer takes about 9 seconds, GLO scoring is 18 to 34 times faster, the GLO optimization is 4 to 7 times faster, and the score data takes 128 MB of memory instead of 1.3 GB.
* Restriction-site removal is much faster (about 30 times for the "High expression" method with many enzymes, 15 times for the other methods), with identical results. Together with the piRNA and GLO speed-ups, a whole design of a 3 kb gene now takes about 1 second instead of 13 to 21 seconds.
* The restriction-enzyme list is now New England Biolabs' full catalogue of 234 enzymes (it had 21), including enzymes whose sites contain ambiguity letters or runs of any base (for example AhdI, SfiI, XcmI). DNA methylation (Dam/Dcm) is no longer recorded or used: a synthesized transgene is not methylated, so a site is treated as a site whatever the cloning host later does.
* The restriction-enzyme menu is grouped: "Commonly used" (22 enzymes), "Golden Gate" (8, including BsaI and BsmBI, which are listed only there) and "All other enzymes". There is deliberately no "select all", because removing many sites lowers how far a sequence can be optimized. No enzyme is pre-selected any more (BsaI used to be).
* The promoter list is expanded to 38 validated promoters (300 bp) in two groups, "Validated single neuron promoters" and "Validated tissue-specific promoters". The menu shows the name and the gene (for example "pADF (WBGene00005359, srh-142)"); the sequence viewer and the GenBank file use the name only ("pADF"). Group headings cannot be chosen. The earlier five promoters were replaced by this list (the old "DD" promoter, Y48G10A.6, is not in the new list).
* Reproducible build: R packages are pinned with `renv.lock`, and the new Docker image is built from R 4.5.3 on Ubuntu 24.04 (pinned by digest) with ViennaRNA 2.7.2 compiled from checksum-verified source. The About tab and the start-up message show the app, R and ViennaRNA versions. A fixed run gives identical output on the Mac and in the arm64 and amd64 images.
* "Force in reading frame" for introns now works. Ticked (the default, and what the app always did before), an intron goes only between codons; unticked, it may go after any base, using the same AAG/CAG preference.
* The app no longer depends on Bioconductor (Biostrings) or stringr; it needs only CRAN packages and the external program ViennaRNA. One small difference: the optional piRNA highlighting in the sequence viewer no longer marks matches that would hang off the end of a sequence.
* Added an automated test suite and benchmark scripts.

# WormTransgeneBuilder 2.0.0
* Major architectural overhaul using the Golem framework.
* Migrated to Docker for universal deployment.
* Integrated modern `bslib` UI components.
* Optimized piRNA shielding and GLO algorithms.