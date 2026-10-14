# piRNA databases for other *Caenorhabditis* species

*Status: these databases are used by the app (`inst/extdata/piRNA_Database_<species>.csv.gz`) when the user selects the species under "Minimize piRNA homology". They were promoted from `data-raw/staging/other_species/` by `dev/make_pirna_data.R`. The species databases are independent of each other and of the *C. elegans* database.*

**In plain words.** The app can avoid sequences that look like the piRNAs of *C. briggsae*, *C. remanei* and *C. brenneri* as well as those of *C. elegans*, if the transgene will be used in those species. The *C. briggsae* list is high quality: it comes from a published experiment that pulled down the piRNA-binding protein PRG-1 and sequenced what was bound. For *C. remanei* and *C. brenneri* no such experiment exists, so we identified piRNAs ourselves in a published small-RNA sequencing data set, with strict rules (the piRNA must be 21 bases starting with T, carry the typical upstream motif and be seen at a reasonable level). Tested on *C. elegans*, where the answer is known, this approach finds roughly one in five to one in seven of the real piRNAs, and about nine in ten of what it finds are real. The two lists are therefore **short and cautious**: using them removes more homology to known piRNAs, but a transgene is not guaranteed to be free of homology to piRNAs these lists do not contain. The technical details follow, for readers who want to check them.

These databases follow the format of the *C. elegans* database ([PIRNA_DATABASE.md](PIRNA_DATABASE.md)): `Target_20mer` is the
reverse complement of `Full_21mer` without its first base, so that sense-strand windows of a transgene can be compared with it
directly. Each entry also records species, strain, assembly, confidence, evidence and (where available) tier and expression.

## Contents

| Species (strain) | Entries | Confidence | Source |
|---|---|---|---|
| *C. briggsae* (AF16) | 25,883 | **high** | PRG-1 immunoprecipitation, Pastore, Hertz & Tang 2022 (GEO GSE203297), Table S1 |
| *C. remanei* (PB4641) | 11,801 | **low** | called from Shi et al. 2013 reads (GEO GSE41461) |
| *C. brenneri* (PB2801) | 17,991 | **low** | called from Shi et al. 2013 reads (GEO GSE41461) |

All names are internal labels (`Cbr-21UR-N` as in the briggsae paper, `Cre_21U_<scaffold>_<start>_<plus|minus>`,
`Cbn_21U_...`); none is a formal gene name.

## *C. briggsae*

The Pastore et al. table lists 25,883 piRNA loci defined by their enrichment in a *C. briggsae* PRG-1 immunoprecipitation (chromosome,
start, end, name, 21-base sequence, strand, mean IP RPM). We took it as it is, after checking every sequence against the CB4 genome
(WormBase ParaSite WBPS19, BioProject PRJNA10731): all 25,883 are 21 bases, start with T, are unique, and equal the genome at their
coordinates (the table's `Start` is 0-based). 17,030 are on chromosome IV and 2,380 on unplaced contigs. The Ruby-motif flag
(below) is set for 49% of them; for briggsae it is information only and was not used to select entries.

## *C. remanei* and *C. brenneri*: how the low-confidence sets were made

No PRG-1 IP data exist for these two species, so the piRNAs were called from the small-RNA sequencing of Shi, Montgomery, Qi & Ruvkun
(2013, Genome Res 23:497-508). We used the adult (males and females pooled), embryo and adult-male libraries of each species (one library
per condition, no replicates; runs SRR585541-SRR585546), and the *C. elegans* adult library (SRR585535) as a control.

1. **Download.** Reads came as NCBI SRA files (public AWS mirror), checked with `vdb-validate`.
2. **Trimming.** The 3' adapter was determined from the reads themselves: 98% contain the Illumina small-RNA adapter
   `TGGAATTCTCGGGTGCCAAGG`. It was removed with cutadapt, keeping reads of 18 to 28 nt (97% of reads passed in the first library tested).
   Reads are collapsed to unique sequences with counts.
3. **Mapping.** bowtie 1.3.1, **0 mismatches**, all alignments, to the genome of the matching strain: PB4641 (BioProject PRJNA53967, 3,670
   scaffolds, 145 Mb) and PB2801 (PRJNA20035, 3,305 scaffolds, 190 Mb), the assemblies Shi et al. used (served by WormBase ParaSite WBPS19).
   Between 51% and 79% of the trimmed reads map.
4. **Candidates.** A candidate is a 21-nt read that starts with T, maps exactly, and is not a miRNA. A read is treated as a miRNA if it shares
   a stretch of 19 bases with a mature miRNA of the species (miRBase 22 entries for *C. remanei* and *C. brenneri*, plus the new miRNAs of Shi
   et al. 2013, Supplementary Table S3). RPM is the read count divided by the number of loci the sequence maps to, per million mapped reads.
5. **Ruby motif.** An upstream core `GTTTC` that lies completely within positions -50 to -35 relative to the first base of the piRNA
   (the core motif is about 35 to 42 bases upstream; Billi et al. 2013).
6. **Selection.** An entry is kept if the motif is present and the sequence has at least 1 RPM in some library ("tier A"). Candidates without the
   motif ("tier B": 4,633 in remanei and 8,151 in brenneri at 1 RPM) and candidates below 1 RPM are not included: in the control, motif-less
   calls had only 77% precision, against 94% with the motif.

### The control, and why 1 RPM

The same pipeline was run on the *C. elegans* adult library and compared with the *C. elegans* database ([PIRNA_DATABASE.md](PIRNA_DATABASE.md)):

| Rule | Calls | Precision (known type 1) | Recall of known type 1 |
|---|---|---|---|
| any read | 65,921 | 14% | 63% |
| at least 2 reads | 23,754 | 34% | 54% |
| at least 2 reads and motif | 7,354 | 79% | 39% |
| **at least 1 RPM** | 3,550 | **88%** | 21% |
| at least 1 RPM and motif (tier A) | 2,307 | 94% | 15% |
| at least 1 RPM, no motif (tier B) | 1,243 | 77% | 6% |
| at least 5 RPM | 1,320 | 93% | 6% |

Only 63% of the known type 1 piRNAs have any read in a single adult library, which caps the recall of any rule. The 1 RPM rule keeps the
precision high at the price of recall, so these databases are *conservative lists of piRNAs that are well expressed in these libraries*.
That is also the set most likely to matter for silencing a transgene. Precision was measured only for *C. elegans*; for the other two species it
is unknown. The Ruby motif rule fires at 2.3% (remanei) and 2.1% (brenneri) of random genome positions, against 2.7% in *C. elegans*.

The totals (11,801 and 17,991) are lower than the numbers Shi et al. report as observed (23.7k and 29.8k) and than their predictions
(33.8k and 54.8k); the differences are the 1 RPM cut and the requirement for the motif. They are conservative lists.

## Limitations

* One library per condition and no replicates; no PRG-1 IP for remanei and brenneri.
* Both assemblies are scaffold-level (contig N50 about 30 kb): a piRNA at a scaffold end can be missed, and the same piRNA can be split.
* Reads were mapped to the reference strains with no mismatches, which is correct only if the sequenced animals match the assembly. Barrière et al. (2009) show that the *C. remanei* and *C. brenneri* assemblies contain both alleles of many loci (about 10% and 30% of the genome), and that the strains have changed since they were sequenced (allele frequencies shifted and some alleles were lost). A piRNA from a locus that differs between the assembly and the strain now used may therefore be missed or listed twice (a second allele is a second, nearly identical, locus; reads that map to both are counted once per locus when RPM is computed).
* Requiring the Ruby motif leaves out piRNAs that lack it. In *C. elegans* such piRNAs are a real class (type 2, Gu et al. 2012), and in *C. briggsae*
  Pastore et al. found that the motif is absent from the piRNA loci outside the large clusters, so type 2-like piRNAs of remanei and brenneri
  are probably missing from these lists.
* The embryo libraries carry less piRNA signal (about 40% of their 21-nt reads start with T) but were kept because they can contribute
  entries not seen in adults.
* A sex-enrichment comparison (as in Shi et al. Fig. 4) was not done: with one library per condition it would be only exploratory.
* The sequence-based prediction of piRNA expression that is being developed for *C. elegans* could later replace the read-based rule.

## Reproducing

The scripts are in `data-raw/other_species/`: `fetch_sra.sh` (downloads), `trim_and_collapse.sh`, `map_library.sh`, `call_piRNAs.R`,
`evaluate_control.R`, `build_species.R` (remanei, brenneri), `build_briggsae.R`, `species_sensitivity.R` and `calibrate_motif.R`. Raw reads, genomes and
bowtie indices (about 8 GB) are deliberately not part of the repository; the scripts regenerate them from public sources (SRA run accessions, WormBase ParaSite
files listed above, miRBase, Shi et al. Table S3). The working folder is set at the top of `call_piRNAs.R`, `build_species.R` and
`run_remanei_brenneri.sh`; in short, run `fetch_sra.sh` for each run, `trim_and_collapse.sh`, `map_library.sh`, `call_piRNAs.R`, then `build_species.R`. Tools: cutadapt 5.2, bowtie 1.3.1, samtools 1.24, sra-tools 3.4.1.
Tests: `tests/testthat/test-pirna_other_species.R`.

## References

* Pastore B, Hertz HL, Tang W (2022). Comparative analysis of piRNA sequences, targets and functions in nematodes. RNA Biol 19(1):1276-1292.
  doi:10.1080/15476286.2022.2149170 (open access, CC BY 4.0; data: GEO GSE203297; supplementary Table S1 lists the 25,883 piRNA loci).
  In that paper, piRNA loci were defined from two PRG-1 IP replicates aligned to WormBase WS279 (at least 4-fold enriched, p < 0.05, 21 nt, 5' U),
  and 302 of the loci lie on contigs without chromosome coordinates.
* Shi Z, Montgomery TA, Qi Y, Ruvkun G (2013). High-throughput sequencing reveals extraordinary fluidity of miRNA, piRNA, and siRNA
  pathways in nematodes. Genome Res 23:497-508. (data: GEO GSE41461)
* Billi AC, et al. (2013). PLoS Genet 9:e1003392 (Ruby upstream motif); Gu W, et al. (2012). Cell 151:1488-1500 (type 1 and type 2).
* Genome sequences (from WormBase ParaSite; see also [DATA.md](DATA.md)): *C. briggsae*: Stein LD, Bao Z, Blasiar D, et al. (2003), The genome sequence of *Caenorhabditis briggsae*: a platform for comparative genomics, PLoS Biol 1:166-192, doi:10.1371/journal.pbio.0000045 (the CB4 assembly used here is the later, chromosome-level version of that genome). *C. remanei* PB4641 and *C. brenneri* PB2801: the assemblies were made at the Genome Sequencing Center of Washington University and **have no dedicated genome paper**; Barrière et al. (2009) say the genomes "will be described in detail elsewhere". The paper to cite for their properties is Barrière A, Yang SP, Pekarek E, Thomas CG, Haag ES, Ruvinsky I (2009), Detecting heterozygosity in shotgun genome assemblies: lessons from obligately outcrossing nematodes, Genome Res 19:470-480, doi:10.1101/gr.081851.108 (read 2026-10-14). It shows that, in these assemblies, about 10% (*C. remanei*) and 30% (*C. brenneri*) of the genome is present as two alleles, almost all on autosomes (on chromosome IV in *C. remanei*), because the sequenced strains were not fully homozygous. See the limitation below.
