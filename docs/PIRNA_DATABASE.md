# The piRNA database

*Status: this is the database the app uses (`inst/extdata/piRNA_Database.csv`, built from the staging file
`data-raw/staging/piRNA_Database_staging.csv` by `dev/make_pirna_data.R`). It replaces the 17,849-entry database used before.*

The app screens a transgene for 20-base stretches that are close to a piRNA target (see [METHODS.md](METHODS.md)). This page
says what is in the database, where each part comes from, and why.

## Composition (18,207 entries)

| Type | Entries | What it is | Source |
|---|---|---|---|
| **type 1** | 15,363 | All piRNAs (21U-RNAs) annotated in WormBase release **WS298** | WormBase annotation (`piRNA` features); the 21-mer is cut from the genome sequence of the same release |
| **type 2** | 2,485 | 21U-RNAs that are processed from capped small RNAs at Pol II promoters and need not have the type 1 upstream motif; PRG-1-bound | Mello lab (Gu et al. 2012); the subset used by pirScan (Wu et al. 2018), see below |
| **PRG-1 IP new** | 359 | 359 distinct sequences (at 436 genomic loci) that are enriched in PRG-1 immunoprecipitation but not annotated as piRNAs | Claycomb lab, Seroussi et al. 2023 (supplementary table 6) |

Type 1 entries keep their WormBase 21ur name. Type 2 entries are named "21ur-N (type 2)" as before. The new entries have
**internal labels only** (`new21U_<chromosome>_<leftmost base>_<plus|minus>`); they are not approved gene names and no
21ur- name has been assigned. A sequence found at several loci is one entry that lists all of them.

## What type 2 piRNAs are, and why the set has 2,485 entries

Gu et al. (2012) divided *C. elegans* 21U-RNAs into two classes. **Type 1** 21U-RNAs are the classical ones: they share an
8-nucleotide upstream motif (CTGTTTCA) and are clustered on chromosome IV. **Type 2** 21U-RNAs "need not have an 8 nt motif
and are processed from csRNAs [capped small RNAs] derived from the promoters of Pol II genes throughout the genome". By
deep-sequencing small RNAs immunoprecipitated with PRG-1, Gu et al. identified 12,183 new 21U-RNA species, "of which
~10,000 exhibit a poor motif score and are associated with TS [transcription start] sites throughout the genome".
Type 2 21U-RNAs are germline-specific and 3'-end modified like type 1, and most are expressed at low levels (about 5% of
total 21U-RNA), although the single most abundant 21U-RNA is a type 2 RNA expressed from the X chromosome.

The 2,485 entries used here are a much smaller subset of that candidate list. They are the set of type II piRNAs used by
**pirScan**, the web server for predicting piRNA targeting sites and avoiding transgene silencing (Wu et al. 2018, Lee
laboratory). The paper states: "The *C. elegans* piRNA sequence data, which includes 15364 type I piRNAs that can be
downloaded from WormMine and 2485 type II piRNAs, which was collected from PRG-1 IP data (26) with at least 1 read per
million reads in (27)", where (26) is Gu et al. (2012) and (27) is Tang et al. (2016). The criterion is therefore a
minimum PRG-1 immunoprecipitation signal of 1 read per million, and this is the reason the set is smaller than the full
candidate list: it keeps the type 2 piRNAs for which there is direct evidence that PRG-1 binds them at a detectable level,
and leaves out candidates that were only inferred from their start-site context. We have not re-derived the 2,485 from the
source data (their per-entry read counts are in supplementary tables we did not use); the wording of the pirScan paper is
ambiguous about which of the two studies the 1 RPM threshold was applied to.

Type 1 (15,364 in pirScan; 15,363 in WormBase WS298) and type 2 together make 17,849 entries, which is the size and
composition (15,364 + 2,485) of the database this app used before the rebuild, so that database was evidently the pirScan
set. The app offers type 2 piRNAs as an option ("Include C. elegans type 2 piRNAs", on by default, also when other species are selected), so that users can choose whether to avoid them.

Per-entry read counts for type 2 are not available to us, so type 2 entries carry no abundance of their own. Thirty of them
are also among the Claycomb PRG-1-enriched loci (next section) and have an abundance from that table.

## The 466 new PRG-1-IP loci (Claycomb lab)

Seroussi et al. (2023) detected reads for 14,568 of the 15,363 annotated piRNAs (WS276) and found 5,943 enriched at
least 2-fold, with at least 5 RPM, in PRG-1 immunoprecipitation (5,932 of them with PRG-1). They also report 466 new 21U small
RNAs that are enriched at least 2-fold in both PRG-1 immunoprecipitation replicates, are not perfect matches to annotated
piRNAs (some are shifted by a few bases, for example 3 bases from 21ur-2789), are depleted in *prg-1* mutants (222
significantly) and have Ruby-like upstream motifs; 375 of the 466 lie in the chromosome IV cluster at 4.5-7.0 Mb. The
paper concludes that they are previously uncharacterized piRNAs. The table gives positions only. We cut the sequences from the genome after
checking how the position is defined: it is the **leftmost base** of the 21-base stretch (for a minus-strand RNA the 5'
end is the right end). With this reading all 466 sequences begin with T; every other reading gave 19% to 78%. Of the 466,
30 are identical to old type 2 entries, 186 overlap (are shifted from) an annotated piRNA, and 250 are independent.

## Abundance and enrichment columns

* **Abundance**: the PRG-1 immunoprecipitation signal from Claycomb (RPM, average of the two replicates), written as
  "5.96 RPM (PRG-1 IP)". "not detected in PRG-1 IP" for type 1 piRNAs that are not in their table or have no reads;
  "not in Claycomb PRG-1 IP data" for type 2 entries that are not in it.
* **Enrichment**: male/female enrichment from Billi et al. (2013), kept as before. Billi et al. scored each 21U RNA over 23
  male-versus-female small-RNA library comparisons (+1 if more than 5-fold higher in the female library, -1 if more than
  5-fold higher in the male library); a score of -3 or below is "male", +3 or above "female", otherwise "not enriched"
  (their analysis covered 13,711 21U RNAs aligned to WS220). The table shows it as "male (score: -4)"; "No Data" is used for piRNAs that Billi et al. do not
  list or list as not determined, and for type 2 and new entries.
* Further columns record the source, the genomic position (WS298), the Claycomb IP values and the depletion in *prg-1*
  mutants (Claycomb supplementary table 4; where a gene has two rows the most significant one is kept).

## Known features

* 579 of the 15,363 type 1 piRNAs (3.8%) do not begin with T although the first base of a 21U-RNA is usually U. These
  are the WormBase annotation as it is, and a known class: Gu et al. (2012) report that about 400 previously annotated
  21U-RNAs (3% of the annotated ones) start with a nucleotide other than U ("21A/G/C" piRNAs), that these loci make
  normal levels of capped small RNA precursors, and that the mature 21-nt RNAs are on average 10-fold less abundant relative
  to their precursors, suggesting that they are processed inefficiently or are unstable. This agrees with what we see in
  the Claycomb data: they are detected in PRG-1 immunoprecipitation less often than the others (84% against 98%) and
  enriched less often (15% against 44%). They are kept, because they are annotated piRNAs.
* The first base does not take part in the matching: the app compares the 20 bases that follow it. This follows the
  statement in pirScan that "the first nucleotide does not contribute to piRNA targeting" (Wu et al. 2018), so
  these entries are screened like all others.
* A few type 2 sequences map to several places in the genome (99 of the 2,485). Gu et al. (2012) describe a family of six
  homologous loci on the X chromosome that includes the most abundant 21U-RNA.

Databases for *C. briggsae*, *C. remanei* and *C. brenneri* are described in [PIRNA_OTHER_SPECIES.md](PIRNA_OTHER_SPECIES.md).

## How it was made

`data-raw/` holds the scripts: `extract_pirna_from_gff3.R` (cuts the 21-mers from the genome), `validate_new21U_and_type2.R`
(locates type 2 sequences, tests the position convention), `build_pirna_staging.R` (joins everything), with the checks in
`tests/testthat/test-pirna_staging.R`. The 7 GB annotation file and the supplementary tables are not part of the repository.

## References

* Wu WS, Huang WC, Brown JS, Zhang D, Song X, Chen H, Tu S, Weng Z, Lee HC (2018). pirScan: a webserver to predict piRNA
  targeting sites and to avoid transgene silencing in *C. elegans*. Nucleic Acids Res 46:W43-W48. doi:10.1093/nar/gky277
* Gu W, Lee HC, Chaves D, Youngman EM, Pazour GJ, Conte D Jr, Mello CC (2012). CapSeq and CIP-TAP identify Pol II start sites
  and reveal capped small RNAs as *C. elegans* piRNA precursors. Cell 151:1488-1500. doi:10.1016/j.cell.2012.11.023
  (data: GEO GSE40053)
* Tang W, Tu S, Lee HC, Weng Z, Mello CC (2016). The RNase PARN-1 trims piRNA 3' ends to promote transcriptome
  surveillance in *C. elegans*. Cell 164:974-984. doi:10.1016/j.cell.2016.02.008 (data: SRA SRS1021265)
* Zhang D, Tu S, Stubna M, Wu WS, Huang WC, Weng Z, Lee HC (2018). The piRNA targeting rules and the resistance to piRNA
  silencing in endogenous genes. Science 359:587-592.
* Priyadarshini M, Ni JZ, Vargas-Velazquez AM, Gu SG, Frøkjær-Jensen C (2022). Reprogramming the piRNA pathway for
  multiplexed and transgenerational gene silencing in *C. elegans*. Nature Methods 19:187-194. doi:10.1038/s41592-021-01369-z
  (mismatch tolerance of synthetic guide piRNAs; data: GEO GSE165210)
* Seroussi U, Lugowski A, Wadi L, Lao RX, Willis AR, Zhao W, Sundby AE, Charlesworth AG, Reinke AW, Claycomb JM (2023). A
  comprehensive survey of *C. elegans* argonaute proteins reveals organism-wide gene regulatory networks and functions.
  eLife 12:e83853. doi:10.7554/eLife.83853 (data: GEO GSE208702; reads mapped to WormBase WS276; open access, CC BY)
* Billi AC, Freeberg MA, Day AM, Chun SY, Khivansara V, Kim JK (2013). A conserved upstream motif orchestrates autonomous,
  germline-enriched expression of *Caenorhabditis elegans* piRNAs. PLoS Genet 9(3):e1003392.
  doi:10.1371/journal.pgen.1003392 (male/female enrichment; open access, CC BY)
* Bagijn MP, et al. (2012). Function, targets, and evolution of *Caenorhabditis elegans* piRNAs. Science 337:574-578.
