#!/usr/bin/env bash
# Step 5 driver: indexes the PB4641 (remanei) and PB2801 (brenneri) genomes, maps the six libraries with 0 mismatches and calls
# candidates. Paths are the working folder on the external drive.
set -eu
S=/Volumes/CFJ_backup1S1TB/pirna_other_species
env=$HOME/miniconda3/envs/pirna_other_species/bin
for sp in remanei brenneri; do
  case $sp in remanei) gz=$S/genomes/caenorhabditis_remanei.PRJNA53967.WBPS19.genomic.fa.gz; libs="SRR585541 SRR585542 SRR585543"; mir=crm; s3=Cre;;
              brenneri) gz=$S/genomes/caenorhabditis_brenneri.PRJNA20035.WBPS19.genomic.fa.gz; libs="SRR585544 SRR585545 SRR585546"; mir=cbn; s3=Cbn;; esac
  fa=$S/genomes/$sp.fa
  [ -s "$fa" ] || gunzip -c "$gz" > "$fa"
  [ -s "$S/index/$sp.1.ebwt" ] || "$env/bowtie-build" --threads 4 "$fa" "$S/index/$sp" > "$S/index/${sp}_build.log" 2>&1
  for r in $libs; do
    bash data-raw/other_species/map_library.sh "$S/collapsed/$r.collapsed.tsv" "$S/index/$sp" "$S/collapsed/$r"
    Rscript data-raw/other_species/call_piRNAs.R "$S/collapsed/$r" "$fa" "$mir" "$s3" > "$S/collapsed/$r.call.log" 2>&1
    echo "$r done"
  done
done
echo ALLDONE
