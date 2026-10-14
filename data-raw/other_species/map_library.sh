#!/usr/bin/env bash
# Step 4/5 of data-raw/pirna_other_species_brief.md: maps the collapsed reads of one library to a genome with bowtie,
# with 0 mismatches (as Shi et al. 2013 did), and writes
#   <prefix>.mapped_any.tsv   sequence, count     every unique 18-28 nt sequence that maps at least once (for the
#                                                   total number of mapped reads, the denominator of RPM)
#   <prefix>.t21.aln          bowtie alignments (all of them) of the 21-nt sequences that start with T: the piRNA candidates
# Usage: bash data-raw/other_species/map_library.sh <collapsed.tsv> <bowtie index prefix> <output prefix>
set -eu
col=$1; index=$2; prefix=$3
env=$HOME/miniconda3/envs/pirna_other_species/bin
tmp=$(mktemp -d "$(dirname "$prefix")/map.XXXXXX"); trap 'rm -rf "$tmp"' EXIT
# all unique sequences, as FASTA with the sequence itself as the name
awk -F'\t' '{print ">" $1 "\n" $1}' "$col" > "$tmp/all.fa"
"$env/bowtie" -f -v 0 -k 1 -p 4 --suppress 2,3,4,5,6,7,8 "$index" "$tmp/all.fa" "$tmp/any.txt" 2> "$prefix.bowtie_all.log"
sort -u "$tmp/any.txt" > "$tmp/any.sorted"
LC_ALL=C sort -t "$(printf '\t')" -k1,1 "$col" > "$tmp/col.sorted"
LC_ALL=C join -t "$(printf '\t')" "$tmp/any.sorted" "$tmp/col.sorted" > "$prefix.mapped_any.tsv"
# candidates: 21 nt, first base T
awk -F'\t' 'length($1) == 21 && substr($1, 1, 1) == "T" {print ">" $1 "\n" $1}' "$col" > "$tmp/t21.fa"
"$env/bowtie" -f -v 0 -a -p 4 --suppress 5,6,7,8 "$index" "$tmp/t21.fa" "$prefix.t21.aln" 2> "$prefix.bowtie_t21.log"
