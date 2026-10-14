#!/usr/bin/env bash
# Step 3 of data-raw/pirna_other_species_brief.md: trims the 3' adapter off one small-RNA library, keeps 18-28 nt reads
# and collapses them to unique sequences with read counts.
#   bash data-raw/other_species/trim_and_collapse.sh <reads.fastq.gz | run.sra> <output prefix> [max reads]
# (an .sra file is read with fastq-dump, as FASTA, since the base qualities are not used)
# The adapter was determined from the reads themselves (98% of the reads of SRR585541 contain it): the standard Illumina
# small-RNA 3' adapter TGGAATTCTCGGGTGCCAAGG. Writes <prefix>.collapsed.tsv (sequence, count), <prefix>.cutadapt.txt and
# <prefix>.stats.txt (length and 5'-nucleotide distribution).
set -eu   # no pipefail: "head" closing the pipe in the test mode is expected
fq=$1; prefix=$2; max=${3:-0}
env=$HOME/miniconda3/envs/pirna_other_species/bin
adapter=TGGAATTCTCGGGTGCCAAGG
case "$fq" in
  *.sra) reader="$env/fastq-dump --fasta 0 -Z $fq"; per=2 ;;
  *)     reader="gzip -dc $fq"; per=4 ;;
esac
if [ "$max" != 0 ]; then reader="$reader | head -n $((max * per))"; fi
eval "$reader" | "$env/cutadapt" -a "$adapter" -m 18 -M 28 --discard-untrimmed --fasta -j 4 -o - - 2> "$prefix.cutadapt.txt" \
  | awk '!/^>/' | LC_ALL=C sort -S 3G -T "$(dirname "$prefix")" | LC_ALL=C uniq -c \
  | awk '{print $2 "\t" $1}' > "$prefix.collapsed.tsv"
awk -F'\t' '{l[length($1)] += $2; if (length($1) == 21) f[substr($1, 1, 1)] += $2; n += $2}
  END {printf "reads kept (18-28 nt, adapter found): %d\nunique sequences: %d\n", n, NR; for (k = 18; k <= 28; k++) printf "length %d\t%d\t%.1f%%\n", k, l[k], 100 * l[k] / n;
       for (b in f) printf "21-nt reads, 5p %s\t%d\n", b, f[b]}' "$prefix.collapsed.tsv" | sort > "$prefix.stats.txt"
