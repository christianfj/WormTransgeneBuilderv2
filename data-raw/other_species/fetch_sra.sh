#!/usr/bin/env bash
# Downloads one run in NCBI's SRA format from the public AWS mirror (about 3x faster than the ENA FASTQ files from this
# machine) and validates it with vdb-validate. Resumes interrupted downloads; re-downloads if validation fails.
# Usage: fetch_sra.sh SRRxxxxxx
set -u
cd "$(dirname "$0")/sra"
r=$1
env=$HOME/miniconda3/envs/pirna_other_species/bin
for attempt in 1 2 3; do
  until curl -sS -f -C - --retry 5 -o "$r.sra" "https://sra-pub-run-odp.s3.amazonaws.com/sra/$r/$r"; do echo "curl retry $r"; sleep 5; done
  if "$env/vdb-validate" "$r.sra" > "$r.validate.txt" 2>&1 && ! grep -qiE "corrupt|invalid|error" "$r.validate.txt"; then
    echo "$r OK attempt $attempt $(stat -f %z $r.sra) bytes"; exit 0
  fi
  echo "$r validation failed on attempt $attempt:"; tail -2 "$r.validate.txt"; rm -f "$r.sra"
done
echo "$r FAILED"
