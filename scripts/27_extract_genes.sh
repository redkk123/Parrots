#!/usr/bin/env bash
# Extract the 101 H8 panel genes (data/pilot_genes/canonical_proteins.faa, Melopsittacus RefSeq) from one assembly.
# Usage: scripts/27_extract_genes.sh <GCA_accession> "<species>" [threads]
# Download -> (chromosome-scale: 5 Mb windows overlapping by 500 kb, as in 10_audit_assembly.sh) -> miniprot
# (vertebrate splice model) -> scripts/28_parse_miniprot.py. The genome is deleted afterwards.
# Outputs: results/phase2_genes/{hits/<acc>.tsv, cds/<acc>.fna, gff/<acc>.gff.gz}
set -euo pipefail
cd "$(dirname "$0")/.."
ACC=$1; SP=$2; T=${3:-4}
BIN=$PWD/.bioenv/bin; export PATH=$BIN:$PATH
OUT=results/phase2_genes; TMP=data/external/genes_tmp/$ACC
if [ -s "$OUT/hits/$ACC.tsv" ]; then echo "$ACC done"; exit 0; fi
mkdir -p "$OUT/gff" "$TMP"
if [ ! -s "$TMP/genome.fna" ]; then
  datasets download genome accession "$ACC" --include genome --filename "$TMP/pkg.zip" --no-progressbar
  unzip -p "$TMP/pkg.zip" 'ncbi_dataset/data/*/*.fna' > "$TMP/genome.fna.part"
  mv "$TMP/genome.fna.part" "$TMP/genome.fna"; rm "$TMP/pkg.zip"
fi
GEN=$TMP/genome.fna
if [ "$(seqkit stats -T "$GEN" | awk -F'\t' 'NR==2{print $8}')" -gt 20000000 ]; then
  seqkit sliding -g -W 5000000 -s 4500000 -j "$T" "$GEN" > "$TMP/windows.fna"; GEN=$TMP/windows.fna
fi
/usr/bin/time -v miniprot -t "$T" -j 2 --gff --outn=5 "$GEN" data/pilot_genes/canonical_proteins.faa \
  > "$TMP/miniprot.gff" 2> "$TMP/time.txt"
.venv/bin/python scripts/28_parse_miniprot.py "$ACC" "$SP" "$GEN" "$TMP/miniprot.gff"
gzip -c "$TMP/miniprot.gff" > "$OUT/gff/$ACC.gff.gz"
grep -E "Maximum resident|Elapsed" "$TMP/time.txt" | sed "s/^/$ACC /" >> logs/27_extract_genes.time.log
rm -rf "$TMP"
