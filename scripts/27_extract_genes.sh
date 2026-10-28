#!/usr/bin/env bash
# Extract the Melopsittacus canonical proteome (data/genome_genes, 16,045 genes incl. the 101 H8 panel genes) from one
# assembly in a single miniprot run (scripts/33_prepare_proteome.py builds the query set).
# Usage: scripts/27_extract_genes.sh <GCA_accession> "<species>" [threads]
# Download -> (chromosome-scale: 5 Mb windows overlapping by 500 kb, as in 10_audit_assembly.sh) -> miniprot
# (vertebrate splice model) -> scripts/28_parse_miniprot.py. The genome is deleted afterwards.
# Outputs: panel results/phase2_genes/{hits/<acc>.tsv, cds/<acc>.fna}; whole proteome
# results/phase3_genome/calls/<acc>.tsv.gz and data/genome_cds/<acc>.fna.gz.
set -euo pipefail
cd "$(dirname "$0")/.."
ACC=$1; SP=$2; T=${3:-8}
Q=data/genome_genes/canonical_proteins.faa
[ -s "$Q" ] || .venv/bin/python scripts/33_prepare_proteome.py
BIN=$PWD/.bioenv/bin; export PATH=$BIN:$PATH
OUT=results/phase2_genes; TMP=data/external/genes_tmp/$ACC
if [ -s "$OUT/hits/$ACC.tsv" ] && [ -s "results/phase3_genome/calls/$ACC.tsv.gz" ]; then echo "$ACC done"; exit 0; fi
mkdir -p "$TMP"
if [ ! -s "$TMP/genome.fna" ]; then
  datasets download genome accession "$ACC" --include genome --filename "$TMP/pkg.zip" --no-progressbar
  unzip -p "$TMP/pkg.zip" 'ncbi_dataset/data/*/*.fna' > "$TMP/genome.fna.part"
  mv "$TMP/genome.fna.part" "$TMP/genome.fna"; rm "$TMP/pkg.zip"
fi
GEN=$TMP/genome.fna
if [ "$(seqkit stats -T "$GEN" | awk -F'\t' 'NR==2{print $8}')" -gt 20000000 ]; then
  seqkit sliding -g -W 5000000 -s 4500000 -j "$T" "$GEN" > "$TMP/windows.fna"; GEN=$TMP/windows.fna
fi
/usr/bin/time -v miniprot -t "$T" -K 500k -j 2 --gff --outn=3 "$GEN" "$Q" \
  > "$TMP/miniprot.gff" 2> "$TMP/time.txt"
.venv/bin/python scripts/28_parse_miniprot.py "$ACC" "$SP" "$GEN" "$TMP/miniprot.gff"
grep -E "Maximum resident|Elapsed" "$TMP/time.txt" | sed "s/^/$ACC /" >> logs/27_extract_genes.time.log
rm -rf "$TMP"
