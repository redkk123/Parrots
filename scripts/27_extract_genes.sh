#!/usr/bin/env bash
# Extract the Melopsittacus canonical proteome (data/genome_genes, 16,045 genes incl. the 101 H8 panel genes) from one
# assembly in a single miniprot run (scripts/33_prepare_proteome.py builds the query set).
# Usage: scripts/27_extract_genes.sh <GCA_accession> "<species>" [threads] [--download-only]   (prefetch)
# Download -> (chromosome-scale: 5 Mb windows overlapping by 500 kb, as in 10_audit_assembly.sh) -> split into parts of
# ~150 Mb, miniprot (vertebrate splice model) on one part at a time (a whole genome needs ~7 GB; a 142 Mb part used
# 2.0 GB in a test), IDs prefixed per part -> scripts/28_parse_miniprot.py on all parts. The genome is deleted after.
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
  for try in 1 2 3 4 5 6; do   # NCBI API outages are common; wait instead of failing the genome
    datasets download genome accession "$ACC" --include genome --filename "$TMP/pkg.zip" --no-progressbar && break
    rm -f "$TMP/pkg.zip"; [ "$try" -lt 6 ] && sleep 300
  done
  [ -s "$TMP/pkg.zip" ] || { echo "$ACC download failed"; exit 1; }
  unzip -p "$TMP/pkg.zip" 'ncbi_dataset/data/*/*.fna' > "$TMP/genome.fna.part"
  mv "$TMP/genome.fna.part" "$TMP/genome.fna"; rm "$TMP/pkg.zip"
fi
[ "${4:-}" = "--download-only" ] && exit 0
GEN=$TMP/genome.fna
if [ "$(seqkit stats -T "$GEN" | awk -F'\t' 'NR==2{print $8}')" -gt 20000000 ]; then
  seqkit sliding -g -W 5000000 -s 4500000 -j "$T" "$GEN" > "$TMP/windows.fna"; GEN=$TMP/windows.fna
fi
P=$(seqkit stats -T "$GEN" | awk -F'\t' 'NR==2{print int($5 / 150000000) + 1}')
seqkit split2 -p "$P" -O "$TMP/parts" "$GEN" 2> /dev/null
: > "$TMP/miniprot.gff"; : > "$TMP/time.txt"
i=0
for part in "$TMP"/parts/*; do
  i=$((i + 1))
  /usr/bin/time -v miniprot -t "$T" -K 500k -j 2 --gff --outn=3 -P "P${i}_" "$part" "$Q" \
    >> "$TMP/miniprot.gff" 2>> "$TMP/time.txt"
done
.venv/bin/python scripts/28_parse_miniprot.py "$ACC" "$SP" "$GEN" "$TMP/miniprot.gff"
grep -E "Maximum resident|Elapsed" "$TMP/time.txt" | sed "s/^/$ACC ($P parts) /" >> logs/27_extract_genes.time.log
rm -rf "$TMP"
