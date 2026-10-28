#!/usr/bin/env bash
# Audit one NCBI assembly: download -> sequence stats -> mitochondrial identity -> BUSCO-equivalent completeness.
# Usage: scripts/10_audit_assembly.sh <GCA_accession> "<claimed species>" [threads]
# Outputs (versioned, small):
#   results/phase1_audit/busco/<acc>/{seqkit_stats.tsv,summary.txt,full_table.tsv,time.txt}
#   results/phase1_audit/mito/<acc>.json
# The FASTA is downloaded to data/external/assemblies/<acc>/ and deleted afterwards (disk budget).
# Completeness: compleasm (miniprot + hmmsearch) with BUSCO aves_odb10, as fixed in the protocol
# (compleasm >= 0.2.7 defaults to odb12, so --odb odb10 is passed explicitly).
# Chromosome-scale assemblies (longest sequence > 20 Mb) are cut into 5 Mb windows overlapping by 500 kb before
# compleasm: miniprot memory grows with sequence length and passed 10 GB (OOM on 11 GB RAM) on 100-200 Mb chromosomes.
# Genes in an overlap are found twice, so S moves to D; C = S + D is unaffected.
set -euo pipefail
cd "$(dirname "$0")/.."
ACC=$1; SP=$2; T=${3:-8}
BIN=$PWD/.bioenv/bin; export PATH=$BIN:$PATH
LIB=data/external/busco_lineages; OUT=results/phase1_audit/busco/$ACC; TMP=data/external/assemblies/$ACC
MITO=results/phase1_audit/mito/$ACC.json
if [ -s "$OUT/summary.txt" ] && [ -s "$MITO" ]; then echo "$ACC done"; exit 0; fi
mkdir -p "$LIB" "$OUT" "$TMP"
[ -d "$LIB/aves_odb10" ] || compleasm download aves --odb odb10 -L "$LIB"
if [ ! -s "$TMP/genome.fna" ]; then
  datasets download genome accession "$ACC" --include genome --filename "$TMP/pkg.zip" --no-progressbar
  unzip -p "$TMP/pkg.zip" 'ncbi_dataset/data/*/*.fna' > "$TMP/genome.fna"
  rm "$TMP/pkg.zip"
fi
[ -s "$OUT/seqkit_stats.tsv" ] || seqkit stats -a -T -j "$T" "$TMP/genome.fna" > "$OUT/seqkit_stats.tsv"
[ -s "$MITO" ] || .venv/bin/python scripts/11_mito_identity.py "$ACC" "$TMP/genome.fna" "$SP" > "$OUT/mito.log"
if [ ! -s "$OUT/summary.txt" ]; then
  rm -rf "$OUT/run"   # a killed run leaves partial miniprot output that compleasm would reuse
  GEN=$TMP/genome.fna
  if [ "$(awk -F'\t' 'NR==2{print $8}' "$OUT/seqkit_stats.tsv")" -gt 20000000 ]; then
    seqkit sliding -g -W 5000000 -s 4500000 -j "$T" "$GEN" > "$TMP/windows.fna"
    GEN=$TMP/windows.fna; echo "windows 5 Mb / overlap 500 kb" > "$OUT/chunked.txt"
  fi
  /usr/bin/time -v compleasm run -a "$GEN" -o "$OUT/run" -l aves --odb odb10 -L "$LIB" -t "$T" \
    2> "$OUT/time.txt" > "$OUT/compleasm.log"
  cp "$OUT/run/summary.txt" "$OUT/summary.txt"
  cp "$OUT/run/aves_odb10/full_table.tsv" "$OUT/full_table.tsv" 2>/dev/null || true
  rm -rf "$OUT/run"
fi
rm -rf "$TMP"
echo "$ACC ok"
