#!/usr/bin/env bash
# Fast pass: mitochondrial identity / foreign-mitogenome check (scripts/11_mito_identity.py) for every Iridian
# assembly, independent of the slow BUSCO batch (scripts/12_audit_batch.sh), which skips this step when the
# JSON already exists. Uses its own download directory so the two batches never touch the same files.
set -uo pipefail
cd "$(dirname "$0")/.."
export PATH=$PWD/.bioenv/bin:$PATH
LOG=logs/15_identity_batch.log
[ -s data/external/audit_queue.tsv ] || { echo "run scripts/12_audit_batch.sh first (creates the queue)"; exit 1; }
while IFS=$'\t' read -r acc sp; do
  [ -s "results/phase1_audit/mito/$acc.json" ] && continue
  d=data/external/assemblies_id/$acc; mkdir -p "$d"; start=$(date +%s)
  if datasets download genome accession "$acc" --include genome --filename "$d/pkg.zip" --no-progressbar < /dev/null \
     && unzip -p "$d/pkg.zip" 'ncbi_dataset/data/*/*.fna' > "$d/genome.fna" \
     && nice -n 5 .venv/bin/python scripts/11_mito_identity.py "$acc" "$d/genome.fna" "$sp" < /dev/null >> "$LOG.detail" 2>&1
  then st=ok; else st=FAIL; fi
  rm -rf "$d"
  echo -e "$(date -Is)\t$acc\t$sp\t$st\t$(( $(date +%s) - start ))s" | tee -a "$LOG"
done < data/external/audit_queue.tsv
