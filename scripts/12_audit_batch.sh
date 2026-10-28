#!/usr/bin/env bash
# Run scripts/10_audit_assembly.sh over every Iridian assembly in results/phase1_audit/assembly_inventory.tsv.
# Sequential (one genome at a time: 11 GB RAM, ~37 GB free disk). Resumable: finished accessions are skipped.
# Order is shuffled with a fixed seed so that a partial run is an unbiased sample of the batch.
# With an argument, audits that queue instead (accession<TAB>species; e.g. results/phase1_audit/audit_queue_other.tsv
# from scripts/26_select_genomes.py), read from a snapshot so the source can be regenerated meanwhile.
set -uo pipefail
cd "$(dirname "$0")/.."
LOG=logs/12_audit_batch.log
Q=data/external/audit_queue.tsv
if [ -n "${1:-}" ]; then Q=data/external/$(basename "$1"); cp "$1" "$Q"; else
.venv/bin/python - <<'EOF' > "$Q"
import pandas as pd
t = pd.read_csv("results/phase1_audit/assembly_inventory.tsv", sep="\t")
t = t[t.iridian_n > 0].sample(frac=1, random_state=20261014)
for r in t.itertuples():
    print(f"{r.iridian_accession}\t{r.species}")
EOF
fi
while IFS=$'\t' read -r acc sp; do
  start=$(date +%s)
  if nice -n 19 bash scripts/10_audit_assembly.sh "$acc" "$sp" 8 < /dev/null >> "$LOG.detail" 2>&1; then st=ok; else st=FAIL; rm -rf "data/external/assemblies/$acc"; fi
  echo -e "$(date -Is)\t$acc\t$sp\t$st\t$(( $(date +%s) - start ))s" | tee -a "$LOG"
done < "$Q"
