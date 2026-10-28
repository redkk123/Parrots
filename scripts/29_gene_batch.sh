#!/usr/bin/env bash
# Gene stage (analysis plan §10, Nível 1), started by scripts/keep_audit_running.sh once the genome audit is done:
#   27_extract_genes.sh on every selected genome (results/phase1_audit/species_genomes.tsv, status included*),
#   then 30_align_genes.py, then 31_selection_tests.sh (HyPhy RELAX / aBSREL), then 32_h8_summary.py.
# Resumable: each step skips what is done. An accession that fails 3 times is skipped. Log: logs/29_gene_batch.log
set -uo pipefail
cd "$(dirname "$0")/.."
LOG=logs/29_gene_batch.log
Q=data/external/gene_queue.tsv
.venv/bin/python scripts/26_select_genomes.py > /dev/null
.venv/bin/python - <<'EOF' > "$Q"
import pandas as pd
t = pd.read_csv("results/phase1_audit/species_genomes.tsv", sep="\t")
for r in t[t.status.str.startswith("included")].itertuples():
    print(f"{r.accession}\t{r.species}")
EOF
while IFS=$'\t' read -r acc sp; do
  if [ "$(grep -c "	$acc	.*	FAIL	" "$LOG" 2>/dev/null)" -ge 3 ]; then continue; fi
  [ -s "results/phase2_genes/hits/$acc.tsv" ] && continue
  start=$(date +%s)
  if nice -n 19 bash scripts/27_extract_genes.sh "$acc" "$sp" 4 < /dev/null >> "$LOG.detail" 2>&1; then st=ok; else st=FAIL; rm -rf "data/external/genes_tmp/$acc"; fi
  echo -e "$(date -Is)\t$acc\t$sp\t$st\t$(( $(date +%s) - start ))s" | tee -a "$LOG"
done < "$Q"
.venv/bin/python scripts/30_align_genes.py >> "$LOG" 2>&1 || { echo "$(date -Is) alignment FAILED" >> "$LOG"; exit 1; }
if [ -x scripts/31_selection_tests.sh ]; then
  bash scripts/31_selection_tests.sh >> "$LOG.detail" 2>&1 || { echo "$(date -Is) selection tests FAILED" >> "$LOG"; exit 1; }
  .venv/bin/python scripts/32_h8_summary.py >> "$LOG" 2>&1
fi
echo "$(date -Is) gene stage finished" >> "$LOG"
touch results/phase2_genes/DONE
