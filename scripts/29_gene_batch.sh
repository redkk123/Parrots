#!/usr/bin/env bash
# Gene stage (H8 §10 Nível 1 and H6), started by scripts/keep_audit_running.sh (priority over the remaining audit):
#   27_extract_genes.sh on every genome of the frozen list (results/phase1_audit/analysis_species.tsv;
#   whole Melopsittacus proteome, panel included), then 30_align_genes.py (panel), 34_align_genome.py (universe),
#   35_prepare_rer_inputs.py, 31_selection_tests.py (HyPhy RELAX / aBSREL), 32_h8_summary.py and, when the R package
#   is installed, R/rerconverge_h6_h8.R on both trees (H6 and H8 within clades); last, 37_h8_null.py (RELAX on a
#   300-gene pool of the universe, the longest step).
# Species list: results/phase1_audit/analysis_species.tsv (frozen 2026-10-28, plan history); the assemblies audited
# later enter a pre-declared extension (regenerate that file from species_genomes.tsv and rerun).
# Resumable: each step skips what is done. An accession that fails 3 times is skipped. Log: logs/29_gene_batch.log
set -uo pipefail
cd "$(dirname "$0")/.."
LOG=logs/29_gene_batch.log
Q=data/external/gene_queue.tsv
.venv/bin/python - <<'EOF' > "$Q"
import pandas as pd
t = pd.read_csv("results/phase1_audit/analysis_species.tsv", sep="\t")
for r in t.itertuples():
    print(f"{r.accession}\t{r.species}")
EOF
todo=()   # accession<TAB>species still to extract; the next genome is downloaded while the current one runs
while IFS=$'\t' read -r acc sp; do
  if [ "$(grep -c "	$acc	.*	FAIL	" "$LOG" 2>/dev/null)" -ge 3 ]; then continue; fi
  [ -s "results/phase2_genes/hits/$acc.tsv" ] && [ -s "results/phase3_genome/calls/$acc.tsv.gz" ] && continue
  todo+=("$acc	$sp")
done < "$Q"
for i in "${!todo[@]}"; do
  IFS=$'\t' read -r acc sp <<< "${todo[$i]}"
  if [ $((i + 1)) -lt ${#todo[@]} ]; then
    IFS=$'\t' read -r nacc nsp <<< "${todo[$((i + 1))]}"
    nice -n 19 bash scripts/27_extract_genes.sh "$nacc" "$nsp" 2 --download-only < /dev/null >> "$LOG.detail" 2>&1 &
    pre=$!
  else pre=""; fi
  start=$(date +%s)
  if nice -n 19 bash scripts/27_extract_genes.sh "$acc" "$sp" 8 < /dev/null >> "$LOG.detail" 2>&1; then st=ok; else st=FAIL; rm -rf "data/external/genes_tmp/$acc"; fi
  echo -e "$(date -Is)\t$acc\t$sp\t$st\t$(( $(date +%s) - start ))s" | tee -a "$LOG"
  [ -n "$pre" ] && wait "$pre"
done
.venv/bin/python scripts/30_align_genes.py >> "$LOG" 2>&1 || { echo "$(date -Is) alignment FAILED" >> "$LOG"; exit 1; }
.venv/bin/python scripts/34_align_genome.py >> "$LOG" 2>&1 || echo "$(date -Is) genome alignment FAILED" >> "$LOG"
.venv/bin/python scripts/35_prepare_rer_inputs.py >> "$LOG" 2>&1
.venv/bin/python scripts/31_selection_tests.py >> "$LOG.detail" 2>&1 || { echo "$(date -Is) selection tests FAILED" >> "$LOG"; exit 1; }
.venv/bin/python scripts/32_h8_summary.py >> "$LOG" 2>&1
if .renv45/bin/Rscript -e 'quit(status = !requireNamespace("RERconverge", quietly = TRUE))' 2> /dev/null; then
  for tr in burgio2019 smith2024; do
    .renv45/bin/Rscript scripts/R/rerconverge_h6_h8.R "$tr" >> "$LOG.detail" 2>&1 \
      || echo "$(date -Is) RERconverge $tr FAILED (see $LOG.detail)" >> "$LOG"
  done
else
  echo "$(date -Is) RERconverge not installed: H6 and H8-within-clade tests skipped (docs/HANDOFF_GENES.md)" >> "$LOG"
fi
.venv/bin/python scripts/37_h8_null.py >> "$LOG.detail" 2>&1 || echo "$(date -Is) H8 null FAILED" >> "$LOG"
echo "$(date -Is) gene stage finished" >> "$LOG"
touch results/phase2_genes/DONE
