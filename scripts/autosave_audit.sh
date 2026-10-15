#!/usr/bin/env bash
# Commit and push finished audit results every hour while the batches run.
# Only completed outputs are added: BUSCO folders with summary.txt, identity JSONs, the consolidated table,
# batch logs and, once scripts/16_mash_batch.sh has finished, the Mash identity tables.
# Usage: nohup setsid bash scripts/autosave_audit.sh > /dev/null 2>&1 &     (one pass: --once)
set -uo pipefail
cd "$(dirname "$0")/.."
BRANCH=claude/new-session-lwh7xh
LOG=logs/autosave_audit.log

pass() {
  .venv/bin/python scripts/14_audit_summary.py > /dev/null 2>&1
  if ! pgrep -f '[1]6_mash_batch.sh' > /dev/null && ls data/external/mash/*.msh > /dev/null 2>&1; then
    .venv/bin/python scripts/17_mash_identity.py > /dev/null 2>&1 \
      && git add results/phase1_audit/mash_identity.tsv results/phase1_audit/mash_distances.tsv
  fi
  git add results/phase1_audit/audit_table.tsv results/phase1_audit/mito logs/1*_batch.log
  for d in results/phase1_audit/busco/*/; do [ -s "$d/summary.txt" ] && git add "$d"; done
  if git diff --cached --quiet; then
    echo "$(date -Is) nothing new" >> "$LOG"
    return
  fi
  n=$(ls results/phase1_audit/busco/*/summary.txt 2>/dev/null | wc -l)
  git commit -q -m "Audit autosave: BUSCO ${n}/96

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" \
    && git push -q origin "$BRANCH" \
    && echo "$(date -Is) pushed: BUSCO ${n}/96" >> "$LOG" \
    || echo "$(date -Is) commit/push FAILED" >> "$LOG"
}

if [ "${1:-}" = "--once" ]; then pass; exit; fi
while pgrep -f '[1]2_audit_batch.sh' > /dev/null || pgrep -f '[1]6_mash_batch.sh' > /dev/null; do
  pass
  sleep 3600
done
pass   # final save after the batches end
echo "$(date -Is) batches finished; autosave stopped" >> "$LOG"
