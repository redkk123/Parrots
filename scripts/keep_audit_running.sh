#!/usr/bin/env bash
# Keep the genome audit going (user's crontab: @reboot and every 6 h).
#  - Relaunches the BUSCO batch and the hourly autosave if they stopped (reboot, crash, power cut).
#  - Once a week re-downloads the NCBI Psittaciformes report (09_ncbi_assembly_inventory.py --refresh), so assemblies
#    deposited later for the 105 analysis species enter the queues (14_audit_summary.py, 26_select_genomes.py).
#  - Does nothing while a batch runs. When every queued assembly is done (or failed 3 times; to retry one, delete its
#    FAIL lines from logs/12_audit_batch.log), starts the gene stage (scripts/29_gene_batch.sh) until
#    results/phase2_genes/DONE exists.
# Usage: bash scripts/keep_audit_running.sh [--refresh]      Log: logs/keep_audit.log
cd "$(dirname "$0")/.." || exit 1
LOG=logs/keep_audit.log
UIDN=$(id -u)
export XDG_RUNTIME_DIR=/run/user/$UIDN DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/$UIDN/bus
export PATH=$HOME/.local/bin:/usr/local/bin:/usr/bin:/bin
exec 9> logs/.keep_audit.lock; flock -n 9 || exit 0   # @reboot and the 6-hourly run must not overlap
log() { echo "$(date -Is) $*" >> "$LOG"; }
for i in $(seq 1 60); do curl -s -m 10 -o /dev/null https://api.ncbi.nlm.nih.gov && break; sleep 10; done

REPORT=data/external/ncbi/psittaciformes_dataset_report.json
if [ "${1:-}" = "--refresh" ] || [ -z "$(find "$REPORT" -mtime -7 2>/dev/null)" ]; then
  .venv/bin/python scripts/09_ncbi_assembly_inventory.py --refresh >> "$LOG" 2>&1 && log "NCBI inventory refreshed" \
    || log "NCBI refresh failed; keeping the previous inventory"
fi

running() { systemctl --user -q is-active parrots-busco || pgrep -f '^(/usr/bin/)?bash scripts/12_audit_batch[.]sh' > /dev/null; }
# systemd-run returns at once; own units, so an out-of-memory kill in BUSCO does not take the desktop app with it
run() { systemd-run --user --unit="$1" --collect --working-directory="$PWD" --setenv=PATH="$PATH" "${@:2}" >> "$LOG" 2>&1; }
autosave() {
  systemctl --user -q is-active parrots-autosave || pgrep -f '^(/usr/bin/)?bash scripts/autosave_audit[.]sh' > /dev/null \
    || run parrots-autosave bash scripts/autosave_audit.sh || { nohup setsid bash scripts/autosave_audit.sh > /dev/null 2>&1 & }
}
if running; then autosave; log "batch running"; exit 0; fi

.venv/bin/python scripts/14_audit_summary.py > /dev/null 2>&1; .venv/bin/python scripts/26_select_genomes.py > /dev/null 2>&1
Q1=results/phase1_audit/audit_queue.tsv; [ -s "$Q1" ] || Q1=data/external/audit_queue.tsv
Q2=results/phase1_audit/audit_queue_other.tsv   # non-Iridian assemblies still to audit (26_select_genomes.py)
left=$(cut -f1 "$Q1" "$Q2" 2>/dev/null | while read -r a; do
  [ -s "results/phase1_audit/busco/$a/summary.txt" ] || [ "$(grep -c "	$a	.*	FAIL	" logs/12_audit_batch.log)" -ge 3 ] || echo
done | wc -l)
if [ "$left" -eq 0 ]; then
  if [ -e results/phase2_genes/DONE ]; then log "audit and gene stage done"; exit 0; fi
  if systemctl --user -q is-active parrots-genes || pgrep -f '^(/usr/bin/)?bash scripts/29_gene_batch[.]sh' > /dev/null; then
    autosave; log "gene stage running"; exit 0
  fi
  rm -rf data/external/genes_tmp
  run parrots-genes -p OOMScoreAdjust=800 -p OOMPolicy=continue bash scripts/29_gene_batch.sh \
    || { nohup setsid bash scripts/29_gene_batch.sh > /dev/null 2>&1 & }
  autosave; log "audit done; gene stage started"; exit 0
fi
for d in data/external/assemblies/*/; do   # genome interrupted mid-run
  [ -d "$d" ] || continue
  a=$(basename "$d"); rm -rf "$d"
  [ -s "results/phase1_audit/busco/$a/summary.txt" ] || rm -rf "results/phase1_audit/busco/$a"
done
BOTH="bash scripts/12_audit_batch.sh; [ -s $Q2 ] && bash scripts/12_audit_batch.sh $Q2"   # finished ones are skipped
run parrots-busco -p OOMScoreAdjust=800 -p OOMPolicy=continue bash -c "$BOTH" \
  || { nohup setsid bash -c "$BOTH" > /dev/null 2>&1 & }
autosave
log "relaunched ($left assemblies to audit)"
