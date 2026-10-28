#!/usr/bin/env bash
# Relaunch the BUSCO audit batch and the hourly autosave after a reboot (called from the user's crontab, @reboot).
# Waits for the network, cleans the genome that was interrupted, and does nothing if the batch is already running
# or all assemblies in the queue are done. Log: logs/resume_on_boot.log
cd "$(dirname "$0")/.." || exit 1
LOG=logs/resume_on_boot.log
UIDN=$(id -u)
# git push uses the gh token in the desktop keyring: reachable once the user is logged in
export XDG_RUNTIME_DIR=/run/user/$UIDN DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/$UIDN/bus
export PATH=$HOME/.local/bin:/usr/local/bin:/usr/bin:/bin
echo "$(date -Is) boot: start" >> "$LOG"
for i in $(seq 1 60); do curl -s -m 10 -o /dev/null https://api.ncbi.nlm.nih.gov && break; sleep 10; done
q=results/phase1_audit/audit_queue.tsv; [ -s "$q" ] || q=data/external/audit_queue.tsv
Q2=results/phase1_audit/audit_queue_other.tsv   # non-Iridian assemblies (scripts/26_select_genomes.py)
total=$(cat "$q" "$Q2" 2>/dev/null | wc -l)
done_n=$(cut -f1 "$q" "$Q2" 2>/dev/null | while read -r a; do [ -s "results/phase1_audit/busco/$a/summary.txt" ] && echo; done | wc -l)
if [ "$done_n" -ge "$total" ]; then echo "$(date -Is) all $done_n/$total done; nothing to do" >> "$LOG"; exit 0; fi
if pgrep -f '[1]2_audit_batch.sh' > /dev/null; then echo "$(date -Is) batch already running" >> "$LOG"; exit 0; fi
for d in data/external/assemblies/*/; do
  [ -d "$d" ] || continue
  a=$(basename "$d"); rm -rf "$d"
  [ -s "results/phase1_audit/busco/$a/summary.txt" ] || rm -rf "results/phase1_audit/busco/$a"
done
# Own systemd units, so an out-of-memory kill in BUSCO (miniprot peaks at ~9.6 GB) does not take the desktop app with it
run() { systemd-run --user --unit="$1" --collect --working-directory="$PWD" --setenv=PATH="$PATH" "${@:2}" >> "$LOG" 2>&1; }
BOTH="bash scripts/12_audit_batch.sh; [ -s $Q2 ] && bash scripts/12_audit_batch.sh $Q2"   # finished ones are skipped
# systemd-run returns at once; not backgrounded, so it completes even when the caller is a unit that stops on exit
run parrots-busco -p OOMScoreAdjust=800 -p OOMPolicy=continue bash -c "$BOTH" \
  || { nohup setsid bash -c "$BOTH" > /dev/null 2>&1 & }
pgrep -f '[a]utosave_audit.sh' > /dev/null || run parrots-autosave bash scripts/autosave_audit.sh \
  || { nohup setsid bash scripts/autosave_audit.sh > /dev/null 2>&1 & }
echo "$(date -Is) relaunched (BUSCO $done_n/$total)" >> "$LOG"
