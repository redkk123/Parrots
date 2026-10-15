#!/usr/bin/env bash
# Nuclear identity check: Mash sketches (k=21, s=50000) of
#   (a) every Iridian assembly of the analysis set,
#   (b) every non-Iridian current assembly of an analysis-set species (independent references),
#   (c) the Iridian assembly of each species whose mitogenome was found in a 'conflict' assembly (donors).
# Sketches go to data/external/mash/; distances are computed by scripts/17_mash_identity.py.
# Waits for scripts/15_identity_batch.sh to finish (list (c) depends on it). Resumable.
set -uo pipefail
cd "$(dirname "$0")/.."
export PATH=$PWD/.bioenv/bin:$PATH
LOG=logs/16_mash_batch.log
while pgrep -f '[1]5_identity_batch.sh' > /dev/null; do sleep 60; done
mkdir -p data/external/mash
.venv/bin/python - <<'EOF' > data/external/mash_queue.tsv
import glob, json
import pandas as pd
flat = pd.read_csv("data/external/ncbi/psittaciformes_assemblies.tsv", sep="\t")
cur = flat[flat.status == "current"]
inv = pd.read_csv("results/phase1_audit/assembly_inventory.tsv", sep="\t")
q = [(a, "iridian_analysis") for a in inv.iridian_accession.dropna()]
q += [(a, "independent_reference") for a in cur[cur.species.isin(inv.species) & ~cur.iridian].accession]
donors = set()
for p in glob.glob("results/phase1_audit/mito/*.json"):
    r = json.load(open(p))
    if r["overall"] in ("conflict", "confirmed_with_foreign_mito"):
        donors |= {s.replace("_", " ") for s in r.get("foreign_species", [])}
d = cur[cur.iridian & cur.species.isin(donors)].sort_values("contig_n50", ascending=False).drop_duplicates("species")
q += [(a, "mito_donor") for a in d.accession]
seen = set()
for a, k in q:
    if a not in seen and not a.startswith("GCF_"):
        seen.add(a)
        print(f"{a}\t{k}")
EOF
while IFS=$'\t' read -r acc kind; do
  [ -s "data/external/mash/$acc.msh" ] && continue
  d=data/external/assemblies_mash/$acc; mkdir -p "$d"; start=$(date +%s)
  if datasets download genome accession "$acc" --include genome --filename "$d/pkg.zip" --no-progressbar < /dev/null \
     && unzip -p "$d/pkg.zip" 'ncbi_dataset/data/*/*.fna' > "$d/genome.fna" \
     && nice -n 5 mash sketch -k 21 -s 50000 -p 4 -o "data/external/mash/$acc" "$d/genome.fna" < /dev/null > /dev/null 2>&1
  then st=ok; else st=FAIL; rm -f "data/external/mash/$acc.msh"; fi
  rm -rf "$d"
  echo -e "$(date -Is)\t$acc\t$kind\t$st\t$(( $(date +%s) - start ))s" | tee -a "$LOG"
done < data/external/mash_queue.tsv
