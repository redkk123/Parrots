#!/usr/bin/env bash
# Fetch external inputs at pinned versions (see config/sources.tsv).
# Everything lands in data/external/, which is not versioned.
set -euo pipefail
cd "$(dirname "$0")/.."
EXT=data/external
mkdir -p "$EXT/genomeark" logs

clone_pinned () {  # url dir commit
  if [ ! -d "$EXT/$2/.git" ]; then git clone --quiet "$1" "$EXT/$2"; fi
  git -C "$EXT/$2" fetch --quiet --depth 1 origin "$3" 2>/dev/null || true
  git -C "$EXT/$2" checkout --quiet "$3"
  echo "$2 $(git -C "$EXT/$2" rev-parse HEAD)"
}

{
  echo "# fetched $(date -u +%FT%TZ)"
  clone_pinned https://github.com/simeonqs/Coevolution_of_relative_brain_size_and_life_expectancy_in_parrots smeele2022 311066c9a68ee2b4866754a2ef24e5da8f77df0e
  clone_pinned https://github.com/trashbirdecology/parrottraits parrottraits 79a8334db7677c9ab4e5d79e6c2e43fb102274b5
  clone_pinned https://github.com/RS-eco/traitdata traitdata 5560f5a609e2c891c23c36ae0cc25c56c81c64ac

  # GenomeArk: list every species folder, keep the parrot genera, then list assembly folders.
  python3 scripts/lib/genomeark_listing.py "$EXT/genomeark"
} | tee logs/00_fetch_external.log

(cd "$EXT" && find smeele2022/ANALYSIS/RESULTS parrottraits/data traitdata/data/an_age.rda traitdata/data/amniota.rda genomeark -type f \
   \( -name '*.RData' -o -name '*.rda' -o -name '*.tsv' -o -name '*.txt' \) -print0 \
   | sort -z | xargs -0 sha256sum) > logs/00_fetch_external.sha256
