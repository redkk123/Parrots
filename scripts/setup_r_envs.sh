#!/usr/bin/env bash
# Create the two R environments used by scripts/R/ (conda-forge via micromamba).
#   .renv   : R 4.4 + ape, phytools, nlme, phylolm, geiger, caper   (validate_phase0.R)
#   .renv45 : R 4.5 + PhylogeneticEM                                 (regime_shifts.R)
set -euo pipefail
cd "$(dirname "$0")/.."
MM=${MICROMAMBA:-./.micromamba/micromamba}
if [ ! -x "$MM" ]; then
  mkdir -p .micromamba
  curl -sL -o .micromamba/micromamba https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64
  chmod +x .micromamba/micromamba
  MM=./.micromamba/micromamba
fi
export MAMBA_ROOT_PREFIX=$PWD/.micromamba/root
"$MM" create -y -q -p "$PWD/.renv" -c conda-forge r-base=4.4 r-ape r-phytools r-nlme r-phylolm r-geiger r-caper r-jsonlite
"$MM" create -y -q -p "$PWD/.renv45" -c conda-forge r-base=4.5 r-phylogeneticem r-ape r-phytools r-nlme r-phylolm r-jsonlite
echo "Rscript: .renv/bin/Rscript and .renv45/bin/Rscript"
