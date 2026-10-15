"""Does assembly quality track relative longevity? (TODO 1.1, item 4)

PGLS (Pagel's lambda, Burgio 2019 supertree) of each quality metric on relative log life expectancy and
relative log max longevity, over the Iridian assemblies of the analysis set. A quality metric that
predicts longevity must enter later gene-level models as a covariate.
Metrics: log10 contig N50, |total length - 1.2 Gb|, GC %, and BUSCO complete % (when available).
Writes results/phase1_audit/quality_vs_longevity.tsv
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from lib.phylo import Tree, fit_pgls  # noqa: E402

A = ROOT / "results/phase1_audit"


def busco_table():
    rows = []
    for p in (A / "busco").glob("*/summary.txt"):
        d = {"accession": p.parent.name}
        for line in p.read_text().splitlines():
            line = line.strip()
            if line[:2] in ("S:", "D:", "F:", "I:", "M:") and "%" in line:
                d[f"busco_{line[0]}"] = float(line.split(":")[1].split("%")[0])
        if "busco_S" in d:
            d["busco_C"] = d["busco_S"] + d.get("busco_D", 0.0)
            rows.append(d)
    return pd.DataFrame(rows)


def main():
    inv = pd.read_csv(A / "assembly_inventory.tsv", sep="\t")
    inv = inv[inv.iridian_n > 0].copy()
    res = pd.read_csv(ROOT / "results/phase0b/species_data_matrix_with_residuals.tsv", sep="\t")
    d = inv.merge(res[["species", "tree_tip", "rel_log_life_expectancy", "rel_log_max_longevity"]], on="species")
    d["log10_contig_n50"] = np.log10(d.iridian_contig_n50)
    d["abs_len_dev_gb"] = (d.iridian_total_len / 1e9 - 1.2).abs()
    d["gc"] = d.iridian_gc
    b = busco_table()
    metrics = ["log10_contig_n50", "abs_len_dev_gb", "gc"]
    if len(b):
        d = d.merge(b, left_on="iridian_accession", right_on="accession", how="left")
        metrics.append("busco_C")
    full = Tree.from_newick((ROOT / "results/phase0/supertree_burgio2019_full.nwk").read_text())
    out = []
    for y in ("rel_log_life_expectancy", "rel_log_max_longevity"):
        for x in metrics:
            s = d.dropna(subset=[x, y])
            if len(s) < 10:
                continue
            t = full.prune(set(s.tree_tip))
            order = t.tip_labels()
            s = s.set_index("tree_tip").loc[order]
            X = np.column_stack([np.ones(len(s)), s[x].to_numpy(float)])
            f = fit_pgls(s[y].to_numpy(float), X, t.vcv(order), "lambda")
            tval = f["beta"][1] / f["se"][1]
            p = 2 * stats.t.sf(abs(tval), len(s) - 2)
            out.append({"outcome": y, "metric": x, "n": len(s), "slope": f["beta"][1], "se": f["se"][1],
                        "t": tval, "p": p, "lambda": f["param"],
                        "spearman_rho": stats.spearmanr(s[x], s[y]).statistic})
    o = pd.DataFrame(out)
    o.to_csv(A / "quality_vs_longevity.tsv", sep="\t", index=False, float_format="%.4g")
    print(o.to_string(index=False, float_format=lambda v: f"{v:.3g}"))


if __name__ == "__main__":
    main()
