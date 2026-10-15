"""Psittaculidae contrast with per-species measurement error (TODO §4).

Phase 0b and the Smith-tree sensitivity treated relative life expectancy as exact. Here the known standard
error of each species' log life expectancy (Smeele et al. 2022 posterior SD) enters the covariance:
V = sigma2 * C_lambda + diag(se^2) (scripts/lib/phylo.py::fit_pgls_me, validated against phytools).

  python scripts/19_measurement_error.py validate   # export data, run R reference, compare
  python scripts/19_measurement_error.py            # contrast on both trees, with and without ME
  python scripts/19_measurement_error.py se95       # matrices without species whose SE > 95th percentile
                                                    # (criterion fixed before looking), for regime_shifts.R

Outcome: rel_log_life_expectancy (the only outcome with per-species SE). Clade: Phase 0b shifted clade
(results/phase0b/regime_shifts.json), expanded to its MRCA on each tree, as in scripts/06_regime_followup.py.
Writes results/phase1_trees/measurement_error_contrast.tsv.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from lib.phylo import Tree, fit_pgls, fit_pgls_me  # noqa: E402

OUT = ROOT / "results/phase1_trees"
VAL = ROOT / "data/external/me_validation"
TREES = {
    "burgio2019": (ROOT / "results/phase0/supertree_burgio2019_full.nwk",
                   ROOT / "results/phase0b/species_data_matrix_with_residuals.tsv"),
    "smith2024": (OUT / "smith2024_species.nwk", OUT / "species_data_matrix_smith2024.tsv"),
}
SUBCLADES = {  # as in scripts/06_regime_followup.py
    "Loriini": ["Chalcopsitta", "Charmosyna", "Eos", "Glossopsitta", "Lorius", "Neopsittacus", "Oreopsittacus",
                "Phigys", "Pseudeos", "Psitteuteles", "Trichoglossus", "Vini"],
    "Psittaculini": ["Psittacula", "Tanygnathus", "Eclectus", "Psittinus", "Polytelis", "Alisterus", "Aprosmictus",
                     "Geoffroyus", "Prioniturus"],
}
Y, SE = "rel_log_life_expectancy", "log_life_expectancy_se"


def load(name):
    tpath, mpath = TREES[name]
    full = Tree.from_newick(Path(tpath).read_text())
    m = pd.read_csv(mpath, sep="\t").dropna(subset=[Y, SE])
    return full, m


def clade_members(full):
    shifts = json.loads((ROOT / "results/phase0b/regime_shifts.json").read_text())["rel_LE"]["selection"]["BGHuni"]["shifts"]
    clade = set(max(shifts, key=lambda s: s["n_tips"])["tips"].split(";"))
    desc = full.descendants_tips()
    labels = set(full.tip_labels())
    clade &= labels
    mrca = min((n for n in range(len(full.parent)) if clade <= {full.label[t] for t in desc[n]}),
               key=lambda n: len(desc[n]))
    return {full.label[t] for t in desc[mrca]}


def validate():
    full, m = load("smith2024")
    m["x"] = m["tree_tip"].isin(clade_members(full)).astype(float)
    t = full.prune(set(m.tree_tip))
    o = t.tip_labels()
    d = m.set_index("tree_tip").loc[o]
    VAL.mkdir(parents=True, exist_ok=True)
    (VAL / "tree.nwk").write_text(t.to_newick(digits=10) + "\n")
    pd.DataFrame({"tip": o, "y": d[Y].values, "x": d.x.values, "se": d[SE].values}).to_csv(VAL / "data.tsv", sep="\t",
                                                                                          index=False)
    subprocess.run([str(ROOT / ".renv45/bin/Rscript"), "scripts/R/validate_pgls_me.R"], cwd=ROOT, check=True,
                   capture_output=True)
    r = json.loads((VAL / "r_results.json").read_text())
    C = t.vcv(o)
    X = np.column_stack([np.ones(len(o)), d.x.values])
    a = fit_pgls_me(d[Y].values, X, C, d[SE].values, model="BM")
    c = fit_pgls_me(d[Y].values, X, C, d[SE].values, model="lambda")
    b = fit_pgls_me(d[Y].values, np.ones((len(o), 1)), C, d[SE].values, model="lambda")
    py = lambda f: {"beta": f["beta"].tolist(), "se": f["se"].tolist(), "lambda": f["param"], "loglik": f["loglik"]}
    res = {"n": len(o),
           "bm_me": {"python": py(a), "R_direct_ML": r["bm_me"]},
           "lambda_me_regression": {"python": py(c), "R_direct_ML": r["lambda_me_regression"]},
           "lambda_me_intercept": {"python": {"lambda": b["param"], "sigma2": b["sigma2"], "loglik": b["loglik"]},
                                   "R_phytools_phylosig": r["lambda_me"]},
           "note": "phytools::pgls.SEy re-estimates a residual scale (different model); not used as reference",
           "R_pgls_SEy": r["pgls_SEy_different_model"]}
    diff = lambda u, v: float(np.max(np.abs(np.atleast_1d(np.array(u, float)) - np.atleast_1d(np.array(v, float)))))
    res["max_abs_diff"] = {
        "bm_me_beta": diff(r["bm_me"]["beta"], a["beta"]), "bm_me_se": diff(r["bm_me"]["se"], a["se"]),
        "bm_me_loglik": diff(r["bm_me"]["loglik"], a["loglik"]),
        "lambda_me_reg_beta": diff(r["lambda_me_regression"]["beta"], c["beta"]),
        "lambda_me_reg_se": diff(r["lambda_me_regression"]["se"], c["se"]),
        "lambda_me_reg_lambda": diff(r["lambda_me_regression"]["lambda"], c["param"]),
        "lambda_me_reg_loglik": diff(r["lambda_me_regression"]["loglik"], c["loglik"]),
        "phylosig_lambda": diff(r["lambda_me"]["lambda"], b["param"]),
        "phylosig_loglik": diff(r["lambda_me"]["loglik"], b["loglik"])}
    (OUT / "validation_pgls_me.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


def contrast():
    rows = []
    for name in TREES:
        full, m = load(name)
        members = clade_members(full)
        m["in_clade"] = m.tree_tip.isin(members).astype(float)
        for drop, gen in [("none", [])] + list(SUBCLADES.items()):
            d0 = m[~m.genus.isin(gen)]
            t = full.prune(set(d0.tree_tip))
            o = t.tip_labels()
            d = d0.set_index("tree_tip").loc[o]
            X = np.column_stack([np.ones(len(o)), d.in_clade.values])
            C = t.vcv(o)
            for me, f in (("no", fit_pgls(d[Y].values, X, C, "lambda")),
                          ("yes", fit_pgls_me(d[Y].values, X, C, d[SE].values, "lambda"))):
                tv = f["beta"][1] / f["se"][1]
                rows.append({"tree": name, "measurement_error": me, "dropped": drop, "n": len(o),
                             "n_in_clade": int(d.in_clade.sum()), "lambda": f["param"], "effect_log": f["beta"][1],
                             "se": f["se"][1], "pct_effect": (np.exp(f["beta"][1]) - 1) * 100,
                             "p": 2 * stats.t.sf(abs(tv), len(o) - 2), "aic": f["aic"]})
    r = pd.DataFrame(rows)
    r.to_csv(OUT / "measurement_error_contrast.tsv", sep="\t", index=False, float_format="%.4g")
    print(r.to_string(index=False, float_format=lambda v: f"{v:.4g}"))


def se95():
    for name, (_, mpath) in TREES.items():
        m = pd.read_csv(mpath, sep="\t")
        drop = m[SE] > m[SE].quantile(0.95)
        m.loc[drop, Y] = np.nan
        m.to_csv(ROOT / f"data/external/matrix_{name}_se95.tsv", sep="\t", index=False)
        print(name, int(drop.sum()), "species dropped:", sorted(pd.read_csv(mpath, sep="\t")[drop].species))


if __name__ == "__main__":
    {"validate": validate, "se95": se95}.get(sys.argv[1] if len(sys.argv) > 1 else "", contrast)()
