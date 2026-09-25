"""Follow-up of the regime-shift analysis: test the main shifted clade directly.

Reads results/phase0b/regime_shifts.json (scripts/R/regime_shifts.R), takes the largest shifted
clade, expands it to its MRCA in the full supertree and tests it with a PGLS indicator under
three outcomes, plus a leave-one-subclade-out check. Because the contrast is a single
evolutionary event, anything else that changed on the same branch is confounded with it.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent / "lib"))
from phylo import Tree, fit_pgls  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase0b"
OUTCOMES = ["rel_log_life_expectancy", "rel_log_max_longevity", "rel_log_max_longevity_ssadj"]
SUBCLADES = {
    "Loriini": ["Chalcopsitta", "Charmosyna", "Eos", "Glossopsitta", "Lorius", "Neopsittacus", "Oreopsittacus",
                "Phigys", "Pseudeos", "Psitteuteles", "Trichoglossus", "Vini"],
    "Platycercini + Pezoporini": ["Barnardius", "Cyanoramphus", "Eunymphicus", "Lathamus", "Neophema", "Neopsephotus",
                                  "Northiella", "Platycercus", "Psephotellus", "Psephotus", "Purpureicephalus",
                                  "Prosopeia", "Pezoporus"],
    "Psittaculini": ["Psittacula", "Tanygnathus", "Eclectus", "Psittinus", "Polytelis", "Alisterus", "Aprosmictus",
                     "Geoffroyus", "Prioniturus"],
    "Agapornis + Loriculus": ["Agapornis", "Loriculus"],
    "Cyclopsittini": ["Cyclopsitta", "Psittaculirostris", "Bolbopsittacus"],
    "Melopsittacus": ["Melopsittacus"],
}


def main():
    shifts = json.loads((OUT / "regime_shifts.json").read_text())["rel_LE"]["selection"]["BGHuni"]["shifts"]
    main_shift = max(shifts, key=lambda s: s["n_tips"])
    clade = set(main_shift["tips"].split(";"))
    full = Tree.from_newick((ROOT / "results/phase0/supertree_burgio2019_full.nwk").read_text())
    desc = full.descendants_tips()
    mrca = min((n for n in range(len(full.parent)) if clade <= {full.label[t] for t in desc[n]}),
               key=lambda n: len(desc[n]))
    members = {full.label[t] for t in desc[mrca]}
    m = pd.read_csv(OUT / "species_data_matrix_with_residuals.tsv", sep="\t")
    m["in_shift_clade"] = m["tree_tip"].isin(members).astype(float)

    def test(col, drop=()):
        d = m.dropna(subset=[col])
        d = d[~d["genus"].isin(drop)]
        t = full.prune(set(d["tree_tip"]))
        o = t.tip_labels()
        dd = d.set_index("tree_tip").loc[o]
        X = np.column_stack([np.ones(len(o)), dd["in_shift_clade"]])
        f = fit_pgls(dd[col].to_numpy(float), X, t.vcv(o), "lambda")
        tval = f["beta"][1] / f["se"][1]
        return {"n": len(o), "n_in_clade": int(dd["in_shift_clade"].sum()), "lambda": f["param"],
                "effect_log": f["beta"][1], "se": f["se"][1], "pct_effect": (np.exp(f["beta"][1]) - 1) * 100,
                "p": 2 * stats.t.sf(abs(tval), len(o) - 2)}

    rows = [{"outcome": c, "dropped": "none", **test(c)} for c in OUTCOMES]
    rows += [{"outcome": "rel_log_life_expectancy", "dropped": g, **test("rel_log_life_expectancy", gen)}
             for g, gen in SUBCLADES.items()]
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "shift_clade_contrast.tsv", sep="\t", index=False, float_format="%.4g")
    info = {"shift_value_phyloem": main_shift["shift_value"], "n_tips_in_fit": main_shift["n_tips"],
            "n_tips_mrca_full_tree": len(members),
            "genera": sorted({x.split("_")[0] for x in members})}
    (OUT / "shift_clade_members.json").write_text(json.dumps(info, indent=2))
    print(json.dumps(info))
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
