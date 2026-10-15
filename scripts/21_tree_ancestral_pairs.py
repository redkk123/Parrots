"""Tree sensitivity, part 2: BM ancestral states and sister-pair contrasts on the Burgio 2019 supertree vs the
Smith et al. 2024 UCE tree (each with its own allometric residuals; see scripts/18_smith_tree.py).

Ancestral states: root, crown Psittaculidae (MRCA of the Phase 0b shifted clade) and the change along its stem.
Sister pairs: cherries among the gene-level usable species with life expectancy (as scripts/05_phase0b_analysis.py).
Writes results/phase1_trees/ancestral_pairs_by_tree.json and results/phase1_trees/sister_pairs_smith2024.tsv.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from lib.phylo import Tree, bm_ancestral  # noqa: E402

OUT = ROOT / "results/phase1_trees"
TREES = {
    "burgio2019": (ROOT / "results/phase0/supertree_burgio2019_full.nwk",
                   ROOT / "results/phase0b/species_data_matrix_with_residuals.tsv"),
    "smith2024": (OUT / "smith2024_species.nwk", OUT / "species_data_matrix_smith2024.tsv"),
}
Y = "rel_log_life_expectancy"


def mrca(t, labels):
    desc = t.descendants_tips()
    return min((n for n in range(len(t.parent)) if labels <= {t.label[x] for x in desc[n]}), key=lambda n: len(desc[n]))


def main():
    shifts = json.loads((ROOT / "results/phase0b/regime_shifts.json").read_text())["rel_LE"]["selection"]["BGHuni"]["shifts"]
    clade = set(max(shifts, key=lambda s: s["n_tips"])["tips"].split(";"))
    summary, pair_sets = {}, {}
    for name, (tp, mp) in TREES.items():
        full = Tree.from_newick(Path(tp).read_text())
        m = pd.read_csv(mp, sep="\t")
        d = m.dropna(subset=[Y])
        t = full.prune(set(d.tree_tip))
        est, se, mu, s2 = bm_ancestral(t, dict(zip(d.tree_tip, d[Y])))
        c = mrca(t, clade & set(t.tip_labels()))
        p = t.parent[c]
        summary[name] = {"n_species": len(d), "root_est": float(est[t.root]), "root_se": float(se[t.root]),
                         "crown_psittaculidae_est": float(est[c]), "crown_psittaculidae_se": float(se[c]),
                         "stem_parent_est": float(est[p]), "stem_change": float(est[c] - est[p]),
                         "stem_branch_myr": float(t.length[c])}
        use = d[d.usable_for_gene_level]
        tu = full.prune(set(use.tree_tip))
        depth = tu.depths()
        H = depth[tu.tips()].max()
        rel = dict(zip(use.tree_tip, use[Y]))
        pairs = []
        for n in range(len(tu.parent)):
            kids = tu.children[n]
            if len(kids) == 2 and all(tu.is_tip(k) for k in kids):
                a, b = sorted((tu.label[k] for k in kids), key=lambda x: -rel[x])
                pairs.append({"longer_lived": a.replace("_", " "), "shorter_lived": b.replace("_", " "),
                              "divergence_myr": H - depth[n], "delta_rel_logLE": rel[a] - rel[b]})
        pairs = pd.DataFrame(pairs).sort_values("delta_rel_logLE", ascending=False)
        pair_sets[name] = {frozenset((r.longer_lived, r.shorter_lived)) for r in pairs.itertuples()}
        summary[name]["n_sister_pairs"] = len(pairs)
        summary[name]["n_pairs_delta_ge_0.3"] = int((pairs.delta_rel_logLE >= 0.3).sum())
        if name == "smith2024":
            pairs.to_csv(OUT / "sister_pairs_smith2024.tsv", sep="\t", index=False, float_format="%.3g")
    shared = pair_sets["burgio2019"] & pair_sets["smith2024"]
    summary["sister_pairs_shared"] = len(shared)
    summary["sister_pairs_only_burgio"] = sorted(" / ".join(sorted(p)) for p in pair_sets["burgio2019"] - shared)
    summary["sister_pairs_only_smith"] = sorted(" / ".join(sorted(p)) for p in pair_sets["smith2024"] - shared)
    (OUT / "ancestral_pairs_by_tree.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
