"""Smith, Thom & Joseph 2024 dated UCE tree as the alternative tree required by the analysis plan (§3).

Input: data/external/smith2024/tree/Time_calibrated_tree/MCC-...-revised-names.tre (Zenodo 10867437,
mirror of Dryad doi:10.5061/dryad.h70rxwdrx; treePL MCC tree, CC0).
1. NEXUS (Translate + [&...] annotations) -> Newick; drop outgroups; one tip per species.
2. Map Smith 2024 names to Burgio 2019 supertree species (NameResolver, then unique-epithet match, then the
   manual table below). Every mapping is written to results/phase1_trees/smith2024_name_map.tsv.
3. Recompute the Phase 0b allometries (PGLS, Pagel's lambda) on this tree and write
   results/phase1_trees/species_data_matrix_smith2024.tsv with the same residual columns.
Tips of results/phase1_trees/smith2024_species.nwk carry Burgio tree_tip labels so downstream scripts are unchanged.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from lib.names import NameResolver  # noqa: E402
from lib.phylo import Tree, fit_pgls  # noqa: E402

SRC = next((ROOT / "data/external/smith2024/tree/Time_calibrated_tree").glob("MCC-*.tre"))
OUT = ROOT / "results/phase1_trees"
OUTGROUPS = {"Caracara", "Calyptomena", "Icterus"}
SIZE_ORDER = {"tiny": 0, "small": 1, "medium": 2, "large": 3, "huge": 4}
# Smith 2024 binomial -> Burgio 2019 species, where neither the resolver nor a unique epithet decides.
# Checked against Supplementary_Table_S1.csv (column "Organism") on 2026-10-15.
MANUAL = {
    "Thectocercus acuticaudata": "Psittacara acuticaudatus",
    "Eunymphicus cornotus": "Eunymphicus cornutus",        # typo in the tree file
    "Eunymphicus uvaensis": "Eunymphicus uvaeensis",
    "Palaeornis calthorpae": "Psittacula calthrapae",
    "Himalayapsitta finschii": "Psittacula finschii",
    "Charmosynopsis pulchella": "Charmosyna pulchella",
    "Vini amabilis": "Charmosyna amabilis",
    "Vini meeki": "Charmosyna meeki",
    "Clarkona varia": "Psephotellus varius",               # S1: Psephotellus_varius (mulga parrot), not a lorikeet
    "Myiopsitta fuchsi": "Myiopsitta luchsi",              # typo in the tree file; S1: Myiopsitta_luchsi
    "Touit melanotus": "Touit melanonotus",
    "Psittacara strenuus": "Psittacara holochlorus",       # lumped in the supertree
}
# When several Smith species map to one Burgio species, keep the best-matching one.
PRIORITY = {"exact": 0, "subspecies_collapsed": 0, "spelling": 1, "burgio_synonym": 1, "genus_change": 1,
            "manual": 1, "unique_epithet": 2, "split_lumped_in_tree": 3}


def nexus_to_newick(text):
    tr = dict(re.findall(r"^\s*(\d+)\s+([^,;\s]+)[,;]?\s*$", text.split("Translate", 1)[1].split(";", 1)[0] + ";",
                         re.M))
    nwk = re.search(r"tree\s+\S+\s*=\s*(?:\[&[RU]\]\s*)?(.*;)", text, re.S).group(1)
    nwk = re.sub(r"\[&[^\]]*\]", "", nwk)
    nwk = re.sub(r"(?<=[(,])(\d+)(?=:)", lambda m: tr[m.group(1)], nwk)
    return nwk


def map_names(tips, R, tree_species):
    by_epithet = {}
    for s in tree_species:
        by_epithet.setdefault(s.split()[1], []).append(s)
    rows = []
    for tip in tips:
        parts = tip.split("_")
        binom = f"{parts[0]} {parts[1]}"
        if parts[0] in OUTGROUPS:
            continue
        sp, how = R.resolve(binom)
        if sp is None and binom in MANUAL:
            sp, how = MANUAL[binom], "manual"
        if binom == "Psittacara strenuus":
            how = "split_lumped_in_tree"
        if sp is None:
            c = by_epithet.get(parts[1], [])
            if len(c) == 1:
                sp, how = c[0], "unique_epithet"
            else:
                how = f"unresolved ({len(c)} epithet matches)"
        rows.append({"smith_tip": tip, "smith_binomial": binom, "burgio_species": sp, "map_type": how})
    d = pd.DataFrame(rows)
    d["priority"] = d.map_type.map(PRIORITY).fillna(9)
    best = d.dropna(subset=["burgio_species"]).sort_values(["priority", "smith_tip"]).drop_duplicates("burgio_species")
    d["kept"] = d.index.isin(best.index)
    return d.drop(columns="priority")


def fit_resid(tree, df, ycol, xcols):
    t = tree.prune(set(df["tree_tip"]))
    order = t.tip_labels()
    d = df.set_index("tree_tip").loc[order]
    X = np.column_stack([np.ones(len(d))] + [d[c].to_numpy(float) for c in xcols])
    f = fit_pgls(d[ycol].to_numpy(float), X, t.vcv(order), "lambda")
    return f, pd.Series(d[ycol].to_numpy(float) - X @ f["beta"], index=order)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    nwk = nexus_to_newick(SRC.read_text())
    tree = Tree.from_newick(nwk)
    R = NameResolver()
    burgio = pd.read_csv(ROOT / "results/phase0/species.tsv", sep="\t")
    nm = map_names(tree.tip_labels(), R, set(burgio.species))
    nm.to_csv(OUT / "smith2024_name_map.tsv", sep="\t", index=False)
    keep = nm[nm.kept]
    t = tree.prune(set(keep.smith_tip))
    lab = dict(zip(keep.smith_tip, keep.burgio_species.str.replace(" ", "_")))
    t.label = [lab.get(x, x) if x else x for x in t.label]
    (OUT / "smith2024_species.nwk").write_text(t.to_newick() + "\n")

    m = pd.read_csv(ROOT / "results/phase0b/species_data_matrix.tsv", sep="\t")
    m = m[m.tree_tip.isin(set(t.tip_labels()))].copy()
    s = {"source": SRC.name, "n_tips_smith": len(nm) + 3, "n_ingroup_mapped": int(keep.shape[0]),
         "map_types": nm.map_type.value_counts().to_dict()}
    m["log_max_longevity"] = np.log(m["max_longevity_years"])
    mm = m.dropna(subset=["log_max_longevity", "log_body_mass"])
    f_ml, r_ml = fit_resid(t, mm, "log_max_longevity", ["log_body_mass"])
    an = mm[mm["max_longevity_sample_size"].isin(SIZE_ORDER)].copy()
    an["sample_size_rank"] = an["max_longevity_sample_size"].map(SIZE_ORDER)
    _, r_ss = fit_resid(t, an, "log_max_longevity", ["log_body_mass", "sample_size_rank"])
    le = m.dropna(subset=["log_life_expectancy", "log_body_mass"])
    f_le, r_le = fit_resid(t, le, "log_life_expectancy", ["log_body_mass"])
    m["rel_log_max_longevity"] = m.tree_tip.map(r_ml)
    m["rel_log_max_longevity_ssadj"] = m.tree_tip.map(r_ss)
    m["rel_log_life_expectancy"] = m.tree_tip.map(r_le)
    m.to_csv(OUT / "species_data_matrix_smith2024.tsv", sep="\t", index=False, float_format="%.5g")

    old = pd.read_csv(ROOT / "results/phase0b/species_data_matrix_with_residuals.tsv", sep="\t").set_index("tree_tip")
    for name, f, col in (("life_expectancy", f_le, "rel_log_life_expectancy"),
                         ("max_longevity", f_ml, "rel_log_max_longevity")):
        both = m.set_index("tree_tip")[col].dropna()
        both = pd.concat([both, old[col]], axis=1, join="inner").dropna()
        s[f"{name}_allometry"] = {"n": int(f["V"].shape[0]), "lambda": f["param"], "slope": float(f["beta"][1]),
                                  "slope_se": float(f["se"][1]),
                                  "pearson_r_resid_vs_burgio": float(np.corrcoef(both.iloc[:, 0], both.iloc[:, 1])[0, 1])}
    json.dump(s, open(OUT / "smith2024_summary.json", "w"), indent=1, default=float)
    print(json.dumps(s, indent=1, default=float))
    print(nm[~nm.kept & nm.burgio_species.isna()].smith_binomial.tolist())


if __name__ == "__main__":
    main()
