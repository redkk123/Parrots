"""Phenotypes and master trees for RERconverge (analysis plan, 2026-10-28: H6 and H8 within clades).

Traits (natural-log scale, PGLS-lambda residuals on the supertree, as scripts/05_phase0b_analysis.py):
  rel_log_life_expectancy   already in results/phase0b/species_data_matrix_with_residuals.tsv
  rel_log_brain             log brain mass on log body mass, all species with both (results/phase0/phenotypes.tsv)
  rel_brain_given_le        rel_log_brain on rel_log_life_expectancy (brain "beyond longevity"), species with both
Master trees: supertree (main) and Smith 2024 (sensitivity), pruned to the species of results/phase1_audit/analysis_species.tsv (frozen 2026-10-28), tips named Genus_species as in the alignments.
Writes results/phase3_genome/{rer_traits.tsv, master_burgio2019.nwk, master_smith2024.nwk}.
"""
import importlib.util
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from lib.phylo import Tree  # noqa: E402

OUT = ROOT / "results/phase3_genome"
spec = importlib.util.spec_from_file_location("p0b", ROOT / "scripts/05_phase0b_analysis.py")
p0b = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p0b)
spec = importlib.util.spec_from_file_location("sel", ROOT / "scripts/26_select_genomes.py")
sel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sel)


def main():
    full = Tree.from_newick((ROOT / "results/phase0/supertree_burgio2019_full.nwk").read_text())
    ph = pd.read_csv(ROOT / "results/phase0/phenotypes.tsv", sep="\t")
    ph = ph[ph.tree_tip.isin(set(full.tip_labels()))]
    b = ph.dropna(subset=["log_brain_mass", "log_body_mass"])
    f_b, r_b, _ = p0b.subset_fit(full, b, "log_brain_mass", ["log_body_mass"])
    m = pd.read_csv(ROOT / "results/phase0b/species_data_matrix_with_residuals.tsv", sep="\t")
    t = m[["species", "tree_tip", "rel_log_life_expectancy"]].copy()
    t["rel_log_brain"] = t.tree_tip.map(r_b)
    both = t.dropna(subset=["rel_log_brain", "rel_log_life_expectancy"])
    f_c, r_c, _ = p0b.subset_fit(full, both, "rel_log_brain", ["rel_log_life_expectancy"])
    t["rel_brain_given_le"] = t.tree_tip.map(r_c)
    sg = pd.read_csv(ROOT / "results/phase1_audit/analysis_species.tsv", sep="\t")
    keep = set(sg.species[sg.status.str.startswith("included")])
    t = t[t.species.isin(keep)].copy()
    t["species_tip"] = t.species.str.replace(" ", "_")
    t["clade"] = ["Psittaculidae" if sel.group(s.split()[0]) != "other_Psittaciformes" else "rest" for s in t.species]
    OUT.mkdir(parents=True, exist_ok=True)
    t.drop(columns="tree_tip").to_csv(OUT / "rer_traits.tsv", sep="\t", index=False, float_format="%.5g")
    for name, tp, mp in (("burgio2019", ROOT / "results/phase0/supertree_burgio2019_full.nwk",
                          ROOT / "results/phase0b/species_data_matrix_with_residuals.tsv"),
                         ("smith2024", ROOT / "results/phase1_trees/smith2024_species.nwk",
                          ROOT / "results/phase1_trees/species_data_matrix_smith2024.tsv")):
        tr = Tree.from_newick(tp.read_text())
        mm = pd.read_csv(mp, sep="\t")
        tip_of = dict(zip(mm.species, mm.tree_tip))
        tips = {tip_of[s]: s.replace(" ", "_") for s in t.species if s in tip_of}
        pr = tr.prune(set(tips))
        pr.label = [tips.get(x, x) for x in pr.label]
        (OUT / f"master_{name}.nwk").write_text(pr.to_newick())
        print(name, len(pr.tips()), "tips")
    print(f"brain allometry: lambda {f_b['param']:.2f}, slope {f_b['beta'][1]:.3f} (n {len(b)}); "
          f"brain on longevity: slope {f_c['beta'][1]:.3f} ± {f_c['se'][1]:.3f} (n {len(both)})")
    print(f"selected species {len(t)}; with longevity {t.rel_log_life_expectancy.notna().sum()}, "
          f"with brain {t.rel_log_brain.notna().sum()}; Psittaculidae {(t.clade == 'Psittaculidae').sum()}")


if __name__ == "__main__":
    main()
