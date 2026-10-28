"""H8 Nível 1 summary by module (analysis plan §10): gene integrity (item 1), RELAX / aBSREL (items 2 and 4) and
module heterogeneity (item 3), per tree (supertree = main, Smith 2024 = sensitivity).

RELAX: BH-FDR within the panel; relaxed = q < 0.10 and k < 1, intensified = q < 0.10 and k > 1. Per module: genes
tested, proportion relaxed, mean log2 k. Heterogeneity: two-sided permutation of module labels (10,000), statistics
= variance of per-module proportion relaxed and of per-module mean log2 k.
Not done here (see docs/HANDOFF_GENES.md): the null from 1,000 random gene sets of the genome-wide universe,
RERconverge within clades, read-level checks of disrupted genes.
Writes results/phase2_genes/h8_relax_by_gene.tsv, h8_modules.tsv, h8_integrity.tsv.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
G = ROOT / "results/phase2_genes"
Q = 0.10
rng = np.random.default_rng(20261028)


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


def relax_row(js):
    r = json.load(open(js))["test results"]
    return {"relax_p": r["p-value"], "k": r["relaxation or intensification parameter"], "relax_LRT": r["LRT"]}


def absrel_row(js):
    d = json.load(open(js))
    br = [v for v in d["branch attributes"]["0"].values() if "Corrected P-value" in v and v.get("Corrected P-value") is not None]
    p = min((v["Corrected P-value"] for v in br), default=np.nan)
    return {"absrel_stem_p": p}


def heterogeneity(df, col, stat):
    obs = stat(df.groupby("module")[col])
    labels = df.module.to_numpy()
    perm = np.array([stat(df.assign(module=rng.permutation(labels)).groupby("module")[col]) for _ in range(10000)])
    return obs, (1 + (perm >= obs).sum()) / (1 + len(perm))


def main():
    genes = pd.read_csv(ROOT / "data/pilot_genes/canonical_genes.tsv", sep="\t").rename(columns={"gene_symbol": "gene"})
    rows = []
    for tree in ("burgio2019", "smith2024"):
        for js in sorted((G / "hyphy" / tree).glob("*.relax.json")):
            g = js.name.split(".")[0]
            r = {"tree": tree, "gene": g, **relax_row(js)}
            ab = js.with_name(f"{g}.absrel.json")
            if ab.exists():
                r.update(absrel_row(ab))
            rows.append(r)
    t = pd.DataFrame(rows)
    mods = []
    if len(t):
        t = t.merge(genes[["gene", "module"]], on="gene")
        t["relax_q"] = t.groupby("tree").relax_p.transform(bh)
        t["log2_k"] = np.log2(t.k.clip(lower=1e-3))
        t["relaxed"] = (t.relax_q < Q) & (t.k < 1)
        t["intensified"] = (t.relax_q < Q) & (t.k > 1)
        if "absrel_stem_p" in t:
            t["absrel_stem_q"] = t.groupby("tree").absrel_stem_p.transform(lambda p: bh(p.fillna(1)))
        t.to_csv(G / "h8_relax_by_gene.tsv", sep="\t", index=False, float_format="%.4g")
        for tree, d in t.groupby("tree"):
            m = d.groupby("module").agg(n_genes=("gene", "size"), prop_relaxed=("relaxed", "mean"),
                                        prop_intensified=("intensified", "mean"), mean_log2_k=("log2_k", "mean"))
            m["tree"] = tree
            for col in ("relaxed", "log2_k"):
                obs, p = heterogeneity(d, col, lambda gb: gb.mean().var())
                m[f"heterogeneity_{col}_p"] = p
            mods.append(m.reset_index())
        pd.concat(mods).to_csv(G / "h8_modules.tsv", sep="\t", index=False, float_format="%.4g")
    calls = pd.read_csv(G / "gene_calls.tsv", sep="\t", index_col="gene")
    sel = __import__("importlib.util").util
    spec = sel.spec_from_file_location("s", ROOT / "scripts/26_select_genomes.py")
    s = sel.module_from_spec(spec)
    spec.loader.exec_module(s)
    psl = [c for c in calls.columns if s.group(c.split("_")[0]) != "other_Psittaciformes"]
    rest = [c for c in calls.columns if c not in psl]
    integ = pd.DataFrame({
        "disrupted_psittaculidae": (calls[psl] == "disrupted").sum(1), "complete_psittaculidae": (calls[psl] == "complete").sum(1),
        "disrupted_rest": (calls[rest] == "disrupted").sum(1), "complete_rest": (calls[rest] == "complete").sum(1),
        "disrupted_melopsittacus": calls.get("Melopsittacus_undulatus", pd.Series(dtype=str)).eq("disrupted")})
    integ = integ.join(genes.set_index("gene").module)
    integ.to_csv(G / "h8_integrity.tsv", sep="\t")
    print(f"RELAX genes: {t.groupby('tree').size().to_dict() if len(t) else 0}; integrity table: {len(integ)} genes")
    if mods:
        print(pd.concat(mods).to_string(index=False))


if __name__ == "__main__":
    main()
