"""Null for the H8 module results (analysis plan §10, Nível 1, item 2; 2026-10-28 entries), main tree only.

RELAX cannot run on the whole universe on this machine, so the null is built from a pool: N_POOL genes drawn at
random (fixed seed) from the eligible universe outside the panel (results/phase3_genome/genome_alignment_summary.tsv),
run with the same RELAX labelling as the panel (scripts/31_selection_tests.py). Each of 1,000 pseudo-panels replaces
every panel gene by a random pool gene of the same alignment-length and background-rate tertiles (background rate =
total tree length of the RELAX nucleotide fit); BH is applied within each pseudo-panel as in the real one. Per module:
observed proportion relaxed and mean log2 k against the pseudo-panels, two-sided empirical p.
Writes results/phase2_genes/h8_null_modules.tsv; RELAX JSON in results/phase2_genes/hyphy_null/ (gitignored).
"""
import gzip
import importlib.util
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
G = ROOT / "results/phase2_genes"
NULL = G / "hyphy_null" / "burgio2019"
N_POOL, N_SETS, Q = 300, 1000, 0.10
rng = np.random.default_rng(20261028)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


sel = load("sel", ROOT / "scripts/31_selection_tests.py")
summ = load("summ", ROOT / "scripts/32_h8_summary.py")


def tree_length(js):
    d = json.load(open(js))
    return sum(v.get("Nucleotide GTR", 0) or 0 for v in d["branch attributes"]["0"].values())


def run_pool(pool):
    tp, mp = sel.TREES["burgio2019"]
    full = sel.Tree.from_newick(tp.read_text())
    m = pd.read_csv(mp, sep="\t")
    tip_of = {s.replace(" ", "_"): t for s, t in zip(m.species, m.tree_tip)}
    NULL.mkdir(parents=True, exist_ok=True)
    for g in pool:
        js = NULL / f"{g}.relax.json"
        if js.exists() or (NULL / f"{g}.skip").exists():
            continue
        aln = NULL / f"{g}.codon.fasta"
        with gzip.open(ROOT / "data/genome_alignments/codon" / f"{g}.fa.gz", "rt") as f:
            aln.write_text(f.read())
        species = [line[1:].strip() for line in open(aln) if line.startswith(">")]
        nwk, status = sel.tree_for(full, tip_of, species)
        if nwk is None:
            (NULL / f"{g}.skip").write_text(status)
            continue
        (NULL / f"{g}.relax.nwk").write_text(nwk["relax"])
        r = subprocess.run([sel.HYPHY, f"CPU={sel.CPU}", "relax", "--alignment", aln, "--tree", NULL / f"{g}.relax.nwk",
                            "--test", "TEST", "--reference", "REFERENCE", "--models", "Minimal",
                            "--kill-zero-lengths", "No", "--output", js.with_suffix(".tmp")],
                           capture_output=True, text=True, cwd=NULL)
        if r.returncode == 0 and js.with_suffix(".tmp").exists():
            js.with_suffix(".tmp").replace(js)
        print("null", g, "ok" if js.exists() else "FAILED", flush=True)


def table(items):
    rows = []
    for g, js, aln in items:
        n = len(next(s for _, s in sel_fasta(aln))) // 3
        rows.append({"gene": g, **summ.relax_row(js), "aln_codons": n, "tree_length": tree_length(js)})
    return pd.DataFrame(rows)


def sel_fasta(path):
    name, buf = None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                yield name, "".join(buf)
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def module_stats(df):
    q = summ.bh(df.relax_p.to_numpy())
    d = df.assign(relaxed=(q < Q) & (df.k < 1), log2_k=np.log2(df.k.clip(lower=1e-3)))
    return d.groupby("module").agg(prop_relaxed=("relaxed", "mean"), mean_log2_k=("log2_k", "mean"))


def main():
    s = pd.read_csv(ROOT / "results/phase3_genome/genome_alignment_summary.tsv", sep="\t")
    panel_genes = set(pd.read_csv(ROOT / "data/pilot_genes/canonical_genes.tsv", sep="\t").gene_symbol)
    universe = sorted(set(s.gene[s.eligible & (s.error.fillna("") == "")]) - panel_genes)
    pool_file = G / "h8_null_pool.txt"
    if not pool_file.exists():
        pool_file.write_text("\n".join(sorted(rng.choice(universe, size=min(N_POOL, len(universe)), replace=False))))
    pool = pool_file.read_text().split()
    run_pool(pool)
    genes = pd.read_csv(ROOT / "data/pilot_genes/canonical_genes.tsv", sep="\t").rename(columns={"gene_symbol": "gene"})
    obs = table([(g, G / "hyphy/burgio2019" / f"{g}.relax.json", G / "alignments" / f"{g}.codon.fasta")
                 for g in genes.gene if (G / "hyphy/burgio2019" / f"{g}.relax.json").exists()])
    obs = obs.merge(genes[["gene", "module"]])
    nul = table([(g, NULL / f"{g}.relax.json", NULL / f"{g}.codon.fasta") for g in pool
                 if (NULL / f"{g}.relax.json").exists()])
    both = pd.concat([obs, nul])
    for c in ("aln_codons", "tree_length"):
        cuts = np.quantile(both[c], [1 / 3, 2 / 3])
        obs[c + "_bin"], nul[c + "_bin"] = np.digitize(obs[c], cuts), np.digitize(nul[c], cuts)
    strata = {k: v.index.to_numpy() for k, v in nul.groupby(["aln_codons_bin", "tree_length_bin"])}
    real = module_stats(obs)
    sims = []
    for _ in range(N_SETS):
        idx = [rng.choice(strata.get((r.aln_codons_bin, r.tree_length_bin), nul.index.to_numpy())) for r in obs.itertuples()]
        fake = nul.loc[idx, ["relax_p", "k"]].reset_index(drop=True).assign(module=obs.module.to_numpy())
        sims.append(module_stats(fake))
    out = []
    for mod, r in real.iterrows():
        row = {"module": mod, "n_genes": int((obs.module == mod).sum()), "n_pool": len(nul)}
        for c in ("prop_relaxed", "mean_log2_k"):
            null = np.array([s.loc[mod, c] for s in sims])
            row[c] = r[c]
            row[c + "_null_mean"] = null.mean()
            row[c + "_p"] = min(1.0, 2 * min((1 + (null >= r[c]).sum()), (1 + (null <= r[c]).sum())) / (1 + N_SETS))
        out.append(row)
    pd.DataFrame(out).to_csv(G / "h8_null_modules.tsv", sep="\t", index=False, float_format="%.4g")
    print(pd.DataFrame(out).to_string(index=False))


if __name__ == "__main__":
    main()
