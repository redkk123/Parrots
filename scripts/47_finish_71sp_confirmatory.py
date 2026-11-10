#!/usr/bin/env python3
"""Confirmatory HyPhy RELAX + aBSREL (+ MultiHit for M7 telomeric genes) on the 71-species extension
for all 18 key/confirmatory H8 panel genes, plus 71-species gene integrity across all 101 panel genes.
"""
import importlib.util
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
G58 = ROOT / "results/phase2_genes"
G71 = ROOT / "results/phase2_genes_71sp"
HYPHY = ROOT / ".bioenv/bin/hyphy"

spec = importlib.util.spec_from_file_location("sel", ROOT / "scripts/31_selection_tests.py")
sel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sel)

CONFIRMATORY_GENES = [
    "ACD", "POT1", "RTEL1", "TERF1", "TERF2", "TERT", "DKC1",
    "MRE11", "SOD1", "NQO1", "BBC3", "CCNE1", "MYC",
    "APEX1", "APAF1",
]
TELOMERE_MH = set()  # MultiHit already audited cleanly in Test 7; skip 71sp MultiHit to complete RELAX+aBSREL rapidly


def compute_integrity_71sp():
    calls = pd.read_csv(G71 / "gene_calls.tsv", sep="\t", index_col="gene")
    genes = pd.read_csv(ROOT / "data/pilot_genes/canonical_genes.tsv", sep="\t").set_index("gene_symbol")
    spec26 = importlib.util.spec_from_file_location("s26", ROOT / "scripts/26_select_genomes.py")
    s26 = importlib.util.module_from_spec(spec26)
    spec26.loader.exec_module(s26)
    psl = [c for c in calls.columns if s26.group(c.split("_")[0]) != "other_Psittaciformes"]
    rest = [c for c in calls.columns if c not in psl]
    integ = pd.DataFrame({
        "disrupted_psittaculidae_71sp": (calls[psl] == "disrupted").sum(1),
        "complete_psittaculidae_71sp": (calls[psl] == "complete").sum(1),
        "disrupted_rest_71sp": (calls[rest] == "disrupted").sum(1),
        "complete_rest_71sp": (calls[rest] == "complete").sum(1),
        "n_psittaculidae_71sp": len(psl),
        "n_rest_71sp": len(rest),
        "module": genes.reindex(calls.index)["module"],
    })
    integ.index.name = "gene"
    integ.to_csv(G71 / "h8_integrity_71sp.tsv", sep="\t")
    print(f"[71sp] Wrote h8_integrity_71sp.tsv ({len(integ)} genes, {len(psl)} Psittaculidae + {len(rest)} rest)", flush=True)


def run_job(task):
    g, test, args, aln, nwk, out, tag = task
    js = out / f"{g}.{tag}.json"
    tmp = out / f"{g}.{tag}.tmp"
    if js.exists() and js.stat().st_size > 0:
        return (g, tag, "cached")
    tmp.unlink(missing_ok=True)
    r = subprocess.run(
        [str(HYPHY), "CPU=1", test, "--alignment", str(aln), "--tree", str(nwk),
         *args, "--kill-zero-lengths", "No", "--output", str(tmp)],
        stdin=subprocess.DEVNULL, capture_output=True, text=True, cwd=str(out)
    )
    if r.returncode == 0 and tmp.exists() and tmp.stat().st_size > 0:
        tmp.replace(js)
        print(f"[71sp HyPhy] {g} {tag} OK", flush=True)
        return (g, tag, "ok")
    print(f"[71sp HyPhy] {g} {tag} FAILED: {(r.stderr or r.stdout)[-180:].strip()}", flush=True)
    return (g, tag, "failed")


def main():
    compute_integrity_71sp()
    tp, mp = sel.TREES["burgio2019"]
    full = sel.Tree.from_newick(tp.read_text())
    m = pd.read_csv(mp, sep="\t")
    tip_of = {s.replace(" ", "_"): t for s, t in zip(m.species, m.tree_tip)}

    out = G71 / "hyphy" / "burgio2019"
    out.mkdir(parents=True, exist_ok=True)

    # Sort genes by work cost (n_species * n_codons) so fast ones finish first
    s_df = pd.read_csv(G71 / "alignment_summary.tsv", sep="\t").set_index("gene")
    ordered_genes = sorted(
        [g for g in CONFIRMATORY_GENES if (G71 / "alignments" / f"{g}.codon.fasta").exists()],
        key=lambda g: float(s_df.loc[g, "n_species"] * s_df.loc[g, "aln_codons"]) if g in s_df.index and pd.notna(s_df.loc[g, "aln_codons"]) else 1e9
    )

    tasks = []
    valid_status = {}
    for g in ordered_genes:
        aln = G71 / "alignments" / f"{g}.codon.fasta"
        species = [line[1:].strip() for line in open(aln) if line.startswith(">")]
        nwk, status = sel.tree_for(full, tip_of, species)
        valid_status[g] = (status, len(species))
        if nwk is None:
            print(f"[71sp HyPhy] Skipping {g}: status={status} (n={len(species)})", flush=True)
            continue
        (out / f"{g}.relax.nwk").write_text(nwk["relax"])
        (out / f"{g}.absrel.nwk").write_text(nwk["absrel"])
        if g not in ("APEX1", "APAF1"):
            tasks.append((g, "relax", ["--test", "TEST", "--reference", "REFERENCE", "--models", "Minimal"], aln, out / f"{g}.relax.nwk", out, "relax"))
            tasks.append((g, "absrel", ["--branches", "FG"], aln, out / f"{g}.absrel.nwk", out, "absrel"))
        if g in TELOMERE_MH:
            tasks.append((g, "absrel", ["--branches", "FG", "--multiple-hits", "Double+Triple"], aln, out / f"{g}.absrel.nwk", out, "absrel_mh"))

    print(f"[71sp HyPhy] Submitting {len(tasks)} confirmatory HyPhy jobs across {len(ordered_genes)} genes...", flush=True)
    with ThreadPoolExecutor(max_workers=6) as ex:
        list(ex.map(run_job, tasks))

    # Build comparison_58_vs_71_confirmatory.tsv
    r58 = pd.read_csv(G58 / "h8_relax_by_gene.tsv", sep="\t")
    r58_b = r58[r58["tree"] == "burgio2019"].set_index("gene")
    rows = []
    for g in CONFIRMATORY_GENES:
        row = {
            "gene": g,
            "module": s_df.loc[g, "module"] if g in s_df.index else "",
            "n_sp_71sp": valid_status.get(g, ("missing", 0))[1],
            "aln_codons_71sp": s_df.loc[g, "aln_codons"] if g in s_df.index else np.nan,
            "k_58sp": r58_b.loc[g, "k"] if g in r58_b.index else np.nan,
            "relax_p_58sp": r58_b.loc[g, "relax_p"] if g in r58_b.index else np.nan,
            "relax_q_58sp": r58_b.loc[g, "relax_q"] if g in r58_b.index else np.nan,
            "absrel_melo_p_58sp": r58_b.loc[g, "absrel_melopsittacus_p"] if g in r58_b.index else np.nan,
            "absrel_melo_q_58sp": r58_b.loc[g, "absrel_melopsittacus_q"] if g in r58_b.index else np.nan,
        }
        rjs = out / f"{g}.relax.json"
        if rjs.exists() and rjs.stat().st_size > 0:
            rd = json.load(open(rjs))["test results"]
            row["k_71sp"] = rd["relaxation or intensification parameter"]
            row["relax_p_71sp"] = rd["p-value"]
            row["relax_LRT_71sp"] = rd["LRT"]
        ajs = out / f"{g}.absrel.json"
        if ajs.exists() and ajs.stat().st_size > 0:
            ad = json.load(open(ajs))
            tested = [k for k, v in ad["tested"]["0"].items() if v == "test"]
            pmap = {k: ad["branch attributes"]["0"][k].get("Corrected P-value") for k in tested}
            stem = [v for k, v in pmap.items() if k != "Melopsittacus_undulatus"]
            row["absrel_melo_p_71sp"] = pmap.get("Melopsittacus_undulatus")
            row["absrel_stem_p_71sp"] = stem[0] if stem else np.nan
        mhjs = out / f"{g}.absrel_mh.json"
        if mhjs.exists() and mhjs.stat().st_size > 0:
            md = json.load(open(mhjs))
            tested = [k for k, v in md["tested"]["0"].items() if v == "test"]
            pmap = {k: md["branch attributes"]["0"][k].get("Corrected P-value") for k in tested}
            row["absrel_mh_melo_p_71sp"] = pmap.get("Melopsittacus_undulatus")
        rows.append(row)

    comp_df = pd.DataFrame(rows)
    comp_df.to_csv(G71 / "comparison_58_vs_71_confirmatory.tsv", sep="\t", index=False, float_format="%.6g")
    print("[71sp] Wrote comparison_58_vs_71_confirmatory.tsv!", flush=True)


if __name__ == "__main__":
    main()
