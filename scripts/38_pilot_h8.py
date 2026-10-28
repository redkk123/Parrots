"""Fast H8 pilot (analysis plan §9 and Histórico 2026-10-28): the 101-gene panel on the 8 chromosome-scale parrot
genomes already on disk (data/pilot_genomes/fasta, Gemini's pilot panel). Exploratory, not confirmatory: 4
Psittaculidae (Melopsittacus, Lathamus, Psittacula echo, Agapornis nigrigenis) x 4 others (Strigops, Ara ararauna,
Amazona ochrocephala, Myiopsitta). Same rules as the main analysis: miniprot with the Melopsittacus and kakapo panel
proteins (28_parse_miniprot.py), MAFFT L-INS-i + PAL2NAL (30_align_genes.align), RELAX / aBSREL on the supertree
(31_selection_tests.tree_for; Psittacula echo takes the tip of its congener in the tree), BH within the panel.
Writes results/phase2_pilot/{calls.tsv, h8_pilot_by_gene.tsv, h8_pilot_modules.tsv, alignments/, hyphy/ (gitignored)}.
"""
import importlib.util
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FA = ROOT / "data/pilot_genomes/fasta"
OUT = ROOT / "results/phase2_pilot"
TMP = ROOT / "data/external/pilot_tmp"
BIN = ROOT / ".bioenv/bin"
SPECIES = ["Melopsittacus_undulatus", "Lathamus_discolor", "Psittacula_echo", "Agapornis_nigrigenis",
           "Strigops_habroptila", "Ara_ararauna", "Amazona_ochrocephala", "Myiopsitta_monachus"]
T = os.environ.get("PILOT_THREADS", "4")


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


parse = load("parse", "scripts/28_parse_miniprot.py")
align = load("align", "scripts/30_align_genes.py")
sel = load("sel", "scripts/31_selection_tests.py")
summ = load("summ", "scripts/32_h8_summary.py")


def panel_query():
    q = TMP / "panel_query.faa"
    if not q.exists():
        panel = {n.split("|")[0] for n, _ in parse.fasta(ROOT / "data/pilot_genes/canonical_proteins.faa")}
        recs = [(n, s) for n, s in parse.fasta(ROOT / "data/genome_genes/canonical_proteins.faa")
                if n.split("|")[0] in panel or n.startswith("kakapo_")]
        q.write_text("".join(f">{n}\n{s}\n" for n, s in recs))
    return q


def extract(sp, query):
    acc = f"pilot_{sp}"
    if (OUT / "hits" / f"{acc}.tsv").exists():
        return
    d = TMP / sp
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    run = lambda *a, **k: subprocess.run([str(x) for x in a], check=True, **k)  # noqa: E731
    with open(d / "windows.fna", "w") as f:
        run(BIN / "seqkit", "sliding", "-g", "-W", "5000000", "-s", "4500000", "-j", T, FA / f"{sp}.fna.gz", stdout=f)
    total = int(subprocess.run([BIN / "seqkit", "stats", "-T", d / "windows.fna"], capture_output=True, text=True,
                               check=True).stdout.splitlines()[1].split("\t")[4])
    run(BIN / "seqkit", "split2", "-p", total // 150_000_000 + 1, "-O", d / "parts", d / "windows.fna",
        stderr=subprocess.DEVNULL)
    with open(d / "miniprot.gff", "w") as f:
        for i, part in enumerate(sorted((d / "parts").iterdir()), 1):
            run(BIN / "miniprot", "-t", T, "-j", "2", "--gff", "--outn=3", "-P", f"P{i}_", part, query, stdout=f,
                stderr=subprocess.DEVNULL)
    run(ROOT / ".venv/bin/python", ROOT / "scripts/28_parse_miniprot.py", acc, sp.replace("_", " "),
        d / "windows.fna", d / "miniprot.gff")
    for src, dst in ((ROOT / "results/phase2_genes/hits" / f"{acc}.tsv", OUT / "hits"),
                     (ROOT / "results/phase2_genes/cds" / f"{acc}.fna", OUT / "cds")):
        dst.mkdir(parents=True, exist_ok=True)
        shutil.move(src, dst / src.name)
    for junk in (ROOT / "results/phase3_genome/calls" / f"{acc}.tsv.gz", ROOT / "data/genome_cds" / f"{acc}.fna.gz"):
        junk.unlink(missing_ok=True)
    shutil.rmtree(d)
    print("extracted", sp, flush=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    TMP.mkdir(parents=True, exist_ok=True)
    q = panel_query()
    for sp in SPECIES:
        extract(sp, q)
    genes = pd.read_csv(ROOT / "data/pilot_genes/canonical_genes.tsv", sep="\t")
    calls = pd.DataFrame({sp: pd.read_csv(OUT / "hits" / f"pilot_{sp}.tsv", sep="\t").set_index("gene").call
                          for sp in SPECIES})
    calls.to_csv(OUT / "calls.tsv", sep="\t")
    cds = {}
    for sp in SPECIES:
        for name, seq in parse.fasta(OUT / "cds" / f"pilot_{sp}.fna"):
            cds.setdefault(name.split("|")[1], {})[sp] = seq
    (OUT / "alignments").mkdir(exist_ok=True)
    for g in genes.gene_symbol:
        out = OUT / "alignments" / f"{g}.codon.fasta"
        if len(cds.get(g, {})) >= 4 and not out.exists():
            try:
                align.align(cds[g], out)
            except subprocess.CalledProcessError as e:
                print("align failed", g, (e.stderr or "")[:100])
    tp, mp = sel.TREES["burgio2019"]
    full = sel.Tree.from_newick(tp.read_text())
    m = pd.read_csv(mp, sep="\t")
    tip_of = {s.replace(" ", "_"): t for s, t in zip(m.species, m.tree_tip)}
    congener = next(t for t in full.tip_labels() if t.startswith("Psittacula_"))
    tip_of.setdefault("Psittacula_echo", congener)
    hy = OUT / "hyphy"
    hy.mkdir(exist_ok=True)
    log = []
    with ThreadPoolExecutor(int(os.environ.get("HYPHY_JOBS", "2"))) as pool:
        for aln in sorted((OUT / "alignments").glob("*.codon.fasta")):
            g = aln.name.split(".")[0]
            species = [line[1:].strip() for line in open(aln) if line.startswith(">")]
            nwk, status = sel.tree_for(full, tip_of, species)
            log.append({"gene": g, "status": status, "n_species": len(species)})
            if nwk is None:
                continue
            for test, text in nwk.items():
                (hy / f"{g}.{test}.nwk").write_text(text)
            for test, args in (("relax", ["--test", "TEST", "--reference", "REFERENCE", "--models", "Minimal"]),
                               ("absrel", ["--branches", "FG"])):
                pool.submit(sel.run, test, args, aln, hy, g, log[-1])
    rows = []
    for g in genes.gene_symbol:
        js = hy / f"{g}.relax.json"
        r = {"gene": g}
        if js.exists():
            r.update(summ.relax_row(js))
            if (hy / f"{g}.absrel.json").exists():
                r.update(summ.absrel_row(hy / f"{g}.absrel.json"))
        rows.append(r)
    t = pd.DataFrame(rows).merge(genes[["gene_symbol", "module"]].rename(columns={"gene_symbol": "gene"}))
    t = t.merge(pd.DataFrame(log), on="gene", how="left")
    ok = t.relax_p.notna()
    t.loc[ok, "relax_q"] = summ.bh(t.relax_p[ok])
    t["relaxed"] = (t.relax_q < summ.Q) & (t.k < 1)
    t["intensified"] = (t.relax_q < summ.Q) & (t.k > 1)
    t["log2_k"] = np.log2(t.k.clip(lower=1e-3))
    t.to_csv(OUT / "h8_pilot_by_gene.tsv", sep="\t", index=False, float_format="%.4g")
    mods = t[ok].groupby("module").agg(n_genes=("gene", "size"), prop_relaxed=("relaxed", "mean"),
                                       prop_intensified=("intensified", "mean"), mean_log2_k=("log2_k", "mean"))
    mods.to_csv(OUT / "h8_pilot_modules.tsv", sep="\t", float_format="%.4g")
    print(f"genes tested: {ok.sum()}; relaxed {int(t.relaxed.sum())}, intensified {int(t.intensified.sum())}")
    print(mods.to_string())


if __name__ == "__main__":
    main()
