"""Selection tests on the H8 panel codon alignments (analysis plan §10, Nível 1, item 2), HyPhy 2.5.

For each tree (supertree = main, Smith 2024 = sensitivity), each alignment from scripts/30_align_genes.py is run with
  RELAX   TEST = stem of Psittaculidae + its internal branches; REFERENCE = every branch outside Psittaculidae;
          Psittaculidae terminal branches unlabelled (plan wording: "tronco de Psittaculidae + ramos internos");
  aBSREL  the Psittaculidae stem and the Melopsittacus terminal branch (H8c, item 4: stem vs Melopsittacus only;
          Holm correction across the two tested branches, as aBSREL does).
Zero-length internal branches are kept (--kill-zero-lengths No), so a short stem is never dropped from the test.
Psittaculidae = genera in scripts/26_select_genomes.py. A gene is skipped on a tree when it has < 2 Psittaculidae or
< 2 other species, or when the Psittaculidae present are not monophyletic on that tree.
JSON: results/phase2_genes/hyphy/<tree>/<gene>.<test>.json (gitignored). Resumable. Threads: HYPHY_CPU (default 4).
Usage: 31_selection_tests.py [gene ...]     (default: every alignment). HYPHY_JOBS runs in parallel (default 2).
"""
import importlib.util
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from lib.phylo import Tree  # noqa: E402

G = ROOT / "results/phase2_genes"
HYPHY = ROOT / ".bioenv/bin/hyphy"
CPU = os.environ.get("HYPHY_CPU", "4")
TREES = {"burgio2019": (ROOT / "results/phase0/supertree_burgio2019_full.nwk",
                        ROOT / "results/phase0b/species_data_matrix_with_residuals.tsv"),
         "smith2024": (ROOT / "results/phase1_trees/smith2024_species.nwk",
                       ROOT / "results/phase1_trees/species_data_matrix_smith2024.tsv")}
spec = importlib.util.spec_from_file_location("sel", ROOT / "scripts/26_select_genomes.py")
sel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sel)


def labelled_newick(t, names, tags):
    """Newick with one HyPhy tag per branch ({node: tag}), topology only (time-scaled lengths as starting values
    broke the supertree fit in a test: HyPhy 'ComputeBranchCache' internal error). Internal nodes are named N<i>,
    since HyPhy ignores tags on unnamed internal nodes; tips are renamed to the alignment names."""
    def rec(n):
        s = names[t.label[n]] if t.is_tip(n) else "(" + ",".join(rec(c) for c in t.children[n]) + f")N{n}"
        return s + (f"{{{tags[n]}}}" if n in tags else "")
    return rec(t.root) + ";"


def tree_for(full, tip_of, species):
    tips = {tip_of[s]: s for s in species if s in tip_of}
    t = full.prune(set(tips))
    psl = {n for n in t.tips() if sel.group(tips[t.label[n]].split("_")[0]) != "other_Psittaciformes"}
    if len(psl) < 2 or len(t.tips()) - len(psl) < 2:
        return None, "too_few_species"
    desc = t.descendants_tips()
    crown = min((n for n in range(len(t.parent)) if psl <= set(desc[n])), key=lambda n: len(desc[n]))
    if set(desc[crown]) != psl:
        return None, "psittaculidae_not_monophyletic"
    inside = {n for n in range(len(t.parent)) if set(desc[n]) <= psl}
    test = {n for n in inside if not t.is_tip(n)}          # stem (branch above the crown node) + internal branches
    reference = set(range(len(t.parent))) - inside - {t.root}
    names = {t.label[n]: tips[t.label[n]] for n in t.tips()}
    relax = {**{n: "TEST" for n in test}, **{n: "REFERENCE" for n in reference}}
    fg = {crown: "FG", **{n: "FG" for n in t.tips() if names[t.label[n]] == "Melopsittacus_undulatus"}}
    return {"relax": labelled_newick(t, names, relax), "absrel": labelled_newick(t, names, fg)}, "ok"


def run(test, args, aln, out, g, log_row):
    js = out / f"{g}.{test}.json"
    if js.exists():
        return
    r = subprocess.run([HYPHY, f"CPU={CPU}", test, "--alignment", aln, "--tree", out / f"{g}.{test}.nwk",
                        *args, "--kill-zero-lengths", "No",
                        "--output", js.with_suffix(".tmp")], capture_output=True, text=True,
                       cwd=out)  # HyPhy writes errors.log / messages.log in its working directory
    if r.returncode == 0 and js.with_suffix(".tmp").exists():
        js.with_suffix(".tmp").replace(js)
    else:
        log_row[f"{test}_error"] = (r.stderr or r.stdout)[-300:].replace("\n", " ")
    print(out.name, g, test, "ok" if js.exists() else "FAILED", flush=True)


def main():
    genes = sys.argv[1:] or sorted(p.name.split(".")[0] for p in (G / "alignments").glob("*.codon.fasta"))
    log = []
    pool = ThreadPoolExecutor(int(os.environ.get("HYPHY_JOBS", "2")))  # HyPhy runs in parallel, CPU threads each
    for name, (tp, mp) in TREES.items():
        full = Tree.from_newick(tp.read_text())
        m = pd.read_csv(mp, sep="\t")
        tip_of = {s.replace(" ", "_"): t for s, t in zip(m.species, m.tree_tip)}
        out = G / "hyphy" / name
        out.mkdir(parents=True, exist_ok=True)
        for g in genes:
            aln = G / "alignments" / f"{g}.codon.fasta"
            species = [line[1:].strip() for line in open(aln) if line.startswith(">")]
            nwk, status = tree_for(full, tip_of, species)
            log.append({"tree": name, "gene": g, "status": status, "n_species": len(species)})
            if nwk is None:
                continue
            for test, text in nwk.items():
                (out / f"{g}.{test}.nwk").write_text(text)
            for test, args in (("relax", ["--test", "TEST", "--reference", "REFERENCE", "--models", "Minimal"]),
                               ("absrel", ["--branches", "FG"])):
                pool.submit(run, test, args, aln, out, g, log[-1])  # queued in order: main tree first
    pool.shutdown(wait=True)
    pd.DataFrame(log).to_csv(G / "selection_tests_log.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
