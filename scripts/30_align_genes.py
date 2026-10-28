"""Codon alignments of the H8 panel genes across the selected genomes (analysis plan §10, Nível 1).

Inputs: results/phase1_audit/species_genomes.tsv (one genome per species, status included*) and, for each chosen
accession, results/phase2_genes/{hits,cds}/<acc>.tsv/.fna from scripts/27_extract_genes.sh.
Per gene: complete CDS only (missing otherwise), protein alignment with MAFFT L-INS-i, back-translated to codons
with PAL2NAL; genes with fewer than MIN_SPECIES species are skipped. Sequence names are species with underscores.
Writes results/phase2_genes/alignments/<gene>.codon.fasta, results/phase2_genes/gene_calls.tsv (gene x species call
matrix) and results/phase2_genes/alignment_summary.tsv.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
G = ROOT / "results/phase2_genes"
BIN = ROOT / ".bioenv/bin"
MIN_SPECIES = 4
spec = __import__("importlib.util").util.spec_from_file_location("parse", ROOT / "scripts/28_parse_miniprot.py")
parse = __import__("importlib.util").util.module_from_spec(spec)
spec.loader.exec_module(parse)


def align(seqs, out):
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        (d / "cds.fna").write_text("".join(f">{s}\n{x}\n" for s, x in seqs.items()))
        (d / "prot.faa").write_text("".join(f">{s}\n{parse.translate(x)}\n" for s, x in seqs.items()))
        aln = subprocess.run([BIN / "mafft", "--localpair", "--maxiterate", "1000", "--thread", "4",
                              "--quiet", d / "prot.faa"], capture_output=True, text=True, check=True).stdout
        (d / "prot.aln").write_text(aln)
        codon = subprocess.run([BIN / "pal2nal.pl", d / "prot.aln", d / "cds.fna", "-output", "fasta",
                                "-codontable", "1"], capture_output=True, text=True, check=True).stdout
    out.write_text(codon)


def main():
    sg = pd.read_csv(ROOT / "results/phase1_audit/species_genomes.tsv", sep="\t")
    sel = sg[sg.status.str.startswith("included")]
    genes = pd.read_csv(ROOT / "data/pilot_genes/canonical_genes.tsv", sep="\t")
    calls, cds = {}, {}
    for r in sel.itertuples():
        h = G / "hits" / f"{r.accession}.tsv"
        if not h.exists():
            print("missing extraction:", r.species, r.accession)
            continue
        sp = r.species.replace(" ", "_")
        calls[sp] = pd.read_csv(h, sep="\t").set_index("gene").call
        for name, seq in parse.fasta(G / "cds" / f"{r.accession}.fna"):
            cds.setdefault(name.split("|")[1], {})[sp] = seq
    m = pd.DataFrame(calls).reindex(genes.gene_symbol)
    m.index.name = "gene"
    m.to_csv(G / "gene_calls.tsv", sep="\t")
    (G / "alignments").mkdir(exist_ok=True)
    rows = []
    for g in genes.itertuples():
        seqs = cds.get(g.gene_symbol, {})
        r = {"gene": g.gene_symbol, "module": g.module, "n_species": len(seqs),
             "n_disrupted": int((m.loc[g.gene_symbol] == "disrupted").sum()),
             "n_partial": int((m.loc[g.gene_symbol] == "partial").sum()),
             "n_paralog": int((m.loc[g.gene_symbol] == "paralog").sum()),
             "n_absent": int((m.loc[g.gene_symbol] == "absent").sum())}
        out = G / "alignments" / f"{g.gene_symbol}.codon.fasta"
        if len(seqs) >= MIN_SPECIES and not out.exists():
            try:
                align(seqs, out)
            except subprocess.CalledProcessError as e:
                r["error"] = (e.stderr or "")[:200].replace("\n", " ")
        if out.exists():
            a = dict(parse.fasta(out))
            r["aln_codons"] = len(next(iter(a.values()))) // 3 if a else 0
        rows.append(r)
    s = pd.DataFrame(rows)
    s.to_csv(G / "alignment_summary.tsv", sep="\t", index=False)
    print(f"{len(calls)} genomes; {(s.n_species >= MIN_SPECIES).sum()} of {len(s)} genes aligned; "
          f"median species per gene {s.n_species.median():.0f}")


if __name__ == "__main__":
    main()
