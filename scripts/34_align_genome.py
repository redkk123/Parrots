"""Genome-wide alignments for RERconverge (analysis plan, 2026-10-28: universe for H6 and the H8 null).

Inputs: results/phase1_audit/species_genomes.tsv (final list, status included*) and data/genome_cds/<acc>.fna.gz from
scripts/27_extract_genes.sh (complete ORFs only). A gene is eligible with a complete ORF in >= 50 % of the species.
Eligible genes: MAFFT --auto on proteins (L-INS-i is kept for the 101-gene panel only, for cost), PAL2NAL codons.
Writes data/genome_alignments/{aa/<gene>.fa, codon/<gene>.fa.gz} (local, large) and
results/phase3_genome/genome_alignment_summary.tsv. Resumable. Workers: ALIGN_WORKERS (default 4).
"""
import gzip
import importlib.util
import os
import subprocess
import tempfile
from multiprocessing import Pool
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / ".bioenv/bin"
ALN = ROOT / "data/genome_alignments"
OUT = ROOT / "results/phase3_genome"
MIN_FRACTION = 0.5
spec = importlib.util.spec_from_file_location("parse", ROOT / "scripts/28_parse_miniprot.py")
parse = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parse)


def align(job):
    gene, seqs = job
    aa, codon = ALN / "aa" / f"{gene}.fa", ALN / "codon" / f"{gene}.fa.gz"
    if aa.exists() and codon.exists():
        return gene, None
    try:
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "cds.fna").write_text("".join(f">{s}\n{x}\n" for s, x in seqs.items()))
            (d / "prot.faa").write_text("".join(f">{s}\n{parse.translate(x)}\n" for s, x in seqs.items()))
            prot = subprocess.run([BIN / "mafft", "--auto", "--thread", "1", "--quiet", d / "prot.faa"],
                                  capture_output=True, text=True, check=True).stdout
            (d / "prot.aln").write_text(prot)
            cod = subprocess.run([BIN / "pal2nal.pl", d / "prot.aln", d / "cds.fna", "-output", "fasta",
                                  "-codontable", "1"], capture_output=True, text=True, check=True).stdout
        with gzip.open(codon, "wt") as f:
            f.write(cod)
        aa.write_text(prot)
        return gene, None
    except subprocess.CalledProcessError as e:
        return gene, (e.stderr or "")[:200].replace("\n", " ")


def main():
    sg = pd.read_csv(ROOT / "results/phase1_audit/species_genomes.tsv", sep="\t")
    sel = sg[sg.status.str.startswith("included")]
    cds, n_sp = {}, 0
    for r in sel.itertuples():
        f = ROOT / "data/genome_cds" / f"{r.accession}.fna.gz"
        if not f.exists():
            print("missing extraction:", r.species, r.accession)
            continue
        n_sp += 1
        with gzip.open(f, "rt") as fh:
            name, buf = None, []
            for line in fh:
                if line.startswith(">"):
                    if name:
                        cds.setdefault(name.split("|")[1], {})[name.split("|")[0]] = "".join(buf)
                    name, buf = line[1:].strip(), []
                else:
                    buf.append(line.strip())
            if name:
                cds.setdefault(name.split("|")[1], {})[name.split("|")[0]] = "".join(buf)
    need = MIN_FRACTION * n_sp
    jobs = [(g, s) for g, s in cds.items() if len(s) >= need]
    for d in ("aa", "codon"):
        (ALN / d).mkdir(parents=True, exist_ok=True)
    with Pool(int(os.environ.get("ALIGN_WORKERS", 4))) as pool:
        errors = dict(pool.imap_unordered(align, jobs, chunksize=8))
    OUT.mkdir(parents=True, exist_ok=True)
    rows = [{"gene": g, "n_species": len(s), "eligible": len(s) >= need, "error": errors.get(g) or ""}
            for g, s in cds.items()]
    pd.DataFrame(rows).sort_values("gene").to_csv(OUT / "genome_alignment_summary.tsv", sep="\t", index=False)
    print(f"{n_sp} genomes; {len(cds)} genes with >= 1 complete ORF; {len(jobs)} eligible (>= {need:.0f} species); "
          f"{sum(1 for e in errors.values() if e)} alignment errors")


if __name__ == "__main__":
    main()
