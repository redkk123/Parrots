"""Check the budgerigar findings of the H8 pilot on a second, independent individual (exploratory).

The pilot used the VGP reference (PacBio CLR, GCF_012275295.1). The Iridian assembly GCA_025448075.1 is another bird
sequenced with Illumina, so assembly errors are unlikely to be shared:
  - TP53 and MRE11 (1-base indel in the VGP genome, NCBI 'LOW QUALITY PROTEIN'): intact in the Iridian bird means
    the VGP indel is probably an assembly error; disrupted in both points to a real frameshift (still needs reads).
  - Telomere genes with budgerigar-specific amino acids (TERT, RTEL1, TERF2, TERF2IP, POT1): at each site where the
    VGP budgerigar differs from all 16 other pilot species, does the Iridian bird carry the same residue?
Writes results/phase2_pilot_plus/budgie_validation.tsv and prints a summary. Retries the NCBI download.
"""
import importlib.util
import os
import subprocess
import time
from collections import Counter
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
os.environ["PILOT_SET"] = "plus"
spec = importlib.util.spec_from_file_location("p38", ROOT / "scripts/38_pilot_h8.py")
p38 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p38)
OUT = ROOT / "results/phase2_pilot_plus"
SP, ACC = "Melopsittacus_iridian", "GCA_025448075.1"
GENES = ["TP53", "MRE11", "TERT", "RTEL1", "TERF2", "TERF2IP", "POT1"]


def aa(seq):
    return "".join("-" if seq[i:i + 3] == "---" else p38.parse.CODE.get(seq[i:i + 3].upper(), "X")
                   for i in range(0, len(seq), 3))


def main():
    p38.EXTRA[SP] = ACC
    p38.TMP.mkdir(parents=True, exist_ok=True)
    for attempt in range(48):  # NCBI download outages: retry for up to ~8 h
        try:
            p38.extract(SP, p38.panel_query())
            break
        except subprocess.CalledProcessError:
            time.sleep(600)
    calls = pd.read_csv(OUT / "hits" / f"pilot_{SP}.tsv", sep="\t").set_index("gene").call
    vgp = pd.read_csv(OUT / "hits" / "pilot_Melopsittacus_undulatus.tsv", sep="\t").set_index("gene").call
    irid = {n.split("|")[1]: s for n, s in p38.parse.fasta(OUT / "cds" / f"pilot_{SP}.fna")}
    rows = []
    for g in GENES:
        r = {"gene": g, "call_vgp": vgp.get(g), "call_iridian": calls.get(g)}
        aln = OUT / "alignments" / f"{g}.codon.fasta"
        if aln.exists() and g in irid:
            seqs = {k: aa(v) for k, v in p38.parse.fasta(aln)}
            b = seqs.get("Melopsittacus_undulatus")
            # place the Iridian protein on the VGP protein by ungapped coordinates (same species: no indels expected)
            ib = p38.parse.translate(irid[g])
            vgp_ungapped = b.replace("-", "") if b else ""
            same = sum(x == y for x, y in zip(vgp_ungapped, ib)) / max(min(len(vgp_ungapped), len(ib)), 1)
            unique, shared = 0, 0
            if b and abs(len(vgp_ungapped) - len(ib)) <= 3 and same > 0.95:
                pos = -1
                for i, c in enumerate(b):
                    if c == "-":
                        continue
                    pos += 1
                    col = [seqs[k][i] for k in seqs if k != "Melopsittacus_undulatus" and seqs[k][i] not in "-X"]
                    if len(col) >= 10 and len(set(col)) == 1 and c != col[0]:
                        unique += 1
                        shared += pos < len(ib) and ib[pos] == c
            r.update({"identity_vgp_iridian": round(same, 4), "budgie_unique_sites": unique,
                      "shared_by_iridian": shared})
        rows.append(r)
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "budgie_validation.tsv", sep="\t", index=False)
    print(t.to_string(index=False))


if __name__ == "__main__":
    main()
