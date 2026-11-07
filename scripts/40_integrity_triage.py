"""Read-free triage of 'disrupted' gene calls (item 3 of docs/HANDOFF_GENES.md, first pass).

A disrupted call (frameshift / internal stop) is split by what congeneric species show for the same gene:
  isolated   : >= 1 congener has the gene complete and none is disrupted -> probably an assembly indel error
  shared     : >= 1 congener also disrupted                              -> candidate real loss (needs read check)
  unresolved : no congener in the set (or congeners only partial/absent/paralog)
This does NOT call losses; it only ranks what is worth checking against reads. Rules fixed before looking at results.
Writes results/phase2_genes/disrupted_triage.tsv and disrupted_triage_by_gene.tsv.
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
G = ROOT / "results/phase2_genes"

calls = pd.read_csv(G / "gene_calls.tsv", sep="\t", index_col=0)
sp = pd.read_csv(ROOT / "results/phase1_audit/analysis_species.tsv", sep="\t")
sp["col"] = sp["species"].str.replace(" ", "_")
sp = sp[sp["col"].isin(calls.columns)].set_index("col")

rows = []
for gene, r in calls.iterrows():
    for s in r.index[r == "disrupted"]:
        genus = sp.loc[s, "genus"]
        cong = [c for c in calls.columns if c != s and c in sp.index and sp.loc[c, "genus"] == genus]
        cc = r[cong]
        if (cc == "disrupted").any():
            cls = "shared"
        elif (cc == "complete").any():
            cls = "isolated"
        else:
            cls = "unresolved"
        rows.append({"gene": gene, "species": s, "genus": genus, "accession": sp.loc[s, "accession"],
                     "source": sp.loc[s, "source"], "level": sp.loc[s, "level"], "n_congeners": len(cong),
                     "congeners_complete": int((cc == "complete").sum()),
                     "congeners_disrupted": int((cc == "disrupted").sum()), "class": cls})
t = pd.DataFrame(rows)
t.to_csv(G / "disrupted_triage.tsv", sep="\t", index=False)
by = t.groupby("gene")["class"].value_counts().unstack(fill_value=0)
by.to_csv(G / "disrupted_triage_by_gene.tsv", sep="\t")
print(t["class"].value_counts().to_string())
print(t.groupby("source")["class"].value_counts().unstack(fill_value=0).to_string())
print("genes with shared disruptions:", sorted(t.loc[t["class"] == "shared", "gene"].unique()))
