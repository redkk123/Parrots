"""Resolve parrot names from external sources to species of the Burgio et al. 2019 supertree."""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


class NameResolver:
    def __init__(self):
        sp = pd.read_csv(ROOT / "results/phase0/species.tsv", sep="\t")
        self.tree = set(sp["species"])
        tax = pd.read_csv(ROOT / "data/external/parrottraits/data-raw/parrot-traits-database.csv", header=1,
                          encoding="utf-8-sig")
        self.db_syn = {}
        for _, r in tax.dropna(subset=["taxon_2 (matching)"]).iterrows():
            target = r["taxon_2 (matching)"].replace("_", " ")
            raw = str(r["synonyms "]) if pd.notna(r["synonyms "]) else ""
            for x in raw.replace(";", ",").split(","):
                x = x.strip()
                if len(x.split()) == 2:
                    self.db_syn.setdefault(x, target)
        cur = pd.read_csv(ROOT / "data/curated/taxonomy_synonyms.tsv", sep="\t")
        self.curated = {r.name_in_source: (r.tree_species_burgio2019, r.map_type) for r in cur.itertuples()}

    def resolve(self, name):
        """Return (tree_species or None, map_type)."""
        name = " ".join(str(name).split())
        if " x " in name or name.endswith(" sp.") or len(name.split()) < 2:
            return None, "hybrid_or_unidentified"
        binom = " ".join(name.split()[:2])
        subsp = "subspecies_collapsed" if len(name.split()) > 2 else "exact"
        if binom in self.curated:
            return self.curated[binom][0], self.curated[binom][1]
        if binom in self.tree:
            return binom, subsp
        if binom in self.db_syn and self.db_syn[binom] in self.tree:
            return self.db_syn[binom], "burgio_synonym"
        return None, "unresolved"
