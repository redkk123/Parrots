"""Nuclear specimen identity from Mash distances (sketches from scripts/16_mash_batch.sh).

For each Iridian assembly of the analysis set:
  - nearest other assembly (by Mash distance) and whether it belongs to the same genus as the claimed species;
  - Mash distance vs patristic distance on the Burgio 2019 supertree (robust linear fit, all pairs of different
    analysis species); a specimen whose distances depart from the tree in a consistent direction is suspect;
  - for assemblies with a foreign mitogenome: Mash distance to the donor species' assembly vs to the nearest
    assembly of the claimed genus. 'nuclear_matches_claim' if the claimed genus is closer than the donor.
    Not testable when every donor is a congener, or when the assembly is fragmented (contig N50 < 5 kb).
Writes results/phase1_audit/mash_identity.tsv and results/phase1_audit/mash_distances.tsv.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from lib.phylo import Tree  # noqa: E402

A = ROOT / "results/phase1_audit"
FRAGMENTED_N50 = 5000  # contig N50 below which Mash distances are inflated (missing k-mers) and not comparable
MASH = ROOT / "data/external/mash"


def main():
    flat = pd.read_csv(ROOT / "data/external/ncbi/psittaciformes_assemblies.tsv", sep="\t").set_index("accession")
    # one sketch per accession (data/external/mash/<accession>.msh); distances in the order of the file list
    d = []
    sk = sorted(MASH.glob("*.msh"))
    acc = [p.stem for p in sk]
    for i, a in enumerate(sk[:-1]):
        o = subprocess.run([str(ROOT / ".bioenv/bin/mash"), "dist", "-p", "8", str(a), *map(str, sk[i + 1:])],
                           capture_output=True, text=True, check=True).stdout.splitlines()
        for j, line in enumerate(o):
            d.append((acc[i], acc[i + 1 + j], float(line.split("\t")[2])))
    D = pd.DataFrame(d, columns=["a", "b", "mash_dist"])
    D["sp_a"] = D.a.map(flat.species)
    D["sp_b"] = D.b.map(flat.species)
    tree = Tree.from_newick((ROOT / "results/phase0/supertree_burgio2019_full.nwk").read_text())
    tips = set(tree.tip_labels())
    sp_tip = {s: s.replace(" ", "_") for s in set(D.sp_a) | set(D.sp_b) if isinstance(s, str) and s.replace(" ", "_") in tips}
    t = tree.prune(set(sp_tip.values()))
    order = t.tip_labels()
    C = t.vcv(order)
    idx = {k: i for i, k in enumerate(order)}
    def patr(x, y):
        if x not in sp_tip or y not in sp_tip:
            return np.nan
        i, j = idx[sp_tip[x]], idx[sp_tip[y]]
        return C[i, i] + C[j, j] - 2 * C[i, j]
    D["patristic"] = [patr(x, y) for x, y in zip(D.sp_a, D.sp_b)]
    D.to_csv(A / "mash_distances.tsv", sep="\t", index=False, float_format="%.5g")

    # expected Mash distance given patristic distance: monotone fit (median in patristic bins), different species
    ok = D[(D.sp_a != D.sp_b) & D.patristic.notna() & (D.mash_dist < 0.25)]
    bins = np.quantile(ok.patristic, np.linspace(0, 1, 21))
    ok = ok.assign(bin=np.clip(np.digitize(ok.patristic, bins) - 1, 0, 19))
    med = ok.groupby("bin").mash_dist.median()
    def expected(p):
        return np.nan if np.isnan(p) else med.get(int(np.clip(np.digitize(p, bins) - 1, 0, 19)), np.nan)

    inv = pd.read_csv(A / "assembly_inventory.tsv", sep="\t")
    sym = pd.concat([D, D.rename(columns={"a": "b", "b": "a", "sp_a": "sp_b", "sp_b": "sp_a"})])
    rows = []
    for r in inv[inv.iridian_n > 0].itertuples():
        a, sp, genus = r.iridian_accession, r.species, r.species.split()[0]
        x = sym[sym.a == a]
        if x.empty:
            continue
        x = x.assign(genus_b=x.sp_b.str.split().str[0])
        nn = x.sort_values("mash_dist").iloc[0]
        other = x[x.sp_b != sp]
        resid = (other.mash_dist - other.patristic.map(expected)).dropna()
        same_sp = x[x.sp_b == sp]
        cong = x[(x.genus_b == genus) & (x.sp_b != sp)]
        m = json.loads((A / "mito" / f"{a}.json").read_text()) if (A / "mito" / f"{a}.json").exists() else {}
        donors = [s.replace("_", " ") for s in m.get("foreign_species", [])]
        # congeneric donors cannot be told apart from the claimed genus by this test
        far_donors = [d for d in donors if d.split()[0] != genus]
        dd = x[x.sp_b.isin(far_donors)]
        cong = cong[~cong.sp_b.isin(donors)]
        row = {"species": sp, "accession": a, "nearest": nn.b, "nearest_species": nn.sp_b,
               "nearest_dist": nn.mash_dist, "nearest_same_genus": nn.genus_b == genus,
               "dist_independent_same_species": same_sp.mash_dist.min() if len(same_sp) else np.nan,
               "dist_nearest_congener": cong.mash_dist.min() if len(cong) else np.nan,
               "median_resid_vs_tree": resid.median() if len(resid) else np.nan,
               "mito_overall": m.get("overall"), "mito_donors": ";".join(donors),
               "dist_nearest_donor": dd.mash_dist.min() if len(dd) else np.nan}
        refs = [v for v in (row["dist_independent_same_species"], row["dist_nearest_congener"]) if not np.isnan(v)]
        ref = min(refs) if refs else np.nan
        if donors and not far_donors:
            row["nuclear_vs_donor"] = "untestable_congeneric_donor"
        elif r.iridian_contig_n50 < FRAGMENTED_N50 and len(dd):
            row["nuclear_vs_donor"] = "untestable_fragmented_assembly"  # Mash distance inflated by fragmentation
        elif len(dd) and refs:
            row["nuclear_vs_donor"] = "nuclear_matches_claim" if ref < row["dist_nearest_donor"] else "nuclear_closer_to_donor"
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(A / "mash_identity.tsv", sep="\t", index=False, float_format="%.4g")
    print(f"sketches: {len(acc)}; Iridian analysis assemblies with distances: {len(out)}")
    print("nearest assembly in claimed genus:", out.nearest_same_genus.value_counts().to_dict())
    if "nuclear_vs_donor" in out:
        print(out.nuclear_vs_donor.value_counts().to_dict())
    print(out.sort_values("median_resid_vs_tree", ascending=False).head(10)[
        ["species", "nearest_species", "nearest_dist", "dist_independent_same_species", "median_resid_vs_tree",
         "mito_overall"]].to_string(index=False))


if __name__ == "__main__":
    main()
