"""Collect the Iridian assembly audit (inventory + BUSCO + mitochondrial identity) into one table.

Pass criteria (analysis plan, revised 2026-10-16):
  BUSCO complete (S+D, aves_odb10) >= 70 % (main); the original BUSCO >= 85 % and contig N50 >= 20 kb is kept as
  a sensitivity set (column original_criterion_pass)
  species identity 'confirmed' by mitochondrial COX1/CYTB/ND2 against independent GenBank sequences;
  a mitochondrial 'conflict' (foreign mitogenome, own not assembled) is not an exclusion by itself: it is
  resolved by the nuclear Mash check (results/phase1_audit/mash_identity.tsv) when available; assemblies listed
  in results/phase1_audit/nuclear_identity_review.tsv (wrong organism/specimen, mixed genomes) fail identity;
  when the mitogenome cannot decide (no own mitogenome, ambiguous, congeneric donor), an assembly with contig N50
  >= 5 kb whose nuclear Mash distances fit the tree (median residual <= RESID_MAX) passes at clade level only
  ('pass_nuclear_placement'; plan, 2026-10-28): kept in the main set, dropped in the identity sensitivity
Contamination: a second (foreign) mitogenome in the assembly is flagged as cross-sample contamination;
GC or total length outside the expected range is flagged as possible contamination/redundancy (no read-level screen yet).
Runs on partial results; rows without BUSCO/mito are 'pending'.
Writes results/phase1_audit/audit_table.tsv and prints a summary.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
A = ROOT / "results/phase1_audit"
BUSCO_MIN = 70.0  # main criterion; original (85 %, N50 >= 20 kb) reported as sensitivity
RESID_MAX = 0.006  # Mash residual vs tree: mito-confirmed (N50 >= 5 kb) reach 0.0031, the four failures start at 0.012


def busco(acc):
    p = A / "busco" / acc / "summary.txt"
    if not p.exists():
        return {}
    d = {}
    for line in p.read_text().splitlines():
        line = line.strip()
        if line[:2] in ("S:", "D:", "F:", "I:", "M:") and "%" in line:
            d[f"busco_{line[0]}"] = float(line.split(":")[1].split("%")[0])
        if line.startswith("N:"):
            d["busco_N"] = int(line[2:].split()[0])
    if "busco_S" in d:
        d["busco_C"] = round(d["busco_S"] + d.get("busco_D", 0.0), 2)
    return d


def mito(acc):
    p = A / "mito" / f"{acc}.json"
    if not p.exists():
        return {}
    r = json.loads(p.read_text())
    return {"mito_overall": r["overall"], "mito_own_genes": ";".join(r.get("own_genes", [])),
            "mito_best_own_pid": r.get("best_own_pid"), "mito_foreign_species": ";".join(r.get("foreign_species", [])),
            "mito_n_foreign_contigs": len(r.get("foreign_contigs", []))}


def main():
    inv = pd.read_csv(A / "assembly_inventory.tsv", sep="\t")
    inv = inv[inv.iridian_n > 0]
    rows = []
    for r in inv.itertuples():
        d = {"species": r.species, "family": r.family, "accession": r.iridian_accession,
             "release": r.iridian_release, "method": r.iridian_method, "total_len_gb": r.iridian_total_len / 1e9,
             "contig_n50": r.iridian_contig_n50, "scaffold_n50": r.iridian_scaffold_n50, "gc": r.iridian_gc,
             "coverage_x": r.iridian_coverage, "alt_assembly": r.best_other_accession}
        d.update(busco(r.iridian_accession))
        d.update(mito(r.iridian_accession))
        rows.append(d)
    t = pd.DataFrame(rows)
    for c in ("busco_C", "mito_overall"):
        if c not in t:
            t[c] = np.nan
    t["flag_fragmented"] = t.contig_n50 < 5000
    t["flag_length"] = ~t.total_len_gb.between(1.0, 1.35)
    t["flag_gc"] = ~t.gc.between(40, 44.5)
    t["flag_reference_scaffolded"] = t.scaffold_n50 > 1e6  # short-read only: scaffold N50 > 1 Mb implies reference guidance
    t["busco_pass"] = np.where(t.busco_C.isna(), "pending", np.where(t.busco_C >= BUSCO_MIN, "pass", "fail"))
    mash = A / "mash_identity.tsv"
    if mash.exists():
        mi = pd.read_csv(mash, sep="\t").set_index("accession")
        t["nuclear_vs_donor"] = t.accession.map(mi.get("nuclear_vs_donor", pd.Series(dtype=str)))
        t["nearest_same_genus"] = t.accession.map(mi.nearest_same_genus)
        t["mash_resid_vs_tree"] = t.accession.map(mi.median_resid_vs_tree)
    else:
        t["nuclear_vs_donor"] = t["mash_resid_vs_tree"] = np.nan
    t["original_criterion_pass"] = np.where(t.busco_C.isna(), "pending",
                                            np.where((t.busco_C >= 85) & (t.contig_n50 >= 20000), "pass", "fail"))
    t["identity_pass"] = np.select(
        [t.mito_overall.isna(), t.mito_overall.isin(["confirmed", "confirmed_with_foreign_mito"]),
         t.mito_overall == "conflict"], ["pending", "pass", "check_nuclear"], "unresolved")
    rev = A / "nuclear_identity_review.tsv"  # manual review of Mash anomalies (wrong organism/specimen, mixtures)
    if rev.exists():
        rv = pd.read_csv(rev, sep="\t").set_index("accession")
        t["nuclear_review"] = t.accession.map(rv.call)
        t.loc[t.nuclear_review.notna(), "identity_pass"] = "fail"
    t.loc[(t.identity_pass == "check_nuclear") & (t.nuclear_vs_donor == "nuclear_matches_claim"), "identity_pass"] = "pass"
    t.loc[(t.identity_pass == "check_nuclear") & (t.nuclear_vs_donor == "nuclear_closer_to_donor"), "identity_pass"] = "fail"
    placed = t.identity_pass.isin(["check_nuclear", "unresolved"]) & (t.contig_n50 >= 5000) & (t.mash_resid_vs_tree <= RESID_MAX)
    t.loc[placed, "identity_pass"] = "pass_nuclear_placement"
    t["flag_foreign_mito"] = t.mito_overall.isin(["confirmed_with_foreign_mito", "conflict"])
    ok_id = t.identity_pass.isin(["pass", "pass_nuclear_placement"])
    t["audit"] = np.select(
        [(t.busco_pass == "fail") | (t.identity_pass == "fail"),
         (t.busco_pass == "pass") & ok_id & ~t.flag_length & ~t.flag_gc & ~t.flag_foreign_mito,
         (t.busco_pass == "pass") & ok_id],
        ["exclude", "keep", "keep_check_contamination"], "pending")
    t.sort_values(["audit", "species"]).to_csv(A / "audit_table.tsv", sep="\t", index=False, float_format="%.4g")
    print(f"assemblies: {len(t)}; BUSCO done: {t.busco_C.notna().sum()}; identity done: {t.mito_overall.notna().sum()}")
    print(t.audit.value_counts().to_string())
    print(t[["flag_fragmented", "flag_length", "flag_gc", "flag_reference_scaffolded", "flag_foreign_mito"]].sum().to_string())
    print(t.mito_overall.value_counts(dropna=False).to_string())
    done = t[t.busco_C.notna()]
    if len(done):
        print(done[["species", "contig_n50", "busco_C", "busco_F", "busco_M", "mito_overall", "audit"]]
              .sort_values("busco_C").to_string(index=False))
    bad = t[t.mito_overall.notna() & (t.mito_overall != "confirmed")]
    if len(bad):
        print(bad[["species", "accession", "mito_overall", "mito_own_genes", "mito_foreign_species"]].to_string(index=False))


if __name__ == "__main__":
    main()
