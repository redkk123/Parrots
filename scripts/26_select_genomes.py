"""One genome per species for the gene-level analyses (analysis plan §4, revised 2026-10-16 and 2026-10-28).

Candidates per species: the Iridian assembly (audited in results/phase1_audit/audit_table.tsv) and the best
non-Iridian assembly in results/phase1_audit/assembly_inventory.tsv (audited with the same scripts/10_audit_assembly.sh).
A candidate is usable when BUSCO (aves_odb10, compleasm) >= 70 % and its identity does not fail. Among usable ones,
a chromosome-scale non-Iridian assembly is preferred, then the highest BUSCO.
Identity of a non-Iridian assembly: own mitogenome confirmed, or (no mitogenome in the assembly, common for VGP and
Sanger, which deposit it separately) nuclear Mash distance <= 0.01 to an Iridian assembly of the same species whose
identity passed; otherwise accepted on provenance only. Clade-level ('pass_nuclear_placement') and provenance-only
identities stay in the main set and are dropped in the identity sensitivity.

Writes results/phase1_audit/species_genomes.tsv (one row per species) and results/phase1_audit/audit_queue_other.tsv
(non-Iridian assemblies still to audit, Psittaculidae first). Runs on partial results.
"""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
A = ROOT / "results/phase1_audit"
spec = importlib.util.spec_from_file_location("audit_summary", ROOT / "scripts/14_audit_summary.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)

LORIINI = {"Trichoglossus", "Lorius", "Eos", "Vini", "Chalcopsitta", "Pseudeos", "Charmosyna", "Glossopsitta",
           "Psitteuteles", "Neopsittacus", "Oreopsittacus", "Phigys", "Parvipsitta"}
PSITTACULIDAE = LORIINI | {
    "Melopsittacus", "Cyclopsitta", "Psittaculirostris", "Psittacella",                       # Loriinae, Psittacellinae
    "Agapornis", "Loriculus", "Bolbopsittacus",                                               # Agapornithinae
    "Psittacula", "Tanygnathus", "Prioniturus", "Eclectus", "Geoffroyus", "Alisterus", "Aprosmictus", "Polytelis",
    "Psittinus", "Micropsitta",                                                               # Psittaculinae
    "Platycercus", "Barnardius", "Purpureicephalus", "Lathamus", "Prosopeia", "Eunymphicus", "Cyanoramphus",
    "Pezoporus", "Neophema", "Neopsephotus", "Northiella", "Psephotus", "Psephotellus"}       # Platycercinae
SAME_SPECIES_MASH = 0.01


def group(genus):
    if genus == "Melopsittacus":
        return "Melopsittacus"
    if genus in LORIINI:
        return "Loriini"
    return "Psittaculidae_other" if genus in PSITTACULIDAE else "other_Psittaciformes"


def other_identity(mito_overall, iridian_id, mash_same):
    if mito_overall in ("confirmed", "confirmed_with_foreign_mito"):
        return "pass"
    if mito_overall == "conflict":
        return "check"
    if iridian_id in ("pass", "pass_nuclear_placement") and mash_same <= SAME_SPECIES_MASH:
        return "pass_nuclear_same_species"
    return "pass_provenance"


def main():
    inv = pd.read_csv(A / "assembly_inventory.tsv", sep="\t")
    ir = pd.read_csv(A / "audit_table.tsv", sep="\t").set_index("species")
    mash = pd.read_csv(A / "mash_identity.tsv", sep="\t").set_index("accession")
    rows, queue = [], []
    for r in inv.itertuples():
        cands = []
        if r.species in ir.index:
            i = ir.loc[r.species]
            cands.append({"accession": i.accession, "source": "Iridian", "level": r.iridian_level,
                          "contig_n50": i.contig_n50, "busco_C": i.busco_C, "identity": i.identity_pass})
        if isinstance(r.best_other_accession, str):
            b, m = audit.busco(r.best_other_accession), audit.mito(r.best_other_accession)
            iid = cands[0]["identity"] if cands else None
            ms = mash.dist_independent_same_species.get(cands[0]["accession"], np.nan) if cands else np.nan
            cands.append({"accession": r.best_other_accession, "source": r.best_other_submitter,
                          "level": r.best_other_level, "contig_n50": r.best_other_contig_n50,
                          "busco_C": b.get("busco_C", np.nan),
                          "identity": other_identity(m.get("mito_overall"), iid, ms) if m else "pending"})
        for c in cands:
            c["busco"] = "pending" if pd.isna(c["busco_C"]) else ("pass" if c["busco_C"] >= audit.BUSCO_MIN else "fail")
            c["usable"] = c["busco"] == "pass" and c["identity"] not in ("fail", "check", "pending", "unresolved",
                                                                          "check_nuclear")
            c["preferred"] = c["source"] != "Iridian" and c["level"] == "Chromosome"
        use = sorted((c for c in cands if c["usable"]), key=lambda c: (not c["preferred"], -c["busco_C"]))
        waiting = [c for c in cands if c["busco"] == "pending" or c["identity"] == "pending"]
        best = use[0] if use else None
        if best and not any(c["preferred"] for c in waiting if not best["preferred"]):
            status = "included"
        elif best:
            status = "included_provisional"  # a preferred chromosome-scale assembly is still being audited
        elif waiting:
            status = "pending"
        else:
            status = "no_assembly" if not cands else "excluded"
        d = {"species": r.species, "genus": r.genus, "group": group(r.genus),
             "psittaculidae": r.genus in PSITTACULIDAE, "status": status, "n_candidates": len(cands)}
        if best:
            d.update({k: best[k] for k in ("accession", "source", "level", "contig_n50", "busco_C", "identity")})
            d["original_criterion_pass"] = bool(best["busco_C"] >= 85 and best["contig_n50"] >= 20000)
            d["identity_sensitivity_keep"] = best["identity"] == "pass"
        d["candidates"] = "; ".join(f"{c['accession']} ({c['source']}, BUSCO {c['busco_C']}, id {c['identity']})"
                                    for c in cands)
        rows.append(d)
        for c in cands:  # non-Iridian assemblies that could add the species or replace a weaker Iridian one
            if c["source"] != "Iridian" and c["busco"] == "pending":
                ic = next((x for x in cands if x["source"] == "Iridian"), None)
                gain = ic is None or not ic["usable"] or c["preferred"] or c["contig_n50"] > 2 * ic["contig_n50"]
                if gain:
                    queue.append((r.genus not in PSITTACULIDAE, ic is not None and ic["usable"], c["accession"],
                                  r.species))
    t = pd.DataFrame(rows)
    t.to_csv(A / "species_genomes.tsv", sep="\t", index=False, float_format="%.4g")
    q = sorted(queue)
    (A / "audit_queue_other.tsv").write_text("".join(f"{acc}\t{sp}\n" for *_, acc, sp in q))
    print(pd.crosstab(t.status, t.group, margins=True).to_string())
    inc = t[t.status.str.startswith("included")]
    print(f"\nincluded (incl. provisional): {len(inc)}; Psittaculidae {inc.psittaculidae.sum()}; "
          f"original criterion {inc.original_criterion_pass.sum()}; identity sensitivity {inc.identity_sensitivity_keep.sum()}")
    print(f"non-Iridian queue: {len(q)} (Psittaculidae {sum(not x[0] for x in q)})")


if __name__ == "__main__":
    main()
