"""NCBI assembly inventory for Psittaciformes and per-species choice for the analysis set.

1. Downloads the NCBI Datasets genome report for taxon 9223 (all assembly versions) unless already present.
2. Flattens it to data/external/ncbi/psittaciformes_assemblies.tsv (one row per assembly).
3. For the 105 analysis species (gene-level usable + life expectancy, results/phase0b/species_data_matrix.tsv)
   writes results/phase1_audit/assembly_inventory.tsv: best current Iridian assembly (by contig N50), best
   non-Iridian (by scaffold N50) and best of any source (by contig N50).
"""
import json
import sys
import time
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from lib.names import NameResolver  # noqa: E402

NCBI = ROOT / "data/external/ncbi"
REPORT = NCBI / "psittaciformes_dataset_report.json"
FLAT = NCBI / "psittaciformes_assemblies.tsv"
OUT = ROOT / "results/phase1_audit/assembly_inventory.tsv"
API = ("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/9223/dataset_report"
       "?page_size=1000&filters.assembly_version=all_assemblies")
COLS = ["accession", "organism", "isolate", "submitter", "bioproject", "biosample", "release", "method", "total_len",
        "contig_n50", "scaffold_n50", "n_contigs", "gc", "coverage", "level"]


def fetch():
    out, tok = [], None
    while True:
        d = json.load(urllib.request.urlopen(API + (f"&page_token={tok}" if tok else ""), timeout=120))
        out += d.get("reports", [])
        tok = d.get("next_page_token")
        if not tok:
            return out
        time.sleep(0.4)


def flatten(reports):
    rows = []
    for x in reports:
        a, s, o = x["assembly_info"], x.get("assembly_stats", {}), x["organism"]
        b = x.get("annotation_info", {}).get("busco", {})
        rows.append(dict(
            accession=x["accession"], organism=o["organism_name"], tax_id=o["tax_id"],
            isolate=o.get("infraspecific_names", {}).get("isolate", ""), assembly_name=a.get("assembly_name"),
            level=a.get("assembly_level"), status=a.get("assembly_status"), refseq_cat=a.get("refseq_category", ""),
            submitter=a.get("submitter", ""), bioproject=a.get("bioproject_accession"),
            biosample=a.get("biosample", {}).get("accession"), release=a.get("release_date"),
            tech=a.get("sequencing_tech", ""), method=a.get("assembly_method", ""),
            total_len=s.get("total_sequence_length"), contig_n50=s.get("contig_n50"),
            scaffold_n50=s.get("scaffold_n50"), n_contigs=s.get("number_of_contigs"),
            n_scaffolds=s.get("number_of_scaffolds"), gc=s.get("gc_percent"), coverage=s.get("genome_coverage"),
            busco_ncbi=b.get("complete"), busco_lineage=b.get("busco_lineage"), paired=x.get("paired_accession", "")))
    return pd.DataFrame(rows)


def main():
    NCBI.mkdir(parents=True, exist_ok=True)
    if not REPORT.exists():
        json.dump(fetch(), open(REPORT, "w"))
    flat = flatten(json.load(open(REPORT)))
    R = NameResolver()
    flat["species"], flat["name_map"] = zip(*flat.organism.map(R.resolve))
    flat["iridian"] = flat.submitter.str.contains("Iridian", case=False)
    flat.to_csv(FLAT, sep="\t", index=False)

    d = flat[flat.status == "current"]
    m = pd.read_csv(ROOT / "results/phase0b/species_data_matrix.tsv", sep="\t")
    a = m[m.usable_for_gene_level & m.life_expectancy_years.notna()]
    rows = []
    for sp in a.species:
        x = d[d.species == sp]
        r = {"species": sp, "n_ncbi_assemblies": len(x)}
        for pre, t in (("iridian", x[x.iridian].sort_values("contig_n50", ascending=False)),
                       ("best_other", x[~x.iridian].sort_values("scaffold_n50", ascending=False))):
            if len(t):
                r.update({f"{pre}_{c}": t.iloc[0][c] for c in COLS})
            r[f"{pre}_n"] = len(t)
        if len(x):
            b = x.sort_values("contig_n50", ascending=False).iloc[0]
            r.update(best_any_accession=b.accession, best_any_contig_n50=b.contig_n50,
                     best_any_submitter=b.submitter, best_any_level=b.level)
        rows.append(r)
    t = pd.DataFrame(rows).merge(a[["species", "family", "genus", "genome_tier", "in_hains2022_iridian",
                                    "best_short_cov_x", "log_life_expectancy"]], on="species")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    t.to_csv(OUT, sep="\t", index=False)
    print(f"{len(flat)} assemblies ({flat.iridian.sum()} Iridian); analysis set {len(t)}: "
          f"{(t.iridian_n > 0).sum()} with Iridian assembly, best-any contig N50 >= 20 kb: "
          f"{(t.best_any_contig_n50 >= 20000).sum()}")


if __name__ == "__main__":
    main()
