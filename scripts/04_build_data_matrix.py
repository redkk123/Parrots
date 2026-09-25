"""Phase 0b: per-species data matrix (phenotypes x genomic data) for the whole supertree.

Adds to Phase 0:
  - maximum longevity from AnAge and the Amniote Life History Database (via traitdata)
  - public sequencing data per species from the SRA metadata snapshot (script 03)
Outputs results/phase0b/species_data_matrix.tsv and max_longevity_sources.tsv.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rdata

sys.path.insert(0, str(Path(__file__).parent / "lib"))
from names import NameResolver  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "data/external"
P0 = ROOT / "results/phase0"
OUT = ROOT / "results/phase0b"
OUT.mkdir(parents=True, exist_ok=True)
GENOME_GB = 1.2
TIERS = ["A_chromosome_scale", "B_long_read_ge15x", "C_short_read_ge20x", "D_short_read_5_20x",
         "E_low_pass_lt5x", "F_literature_only", "none"]


def load_rda(name):
    return rdata.read_rda(str(EXT / f"traitdata/data/{name}.rda"))[name]


def max_longevity(res):
    an = load_rda("an_age")
    an = an[an["Order"] == "Psittaciformes"].copy()
    an["name"] = an["Genus"].str.strip() + " " + an["Species"].str.strip()
    an = an.rename(columns={"Maximum.longevity..yrs.": "max_longevity_years", "Specimen.origin": "origin",
                            "Sample.size": "sample_size", "Data.quality": "data_quality"})
    an["source"] = "AnAge"
    am = load_rda("amniota")
    am = am[am["Order"] == "Psittaciformes"].copy()
    am["name"] = am["Genus"].str.strip() + " " + am["Species"].str.strip()
    am = am.rename(columns={"maximum_longevity_y": "max_longevity_years"})
    am["max_longevity_years"] = am["max_longevity_years"].where(am["max_longevity_years"] > 0)
    for c in ["origin", "sample_size", "data_quality"]:
        am[c] = ""
    am["source"] = "Amniote LHD (Myhrvold et al. 2015)"
    cols = ["name", "max_longevity_years", "origin", "sample_size", "data_quality", "source"]
    ml = pd.concat([an[cols], am[cols]], ignore_index=True)
    ml[["tree_species", "map_type"]] = ml["name"].apply(lambda n: pd.Series(res.resolve(n)))
    ml.to_csv(OUT / "max_longevity_sources.tsv", sep="\t", index=False)
    ok = ml.dropna(subset=["max_longevity_years", "tree_species"])
    ok = ok[ok["map_type"].isin(["exact", "subspecies_collapsed", "genus_change", "burgio_synonym", "spelling",
                                 "same_taxon"])]
    ok = ok[ok["data_quality"].fillna("") != "low"]
    # AnAge first (it records origin, sample size and quality); Amniote LHD fills gaps.
    ok = ok.assign(pref=(ok["source"] != "AnAge").astype(int)).sort_values(["tree_species", "pref"])
    best = ok.drop_duplicates("tree_species").set_index("tree_species")
    return best, ml


def sra_tiers(res):
    runs = pd.read_csv(ROOT / "data/processed/sra/sra_parrot_runs.tsv", sep="\t")
    # 'WGS' runs built by sequence capture or restriction digestion are not whole-genome coverage.
    runs = runs[(runs["assay_type"] == "WGS") & (runs["librarysource"] == "GENOMIC")
                & ~runs["libraryselection"].isin(["Hybrid Selection", "Restriction Digest"])].copy()
    runs[["tree_species", "map_type"]] = runs["organism"].apply(lambda n: pd.Series(res.resolve(n)))
    unresolved = sorted(runs.loc[runs["map_type"] == "unresolved", "organism"].unique())
    if unresolved:
        print("unresolved SRA names:", unresolved)
    runs = runs.dropna(subset=["tree_species"])
    runs["long"] = runs["platform"].isin(["PACBIO_SMRT", "OXFORD_NANOPORE"])
    runs["split_only"] = runs["map_type"] == "split_lumped_in_tree"
    runs["hains2022"] = runs["center_name"].fillna("").str.upper().eq("IRIDIAN GENOMES")
    per_ind = runs.groupby(["tree_species", "biosample", "long", "split_only"])["mbases"].sum().reset_index()
    per_ind["cov"] = per_ind["mbases"] / 1000 / GENOME_GB

    def best(df, long, allow_split):
        d = df[(df["long"] == long) & (allow_split | ~df["split_only"])]
        return d.groupby("tree_species")["cov"].max()

    out = pd.DataFrame({
        "best_short_cov_x": best(per_ind, False, False),
        "best_long_cov_x": best(per_ind, True, False),
        "best_short_cov_incl_split_taxa_x": best(per_ind, False, True),
    })
    g = runs.groupby("tree_species")
    out["n_wgs_biosamples"] = g["biosample"].nunique()
    out["wgs_bioprojects"] = g["bioproject"].agg(lambda s: ";".join(sorted(set(s.dropna()))))
    out["in_hains2022_iridian"] = g["hains2022"].any()
    out["sra_names"] = g["organism"].agg(lambda s: ";".join(sorted(set(s))))
    return out.fillna({"best_short_cov_x": 0, "best_long_cov_x": 0, "best_short_cov_incl_split_taxa_x": 0})


def main():
    res = NameResolver()
    sp = pd.read_csv(P0 / "species.tsv", sep="\t")
    ph = pd.read_csv(P0 / "phenotypes.tsv", sep="\t").set_index("species")
    ml, _ = max_longevity(res)
    sra = sra_tiers(res)

    m = sp[["species", "tree_tip", "family", "genus", "genome_best_tier"]].rename(
        columns={"genome_best_tier": "phase0_genome_tier"}).set_index("species")
    for c in ["life_expectancy_years", "log_life_expectancy", "log_life_expectancy_se", "body_mass_g", "log_body_mass"]:
        m[c] = ph[c]
    m["max_longevity_years"] = ml["max_longevity_years"]
    m["max_longevity_source"] = ml["source"]
    m["max_longevity_origin"] = ml["origin"]
    m["max_longevity_sample_size"] = ml["sample_size"]
    m["max_longevity_quality"] = ml["data_quality"]
    m = m.join(sra, how="left")
    for c in ["best_short_cov_x", "best_long_cov_x", "best_short_cov_incl_split_taxa_x"]:
        m[c] = m[c].fillna(0.0)
    m["n_wgs_biosamples"] = m["n_wgs_biosamples"].fillna(0).astype(int)
    m["in_hains2022_iridian"] = m["in_hains2022_iridian"].fillna(False).astype(bool)

    ga = pd.read_csv(P0 / "assemblies.tsv", sep="\t", dtype=str).fillna("")
    chrom = set(ga.loc[(ga["source"] == "GenomeArk/VGP S3 listing") &
                       ga["assembly_level_reported"].str.startswith("chromosome"), "tree_species"])
    lit = set(ga.loc[ga["source"] != "GenomeArk/VGP S3 listing", "tree_species"])
    m["genome_tier"] = np.select(
        [m.index.isin(chrom), m["best_long_cov_x"] >= 15, m["best_short_cov_x"] >= 20,
         m["best_short_cov_x"] >= 5, m["best_short_cov_x"] > 0, m.index.isin(lit)],
        TIERS[:-1], "none")
    m["usable_for_gene_level"] = m["genome_tier"].isin(TIERS[:3])
    m = m.reset_index()
    m.to_csv(OUT / "species_data_matrix.tsv", sep="\t", index=False, float_format="%.4g")

    has_le = m["log_life_expectancy"].notna()
    has_ml = m["max_longevity_years"].notna()
    print(f"max longevity: {has_ml.sum()} tree species ({(has_ml & has_le).sum()} also with life expectancy)")
    print(pd.crosstab(m["genome_tier"], has_le.map({True: "with_LE", False: "no_LE"})).reindex(TIERS).fillna(0)
          .astype(int).to_string())
    print(f"usable (tiers A-C) with LE: {(m['usable_for_gene_level'] & has_le).sum()}; "
          f"with max longevity: {(m['usable_for_gene_level'] & has_ml).sum()}")


if __name__ == "__main__":
    main()
