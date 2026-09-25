"""Phase 0: harmonise phenotypes, taxonomy, tree and genome availability.

Inputs (data/external, fetched by 00_fetch_external.sh):
  - Smeele et al. 2022 processed per-species summaries (RData)
  - Burgio et al. 2019 supertree (ape phylo, via the parrottraits package) and taxonomy
  - GenomeArk parrot listing
  - data/curated/assemblies_literature.tsv (hand-curated, verification status per row)
Outputs:
  - results/phase0/species_tree.nwk (+ full supertree)
  - results/phase0/phenotypes.tsv, species.tsv, assemblies.tsv
"""
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadr
import rdata

sys.path.insert(0, str(Path(__file__).parent / "lib"))
from phylo import Tree  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "data/external"
OUT = ROOT / "results/phase0"
OUT.mkdir(parents=True, exist_ok=True)
SMEELE = EXT / "smeele2022/ANALYSIS/RESULTS"

# Species the protocol names as priorities for the audit (no role assigned a priori).
PRIORITY = ["Nestor notabilis", "Cacatua moluccensis", "Cacatua galerita", "Amazona aestiva", "Ara macao",
            "Anodorhynchus hyacinthinus", "Melopsittacus undulatus", "Nymphicus hollandicus"]
PRIORITY_GENERA = ["Agapornis", "Amazona"]


def load_rdata(sub):
    return pyreadr.read_r(str(SMEELE / sub / "master_dat.RData"))["master_dat"]


def main():
    # ---- tree ---------------------------------------------------------------
    st = rdata.read_rda(str(EXT / "parrottraits/data/supertree.rda"))["supertree"]["supertree"]
    tree = Tree.from_ape(st["edge"], st["edge.length"], list(st["tip.label"]))
    (OUT / "supertree_burgio2019_full.nwk").write_text(tree.to_newick() + "\n")
    tips = tree.tip_labels()
    print(f"supertree: {len(tips)} tips, ultrametric={tree.is_ultrametric()}, "
          f"height={tree.depths()[tree.tips()].max():.2f} Myr")

    # ---- phenotypes (Smeele et al. 2022) ------------------------------------
    le = load_rdata("life expectancy")
    wt = load_rdata("weight")
    br = load_rdata("brain")
    cl = load_rdata("clutch size")
    dv = load_rdata("developmental time")
    df = le.merge(wt, on="species", how="outer").merge(br, on="species", how="outer") \
           .merge(cl, on="species", how="outer").merge(dv, on="species", how="outer")
    afr = pd.read_csv(SMEELE / "AFR/5th percentile results.csv", sep=None, engine="python")
    afr_cols = [c for c in afr.columns if c.lower() != "species"]
    print("AFR file columns:", list(afr.columns))
    df["tree_tip"] = df["species"].str.replace(" ", "_")
    df["in_tree"] = df["tree_tip"].isin(tips)

    ph = pd.DataFrame({
        "species": df["species"],
        "tree_tip": df["tree_tip"],
        "in_tree": df["in_tree"],
        # life expectancy: primary candidate outcome
        "life_expectancy_years": np.exp(df["log_mean_life_exp"]).round(3),
        "log_life_expectancy": df["log_mean_life_exp"],
        "log_life_expectancy_se": df["log_SE_life_exp"],
        "life_expectancy_lo95_years": np.exp(df["log_mean_life_exp"] - 1.96 * df["log_SE_life_exp"]).round(3),
        "life_expectancy_hi95_years": np.exp(df["log_mean_life_exp"] + 1.96 * df["log_SE_life_exp"]).round(3),
        "life_expectancy_definition": np.where(df["log_mean_life_exp"].notna(),
            "mean life expectancy at birth; Gompertz bathtub survival model (BaSTA) fitted to ZIMS zoo records; sexes pooled",
            ""),
        "life_expectancy_setting": np.where(df["log_mean_life_exp"].notna(), "captive (ZIMS member institutions)", ""),
        "life_expectancy_source": np.where(df["log_mean_life_exp"].notna(), "Smeele et al. 2022 (processed; raw records not public)", ""),
        # maximum longevity: kept as a separate, not-yet-compiled measure
        "max_longevity_years": np.nan,
        "max_longevity_source": "not compiled (AnAge/ZIMS max records pending; kept separate from life expectancy)",
        "body_mass_g": np.exp(df["log_mean_body_weight"]).round(2),
        "log_body_mass": df["log_mean_body_weight"],
        "log_body_mass_se": df["log_SE_body_weight"],
        "body_mass_source": np.where(df["log_mean_body_weight"].notna(), "Smeele et al. 2022 (ZIMS + literature)", ""),
        "brain_mass_g": np.exp(df["log_mean_brain_size"]).round(3),
        "log_brain_mass": df["log_mean_brain_size"],
        "log_brain_mass_se": df["log_SE_brain_size"],
        "clutch_size": df["clutch_size_n"],
        "incubation_days": df["incubation_days"],
        "fledging_age_days": df["fledging_age_days"],
    })
    if afr_cols:
        afr = afr.rename(columns={afr_cols[0]: "age_first_repro_p05_years_raw"})
        ph = ph.merge(afr[["species", "age_first_repro_p05_years_raw"]], on="species", how="left")
    ph = ph.sort_values("species").reset_index(drop=True)
    ph.to_csv(OUT / "phenotypes.tsv", sep="\t", index=False, na_rep="NA")

    # ---- taxonomy ------------------------------------------------------------
    tax = pd.read_csv(EXT / "parrottraits/data-raw/parrot-traits-database.csv", header=1, encoding="utf-8-sig")
    tax = tax.rename(columns={"taxon_2 (matching)": "tree_tip", "common name": "common_name"})
    tax = tax[["tree_tip", "family", "genus", "common_name"]].dropna(subset=["tree_tip"]).drop_duplicates("tree_tip")

    # ---- assemblies ----------------------------------------------------------
    ga = pd.read_csv(EXT / "genomeark/genomeark_parrots.tsv", sep="\t", dtype=str).fillna("")
    in_progress = ga["assembly_folders"].str.contains("assembly_") & (ga["latest_curated_fasta"] == "")
    ga = ga[(ga["latest_curated_fasta"] != "") | in_progress].copy()
    ga_rows = pd.DataFrame({
        "species": ga["species"], "assembly_name_or_tolid": ga["tolid"],
        "ncbi_accession_version": "to_verify_in_NCBI",
        "source": "GenomeArk/VGP S3 listing", "file_or_reference": ga["latest_curated_fasta"],
        "assembly_level_reported": np.where(ga["latest_curated_fasta"] != "", "chromosome-scale curated (VGP)",
                                            "VGP assembly in progress (only intermediate files public)"),
        "total_length_bp": ga["total_scaffold_length_bp"], "n_scaffolds": ga["n_scaffolds"],
        "scaffold_n50_bp": ga["scaffold_n50_bp"], "n_contigs": ga["n_contigs"],
        "contig_n50_bp": ga["contig_n50_bp"], "busco_complete_pct": "",
        "annotation": "to_verify", "verification_status": "verified_on_genomeark_s3",
        "retrieved": "2026-09-25",
    })
    lit = pd.read_csv(ROOT / "data/curated/assemblies_literature.tsv", sep="\t", dtype=str).fillna("")
    asm = pd.concat([ga_rows, lit], ignore_index=True).fillna("")
    syn = pd.read_csv(ROOT / "data/curated/taxonomy_synonyms.tsv", sep="\t")
    asm.insert(1, "tree_species", asm["species"].replace(dict(zip(syn["name_in_source"], syn["tree_species_burgio2019"]))))
    asm.to_csv(OUT / "assemblies.tsv", sep="\t", index=False)

    best = asm.assign(tier=np.select(
        [asm["assembly_level_reported"].str.startswith("chromosome"),
         asm["assembly_level_reported"].str.contains("short-read|draft|scaffold|in progress", regex=True)],
        [2, 1], 0)).groupby("tree_species")["tier"].max()

    # ---- species table --------------------------------------------------------
    sp = pd.DataFrame({"tree_tip": tips})
    sp["species"] = sp["tree_tip"].str.replace("_", " ")
    sp = sp.merge(tax, on="tree_tip", how="left")
    sp["genus"] = sp["genus"].fillna(sp["species"].str.split().str[0])
    have = ph.set_index("tree_tip")
    sp["has_life_expectancy"] = sp["tree_tip"].map(have["log_life_expectancy"].notna()).fillna(False)
    sp["has_body_mass"] = sp["tree_tip"].map(have["log_body_mass"].notna()).fillna(False)
    sp["has_brain_mass"] = sp["tree_tip"].map(have["log_brain_mass"].notna()).fillna(False)
    sp["has_max_longevity"] = False
    sp["genome_best_tier"] = sp["species"].map(best).fillna(-1).astype(int).map(
        {2: "chromosome_scale", 1: "draft_or_in_progress", 0: "reported_level_unknown", -1: "none_found_in_this_audit"})
    sp["protocol_priority"] = sp["species"].isin(PRIORITY) | sp["genus"].isin(PRIORITY_GENERA)
    sp["role_in_design"] = "undetermined (to be set from adjusted phenotype and tree position)"
    sp.to_csv(OUT / "species.tsv", sep="\t", index=False)

    # ---- pruned analysis tree --------------------------------------------------
    keep = ph.loc[ph["in_tree"] & ph["log_life_expectancy"].notna() & ph["log_body_mass"].notna(), "tree_tip"]
    pruned = tree.prune(set(keep))
    (OUT / "species_tree.nwk").write_text(pruned.to_newick() + "\n")

    missing = ph.loc[~ph["in_tree"] & ph["log_life_expectancy"].notna(), "species"].tolist()
    print(f"phenotype rows: {len(ph)}; with life expectancy: {ph['log_life_expectancy'].notna().sum()}; "
          f"with LE+mass in tree: {len(keep)}; LE species not in tree: {missing}")
    print(sp.groupby(["genome_best_tier", "has_life_expectancy"]).size().to_string())


if __name__ == "__main__":
    main()
