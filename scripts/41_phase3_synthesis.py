"""Synthesize Phase 3 RERconverge genome-wide results (H6, H8 within-clade, GO enrichment, and exploratory axes).

Reads:
  results/phase3_genome/rer_<tree>_<trait>.tsv
  results/phase3_genome/rer_<tree>_enrichment_<trait>.tsv
  results/phase3_genome/rer_<tree>_h8_<clade>_panel.tsv
  data/pilot_genes/canonical_genes.tsv

Writes:
  results/phase3_genome/h6_genome_summary.tsv
  results/phase3_genome/h6_top_genes.tsv
  results/phase3_genome/h6_shared_longevity_brain.tsv
  results/phase3_genome/h8_within_clade_summary.tsv
  results/phase3_genome/exploratory_axes_genes.tsv
  results/phase3_genome/exploratory_axes_summary.tsv
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase3_genome"

AXES = {
    "Neuroprotection_Proteostasis": [
        "PARK7", "SNCAIP", "SOD1", "SQSTM1", "BECN1", "TFEB", "ULK1", "ULK2",
        "MAP1LC3B", "LAMP2", "CTSD", "GBA", "VPS35", "LRRK2", "DNAJB1", "DNAJB6",
        "DNAJC5", "DNAJC6", "HSPA1A", "HSPA5", "HSPA8", "HSPA9", "HSP90AA1", "HSP90AB1",
        "HSP90B1", "VCP", "UBQLN2", "OPTN", "TBK1", "ATXN1", "ATXN2", "ATXN3",
    ],
    "VocalLearning_Language_Plasticity": [
        "FOXP1", "FOXP2", "FOXP4", "CNTNAP2", "ROBO1", "ROBO2", "SLIT1", "SLIT2",
        "BDNF", "NTRK2", "EGR1", "DRD1", "DRD2", "GAP43", "GRIN1", "GRIN2A", "GRIN2B",
        "NEUROD1", "NEUROD2", "NEUROD6", "MEF2C", "ARC", "FOS", "IGF2", "SRPX2",
    ],
    "Antioxidation_Redox": [
        "SOD1", "SOD2", "SOD3", "CAT", "GPX1", "GPX3", "GSR", "TXN", "TXN2",
        "TXNRD1", "TXNRD2", "PRDX1", "PRDX2", "PRDX4", "PRDX5", "PRDX6",
        "NFE2L2", "NQO1", "HMOX1", "GCLM", "SRXN1", "MSRA", "SELENOP", "GLRX", "GLRX2",
    ],
    "ExtremeLongevity_NutrientSensing": [
        "IGF1R", "INSR", "IRS1", "IRS2", "PIK3CA", "AKT1", "FOXO1", "FOXO3",
        "MTOR", "RPTOR", "RICTOR", "RPS6KB1", "EIF4EBP1", "SIRT1", "SIRT2",
        "SIRT3", "SIRT5", "SIRT6", "SIRT7", "PRKAA1", "PRKAA2", "KL", "CISD2", "NAMPT",
    ],
    "Cognition_Synapse_Pallium": [
        "NRXN1", "NRXN2", "NRXN3", "NLGN1", "NLGN2", "NLGN3", "SHANK1", "SHANK2",
        "SHANK3", "DLG4", "SYP", "SYN1", "CAMK2A", "CAMK2B", "CREB1", "RELN",
        "NCAM1", "SNAP25", "STX1A", "GRIA1", "GRIA2", "SYT1", "CASK", "HOMER1",
    ],
    "AntiCancer_DNARepair_Telomere": [
        "PTEN", "BRCA1", "BARD1", "PALB2", "CHEK2", "MDM2", "MDM4", "CDKN1A",
        "CDKN1B", "PMAIP1", "BBC3", "BAX", "BAK1", "CASP3", "CASP8", "MRE11",
        "RAD50", "RAD51", "PARP1", "POLB", "NUDT1", "DKC1", "TERF2IP", "TSC1", "TSC2",
    ],
}


def load_rer(tree, label):
    p = OUT / f"rer_{tree}_{label}.tsv"
    return pd.read_csv(p, sep="\t") if p.exists() else pd.DataFrame()


def main():
    panel_df = pd.read_csv(ROOT / "data/pilot_genes/canonical_genes.tsv", sep="\t")
    mod_map = dict(zip(panel_df.gene_symbol, panel_df.module))

    summary_rows = []
    top_rows = []
    for tree in ("burgio2019", "smith2024"):
        for trait in ("rel_log_life_expectancy", "rel_log_brain", "rel_brain_given_le",
                      "h8_Psittaculidae", "h8_rest"):
            df = load_rer(tree, trait)
            if df.empty:
                continue
            v = df.dropna(subset=["Rho", "P"])
            if v.empty:
                continue
            summary_rows.append({
                "tree": tree,
                "trait": trait,
                "n_genes": len(v),
                "median_N_branches": float(v["N"].median()),
                "n_padj_lt_0.10": int((v["p.adj"] < 0.10).sum()),
                "n_permQ_lt_0.10": int((v["permQ"] < 0.10).sum()),
                "n_P_lt_0.01": int((v["P"] < 0.01).sum()),
                "n_permP_lt_0.01": int((v["permP"] < 0.01).sum()),
                "n_both_P01_permP05": int(((v["P"] < 0.01) & (v["permP"] < 0.05)).sum()),
                "top_pos_gene": v.sort_values("Rho", ascending=False).iloc[0]["gene"],
                "top_pos_Rho": float(v["Rho"].max()),
                "top_neg_gene": v.sort_values("Rho", ascending=True).iloc[0]["gene"],
                "top_neg_Rho": float(v["Rho"].min()),
            })
            if trait in ("rel_log_life_expectancy", "rel_log_brain", "rel_brain_given_le"):
                sub = v[(v["P"] < 0.005) & (v["permP"] < 0.05)].sort_values("P").head(25).copy()
                sub["tree"] = tree
                sub["trait"] = trait
                top_rows.append(sub)

    pd.DataFrame(summary_rows).to_csv(OUT / "h6_genome_summary.tsv", sep="\t", index=False, float_format="%.5g")
    if top_rows:
        pd.concat(top_rows, ignore_index=True).to_csv(OUT / "h6_top_genes.tsv", sep="\t", index=False, float_format="%.5g")

    # Cross-trait correlation and shared genes (H6 core test)
    shared_rows = []
    for tree in ("burgio2019", "smith2024"):
        le = load_rer(tree, "rel_log_life_expectancy").set_index("gene")
        br = load_rer(tree, "rel_log_brain").set_index("gene")
        bg = load_rer(tree, "rel_brain_given_le").set_index("gene")
        m = le.join(br, lsuffix="_le", rsuffix="_brain").join(bg.add_suffix("_brain_given_le"))
        valid = m.dropna(subset=["Rho_le", "Rho_brain", "Rho_brain_given_le"])
        r_le_br, p_le_br = stats.pearsonr(valid["Rho_le"], valid["Rho_brain"])
        r_le_bg, p_le_bg = stats.pearsonr(valid["Rho_le"], valid["Rho_brain_given_le"])
        r_br_bg, p_br_bg = stats.pearsonr(valid["Rho_brain"], valid["Rho_brain_given_le"])
        print(f"[{tree}] Genome-wide Rho correlation (n={len(valid)}): "
              f"cor(LE, Brain) = {r_le_br:.3f} (p={p_le_br:.2e}); "
              f"cor(LE, Brain|LE) = {r_le_bg:.3f}; cor(Brain, Brain|LE) = {r_br_bg:.3f}")

        # Shared genes with nominal P < 0.05 and permP < 0.05 in both longevity and brain with same sign
        cond = (
            (valid["P_le"] < 0.05) & (valid["permP_le"] < 0.05) &
            (valid["P_brain"] < 0.05) & (valid["permP_brain"] < 0.05) &
            (np.sign(valid["Rho_le"]) == np.sign(valid["Rho_brain"]))
        )
        sh = valid[cond].copy()
        sh["tree"] = tree
        sh["gene"] = sh.index
        sh["direction"] = np.where(sh["Rho_le"] > 0, "accelerated_in_long_lived_large_brained",
                                   "constrained_in_long_lived_large_brained")
        sh["persists_brain_beyond_le"] = (
            (sh["P_brain_given_le"] < 0.05) &
            (np.sign(sh["Rho_brain_given_le"]) == np.sign(sh["Rho_brain"]))
        )
        shared_rows.append(sh.reset_index(drop=True))

    if shared_rows:
        sh_df = pd.concat(shared_rows, ignore_index=True)
        cols = ["tree", "gene", "direction", "persists_brain_beyond_le",
                "Rho_le", "P_le", "permP_le", "permQ_le",
                "Rho_brain", "P_brain", "permP_brain", "permQ_brain",
                "Rho_brain_given_le", "P_brain_given_le", "permP_brain_given_le", "permQ_brain_given_le"]
        sh_df[cols].sort_values(["tree", "P_le"]).to_csv(
            OUT / "h6_shared_longevity_brain.tsv", sep="\t", index=False, float_format="%.5g"
        )

    # H8 within-clade module summary
    h8_rows = []
    for tree in ("burgio2019", "smith2024"):
        for cl in ("Psittaculidae", "rest"):
            p = OUT / f"rer_{tree}_h8_{cl}_panel.tsv"
            if not p.exists():
                continue
            df = pd.read_csv(p, sep="\t")
            df["module"] = df["gene"].map(mod_map)
            v = df.dropna(subset=["Rho", "P"])
            if v.empty:
                continue
            for mod, gdf in v.groupby("module"):
                h8_rows.append({
                    "tree": tree,
                    "clade": cl,
                    "module": mod,
                    "n_genes": len(gdf),
                    "mean_Rho": float(gdf["Rho"].mean()),
                    "median_Rho": float(gdf["Rho"].median()),
                    "n_P_lt_0.05": int((gdf["P"] < 0.05).sum()),
                    "n_permP_lt_0.05": int((gdf["permP"] < 0.05).sum()),
                    "n_panelQ_lt_0.10": int((gdf["panelQ"] < 0.10).sum()),
                    "sig_genes_nominal": ",".join(
                        f"{r.gene}({r.Rho:+.2f},p={r.P:.3f},permP={r.permP:.3f})"
                        for r in gdf[gdf["P"] < 0.05].itertuples()
                    ),
                })
    pd.DataFrame(h8_rows).to_csv(OUT / "h8_within_clade_summary.tsv", sep="\t", index=False, float_format="%.5g")

    # Exploratory functional axes
    axis_gene_rows = []
    axis_sum_rows = []
    for tree in ("burgio2019", "smith2024"):
        le = load_rer(tree, "rel_log_life_expectancy").set_index("gene")
        br = load_rer(tree, "rel_log_brain").set_index("gene")
        bg = load_rer(tree, "rel_brain_given_le").set_index("gene")
        m = le.join(br, lsuffix="_le", rsuffix="_brain").join(bg.add_suffix("_brain_given_le"))
        for axis, genes in AXES.items():
            present = [g for g in genes if g in m.index]
            sub = m.loc[present].copy()
            for g, r in sub.iterrows():
                axis_gene_rows.append({
                    "tree": tree,
                    "axis": axis,
                    "gene": g,
                    "N_le": r["N_le"],
                    "Rho_le": r["Rho_le"],
                    "P_le": r["P_le"],
                    "permP_le": r["permP_le"],
                    "Rho_brain": r["Rho_brain"],
                    "P_brain": r["P_brain"],
                    "permP_brain": r["permP_brain"],
                    "Rho_brain_given_le": r["Rho_brain_given_le"],
                    "P_brain_given_le": r["P_brain_given_le"],
                    "permP_brain_given_le": r["permP_brain_given_le"],
                })
            for trait_sfx, label in (("_le", "rel_log_life_expectancy"),
                                     ("_brain", "rel_log_brain"),
                                     ("_brain_given_le", "rel_brain_given_le")):
                rho_col = f"Rho{trait_sfx}"
                p_col = f"P{trait_sfx}"
                pp_col = f"permP{trait_sfx}"
                fg = sub[rho_col].dropna()
                bg_vals = m.loc[~m.index.isin(present), rho_col].dropna()
                if len(fg) >= 3:
                    u_stat, mw_p = stats.mannwhitneyu(fg, bg_vals, alternative="two-sided")
                else:
                    mw_p = np.nan
                notable = sub[(sub[p_col] < 0.05) | (sub[pp_col] < 0.05)]
                axis_sum_rows.append({
                    "tree": tree,
                    "axis": axis,
                    "trait": label,
                    "n_requested": len(genes),
                    "n_present_in_7756": len(fg),
                    "mean_Rho": float(fg.mean()) if len(fg) else np.nan,
                    "median_Rho": float(fg.median()) if len(fg) else np.nan,
                    "mw_p_vs_genome": float(mw_p),
                    "n_nominal_sig": len(notable),
                    "notable_genes": "; ".join(
                        f"{g}(Rho={r[rho_col]:+.2f},P={r[p_col]:.3f},permP={r[pp_col]:.3f})"
                        for g, r in notable.sort_values(p_col).iterrows()
                    ),
                })

    pd.DataFrame(axis_gene_rows).to_csv(OUT / "exploratory_axes_genes.tsv", sep="\t", index=False, float_format="%.5g")
    pd.DataFrame(axis_sum_rows).to_csv(OUT / "exploratory_axes_summary.tsv", sep="\t", index=False, float_format="%.5g")
    print("Phase 3 synthesis tables written.")


if __name__ == "__main__":
    main()
