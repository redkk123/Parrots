"""Phase 0b analyses: outcome sensitivity (max longevity), realised power and sister-pair contrasts.

Uses the species data matrix from script 04. Phylogenetic regressions use the Python
implementation in lib/phylo.py; scripts/R/validate_phase0.R cross-checks them with ape/nlme.
"""
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent / "lib"))
from phylo import Tree, fit_pgls, gls_loglik  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
P0 = ROOT / "results/phase0"
OUT = ROOT / "results/phase0b"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(20260926)
NAMED = ["Nestor notabilis", "Cacatua moluccensis", "Cacatua galerita", "Amazona aestiva", "Ara macao",
         "Anodorhynchus hyacinthinus", "Melopsittacus undulatus", "Nymphicus hollandicus", "Psittacus erithacus"]
SIZE_ORDER = {"tiny": 0, "small": 1, "medium": 2, "large": 3, "huge": 4}


def design(*cols):
    return np.column_stack([np.ones(len(cols[0]))] + list(cols))


def subset_fit(full, df, ycol, xcols, model="lambda"):
    t = full.prune(set(df["tree_tip"]))
    order = t.tip_labels()
    d = df.set_index("tree_tip").loc[order]
    X = design(*[d[c].to_numpy(float) for c in xcols])
    f = fit_pgls(d[ycol].to_numpy(float), X, t.vcv(order), model)
    resid = d[ycol].to_numpy(float) - X @ f["beta"]
    return f, pd.Series(resid, index=order), t


def main():
    m = pd.read_csv(OUT / "species_data_matrix.tsv", sep="\t")
    full = Tree.from_newick((P0 / "supertree_burgio2019_full.nwk").read_text())
    summary = {}

    # ---- 1. maximum longevity: allometry, sample-size effect, agreement with life expectancy ----
    m["log_max_longevity"] = np.log(m["max_longevity_years"])
    mm = m.dropna(subset=["log_max_longevity", "log_body_mass"])
    f_ml, r_ml, _ = subset_fit(full, mm, "log_max_longevity", ["log_body_mass"])
    summary["max_longevity_allometry"] = {"n": len(mm), "lambda": f_ml["param"], "slope": float(f_ml["beta"][1]),
                                          "slope_se": float(f_ml["se"][1])}
    anage = mm[mm["max_longevity_sample_size"].isin(SIZE_ORDER)].copy()
    anage["sample_size_rank"] = anage["max_longevity_sample_size"].map(SIZE_ORDER)
    f_ss, _, _ = subset_fit(full, anage, "log_max_longevity", ["log_body_mass", "sample_size_rank"])
    summary["max_longevity_sample_size_effect"] = {"n": len(anage), "coef_per_size_class": float(f_ss["beta"][2]),
                                                   "se": float(f_ss["se"][2])}
    m["rel_log_max_longevity"] = m["tree_tip"].map(r_ml)
    # Records grow with the number of animals observed (extreme-value sampling). Residuals from the
    # model with AnAge sample-size class are the sample-size-adjusted outcome.
    _, r_ss, _ = subset_fit(full, anage, "log_max_longevity", ["log_body_mass", "sample_size_rank"])
    m["rel_log_max_longevity_ssadj"] = m["tree_tip"].map(r_ss)

    le = m.dropna(subset=["log_life_expectancy", "log_body_mass"])
    f_le, r_le, _ = subset_fit(full, le, "log_life_expectancy", ["log_body_mass"])
    m["rel_log_life_expectancy"] = m["tree_tip"].map(r_le)
    both = m.dropna(subset=["rel_log_life_expectancy", "rel_log_max_longevity"])
    f_ag, _, _ = subset_fit(full, both, "rel_log_max_longevity", ["rel_log_life_expectancy"])
    rho = stats.spearmanr(both["rel_log_life_expectancy"], both["rel_log_max_longevity"])
    both_ss = m.dropna(subset=["rel_log_life_expectancy", "rel_log_max_longevity_ssadj"])
    rho_ss = stats.spearmanr(both_ss["rel_log_life_expectancy"], both_ss["rel_log_max_longevity_ssadj"])
    summary["agreement_rel_LE_vs_rel_maxL"] = {
        "n": len(both), "pgls_slope": float(f_ag["beta"][1]), "pgls_slope_se": float(f_ag["se"][1]),
        "spearman_rho": float(rho.statistic), "spearman_p": float(rho.pvalue),
        "n_ssadj": len(both_ss), "spearman_rho_ssadj": float(rho_ss.statistic)}

    pr = m[m["species"].isin(NAMED)].copy()
    pr["rank_rel_LE"] = m["rel_log_life_expectancy"].rank(ascending=False)[pr.index]
    pr["rank_rel_maxL"] = m["rel_log_max_longevity"].rank(ascending=False)[pr.index]
    pr["pct_rel_LE"] = (np.exp(pr["rel_log_life_expectancy"]) - 1) * 100
    pr["pct_rel_maxL"] = (np.exp(pr["rel_log_max_longevity"]) - 1) * 100
    pr["rank_rel_maxL_ssadj"] = m["rel_log_max_longevity_ssadj"].rank(ascending=False)[pr.index]
    pr["pct_rel_maxL_ssadj"] = (np.exp(pr["rel_log_max_longevity_ssadj"]) - 1) * 100
    pr[["species", "body_mass_g", "life_expectancy_years", "max_longevity_years", "max_longevity_sample_size",
        "pct_rel_LE", "rank_rel_LE", "pct_rel_maxL", "rank_rel_maxL", "pct_rel_maxL_ssadj", "rank_rel_maxL_ssadj",
        "genome_tier"]] \
        .sort_values("pct_rel_maxL", ascending=False) \
        .to_csv(OUT / "priority_species_two_outcomes.tsv", sep="\t", index=False, float_format="%.3g")
    summary["n_ranked"] = {"rel_LE": int(m["rel_log_life_expectancy"].notna().sum()),
                           "rel_maxL": int(m["rel_log_max_longevity"].notna().sum()),
                           "rel_maxL_ssadj": int(m["rel_log_max_longevity_ssadj"].notna().sum())}

    # ---- 2. realised overlap and power ------------------------------------------------
    use = m[m["usable_for_gene_level"] & m["rel_log_life_expectancy"].notna()]
    summary["usable_with_LE"] = int(len(use))
    summary["usable_with_LE_by_family"] = use["family"].value_counts().to_dict()
    summary["usable_with_LE_by_tier"] = use["genome_tier"].value_counts().to_dict()
    summary["usable_with_maxL"] = int((m["usable_for_gene_level"] & m["rel_log_max_longevity"].notna()).sum())
    summary["rel_LE_spread"] = {"all_218_sd": float(m["rel_log_life_expectancy"].std()),
                                "usable_sd": float(use["rel_log_life_expectancy"].std())}
    scen = {"A_B_C_usable_with_LE": set(use["tree_tip"]),
            "A_B_C_D_with_LE": set(m.loc[m["genome_tier"].isin(["A_chromosome_scale", "B_long_read_ge15x",
                                                               "C_short_read_ge20x", "D_short_read_5_20x"])
                                         & m["rel_log_life_expectancy"].notna(), "tree_tip"])}
    rows = []
    for name, keep in scen.items():
        t = full.prune(keep)
        order = t.tip_labels()
        C = t.vcv(order)
        L = np.linalg.cholesky(C + 1e-10 * np.eye(len(order)))
        for r in (0.2, 0.3, 0.4, 0.6):
            for alpha in (0.05, 0.05 / 1e4):
                hits, nsim = 0, 500
                for _ in range(nsim):
                    z1, z2 = L @ RNG.standard_normal(len(order)), L @ RNG.standard_normal(len(order))
                    ll, b, cb, _ = gls_loglik(r * z1 + np.sqrt(1 - r ** 2) * z2, design(z1), C)
                    hits += 2 * stats.t.sf(abs(b[1] / np.sqrt(cb[1, 1])), len(order) - 2) < alpha
                rows.append({"scenario": name, "n_species": len(order), "evolutionary_r": r, "alpha": alpha,
                             "power": hits / nsim, "n_sim": nsim})
    pd.DataFrame(rows).to_csv(OUT / "power_realised_sets.tsv", sep="\t", index=False, float_format="%.4g")

    # ---- 3. sister-pair contrasts among usable species (independent by construction) ----------
    t = full.prune(set(use["tree_tip"]))
    depth = t.depths()
    H = depth[t.tips()].max()
    rel = dict(zip(use["tree_tip"], use["rel_log_life_expectancy"]))
    relm = dict(zip(m["tree_tip"], m["rel_log_max_longevity"]))
    tier = dict(zip(m["tree_tip"], m["genome_tier"]))
    pairs = []
    for n in range(len(t.parent)):
        kids = t.children[n]
        if len(kids) == 2 and all(t.is_tip(k) for k in kids):
            a, b = sorted((t.label[k] for k in kids), key=lambda x: -rel[x])
            pairs.append({"longer_lived": a.replace("_", " "), "shorter_lived": b.replace("_", " "),
                          "divergence_myr": H - depth[n], "delta_rel_logLE": rel[a] - rel[b],
                          "delta_rel_logMaxL": relm.get(a, np.nan) - relm.get(b, np.nan),
                          "tiers": f"{tier[a]} | {tier[b]}"})
    pairs = pd.DataFrame(pairs).sort_values("delta_rel_logLE", ascending=False)
    pairs["same_sign_maxL"] = (pairs["delta_rel_logMaxL"] > 0).where(pairs["delta_rel_logMaxL"].notna())
    pairs.to_csv(OUT / "sister_pair_contrasts.tsv", sep="\t", index=False, float_format="%.3g")
    summary["sister_pairs"] = {"n": len(pairs), "n_delta_ge_0.3": int((pairs["delta_rel_logLE"] >= 0.3).sum()),
                               "n_with_maxL_both": int(pairs["delta_rel_logMaxL"].notna().sum()),
                               "n_maxL_agrees": int(pairs["same_sign_maxL"].sum())}

    # ---- 4. figure: the two outcomes ------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.axhline(0, color="grey", lw=.5)
    ax.axvline(0, color="grey", lw=.5)
    ax.scatter(both["rel_log_life_expectancy"], both["rel_log_max_longevity"], s=12, alpha=.6, c="#3b6ea5")
    for _, r_ in both[both["species"].isin(NAMED)].iterrows():
        ax.annotate(r_["species"], (r_["rel_log_life_expectancy"], r_["rel_log_max_longevity"]), fontsize=6,
                    xytext=(3, -3), textcoords="offset points")
    ax.set_xlabel("relative log life expectancy (captive, at birth)")
    ax.set_ylabel("relative log maximum longevity (AnAge / Amniote LHD)")
    ax.set_title(f"Two longevity outcomes, {len(both)} species (Spearman ρ = {rho.statistic:.2f})", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "fig3_two_outcomes.png", dpi=200)
    plt.close(fig)

    # QC: a record below ~1.3x the mean life expectancy is implausible for the same population
    m["maxL_over_LE"] = m["max_longevity_years"] / m["life_expectancy_years"]
    m["qc_outcomes_inconsistent"] = m["maxL_over_LE"] < 1.3
    summary["qc_outcomes_inconsistent"] = m.loc[m["qc_outcomes_inconsistent"], "species"].tolist()
    m.to_csv(OUT / "species_data_matrix_with_residuals.tsv", sep="\t", index=False, float_format="%.5g")
    (OUT / "phase0b_summary.json").write_text(json.dumps(summary, indent=2, default=float))
    print(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    main()
