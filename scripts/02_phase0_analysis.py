"""Phase 0 analyses: allometry, phylogenetic signal, ancestral estimates, sensitivity, power.

All results are exploratory feasibility checks. They define candidate contrasts from the
phenotype alone, before any molecular data are examined (protocol, "Fase 0").
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
from phylo import Tree, bm_ancestral, fit_pgls, gls_loglik, lambda_transform  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase0"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(20260925)
SE_CUTOFF = 0.2  # log-scale SE above which a life-expectancy estimate is treated as imprecise


def design(x):
    return np.column_stack([np.ones(len(x)), x])


def main():
    ph = pd.read_csv(OUT / "phenotypes.tsv", sep="\t")
    sp = pd.read_csv(OUT / "species.tsv", sep="\t")
    tree = Tree.from_newick((OUT / "species_tree.nwk").read_text())
    tips = tree.tip_labels()
    d = ph.set_index("tree_tip").loc[tips].reset_index()
    d = d.merge(sp[["tree_tip", "family", "genus", "genome_best_tier", "protocol_priority"]], on="tree_tip", how="left")
    y, x = d["log_life_expectancy"].to_numpy(), d["log_body_mass"].to_numpy()
    C = tree.vcv(tips)
    summary = {"n_species": len(tips), "tree_height_myr": float(np.max(np.diag(C)))}

    # ---- 1. allometry under alternative covariance models -----------------------
    fits = {m: fit_pgls(y, design(x), C, m) for m in ["BM", "lambda", "OU"]}
    ols = stats.linregress(x, y)
    rows = [{"model": "OLS (no phylogeny)", "param": "", "slope": ols.slope, "slope_se": ols.stderr,
             "intercept": ols.intercept, "loglik": np.nan, "aic": np.nan}]
    for m, f in fits.items():
        rows.append({"model": m, "param": f["param"], "slope": f["beta"][1], "slope_se": f["se"][1],
                     "intercept": f["beta"][0], "loglik": f["loglik"], "aic": f["aic"]})
    mf = pd.DataFrame(rows)
    mf["delta_aic"] = mf["aic"] - mf["aic"].min()
    mf.to_csv(OUT / "model_fits.tsv", sep="\t", index=False, float_format="%.5g")
    best = fits["lambda"]
    summary["allometry"] = {"lambda": best["param"], "slope": float(best["beta"][1]), "slope_se": float(best["se"][1]),
                            "best_model_by_aic": mf.loc[mf["aic"].idxmin(), "model"]}

    # Residual ("relative") life expectancy from the lambda PGLS. The fitted line is used
    # as a fixed reference; its uncertainty is reported in model_fits.tsv.
    resid = y - design(x) @ best["beta"]
    d["rel_log_life_expectancy"] = resid

    # ---- 2. phylogenetic signal of relative life expectancy ------------------------
    one = np.ones((len(y), 1))
    lam_fit = fit_pgls(resid, one, C, "lambda")
    ll0 = gls_loglik(resid, one, lambda_transform(C / C.max(), 0.0))[0]
    lr = 2 * (lam_fit["loglik"] - ll0)
    summary["signal_relative_LE"] = {"lambda": lam_fit["param"], "LR_vs_lambda0": lr,
                                     "p_value": float(stats.chi2.sf(lr, 1) / 2)}  # boundary: half chi2(1)
    lam_raw = fit_pgls(y, one, C, "lambda")
    summary["signal_log_LE"] = {"lambda": lam_raw["param"]}

    # ---- 3. BM ancestral estimates of relative life expectancy ------------------------
    est, se, mu, s2 = bm_ancestral(tree, dict(zip(tips, resid)))
    depths = tree.depths()
    H = depths[tree.tips()].max()
    desc = tree.descendants_tips()
    anc_rows = []
    for n in range(len(tree.parent)):
        p = tree.parent[n]
        clade = [tree.label[t] for t in desc[n]]
        genera = sorted({c.split("_")[0] for c in clade})
        anc_rows.append({
            "node": n, "is_tip": tree.is_tip(n), "label": tree.label[n] or f"node{n}",
            "clade_key": f"{len(clade)}|{min(clade)}|{max(clade)}",
            "age_myr": H - depths[n], "n_desc_tips": len(clade),
            "genera": ";".join(genera) if len(genera) <= 6 else f"{len(genera)} genera",
            "est_rel_logLE": est[n], "se": se[n],
            "parent_est": est[p] if p >= 0 else np.nan,
            "branch_change": est[n] - est[p] if p >= 0 else np.nan,
            "branch_length_myr": tree.length[n],
        })
    anc = pd.DataFrame(anc_rows)
    anc.to_csv(OUT / "ancestral_states_rel_logLE.tsv", sep="\t", index=False, float_format="%.4g")
    top = anc[anc["branch_change"].notna()].sort_values("branch_change", ascending=False)
    cols = ["label", "genera", "n_desc_tips", "age_myr", "branch_length_myr", "parent_est", "est_rel_logLE",
            "branch_change", "se"]
    pd.concat([top.head(20).assign(direction="increase"), top.tail(20).iloc[::-1].assign(direction="decrease")]) \
        [["direction"] + cols].to_csv(OUT / "largest_branch_changes.tsv", sep="\t", index=False, float_format="%.4g")
    summary["ancestral_root_rel_logLE"] = {"estimate": float(est[tree.root]), "se": float(se[tree.root])}

    # ---- 4. genus summaries and priority lineages ------------------------------------
    d["life_expectancy_years"] = np.exp(y)
    gs = d.groupby("genus").agg(n=("tree_tip", "size"), family=("family", "first"),
                                mean_rel_logLE=("rel_log_life_expectancy", "mean"),
                                min_rel=("rel_log_life_expectancy", "min"), max_rel=("rel_log_life_expectancy", "max"),
                                median_LE_years=("life_expectancy_years", "median"),
                                genome_species=("genome_best_tier", lambda s: int((s != "none_found_in_this_audit").sum())))
    gs = gs.sort_values("mean_rel_logLE", ascending=False)
    gs.to_csv(OUT / "genus_relative_LE.tsv", sep="\t", float_format="%.4g")

    pri = d[d["protocol_priority"]].copy()
    allrank = d["rel_log_life_expectancy"].rank(ascending=False)
    pri["rank_rel_LE_of_218"] = allrank[pri.index]
    pri["pct_rel_LE_vs_expected"] = (np.exp(pri["rel_log_life_expectancy"]) - 1) * 100
    pri[["species", "family", "life_expectancy_years", "life_expectancy_lo95_years", "life_expectancy_hi95_years",
         "log_life_expectancy_se", "body_mass_g", "rel_log_life_expectancy", "pct_rel_LE_vs_expected",
         "rank_rel_LE_of_218", "genome_best_tier"]].sort_values("rel_log_life_expectancy", ascending=False) \
        .to_csv(OUT / "priority_species_phenotypes.tsv", sep="\t", index=False, float_format="%.4g")

    # ---- 5. sensitivity of the allometric slope ------------------------------------
    sens = []

    def refit(mask, label):
        t2 = tree.prune(set(np.array(tips)[mask]))
        order = t2.tip_labels()
        dd = d.set_index("tree_tip").loc[order]
        f = fit_pgls(dd["log_life_expectancy"].to_numpy(), design(dd["log_body_mass"].to_numpy()), t2.vcv(order), "lambda")
        sens.append({"analysis": label, "n": len(order), "lambda": f["param"], "slope": f["beta"][1], "slope_se": f["se"][1]})

    refit(np.ones(len(tips), bool), "all species")
    refit(d["log_life_expectancy_se"].to_numpy() <= SE_CUTOFF, f"exclude log-SE > {SE_CUTOFF}")
    for fam in d["family"].dropna().unique():
        if (d["family"] == fam).sum() >= 3 and (d["family"] != fam).sum() >= 10:
            refit((d["family"] != fam).to_numpy(), f"drop family {fam}")
    for g in gs.index[gs["n"] >= 5]:
        refit((d["genus"] != g).to_numpy(), f"drop genus {g}")
    pd.DataFrame(sens).to_csv(OUT / "sensitivity_allometry.tsv", sep="\t", index=False, float_format="%.4g")
    summary["n_imprecise_LE"] = int((d["log_life_expectancy_se"] > SE_CUTOFF).sum())

    # ---- 6. overlap with genomes and power ------------------------------------------
    gen = d[d["genome_best_tier"] != "none_found_in_this_audit"]
    summary["overlap_LE_and_any_genome"] = int(len(gen))
    summary["overlap_LE_and_chromosome_scale"] = int((d["genome_best_tier"] == "chromosome_scale").sum())
    gen[["species", "family", "genome_best_tier", "life_expectancy_years", "rel_log_life_expectancy"]] \
        .sort_values("rel_log_life_expectancy").to_csv(OUT / "genome_phenotype_overlap.tsv", sep="\t", index=False,
                                                       float_format="%.4g")

    power_rows = []
    scenarios = {"current_overlap": set(gen["tree_tip"]), "all_218_with_LE": set(tips)}
    # the Hains et al. 2022 batch (94 spp.) is not yet mapped; emulate its size by random draws
    for n_target in (40, 94):
        scenarios[f"random_{n_target}_species"] = set(RNG.choice(tips, n_target, replace=False))
    for name, keep in scenarios.items():
        t2 = tree.prune(keep)
        order = t2.tip_labels()
        C2 = t2.vcv(order)
        L = np.linalg.cholesky(C2 + 1e-10 * np.eye(len(order)))
        for r in (0.2, 0.4, 0.6):
            for alpha in (0.05, 0.05 / 1e4):
                hits = 0
                nsim = 400
                for _ in range(nsim):
                    z1, z2 = L @ RNG.standard_normal(len(order)), L @ RNG.standard_normal(len(order))
                    xs, ys = z1, r * z1 + np.sqrt(1 - r ** 2) * z2  # correlated BM evolution
                    ll, b, cb, _ = gls_loglik(ys, design(xs), C2)
                    tval = b[1] / np.sqrt(cb[1, 1])
                    hits += 2 * stats.t.sf(abs(tval), len(order) - 2) < alpha
                power_rows.append({"scenario": name, "n_species": len(order), "evolutionary_r": r,
                                   "alpha": alpha, "power": hits / nsim, "n_sim": nsim})
    pd.DataFrame(power_rows).to_csv(OUT / "power_simulation.tsv", sep="\t", index=False, float_format="%.4g")

    # ---- 7. figures -------------------------------------------------------------
    fam_col = {"Psittacidae": "#3b6ea5", "Cacatuidae": "#c8553d", "Strigopidae": "#2a9d8f"}
    fig, ax = plt.subplots(figsize=(7, 5))
    for fam, sub in d.groupby("family"):
        ax.scatter(sub["log_body_mass"], sub["log_life_expectancy"], s=14, alpha=.7, label=fam, c=fam_col.get(fam, "grey"))
    xx = np.linspace(x.min(), x.max(), 50)
    ax.plot(xx, best["beta"][0] + best["beta"][1] * xx, "k-", lw=1, label=f"PGLS λ={best['param']:.2f}")
    named = ["Nestor notabilis", "Cacatua moluccensis", "Cacatua galerita", "Amazona aestiva", "Ara macao",
             "Anodorhynchus hyacinthinus", "Melopsittacus undulatus", "Nymphicus hollandicus", "Agapornis roseicollis"]
    for _, r_ in d[d["species"].isin(named)].iterrows():
        ax.annotate(r_["species"], (r_["log_body_mass"], r_["log_life_expectancy"]), fontsize=6,
                    xytext=(3, -3), textcoords="offset points")
    ax.set_xlabel("log body mass (g)")
    ax.set_ylabel("log life expectancy at birth (years, captive)")
    ax.legend(fontsize=7, frameon=False)
    ax.set_title(f"Life expectancy vs body mass, {len(tips)} parrot species (Smeele et al. 2022 data)", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "fig2_allometry.png", dpi=200)
    plt.close(fig)

    # tree with relative LE on tips and BM ancestral estimates on branches
    ypos = {}
    for i, t in enumerate(tree.tips()):
        ypos[t] = i
    for n in tree.postorder():
        if not tree.is_tip(n):
            ypos[n] = np.mean([ypos[c] for c in tree.children[n]])
    vmax = np.nanpercentile(np.abs(est), 98)
    cmap = plt.get_cmap("RdBu")
    fig, ax = plt.subplots(figsize=(8, 26))
    for n in range(len(tree.parent)):
        p = tree.parent[n]
        c = cmap(0.5 + 0.5 * np.clip(est[n] / vmax, -1, 1))
        if p >= 0:
            ax.plot([depths[p], depths[n]], [ypos[n], ypos[n]], color=c, lw=1.2)
        if not tree.is_tip(n):
            ys_ = [ypos[ch] for ch in tree.children[n]]
            ax.plot([depths[n], depths[n]], [min(ys_), max(ys_)], color=c, lw=1.2)
    gb = dict(zip(d["tree_tip"], d["genome_best_tier"]))
    for t in tree.tips():
        lab = tree.label[t].replace("_", " ")
        mark = {"chromosome_scale": " ●", "draft_or_in_progress": " ○", "reported_level_unknown": " ◌"}.get(gb.get(tree.label[t]), "")
        ax.text(depths[t] + 0.3, ypos[t], lab + mark, fontsize=4.2, va="center")
    ax.set_xlim(0, H * 1.35)
    ax.set_ylim(-1, len(tree.tips()))
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.set_xlabel("Myr from root (Burgio et al. 2019 supertree)")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(-vmax, vmax))
    fig.colorbar(sm, ax=ax, shrink=0.2, label="relative log life expectancy (BM estimate)")
    ax.set_title("Relative life expectancy mapped on the tree (exploratory BM reconstruction)\n"
                 "● chromosome-scale genome, ○ draft/in progress, ◌ reported (level unknown)", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig1_tree_relative_LE.png", dpi=200)
    plt.close(fig)

    # ---- 8. phylogenetic neighbours of Melopsittacus (cancer axis comparators) ------------
    full = Tree.from_newick((OUT / "supertree_burgio2019_full.nwk").read_text())
    fd = full.depths()
    fH = fd[full.tips()].max()
    anc_sets = {}
    for t in full.tips():
        chain, n = [], t
        while n >= 0:
            chain.append(n)
            n = full.parent[n]
        anc_sets[full.label[t]] = chain
    focal = anc_sets["Melopsittacus_undulatus"]
    focal_set = set(focal)
    phl = ph.set_index("tree_tip")
    rel = []
    for lab, chain in anc_sets.items():
        if lab == "Melopsittacus_undulatus":
            continue
        mrca = next(n for n in chain if n in focal_set)
        rel.append({"species": lab.replace("_", " "), "divergence_from_Melopsittacus_myr": fH - fd[mrca],
                    "has_life_expectancy": bool(lab in phl.index and pd.notna(phl.loc[lab, "log_life_expectancy"])),
                    "life_expectancy_years": phl.loc[lab, "life_expectancy_years"] if lab in phl.index else np.nan,
                    "genome_best_tier": sp.set_index("tree_tip").loc[lab, "genome_best_tier"]})
    rel = pd.DataFrame(rel).sort_values(["divergence_from_Melopsittacus_myr", "species"])
    rel.to_csv(OUT / "melopsittacus_phylogenetic_neighbours.tsv", sep="\t", index=False, float_format="%.3g")
    summary["melopsittacus_neighbours_within_20myr"] = int((rel["divergence_from_Melopsittacus_myr"] <= 20).sum())

    d[["species", "tree_tip", "family", "genus", "log_life_expectancy", "log_body_mass", "rel_log_life_expectancy",
       "genome_best_tier"]].to_csv(OUT / "phenotypes_adjusted.tsv", sep="\t", index=False, float_format="%.5g")
    (OUT / "phase0_summary.json").write_text(json.dumps(summary, indent=2, default=float))
    print(json.dumps(summary, indent=2, default=float))


if __name__ == "__main__":
    main()
