"""Generate publication-grade figures for Phase 1 (Audit), Phase 2 (H8 Selection), and Phase 3 (Genome-wide RERconverge).

Outputs:
  results/phase1_audit/figS1_audit_126_assemblies.png
  results/phase2_genes/fig5_phase2_selection_h8.png
  results/phase3_genome/fig6_phase3_genome_rer.png
"""
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]

plt.rcParams.update({
    "font.sans-serif": "DejaVu Sans",
    "font.size": 9.5,
    "axes.titlesize": 10.5,
    "axes.titleweight": "bold",
    "axes.labelsize": 9.5,
    "xtick.labelsize": 8.5,
    "ytick.labelsize": 8.5,
    "legend.fontsize": 8.0,
    "figure.dpi": 300,
})


def make_fig_s1():
    at = pd.read_csv(ROOT / "results/phase1_audit/audit_table.tsv", sep="\t")
    sg = pd.read_csv(ROOT / "results/phase1_audit/species_genomes.tsv", sep="\t")
    m = pd.read_csv(ROOT / "results/phase0b/species_data_matrix_with_residuals.tsv", sep="\t")
    at_m = at.merge(m[["species", "rel_log_life_expectancy"]], on="species", how="left")

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))

    # Panel A: Contig N50 vs BUSCO Completeness
    ax = axes[0]
    inc_acc = set(sg[sg.status == "included"]["accession"].dropna())
    is_inc = at["accession"].isin(inc_acc)
    ax.scatter(at.loc[~is_inc, "contig_n50"] / 1e3, at.loc[~is_inc, "busco_C"],
               c="#c0392b", alpha=0.65, s=32, label=f"Iridian excluída (n={(~is_inc).sum()})")
    ax.scatter(at.loc[is_inc, "contig_n50"] / 1e3, at.loc[is_inc, "busco_C"],
               c="#2980b9", alpha=0.80, s=36, label=f"Iridian incluída (n={is_inc.sum()})")
    other_inc = sg[(sg.status == "included") & (sg.source != "Iridian")]
    ax.scatter(other_inc["contig_n50"] / 1e3, other_inc["busco_C"],
               c="#27ae60", marker="^", alpha=0.90, s=48,
               label=f"Outra fonte incluída (VGP/B10K/Sanger, n={len(other_inc)})")
    ax.axhline(70.0, color="#2c3e50", linestyle="--", linewidth=1.1, label="Limiar BUSCO ≥ 70%")
    ax.axvline(5.0, color="#7f8c8d", linestyle=":", linewidth=1.0, label="N50 = 5 kb")
    ax.set_xscale("log")
    ax.set_xlabel("Contig N50 (kb, escala log)")
    ax.set_ylabel("BUSCO Completo aves_odb10 (%)")
    ax.set_title("A. Auditoria de 126 montagens (71 spp. incluídas)")
    ax.legend(loc="lower right", frameon=True)
    ax.grid(alpha=0.25)

    # Panel B: Mitochondrial & Nuclear Identity Breakdown (96 Iridian)
    ax = axes[1]
    cats = [
        "Mitogenoma próprio\nisolado (n=22)",
        "Próprio + mitogenoma\nestranho (n=39)",
        "Conflito mito resolvido\npor Mash nuclear (n=21)",
        "Sem mito decisivo\n(Mash clado/fragmentado, n=10)",
        "Reprovadas: quimera ou\nespécime errado (n=4)",
    ]
    counts = [22, 39, 21, 10, 4]
    colors = ["#27ae60", "#f39c12", "#2980b9", "#95a5a6", "#c0392b"]
    y_pos = np.arange(len(cats))
    bars = ax.barh(y_pos, counts, color=colors, edgecolor="#2c3e50", height=0.62)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(cats)
    ax.invert_yaxis()
    ax.set_xlabel("Número de montagens Iridian (total = 96)")
    ax.set_title("B. Auditoria de identidade e contaminação")
    for b, c in zip(bars, counts):
        ax.text(b.get_width() + 0.8, b.get_y() + b.get_height() / 2, f"{c} ({c/96*100:.0f}%)",
                va="center", fontsize=8.5, fontweight="bold")
    ax.set_xlim(0, 48)
    ax.grid(axis="x", alpha=0.25)

    # Panel C: Quality vs Relative Longevity
    ax = axes[2]
    v = at_m.dropna(subset=["busco_C", "rel_log_life_expectancy"])
    ax.scatter(v["busco_C"], v["rel_log_life_expectancy"], c="#34495e", alpha=0.7, s=34)
    slope, intercept = np.polyfit(v["busco_C"], v["rel_log_life_expectancy"], 1)
    xs = np.linspace(v["busco_C"].min(), v["busco_C"].max(), 100)
    ax.plot(xs, slope * xs + intercept, color="#e74c3c", linestyle="--", linewidth=1.4,
            label="PGLS λ=0,51: β = -0,0015 (p = 0,11)")
    ax.axhline(0, color="#7f8c8d", linestyle=":", linewidth=0.9)
    ax.set_xlabel("BUSCO Completo (%)")
    ax.set_ylabel("Longevidade Relativa (resíduo PGLS log LE)")
    ax.set_title("C. Ausência de viés técnico qualidade × longevidade")
    ax.legend(loc="upper right", frameon=True)
    ax.grid(alpha=0.25)

    plt.tight_layout()
    out_path = ROOT / "results/phase1_audit/figS1_audit_126_assemblies.png"
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print("Saved:", out_path)


def make_fig_5():
    rel = pd.read_csv(ROOT / "results/phase2_genes/h8_relax_by_gene.tsv", sep="\t")

    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.8))

    # Panel A: Module-level log2(k) on Burgio 2019
    ax = axes[0]
    b_rel = rel[rel.tree == "burgio2019"].copy()
    b_rel["log2_k_clip"] = np.log2(b_rel["k"].clip(0.05, 20.0))
    mods = [f"M{i}" for i in range(1, 9)]
    data = [b_rel[b_rel.module == m]["log2_k_clip"].dropna().values for m in mods]
    bp = ax.boxplot(data, tick_labels=mods, patch_artist=True, showfliers=False, widths=0.52)
    for patch in bp["boxes"]:
        patch.set_facecolor("#d6eaf8")
        patch.set_edgecolor("#2980b9")
    for i, m in enumerate(mods, start=1):
        vals = b_rel[b_rel.module == m]["log2_k_clip"].dropna().values
        jitter = np.random.default_rng(42 + i).normal(0, 0.07, size=len(vals))
        ax.scatter(np.full_like(vals, i) + jitter, vals, color="#1f618d", alpha=0.75, s=24, zorder=3)
    ax.axhline(0, color="#c0392b", linestyle="--", linewidth=1.2, label="k = 1 (neutro / sem mudança)")
    ax.set_xlabel("Módulos funcionais pré-registrados (M1–M8)")
    ax.set_ylabel("RELAX log₂(k) em Psittaculidae (truncado em [0,05; 20])")
    ax.set_title("A. RELAX por módulo (nulo 298 genes: todos p_Bonf > 0,38)")
    ax.legend(loc="lower left", frameon=True)
    ax.grid(axis="y", alpha=0.25)

    # Panel B: Concordant RELAX genes across Burgio 2019 vs Smith 2024
    ax = axes[1]
    p_b = rel[rel.tree == "burgio2019"].set_index("gene")
    p_s = rel[rel.tree == "smith2024"].set_index("gene")
    both = p_b.join(p_s, lsuffix="_b", rsuffix="_s").dropna(subset=["k_b", "k_s", "relax_q_b", "relax_q_s"])
    both["log2_kb"] = np.log2(both["k_b"].clip(0.05, 20.0))
    both["log2_ks"] = np.log2(both["k_s"].clip(0.05, 20.0))
    sig_both = (both["relax_q_b"] < 0.10) & (both["relax_q_s"] < 0.10) & (both["k_b"] > 0.01) & (both["k_s"] > 0.01)
    ax.scatter(both.loc[~sig_both, "log2_kb"], both.loc[~sig_both, "log2_ks"],
               c="#95a5a6", alpha=0.65, s=28, label="Demais genes do painel (n=89)")
    ax.scatter(both.loc[sig_both, "log2_kb"], both.loc[sig_both, "log2_ks"],
               c="#e74c3c", alpha=0.95, s=55, zorder=4, label="q < 0,10 nas duas árvores (k não-fronteira)")
    for g, r in both[sig_both].iterrows():
        ax.annotate(g, (r["log2_kb"], r["log2_ks"]), xytext=(5, 4), textcoords="offset points",
                    fontsize=8.5, fontweight="bold", color="#922b21")
    ax.axhline(0, color="#2c3e50", linestyle=":", linewidth=0.9)
    ax.axvline(0, color="#2c3e50", linestyle=":", linewidth=0.9)
    ax.set_xlabel("RELAX log₂(k) — Supertree (Burgio 2019)")
    ax.set_ylabel("RELAX log₂(k) — Árvore UCE (Smith 2024)")
    ax.set_title("B. Concordância RELAX entre árvores (TERT, MRE11, SOD1, DKC1, NQO1)")
    ax.legend(loc="upper left", frameon=True)
    ax.grid(alpha=0.25)

    # Panel C: aBSREL on Melopsittacus terminal branch
    ax = axes[2]
    ab_bud = b_rel.dropna(subset=["absrel_melopsittacus_q"]).copy()
    ab_bud["neg_log10_q"] = -np.log10(ab_bud["absrel_melopsittacus_q"].clip(1e-10, 1.0))
    top_ab = ab_bud.sort_values("absrel_melopsittacus_q").head(10).iloc[::-1]
    colors = ["#c0392b" if q < 0.05 else "#7f8c8d" for q in top_ab["absrel_melopsittacus_q"]]
    ax.barh(top_ab["gene"] + " (" + top_ab["module"] + ")", top_ab["neg_log10_q"],
            color=colors, edgecolor="#2c3e50", height=0.62)
    ax.axvline(-np.log10(0.05), color="#2c3e50", linestyle="--", linewidth=1.1, label="FDR q = 0,05")
    ax.set_xlabel("-log₁₀(q) aBSREL no ramo terminal de M. undulatus")
    ax.set_title("C. Seleção episódica no complexo telomérico (M7) e MYC")
    ax.legend(loc="lower right", frameon=True)
    ax.grid(axis="x", alpha=0.25)

    plt.tight_layout()
    out_path = ROOT / "results/phase2_genes/fig5_phase2_selection_h8.png"
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print("Saved:", out_path)


def make_fig_6():
    le = pd.read_csv(ROOT / "results/phase3_genome/rer_burgio2019_rel_log_life_expectancy.tsv", sep="\t").set_index("gene")
    br = pd.read_csv(ROOT / "results/phase3_genome/rer_burgio2019_rel_log_brain.tsv", sep="\t").set_index("gene")
    bg = pd.read_csv(ROOT / "results/phase3_genome/rer_burgio2019_rel_brain_given_le.tsv", sep="\t").set_index("gene")
    m = le.join(br, lsuffix="_le", rsuffix="_brain").join(bg.add_suffix("_bg")).dropna(subset=["Rho_le", "Rho_brain"])

    fig, axes = plt.subplots(1, 3, figsize=(16.2, 4.9))

    # Panel A: Genome-wide Rho_le vs Rho_brain (7,756 genes)
    ax = axes[0]
    shared = (
        (m["P_le"] < 0.05) & (m["permP_le"] < 0.05) &
        (m["P_brain"] < 0.05) & (m["permP_brain"] < 0.05) &
        (np.sign(m["Rho_le"]) == np.sign(m["Rho_brain"]))
    )
    ax.scatter(m.loc[~shared, "Rho_le"], m.loc[~shared, "Rho_brain"],
               c="#bdc3c7", alpha=0.35, s=10, rasterized=True, label=f"Universo elegível (n={len(m)})")
    pos_sh = shared & (m["Rho_le"] > 0)
    neg_sh = shared & (m["Rho_le"] < 0)
    ax.scatter(m.loc[pos_sh, "Rho_le"], m.loc[pos_sh, "Rho_brain"],
               c="#e67e22", alpha=0.85, s=26, label=f"Acelerados em ambos (n={pos_sh.sum()})")
    ax.scatter(m.loc[neg_sh, "Rho_le"], m.loc[neg_sh, "Rho_brain"],
               c="#2980b9", alpha=0.85, s=26, label=f"Conservados em ambos (n={neg_sh.sum()})")

    highlight = ["PTEN", "CTNNB1", "PSMA4", "SSU72", "MAD2L2", "HSPA8", "SOX14", "RTN1", "NCALD", "YWHAE"]
    for g in highlight:
        if g in m.index:
            r = m.loc[g]
            ax.scatter([r["Rho_le"]], [r["Rho_brain"]], c="#c0392b", s=46, zorder=5)
            ax.annotate(g, (r["Rho_le"], r["Rho_brain"]), xytext=(4, 3), textcoords="offset points",
                        fontsize=8.0, fontweight="bold", color="#1a252f")
    ax.axhline(0, color="#7f8c8d", linestyle=":", linewidth=0.8)
    ax.axvline(0, color="#7f8c8d", linestyle=":", linewidth=0.8)
    r_val, _ = stats.pearsonr(m["Rho_le"], m["Rho_brain"])
    ax.set_xlabel("RERconverge Rho — Longevidade Relativa")
    ax.set_ylabel("RERconverge Rho — Tamanho Cerebral Relativo")
    ax.set_title(f"A. Co-evolução genômica Longevidade × Cérebro (r = {r_val:.3f}, 7.756 genes)")
    ax.legend(loc="upper left", frameon=True)
    ax.grid(alpha=0.25)

    # Panel B: Exploratory Functional Axes
    ax = axes[1]
    ax_sum = pd.read_csv(ROOT / "results/phase3_genome/exploratory_axes_genes.tsv", sep="\t")
    ax_b = ax_sum[ax_sum.tree == "burgio2019"]
    axis_order = [
        ("Neuroprotection_Proteostasis", "Neuroproteção &\nProteostase (n=22)"),
        ("Antioxidation_Redox", "Antioxidação &\nRedox (n=22)"),
        ("ExtremeLongevity_NutrientSensing", "Longevidade Extrema &\nmTOR/IIS/Sirtuínas (n=9)"),
        ("VocalLearning_Language_Plasticity", "Aprendizado Vocal &\nPlasticidade (n=13)"),
        ("AntiCancer_DNARepair_Telomere", "Anticâncer, Reparo DNA\n& Telômeros (n=24)"),
    ]
    y_pos = np.arange(len(axis_order))
    le_means = [ax_b[ax_b.axis == k]["Rho_le"].mean() for k, _ in axis_order]
    br_means = [ax_b[ax_b.axis == k]["Rho_brain"].mean() for k, _ in axis_order]
    le_se = [ax_b[ax_b.axis == k]["Rho_le"].sem() for k, _ in axis_order]
    br_se = [ax_b[ax_b.axis == k]["Rho_brain"].sem() for k, _ in axis_order]

    w = 0.34
    ax.barh(y_pos - w/2, le_means, height=w, xerr=le_se, color="#2980b9", alpha=0.85, label="Longevidade Relativa")
    ax.barh(y_pos + w/2, br_means, height=w, xerr=br_se, color="#e67e22", alpha=0.85, label="Cérebro Relativo")
    ax.axvline(0, color="#2c3e50", linestyle="--", linewidth=1.0)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([lbl for _, lbl in axis_order])
    ax.invert_yaxis()
    ax.set_xlabel("Média de RER Rho (± EP)")
    ax.set_title("B. Eixos funcionais exploratórios (Neuroproteção p = 0,027)")
    ax.legend(loc="lower right", frameon=True)
    ax.grid(axis="x", alpha=0.25)

    # Panel C: Top individual candidate genes across the 3 traits
    ax = axes[2]
    cands = ["MAD2L2", "RTN1", "SOX14", "NCALD", "HSPA8", "PSMA2", "SOD1", "TERF2IP", "RAD50", "PTEN", "PSMA4", "CTNNB1"]
    cands = [g for g in cands if g in m.index]
    sub = m.loc[cands]
    y = np.arange(len(cands))
    ax.scatter(sub["Rho_le"], y - 0.18, color="#2980b9", s=42, label="Longevidade (rel_LE)")
    ax.scatter(sub["Rho_brain"], y, color="#e67e22", s=42, marker="s", label="Cérebro (rel_Brain)")
    ax.scatter(sub["Rho_bg"], y + 0.18, color="#27ae60", s=42, marker="^", label="Cérebro | Longevidade")
    for idx_y, g in enumerate(cands):
        r = sub.loc[g]
        ax.plot([min(r["Rho_le"], r["Rho_brain"], r["Rho_bg"]),
                 max(r["Rho_le"], r["Rho_brain"], r["Rho_bg"])],
                [idx_y, idx_y], color="#7f8c8d", alpha=0.4, linewidth=1.2)
    ax.axvline(0, color="#c0392b", linestyle="--", linewidth=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels(cands, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel("Correlação de taxa evolutiva (RER Rho)")
    ax.set_title("C. Genes-chave compartilhados (Conservados < 0 vs. Acelerados > 0)")
    ax.legend(loc="lower right", frameon=True)
    ax.grid(axis="x", alpha=0.25)

    plt.tight_layout()
    out_path = ROOT / "results/phase3_genome/fig6_phase3_genome_rer.png"
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print("Saved:", out_path)


if __name__ == "__main__":
    make_fig_s1()
    make_fig_5()
    make_fig_6()
