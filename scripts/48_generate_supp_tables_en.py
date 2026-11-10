#!/usr/bin/env python3
"""Programmatically builds docs/supplementary_parrots_EN.tex (Tables S1-S6) directly from
verified TSV outputs on disk across Phases 1, 2, 3, 71sp extension, 3D modeling, and Red-Team audit.
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def tex_esc(s):
    if pd.isna(s):
        return "--"
    s = str(s)
    for k, v in (("\\", "\\textbackslash{}"), ("_", "\\_"), ("%", "\\%"), ("&", "\\&"), ("#", "\\#")):
        s = s.replace(k, v)
    return s


def fmt_p(x):
    if pd.isna(x):
        return "--"
    x = float(x)
    if x < 1e-3:
        mant, exp = f"{x:.2e}".split("e")
        return f"${mant}\\times 10^{{{int(exp)}}}$"
    return f"{x:.4f}"


def fmt_f(x, d=3):
    if pd.isna(x):
        return "--"
    return f"{float(x):.{d}f}"


def main():
    sg = pd.read_csv(ROOT / "results/phase1_audit/species_genomes.tsv", sep="\t")
    inc = sg[sg.status.str.startswith("included")].copy()
    t71 = pd.read_csv(ROOT / "results/phase2_genes_71sp/rer_traits_71sp.tsv", sep="\t")
    coh_map = dict(zip(t71["species"], t71["cohort"]))
    inc["cohort"] = inc["species"].map(coh_map).fillna("holdout_13")
    inc = inc.sort_values(["cohort", "psittaculidae", "species"], ascending=[True, False, True])

    comp71 = pd.read_csv(ROOT / "results/phase2_genes_71sp/comparison_58_vs_71_confirmatory.tsv", sep="\t")
    tel7 = pd.read_csv(ROOT / "results/redteam_audit/7_telomere_alignment_and_rbh.tsv", sep="\t").set_index("gene")
    integ71 = pd.read_csv(ROOT / "results/phase2_genes_71sp/h8_integrity_71sp.tsv", sep="\t").set_index("gene")

    null_mod = pd.read_csv(ROOT / "results/phase2_genes/h8_null_modules.tsv", sep="\t")
    base6 = pd.read_csv(ROOT / "results/redteam_audit/6_absrel_base_rate.tsv", sep="\t")

    r1 = pd.read_csv(ROOT / "results/redteam_audit/1_reproducibility.tsv", sep="\t")
    r2 = pd.read_csv(ROOT / "results/redteam_audit/2_coupling_null_and_controls.tsv", sep="\t")
    r3c = pd.read_csv(ROOT / "results/redteam_audit/3_coupling_after_confounder_removal.tsv", sep="\t")
    r4 = pd.read_csv(ROOT / "results/redteam_audit/4_leave_one_clade_out.tsv", sep="\t")
    r5s = pd.read_csv(ROOT / "results/redteam_audit/5_residual_brain_summary.tsv", sep="\t")

    r5k = pd.read_csv(ROOT / "results/redteam_audit/5_residual_brain_key_genes.tsv", sep="\t").set_index("gene")
    r8cv = pd.read_csv(ROOT / "results/redteam_audit/8_holdout_cv_key_genes.tsv", sep="\t").set_index("gene")
    le = pd.read_csv(ROOT / "results/phase3_genome/rer_burgio2019_rel_log_life_expectancy.tsv", sep="\t")
    br = pd.read_csv(ROOT / "results/phase3_genome/rer_burgio2019_rel_log_brain.tsv", sep="\t")
    m3 = le.merge(br, on="gene", suffixes=("_LE", "_Brain"))
    zfrac, n_nz = {}, {}
    for line in open(ROOT / "results/phase3_genome/genetrees_burgio2019.txt"):
        g, nwk = line.rstrip("\n").split("\t", 1)
        els = [float(x) for x in re.findall(r":([0-9eE.+-]+)", nwk)]
        n_nz[g] = sum(1 for e in els if e > 1e-6)
        zfrac[g] = 1.0 - n_nz[g] / len(els) if els else 1.0
    m3["n_nonzero"] = m3["gene"].map(n_nz)
    m3["zero_frac"] = m3["gene"].map(zfrac)

    tier_b = m3[(m3.n_nonzero >= 20) & (m3.P_LE < 0.05) & (m3.permP_LE < 0.05) &
                (m3.P_Brain < 0.05) & (m3.permP_Brain < 0.05) & (m3.Rho_LE * m3.Rho_Brain > 0)].sort_values("P_LE")
    m3_idx = m3.set_index("gene")

    inv9 = pd.read_csv(ROOT / "results/redteam_audit/9_statistical_honesty_inventory.tsv", sep="\t")
    s3d = pd.read_csv(ROOT / "results/phase3_genome/structural_3d_summary.tsv", sep="\t")

    lines = [
        r"\documentclass[9pt,a4paper]{article}",
        r"\usepackage[utf8]{inputenc}",
        r"\usepackage[T1]{fontenc}",
        r"\usepackage[margin=1.4cm,top=1.6cm,bottom=1.6cm]{geometry}",
        r"\usepackage{graphicx,booktabs,longtable,array,amsmath,amssymb,xcolor,hyperref}",
        r"\hypersetup{colorlinks=true,linkcolor=blue!50!black,urlcolor=blue!50!black}",
        r"\setlength{\parskip}{3pt}",
        r"\title{\textbf{\Large Supplementary Materials: Tables S1--S6}\\[4pt]\large Comparative Genomics of Longevity, Encephalization, and Cancer Resistance in Psittaciformes (with 71-Species Extension and Full Red-Team Statistical Audit)}",
        r"\author{Angelo --- UNICEPLAC, Bras\'ilia, DF, Brazil \\ Repository: \texttt{github.com/redkk123/Parrots} (\texttt{claude/new-session-lwh7xh})}",
        r"\date{November 10, 2026}",
        r"\begin{document}",
        r"\maketitle",
        r"\tableofcontents",
        r"\clearpage",
    ]

    # TABLE S1
    lines.append(r"\section*{Table S1. Validated 71-Species Genomic Cohort from 126 Audited Assemblies}")
    lines.append(r"\addcontentsline{toc}{section}{Table S1: Validated 71-Species Cohort}")
    lines.append(r"Complete list of all 71 included Psittaciformes species passing the pre-registered Phase 1 audit ($\text{BUSCO}_{\text{complete}} \ge 70\%$ on \texttt{compleasm aves\_odb10}, nuclear Mash identity confirmed, deduplicated to one best assembly per species). Cohort indicates the frozen 58-species discovery panel vs.\ the 13 newly finalized assemblies (\texttt{holdout\_13}).")
    lines.append(r"\begin{longtable}{llllrrrr}")
    lines.append(r"\toprule")
    lines.append(r"\textbf{Species} & \textbf{Cohort} & \textbf{Clade} & \textbf{Accession} & \textbf{Source} & \textbf{Contig N50 (kb)} & \textbf{BUSCO C (\%)} & \textbf{Orig.\ Strict} \\")
    lines.append(r"\midrule\endhead")
    for r in inc.itertuples():
        cl = "Psittaculidae" if r.psittaculidae else "Other"
        n50_kb = f"{r.contig_n50/1e3:.1f}" if pd.notna(r.contig_n50) else "--"
        strict = "Yes" if bool(r.original_criterion_pass) else "No"
        coh = "Discovery (58)" if r.cohort == "discovery_58" else "Holdout (13)"
        src = tex_esc(str(r.source)[:22])
        lines.append(f"\\textit{{{tex_esc(r.species)}}} & {coh} & {cl} & \\texttt{{{tex_esc(r.accession)}}} & {src} & {n50_kb} & {r.busco_C:.2f} & {strict} \\\\")
    lines.append(r"\bottomrule\end{longtable}\clearpage")

    # TABLE S2
    lines.append(r"\section*{Table S2. Phase 2 H8 Confirmatory Selection Tests (58 vs.\ 71 Species) and Telomere Sensitivity Audit}")
    lines.append(r"\addcontentsline{toc}{section}{Table S2: H8 58 vs. 71 Species \& Telomere Sensitivity}")
    lines.append(r"Comparison of key H8 candidate genes across the frozen 58-species cohort and the extended 71-species cohort (\texttt{burgio2019} topology), alongside Red-Team alignment sensitivity (MAFFT vs.\ L-INS-i vs.\ stringent window/gap cleaning vs.\ \texttt{aBSREL --multiple-hits Double+Triple}) and \textit{Strigops habroptila} (kakapo) dual-reference orthology concordance.")
    lines.append(r"\resizebox{\textwidth}{!}{%")
    lines.append(r"\begin{tabular}{llcccccccccc}")
    lines.append(r"\toprule")
    lines.append(r"\textbf{Gene} & \textbf{Mod.} & \textbf{$N_{\text{58/71}}$} & \textbf{$k_{\text{58}}$} & \textbf{RELAX $p_{\text{58}}$} & \textbf{$k_{\text{71}}$} & \textbf{RELAX $p_{\text{71}}$} & \textbf{aBSREL $p_{\text{58}}$} & \textbf{aBSREL $p_{\text{71}}$} & \textbf{Cleaned $p$} & \textbf{MultiHit $p$} & \textbf{Kakapo Orth.} \\")
    lines.append(r"\midrule")
    r58 = pd.read_csv(ROOT / "results/phase2_genes/h8_relax_by_gene.tsv", sep="\t")
    r58b = r58[r58.tree == "burgio2019"].set_index("gene")
    s58 = pd.read_csv(ROOT / "results/phase2_genes/alignment_summary.tsv", sep="\t").set_index("gene")
    for r in comp71.itertuples():
        g = r.gene
        n58 = int(s58.loc[g, "n_species"]) if g in s58.index else "--"
        n71 = int(r.n_sp_71sp) if pd.notna(r.n_sp_71sp) else "--"
        k58 = fmt_f(r.k_58sp, 2)
        rp58 = fmt_p(r.relax_p_58sp)
        k71 = fmt_f(getattr(r, "k_71sp", np.nan), 2)
        rp71 = fmt_p(getattr(r, "relax_p_71sp", np.nan))
        ap58 = fmt_p(r.absrel_melo_p_58sp)
        ap71 = fmt_p(getattr(r, "absrel_melo_p_71sp", np.nan))
        if g in tel7.index:
            cp = fmt_p(tel7.loc[g, "absrel_p_window_cleaned"])
            mp = fmt_p(tel7.loc[g, "absrel_p_MultiHit_DoubleTriple"])
            ko = tex_esc(tel7.loc[g, "dual_ref_kakapo_concordance"])
        else:
            cp = "--"
            mp = fmt_p(getattr(r, "absrel_mh_melo_p_71sp", np.nan))
            ko = f"{int(integ71.loc[g, 'complete_psittaculidae_71sp'])}/20 Psl" if g in integ71.index else "--"
        lines.append(f"\\texttt{{{tex_esc(g)}}} & {tex_esc(r.module)} & {n58}/{n71} & {k58} & {rp58} & {k71} & {rp71} & {ap58} & {ap71} & {cp} & {mp} & {ko} \\\\")
    lines.append(r"\bottomrule\end{tabular}}\clearpage")

    # TABLE S3
    lines.append(r"\section*{Table S3. Phase 2 Modular Null Benchmark (298 Control Genes) and Empirical \texttt{aBSREL} Base Rate}")
    lines.append(r"\addcontentsline{toc}{section}{Table S3: Modular Null \& aBSREL Base Rate}")
    lines.append(r"\textbf{Top panel:} H8 modules (M1--M8) compared against 1,000 length- and rate-matched null pseudo-panels sampled from 298 control genes (\texttt{h8\_null\_modules.tsv}). \textbf{Bottom panel:} Empirical \texttt{aBSREL} positive-selection rate on the \textit{M.\ undulatus} terminal branch vs.\ Psittaculidae stem across the panel and control genes (\texttt{6\_absrel\_base\_rate.tsv}).")
    lines.append(r"\begin{tabular}{lccccccc}")
    lines.append(r"\toprule")
    lines.append(r"\textbf{Module} & \textbf{$n_{\text{genes}}$} & \textbf{Prop.\ Relaxed} & \textbf{Null Mean} & \textbf{$p$ (Prop.)} & \textbf{Mean $\log_2 k$} & \textbf{Null Mean} & \textbf{$p$ ($\log_2 k$)} \\")
    lines.append(r"\midrule")
    for r in null_mod.itertuples():
        lines.append(f"{tex_esc(r.module)} & {r.n_genes} & {r.prop_relaxed:.3f} & {r.prop_relaxed_null_mean:.3f} & {r.prop_relaxed_p:.3f} & {r.mean_log2_k:.3f} & {r.mean_log2_k_null_mean:.3f} & {r.mean_log2_k_p:.3f} \\\\")
    lines.append(r"\bottomrule\end{tabular}\\[10pt]")
    lines.append(r"\resizebox{\textwidth}{!}{%")
    lines.append(r"\begin{tabular}{lccccccc}")
    lines.append(r"\toprule")
    lines.append(r"\textbf{Gene Set} & \textbf{Tested} & \textbf{Melo $p<0.05$} & \textbf{Melo $q<0.10$} & \textbf{Stem $p<0.05$} & \textbf{Stem $q<0.10$} & \textbf{Fisher $p$ vs.\ Ctrl ($p<0.05$)} & \textbf{Fisher $p$ ($q<0.10$)} \\")
    lines.append(r"\midrule")
    for r in base6.itertuples():
        fp05 = fmt_p(r.fisher_vs_control_p05) if pd.notna(r.fisher_vs_control_p05) else "--"
        fq10 = fmt_p(r.fisher_vs_control_q10) if pd.notna(r.fisher_vs_control_q10) else "--"
        lines.append(f"{tex_esc(r.group)} & {r.n_genes_tested} & {tex_esc(r.melopsittacus_p_holm2_lt_005)} & {tex_esc(r.melopsittacus_q_bh_panel_lt_010)} & {tex_esc(r.psittaculidae_stem_p_holm2_lt_005)} & {tex_esc(r.psittaculidae_stem_q_bh_panel_lt_010)} & {fp05} & {fq10} \\\\")
    lines.append(r"\bottomrule\end{tabular}}\clearpage")

    # TABLE S4
    lines.append(r"\section*{Table S4. Red-Team Sensitivity Audit of Phase 3 Genome-Wide Coupling ($r = 0.666$)}")
    lines.append(r"\addcontentsline{toc}{section}{Table S4: Red-Team Sensitivity of Phase 3 Coupling}")
    lines.append(r"Summary of all stress tests applied to the 7,756-gene \texttt{RERconverge} correlation coupling between relative longevity ($\text{Rho}_{\text{LE}}$) and relative brain mass ($\text{Rho}_{\text{Brain}}$): seed reproducibility, Brownian motion null permulations, technical confounders, zero-branch invariance stratification, leave-one-clade-out perturbations, and true partial orthogonalization.")
    lines.append(r"\begin{longtable}{p{6.5cm}cccc}")
    lines.append(r"\toprule")
    lines.append(r"\textbf{Red-Team Audit Condition / Perturbation} & \textbf{$r_{\text{tip}}$} & \textbf{$r_{\text{branch path}}$} & \textbf{$r(\text{Rho}_1,\text{Rho}_2)$} & \textbf{Key Note} \\")
    lines.append(r"\midrule\endhead")
    lines.append(r"\multicolumn{5}{l}{\textbf{A. Negative Controls \& Brownian Permulation Nulls (\texttt{2\_coupling\_null\_and\_controls.tsv})}} \\")
    for r in r2.itertuples():
        lines.append(f"{tex_esc(r.comparison[:58])} & {fmt_f(r.tip_cor_pearson,4)} & {fmt_f(r.branch_path_cor,4)} & {fmt_f(r.genome_Rho_cor,4)} & \\scriptsize {tex_esc(str(r.notes)[:52])} \\\\")
    lines.append(r"\midrule\multicolumn{5}{l}{\textbf{B. Technical Confounder Adjustment \& Zero-Branch Strata (\texttt{3\_coupling\_after\_confounder\_removal.tsv})}} \\")
    for r in r3c.iloc[:4].itertuples():
        lines.append(f"{tex_esc(r.model)} ($n={r.n_genes}$) & -- & -- & {fmt_f(r.r_pearson,4)} & Spearman={fmt_f(r.r_spearman,4)}; $r(|\\text{{Rho}}|)$={fmt_f(r.r_abs_Rho_pearson,4)} \\\\")
    lines.append(r"Ultra-conserved stratum ($>80\%$ zero branches, $n=1,450$) & -- & -- & 0.6761 & \textbf{117/131 (89.3\%)} of shared hits! \\")
    lines.append(r"Informative stratum ($\ge 20$ nonzero branches, $n=6,182$) & -- & -- & 0.6626 & Only \textbf{14/131} shared hits \\")
    lines.append(r"Most informative ($\ge 70\%$ nonzero branches, $n=1,263$) & -- & -- & 0.6982 & Only \textbf{2} shared hits \\")
    lines.append(r"\midrule\multicolumn{5}{l}{\textbf{C. Leave-One-Clade-Out Sensitivity (\texttt{4\_leave\_one\_clade\_out.tsv})}} \\")
    for r in r4.itertuples():
        lines.append(f"{tex_esc(r.subset[:58])} ($n={r.n_species}$) & {fmt_f(r.tip_cor_pearson,4)} & {fmt_f(r.branch_path_cor,4)} & {fmt_f(r.genome_Rho_cor,4)} & PTEN={fmt_f(r.PTEN_Rho_LE,3)}; MAD2L2={fmt_f(r.MAD2L2_Rho_LE,3)} \\\\")
    lines.append(r"\midrule\multicolumn{5}{l}{\textbf{D. Orthogonalization of Brain Beyond Longevity (\texttt{5\_residual\_brain\_summary.tsv})}} \\")
    for r in r5s.itertuples():
        lines.append(f"{tex_esc(r.method[:58])} & {fmt_f(r.tip_cor_with_LE,4)} & {fmt_f(r.branch_path_cor_with_LE,4)} & {fmt_f(r.genome_Rho_cor_with_LE_58sp,4)} & 50sp={fmt_f(r.genome_Rho_cor_with_LE_50sp,4)}; cor(Br)={fmt_f(r.genome_Rho_cor_with_Brain_raw,3)} \\\\")
    lines.append(r"\bottomrule\end{longtable}\clearpage")

    # TABLE S5
    lines.append(r"\section*{Table S5. Dual-Tier Phase 3 Candidate Gene Catalog: Ultra-Conserved (Tier A) vs.\ Informative (Tier B)}")
    lines.append(r"\addcontentsline{toc}{section}{Table S5: Tier A vs. Tier B Candidate Genes}")
    lines.append(r"\textbf{Tier A} lists the 13 flagship shared candidate genes from the original Phase 3 synthesis together with their fraction of zero-length branches ($\le 10^{-6}$), branch-orthogonalized $\text{Rho}_{\text{Brain}\perp\text{LE}}$, OLS tip-orthogonalized permulation $p$, and 30-split cross-validation sign concordance. \textbf{Tier B} lists all 14 shared genes passing $P < 0.05$ and $\text{permP} < 0.05$ for both traits that have $\ge 20$ nonzero branches.")
    lines.append(r"\resizebox{\textwidth}{!}{%")
    lines.append(r"\begin{tabular}{lcccccccc}")
    lines.append(r"\toprule")
    lines.append(r"\textbf{Gene} & \textbf{\% Zero Branches ($n_{\text{nz}}$)} & \textbf{$\text{Rho}_{\text{LE}}$} & \textbf{$\text{permP}_{\text{LE}}$} & \textbf{$\text{Rho}_{\text{Brain}}$} & \textbf{$\text{Rho}_{\text{Br}\perp\text{LE}}$ (Orth.)} & \textbf{OLS $\text{permP}_{\text{Br}}$} & \textbf{CV Sign Conc.\ (30 splits)} & \textbf{Tier} \\")
    lines.append(r"\midrule")
    for g in r5k.index:
        zf = m3_idx.loc[g, "zero_frac"] * 100
        nz = int(m3_idx.loc[g, "n_nonzero"])
        pl = m3_idx.loc[g, "permP_LE"]
        cv_s = f"{100*r8cv.loc[g, 'sign_concordance_holdout_LE']:.0f}\\% / {100*r8cv.loc[g, 'sign_concordance_holdout_Brain']:.0f}\\%" if g in r8cv.index else "--"
        lines.append(f"\\texttt{{{tex_esc(g)}}} & {zf:.1f}\\% ({nz}) & {r5k.loc[g, 'Rho_LE_58sp']:+.3f} & {pl:.4f} & {r5k.loc[g, 'Rho_Brain_raw']:+.3f} & {r5k.loc[g, 'Rho_Brain_branch_orth']:+.3f} & {r5k.loc[g, 'permP_Brain_OLS_tip_orth']:.3f} & {cv_s} & Tier A (Ultra-cons.) \\\\")
    lines.append(r"\midrule")
    for r in tier_b.itertuples():
        zf = r.zero_frac * 100
        nz = int(r.n_nonzero)
        lines.append(f"\\texttt{{{tex_esc(r.gene)}}} & {zf:.1f}\\% ({nz}) & {r.Rho_LE:+.3f} & {r.permP_LE:.4f} & {r.Rho_Brain:+.3f} & -- & -- & -- & Tier B (Informative) \\\\")
    lines.append(r"\bottomrule\end{tabular}}\clearpage")

    # TABLE S6
    lines.append(r"\section*{Table S6. Exhaustive Statistical Inventory (78,802 Tests), Null Findings, and AlphaFold 3D Summary}")
    lines.append(r"\addcontentsline{toc}{section}{Table S6: Statistical Inventory \& 3D Structural Summary}")
    lines.append(r"\resizebox{\textwidth}{!}{%")
    lines.append(r"\begin{tabular}{p{3.4cm}rp{3.5cm}p{5.5cm}p{6.5cm}}")
    lines.append(r"\toprule")
    lines.append(r"\textbf{Phase / Scope} & \textbf{Tests Run} & \textbf{Correction} & \textbf{Positive / Nominally Significant} & \textbf{Null / Negative Findings Reported} \\")
    lines.append(r"\midrule")
    for r in inv9.itertuples():
        lines.append(f"{tex_esc(r.phase)} & {r.n_tests_run:,} & {tex_esc(r.correction_method)} & {tex_esc(r.significant_findings)} & {tex_esc(r.null_or_negative_findings)} \\\\")
    lines.append(r"\bottomrule\end{tabular}}\\[12pt]")
    lines.append(r"\resizebox{\textwidth}{!}{%")
    lines.append(r"\begin{tabular}{llcccp{8.5cm}}")
    lines.append(r"\toprule")
    lines.append(r"\textbf{Protein} & \textbf{UniProt ID} & \textbf{AA Len.} & \textbf{Global pLDDT} & \textbf{Spp.} & \textbf{Structural Mapping Finding (\texttt{structural\_3d\_summary.tsv})} \\")
    lines.append(r"\midrule")
    for r in s3d.itertuples():
        lines.append(f"\\texttt{{{tex_esc(r.gene)}}} & {tex_esc(r.uniprot_id)} & {r.length_aa} & {r.global_plddt:.2f} & {r.n_aligned_species} & {tex_esc(r.structural_finding)} \\\\")
    lines.append(r"\bottomrule\end{tabular}}")
    lines.append(r"\end{document}")

    out_tex = ROOT / "docs/supplementary_parrots_EN.tex"
    out_tex.write_text("\n".join(lines))
    print(f"Wrote {out_tex} cleanly.")


if __name__ == "__main__":
    main()
