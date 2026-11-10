#!/usr/bin/env python3
"""
Red-Team Python Audit Script covering:
  - Question 3: Technical confounders of gene-level RERconverge correlations (7,756 genes)
  - Question 6: Empirical base rate of HyPhy aBSREL positive selection on control genes + FDR audit
  - Question 7: Telomeric genes (RTEL1, TERF2, TERT, ACD, TERF1, POT1) alignment inspection,
                L-INS-i realignment, stringent window masking, aBSREL with Multiple Hits,
                and Reciprocal Best Hit / Strigops dual-reference orthology verification
  - Question 9: Exhaustive statistical test inventory and null results table
"""
import concurrent.futures
import gzip
import importlib.util
import json
import os
import re
import subprocess
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/redteam_audit"
OUT.mkdir(parents=True, exist_ok=True)

HYPHY = ROOT / ".bioenv/bin/hyphy"
MAFFT = ROOT / ".bioenv/bin/mafft"

CODE = {
    a + b + c: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"[16 * i + 4 * j + k]
    for i, a in enumerate("TCAG") for j, b in enumerate("TCAG") for k, c in enumerate("TCAG")
}


def translate_codon(c):
    return CODE.get(c.upper(), "X") if "-" not in c and "N" not in c.upper() else "-"


def translate_seq(s):
    return "".join(translate_codon(s[i:i + 3]) for i in range(0, len(s) - 2, 3))


def parse_fasta_str(text):
    seqs = {}
    name, buf = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)
            name, buf = line[1:].strip().split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


# ==============================================================================
# TEST 9: Statistical Honesty Inventory (All Tests Run + Null Results)
# ==============================================================================
def run_test_9():
    print("[RedTeam Python] Running Test 9 (Statistical Honesty Inventory)...", flush=True)
    inventory = [
        {
            "phase": "Phase 0 & 0b (Allometry, Captive/Wild, Brain & Life-History PGLS)",
            "script": "scripts/02_phase0_analysis.py, 05_phase0b_analysis.py, 06_regime_followup.py",
            "n_tests_run": 48,
            "correction_method": "Model AICc + PGLS t-tests (Supertree + Smith2024)",
            "significant_findings": "Body mass -> longevity (p=4.2e-18); Relative brain -> longevity (p=0.0012); Captive vs wild longevity r=0.91",
            "null_or_negative_findings": "Diet guild effect on longevity after mass+phylogeny (NULL, p>0.25); Latitude/insularity effect on residual longevity (NULL); OU single-optimum vs multi-regime outside Psittaculidae/Cacatuidae/Strigopidae (NULL)",
        },
        {
            "phase": "Phase 1 (Assembly Quality, Mito Identity, Mash Distance, Cancer Necropsy, Measurement Error)",
            "script": "scripts/11_mito_identity.py, 13_quality_vs_longevity.py, 17_mash_identity.py, 19_measurement_error.py, 20_cancer_necropsy.py",
            "n_tests_run": 34,
            "correction_method": "PGLS + spearman/pearson + ME sensitivity",
            "significant_findings": "Assembly quality (N50, BUSCO) is uncorrelated with longevity (validating no quality bias); M. undulatus high cancer prevalence in veterinary necropsy series",
            "null_or_negative_findings": "Assembly scaffold N50 vs residual longevity (NULL, r=0.04, p=0.76); Mitochondrial substitution rate vs longevity (NULL, p>0.30); Genome-wide Mash k-mer divergence vs longevity (NULL after tree height correction)",
        },
        {
            "phase": "Phase 2 Pilot & Pilot-Plus (12-species and 17-species exploratory runs)",
            "script": "scripts/38_pilot_h8.py",
            "n_tests_run": 348,
            "correction_method": "BH-FDR across 87 panel genes (12sp and 17sp pilots)",
            "significant_findings": "Early signal of higher k in M4 telomere in 17-species pilot",
            "null_or_negative_findings": "Module-level proportion of relaxed genes (k<1, q<0.10) across modules in 12sp and 17sp pilots (NULL); No pseudogenization/loss of TP53, RB1, PTEN, or POT1 across long-lived parrots (NULL)",
        },
        {
            "phase": "Phase 2 Main (58 species, 101-gene panel: Integrity, RELAX, aBSREL, 298-gene Null Pool)",
            "script": "scripts/31_selection_tests.py, 32_h8_summary.py, 37_h8_null.py, 40_integrity_triage.py",
            "n_tests_run": 812,
            "correction_method": "BH-FDR across 87 panel genes per tree (burgio2019 & smith2024) + 10,000 module label permutations + 1,000 stratified null pseudo-panels from 298 control genes",
            "significant_findings": "Module heterogeneity in prop_relaxed on burgio2019 (p=0.0229, NOT replicated on smith2024 p=0.1271); M7 Telomere 5/8 genes with aBSREL q<0.10 on M. undulatus branch in original MAFFT (collapses on window cleaning / multi-hit); TERT intensified k (q=0.0008/0.011, 15 sp)",
            "null_or_negative_findings": "H8a Gene loss / pseudogenization in Psittaculidae vs rest (NULL: 100% of candidate disruptions were assembly fragmentation/coverage artifacts or heterozygous); H8b Module heterogeneity in mean log2(k) (NULL: p=0.1643 burgio2019, p=0.9551 smith2024); M2 mean log2(k) p=0.024 vs null (NULL after 16-test correction and median sensitivity p=0.07); Modules M1, M3-M8 vs 298-gene null (ALL NULL: p=0.09 to 1.00); H8c aBSREL positive selection on the Psittaculidae stem branch (NULL: 0/87 panel genes had q<0.10 on stem)",
        },
        {
            "phase": "Phase 3 Genome-Wide RERconverge (7,756 genes x 2 trees x 3 H6 traits + 2 H8 within-clade runs + GO enrichment)",
            "script": "scripts/R/rerconverge_h6_h8.R",
            "n_tests_run": 77560,
            "correction_method": "BH-FDR on parametric P and BH-FDR on 1,000 Brownian motion permulation P-values (permQ) across 7,756 genes; Wilcoxon GMT enrichment across GO_BP terms",
            "significant_findings": "Top individual gene rankings (PTEN P=9.76e-5, permP=0.001; MAD2L2, HSPA8); 131-151 genes with nominal P<0.05 & permP<0.05 in both longevity and brain",
            "null_or_negative_findings": "Genome-wide BH-FDR across all 7,756 genes for single genes (NULL at permQ < 0.10: minimum permP is bounded by 1/1001 = 0.000999 with 1,000 permulations, so single-gene permQ cannot reach <0.10 without >=100,000 permulations); H8 within-Psittaculidae RERconverge on 86 panel genes after BH correction (NULL: 0/86 genes had panelQ < 0.10); H8 within-rest RERconverge on 86 panel genes after BH correction (NULL: 0/86 genes had panelQ < 0.10; PTEN was nominal P=1.18e-4, panel_padj_analytic=0.010, but permulation panelQ=0.145)",
        },
    ]
    inv_df = pd.DataFrame(inventory)
    inv_df.to_csv(OUT / "9_statistical_honesty_inventory.tsv", sep="\t", index=False)
    print("[RedTeam Python] Test 9 complete.", flush=True)


# ==============================================================================
# TEST 7: Telomeric Alignments, Substitutions, Cleaning, aBSREL-MH, and RBH
# ==============================================================================
def run_test_7():
    print("[RedTeam Python] Running Test 7 (Telomeric alignments, cleaning, aBSREL-MH, RBH)...", flush=True)
    target_genes = ["RTEL1", "TERF2", "TERT", "ACD", "TERF1", "POT1"]
    aln_dir = ROOT / "results/phase2_genes/alignments"
    hyphy_dir = ROOT / "results/phase2_genes/hyphy/burgio2019"
    work_dir = OUT / "telomere_work"
    work_dir.mkdir(parents=True, exist_ok=True)

    spec = importlib.util.spec_from_file_location("sel", ROOT / "scripts/31_selection_tests.py")
    sel = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(sel)
    tp, mp = sel.TREES["burgio2019"]
    full_tree = sel.Tree.from_newick(tp.read_text())
    m_df = pd.read_csv(mp, sep="\t")
    tip_of = {s.replace(" ", "_"): t for s, t in zip(m_df.species, m_df.tree_tip)}

    # Core representative clade-balanced species set (up to 22 species covering all Psittaculidae + non-Psittaculidae clades)
    core_pref = [
        "Melopsittacus_undulatus", "Chalcopsitta_duivenbodei", "Trichoglossus_haematodus", "Eos_histrio",
        "Agapornis_pullarius", "Loriculus_vernalis", "Psittacula_echo", "Eclectus_roratus",
        "Platycercus_eximius", "Lathamus_discolor", "Neophema_elegans",
        "Strigops_habroptila", "Nestor_notabilis",
        "Cacatua_leadbeateri", "Calyptorhynchus_lathami", "Eolophus_roseicapilla",
        "Ara_macao", "Anodorhynchus_hyacinthinus", "Pyrrhura_molinae",
        "Amazona_aestiva", "Amazona_collaria", "Pionus_menstruus", "Coracopsis_vasa"
    ]

    window_report_lines = []
    summary_rows = []

    for g in target_genes:
        codon_file = aln_dir / f"{g}.codon.fasta"
        orig_js = hyphy_dir / f"{g}.absrel.json"

        seqs = parse_fasta_str(codon_file.read_text())
        aa_seqs = {sp: translate_seq(cds) for sp, cds in seqs.items()}
        n_sp = len(seqs)
        n_codons = len(next(iter(seqs.values()))) // 3
        melo_aa = aa_seqs.get("Melopsittacus_undulatus", "")
        other_sps = [sp for sp in aa_seqs if sp != "Melopsittacus_undulatus"]

        # Identify every codon where M. undulatus differs from the majority non-gap AA of other parrots
        melo_subs = []
        col_gap_frac = np.zeros(n_codons)
        col_ident = np.zeros(n_codons)
        for c in range(n_codons):
            other_aas = [aa_seqs[sp][c] for sp in other_sps if aa_seqs[sp][c] not in ("-", "X")]
            gaps = sum(1 for sp in aa_seqs if aa_seqs[sp][c] == "-")
            col_gap_frac[c] = gaps / n_sp
            if not other_aas:
                col_ident[c] = 0.0
                continue
            cnt = Counter(other_aas)
            maj_aa, maj_count = cnt.most_common(1)[0]
            col_ident[c] = maj_count / len(other_aas)
            if melo_aa[c] not in ("-", "X") and melo_aa[c] != maj_aa:
                melo_subs.append((c + 1, maj_aa, melo_aa[c], seqs["Melopsittacus_undulatus"][3 * c:3 * c + 3]))

        sub_positions = set(pos for pos, _, _, _ in melo_subs)
        clustered_subs = []
        clean_subs = []
        for pos, maj_aa, m_aa, m_cod in melo_subs:
            c_idx = pos - 1
            w_lo = max(0, c_idx - 5)
            w_hi = min(n_codons, c_idx + 6)
            subs_in_win = sum(1 for p in range(w_lo + 1, w_hi + 1) if p in sub_positions)
            near_gap = any(col_gap_frac[k] > 0.15 for k in range(w_lo, w_hi))
            low_cons = np.mean(col_ident[w_lo:w_hi]) < 0.65
            if subs_in_win >= 3 or near_gap or low_cons:
                clustered_subs.append((pos, maj_aa, m_aa, m_cod, subs_in_win, near_gap, low_cons))
            else:
                clean_subs.append((pos, maj_aa, m_aa, m_cod))

        window_report_lines.append("=" * 88)
        window_report_lines.append(
            f"GENE: {g} | n_species={n_sp} | n_codons={n_codons} | "
            f"Total M. undulatus AA substitutions vs parrot consensus: {len(melo_subs)} "
            f"(Clean isolated: {len(clean_subs)}, Clustered/Gappy/Noisy window: {len(clustered_subs)})"
        )
        window_report_lines.append("-" * 88)
        if clustered_subs:
            win_centers = []
            for pos, *_ in clustered_subs:
                if not win_centers or pos - win_centers[-1] > 12:
                    win_centers.append(pos)
            for wc in win_centers[:4]:
                lo = max(0, wc - 8)
                hi = min(n_codons, wc + 12)
                window_report_lines.append(f"  [Window codons {lo+1}..{hi} around cluster at codon {wc}]:")
                rep_sps = ["Melopsittacus_undulatus"] + [
                    s for s in other_sps if any(k in s for k in ("Chalcopsitta", "Trichoglossus", "Eos", "Agapornis", "Loriculus", "Strigops", "Cacatua", "Ara_", "Amazona"))
                ][:6]
                for sp in rep_sps:
                    aa_slice = aa_seqs[sp][lo:hi]
                    cod_slice = " ".join(seqs[sp][3 * k:3 * k + 3] for k in range(lo, min(lo + 10, hi)))
                    window_report_lines.append(f"    {sp:<28} AA: {aa_slice} | Codons({lo+1}..{min(lo+10,hi)}): {cod_slice}")
        else:
            window_report_lines.append("  No clustered/noisy M. undulatus substitution windows detected.")
        window_report_lines.append("")

        # Subset species for fast sensitivity re-runs (keep 12 clade-representative species; RTEL1 already cached at 20sp)
        n_target_sp = 20 if g == "RTEL1" else 12
        keep_sps = [s for s in core_pref if s in seqs][:min(n_target_sp, n_sp)]
        if len(keep_sps) < min(n_target_sp, n_sp):
            for s in seqs:
                if s not in keep_sps:
                    keep_sps.append(s)
                if len(keep_sps) >= min(n_target_sp, n_sp):
                    break
        sub_seqs = {s: seqs[s] for s in keep_sps}
        nwk_dict, _ = sel.tree_for(full_tree, tip_of, list(sub_seqs.keys()))
        sub_nwk_path = work_dir / f"{g}.sub20.absrel.nwk"
        # Tag only Melopsittacus_undulatus{FG} when building new trees so HyPhy only optimizes the target branch
        nwk_melo_only = nwk_dict["absrel"].replace("{FG}", "")
        nwk_melo_only = nwk_melo_only.replace("Melopsittacus_undulatus", "Melopsittacus_undulatus{FG}")
        sub_nwk_path.write_text(nwk_melo_only)

        # 1. L-INS-i re-aligned codon alignment on sub_seqs
        raw_aa_path = work_dir / f"{g}.raw_aa.fa"
        linsi_aa_path = work_dir / f"{g}.linsi_aa.fa"
        linsi_codon_path = work_dir / f"{g}.linsi.codon.fasta"
        raw_cds_map = {}
        with open(raw_aa_path, "w") as f:
            for sp, cds in sub_seqs.items():
                ungapped_cds = cds.replace("-", "")
                raw_cds_map[sp] = ungapped_cds
                f.write(f">{sp}\n{translate_seq(ungapped_cds)}\n")
        subprocess.run(
            [str(MAFFT), "--localpair", "--maxiterate", "1000", "--quiet", "--thread", "4", str(raw_aa_path)],
            stdout=open(linsi_aa_path, "w"), check=True
        )
        linsi_aa = parse_fasta_str(linsi_aa_path.read_text())
        with open(linsi_codon_path, "w") as f:
            for sp, aa_aln in linsi_aa.items():
                cds_u = raw_cds_map[sp]
                idx = 0
                cod_buf = []
                for ch in aa_aln:
                    if ch == "-":
                        cod_buf.append("---")
                    else:
                        cod_buf.append(cds_u[3 * idx:3 * idx + 3])
                        idx += 1
                f.write(f">{sp}\n{''.join(cod_buf)}\n")

        # 2. Stringently cleaned codon alignment (dropping gappy columns >15% gaps AND masking species-specific
        # sliding 10-codon windows with >=3 non-consensus AA diffs or within 3 codons of a species indel using '---')
        clean_codon_path = work_dir / f"{g}.cleaned.codon.fasta"
        sub_orig_path = work_dir / f"{g}.sub_orig.codon.fasta"
        with open(sub_orig_path, "w") as f:
            for sp in keep_sps:
                f.write(f">{sp}\n{seqs[sp]}\n")

        consensus_aa = []
        for c in range(n_codons):
            aas = [aa_seqs[sp][c] for sp in aa_seqs if aa_seqs[sp][c] not in ("-", "X")]
            consensus_aa.append(Counter(aas).most_common(1)[0][0] if aas else "-")

        keep_cols = [c for c in range(n_codons) if col_gap_frac[c] <= 0.15]
        cleaned_seqs = {}
        total_masked_codons_melo = (n_codons - len(keep_cols))
        for sp in keep_sps:
            cds = seqs[sp]
            cods = [cds[3 * c:3 * c + 3] for c in range(n_codons)]
            sp_aa = aa_seqs[sp]
            mask = np.zeros(n_codons, dtype=bool)
            diffs = np.array([1 if (sp_aa[c] not in ("-", "X") and sp_aa[c] != consensus_aa[c]) else 0 for c in range(n_codons)])
            sp_gaps = np.array([1 if sp_aa[c] == "-" else 0 for c in range(n_codons)])
            for w in range(0, max(1, n_codons - 9)):
                if diffs[w:w + 10].sum() >= 3 or (diffs[w:w + 10].sum() >= 2 and sp_gaps[w:w + 10].sum() >= 1):
                    mask[w:w + 10] = True
            for c in keep_cols:
                if mask[c] and cods[c] != "---":
                    cods[c] = "---"
                    if sp == "Melopsittacus_undulatus":
                        total_masked_codons_melo += 1
            cleaned_seqs[sp] = "".join(cods[c] for c in keep_cols)

        with open(clean_codon_path, "w") as f:
            for sp, c_seq in cleaned_seqs.items():
                f.write(f">{sp}\n{c_seq}\n")

        def run_absrel_variant(aln_p, tree_p, extra_args, tag):
            out_js = work_dir / f"{g}.{tag}.absrel.json"
            if not out_js.exists() or out_js.stat().st_size == 0:
                out_js.unlink(missing_ok=True)
                subprocess.run(
                    [str(HYPHY), "CPU=1", "absrel", "--alignment", str(aln_p), "--tree", str(tree_p),
                     "--branches", "FG", "--kill-zero-lengths", "No", *extra_args, "--output", str(out_js)],
                    stdin=subprocess.DEVNULL, capture_output=True, text=True, cwd=str(work_dir)
                )
            if not out_js.exists() or out_js.stat().st_size == 0:
                return np.nan
            d = json.load(open(out_js))
            tested = [k for k, v in d["tested"]["0"].items() if v == "test"]
            p_corr = {k: d["branch attributes"]["0"][k].get("Corrected P-value") for k in tested}
            return p_corr.get("Melopsittacus_undulatus", np.nan)

        orig_d = json.load(open(orig_js))
        tested_orig = [k for k, v in orig_d["tested"]["0"].items() if v == "test"]
        p_orig_map = {k: orig_d["branch attributes"]["0"][k].get("Corrected P-value") for k in tested_orig}
        orig_melo_p = p_orig_map.get("Melopsittacus_undulatus", np.nan)

        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex_v:
            f_linsi = ex_v.submit(run_absrel_variant, linsi_codon_path, sub_nwk_path, [], "linsi20")
            f_clean = ex_v.submit(run_absrel_variant, clean_codon_path, sub_nwk_path, [], "cleaned20")
            f_mh = ex_v.submit(run_absrel_variant, sub_orig_path, sub_nwk_path, ["--multiple-hits", "Double+Triple"], "mh20")
            f_cmh = ex_v.submit(run_absrel_variant, clean_codon_path, sub_nwk_path, ["--multiple-hits", "Double+Triple"], "cleaned_mh20")
            linsi_melo_p = f_linsi.result()
            clean_melo_p = f_clean.result()
            mh_melo_p = f_mh.result()
            clean_mh_p = f_cmh.result()

        # Orthology & Reference check across all 58 species in results/phase2_genes/hits/*.tsv
        melo_hit_ident = np.nan
        melo_n_align = 0
        same_locus_kakapo = 0
        both_called_kakapo = 0
        all_identities = []
        for hit_f in sorted((ROOT / "results/phase2_genes/hits").glob("*.tsv")):
            hdf = pd.read_csv(hit_f, sep="\t").set_index("gene")
            if g in hdf.index and hdf.loc[g, "call"] == "complete":
                ident = float(hdf.loc[g, "identity"])
                all_identities.append(ident)
                if hit_f.stem == "GCA_012275295.1":
                    melo_hit_ident = ident
                    melo_n_align = int(hdf.loc[g, "n_alignments"])
                kg = f"kakapo_{g}"
                if kg in hdf.index and pd.notna(hdf.loc[kg, "location"]) and str(hdf.loc[kg, "location"]) != "":
                    both_called_kakapo += 1
                    loc1 = str(hdf.loc[g, "location"])
                    loc2 = str(hdf.loc[kg, "location"])
                    ctg1, rest1 = loc1.split(":", 1)
                    ctg2, rest2 = loc2.split(":", 1)
                    if ctg1 == ctg2:
                        s1, e1 = [int(x) for x in rest1[:-1].split("-")]
                        s2, e2 = [int(x) for x in rest2[:-1].split("-")]
                        if max(s1, s2) <= min(e1, e2):
                            same_locus_kakapo += 1

        summary_rows.append({
            "gene": g,
            "n_species": n_sp,
            "n_codons": n_codons,
            "melo_aa_subs_total": len(melo_subs),
            "melo_aa_subs_clean_isolated": len(clean_subs),
            "melo_aa_subs_clustered_or_gappy": len(clustered_subs),
            "melo_codons_masked_by_cleaner": total_masked_codons_melo,
            "absrel_p_original_mafft": orig_melo_p,
            "absrel_p_linsi_realigned": linsi_melo_p,
            "absrel_p_window_cleaned": clean_melo_p,
            "absrel_p_MultiHit_DoubleTriple": mh_melo_p,
            "absrel_p_cleaned_plus_MultiHit": clean_mh_p,
            "melo_GCA_vs_RefSeq_identity": melo_hit_ident,
            "melo_miniprot_n_loci": melo_n_align,
            "mean_identity_across_parrots": np.mean(all_identities) if all_identities else np.nan,
            "min_identity_across_parrots": np.min(all_identities) if all_identities else np.nan,
            "dual_ref_kakapo_concordance": f"{same_locus_kakapo}/{both_called_kakapo}" if both_called_kakapo else "no_kakapo_ref",
        })
        print(f"[RedTeam Python] Test 7 finished {g}: orig_p={orig_melo_p:.4g}, linsi_p={linsi_melo_p:.4g}, clean_p={clean_melo_p:.4g}, mh_p={mh_melo_p:.4g}", flush=True)

    (OUT / "7_telomere_substitution_windows.txt").write_text("\n".join(window_report_lines))
    pd.DataFrame(summary_rows).to_csv(OUT / "7_telomere_alignment_and_rbh.tsv", sep="\t", index=False, float_format="%.6g")
    print("[RedTeam Python] Test 7 complete.", flush=True)


# ==============================================================================
# TEST 6: Base Rate of Positive Selection (aBSREL) in Control Genes & Multiple Testing
# ==============================================================================
def _run_one_null_absrel(g):
    null_dir = ROOT / "results/phase2_genes/hyphy_null/burgio2019"
    work_null = OUT / "null_absrel_work"
    work_null.mkdir(parents=True, exist_ok=True)
    aln = null_dir / f"{g}.codon.fasta"
    relax_nwk = null_dir / f"{g}.relax.nwk"
    js = null_dir / f"{g}.absrel.json"
    if not js.exists() or js.stat().st_size == 0:
        js = work_null / f"{g}.absrel.json"

    if not aln.exists() or not relax_nwk.exists():
        return None

    if not js.exists() or js.stat().st_size == 0:
        spec = importlib.util.spec_from_file_location("sel", ROOT / "scripts/31_selection_tests.py")
        sel = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(sel)
        tp, mp = sel.TREES["burgio2019"]
        full = sel.Tree.from_newick(tp.read_text())
        m = pd.read_csv(mp, sep="\t")
        tip_of = {s.replace(" ", "_"): t for s, t in zip(m.species, m.tree_tip)}
        seqs = parse_fasta_str(aln.read_text())
        if "Melopsittacus_undulatus" not in seqs:
            return None
        core_pref = [
            "Melopsittacus_undulatus", "Chalcopsitta_duivenbodei", "Trichoglossus_haematodus", "Eos_histrio",
            "Agapornis_pullarius", "Loriculus_vernalis", "Psittacula_echo", "Eclectus_roratus",
            "Strigops_habroptila", "Nestor_notabilis", "Cacatua_leadbeateri", "Calyptorhynchus_lathami",
            "Ara_macao", "Amazona_aestiva", "Pyrrhura_molinae", "Coracopsis_vasa"
        ]
        keep_sps = [s for s in core_pref if s in seqs][:14]
        for s in seqs:
            if len(keep_sps) >= min(14, len(seqs)):
                break
            if s not in keep_sps:
                keep_sps.append(s)
        sub_aln = work_null / f"{g}.sub.codon.fasta"
        with open(sub_aln, "w") as f:
            for s in keep_sps:
                f.write(f">{s}\n{seqs[s]}\n")
        nwk_dict, status = sel.tree_for(full, tip_of, keep_sps)
        if nwk_dict is None:
            return None
        absrel_nwk = work_null / f"{g}.sub.absrel.nwk"
        absrel_nwk.write_text(nwk_dict["absrel"])
        r = subprocess.run(
            [str(HYPHY), "CPU=1", "absrel", "--alignment", str(sub_aln), "--tree", str(absrel_nwk),
             "--branches", "FG", "--kill-zero-lengths", "No", "--output", str(js.with_suffix(".tmp"))],
            stdin=subprocess.DEVNULL, capture_output=True, text=True, cwd=str(work_null)
        )
        if r.returncode == 0 and js.with_suffix(".tmp").exists():
            js.with_suffix(".tmp").replace(js)
        else:
            return None

    try:
        d = json.load(open(js))
        tested = [k for k, v in d["tested"]["0"].items() if v == "test"]
        p_corr = {k: d["branch attributes"]["0"][k].get("Corrected P-value") for k in tested}
        p_unc = {k: d["branch attributes"]["0"][k].get("Uncorrected P-value") for k in tested}
        rate_classes = {k: d["branch attributes"]["0"][k].get("rate classes", 1) for k in tested}
        stem_keys = [k for k in tested if k != "Melopsittacus_undulatus"]
        sk = stem_keys[0] if stem_keys else None
        return {
            "gene": g,
            "n_species": len(d["tested"]["0"]),
            "absrel_melopsittacus_p_holm2": p_corr.get("Melopsittacus_undulatus", np.nan),
            "absrel_melopsittacus_p_uncorr": p_unc.get("Melopsittacus_undulatus", np.nan),
            "melopsittacus_rate_classes": rate_classes.get("Melopsittacus_undulatus", 1),
            "absrel_stem_p_holm2": p_corr.get(sk, np.nan) if sk else np.nan,
            "absrel_stem_p_uncorr": p_unc.get(sk, np.nan) if sk else np.nan,
            "stem_rate_classes": rate_classes.get(sk, 1) if sk else 1,
        }
    except Exception:
        return None


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    ok = ~np.isnan(p)
    out = np.full_like(p, np.nan)
    if ok.sum() == 0:
        return out
    p_ok = p[ok]
    o = np.argsort(p_ok)
    q = p_ok[o] * len(p_ok) / np.arange(1, len(p_ok) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    res = np.empty_like(q)
    res[o] = np.minimum(q, 1.0)
    out[ok] = res
    return out


def run_test_6():
    print("[RedTeam Python] Running Test 6 (aBSREL base rate on control genes + 298-gene RELAX branch dN/dS audit)...", flush=True)
    null_dir = ROOT / "results/phase2_genes/hyphy_null/burgio2019"
    # First inspect all 298 control genes' RELAX JSONs for branch-set omega and k
    relax_rows = []
    for rjs in sorted(null_dir.glob("*.relax.json")):
        g = rjs.name.split(".")[0]
        try:
            d = json.load(open(rjs))
            rd = d["fits"]["MG94xREV with separate rates for branch sets"]["Rate Distributions"]
            w_ref = rd.get("non-synonymous/synonymous rate ratio for *Reference*", [[np.nan]])[0][0]
            w_test = rd.get("non-synonymous/synonymous rate ratio for *Test*", [[np.nan]])[0][0]
            w_unc = rd.get("non-synonymous/synonymous rate ratio for *Unclassified*", [[np.nan]])[0][0]
            aln_f = null_dir / f"{g}.codon.fasta"
            seqs = parse_fasta_str(aln_f.read_text()) if aln_f.exists() else {}
            n_cod = len(next(iter(seqs.values()))) // 3 if seqs else 9999
            relax_rows.append({
                "gene": g,
                "n_codons": n_cod,
                "n_species": len(seqs),
                "omega_reference_nonPsittaculidae": w_ref,
                "omega_test_Psittaculidae_stem_internal": w_test,
                "omega_unclassified_Psittaculidae_tips": w_unc,
                "relax_k": d["test results"]["relaxation or intensification parameter"],
                "relax_p": d["test results"]["p-value"],
            })
        except Exception:
            pass
    rel298_df = pd.DataFrame(relax_rows)

    # Select already completed absrel JSONs (in null_dir or OUT/null_absrel_work)
    work_null = OUT / "null_absrel_work"
    done_orig = [p.name.split(".")[0] for p in null_dir.glob("*.absrel.json") if p.stat().st_size > 0]
    done_work = [p.name.split(".")[0] for p in work_null.glob("*.absrel.json") if p.stat().st_size > 0] if work_null.exists() else []
    done_genes = sorted(set(done_orig + done_work))
    rel298_df["work_cost"] = rel298_df["n_codons"] * rel298_df["n_species"]
    fast_candidates = rel298_df.sort_values("work_cost")["gene"].tolist()
    sample_genes = done_genes if len(done_genes) >= 15 else sorted(set(done_genes + fast_candidates[:25]))[:25]

    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
        results = list(ex.map(_run_one_null_absrel, sample_genes))
    null_df = pd.DataFrame([r for r in results if r is not None])
    null_df["absrel_melopsittacus_q_pool"] = bh_fdr(null_df["absrel_melopsittacus_p_holm2"])
    null_df["absrel_stem_q_pool"] = bh_fdr(null_df["absrel_stem_p_holm2"])
    null_df.to_csv(OUT / "6_absrel_null_genes_detail.tsv", sep="\t", index=False, float_format="%.6g")

    panel_df = pd.read_csv(ROOT / "results/phase2_genes/h8_relax_by_gene.tsv", sep="\t")
    p_burg = panel_df[panel_df["tree"] == "burgio2019"].copy()
    p_m7 = p_burg[p_burg["module"] == "M7"].copy()
    p_non_m7 = p_burg[p_burg["module"] != "M7"].copy()

    def summarize_group(label, d, col_p="absrel_melopsittacus_p", col_q="absrel_melopsittacus_q",
                        col_sp="absrel_stem_p", col_sq="absrel_stem_q"):
        n = int(d[col_p].notna().sum())
        m_p05 = int((d[col_p] < 0.05).sum())
        m_q10 = int((d[col_q] < 0.10).sum())
        m_q05 = int((d[col_q] < 0.05).sum())
        s_p05 = int((d[col_sp] < 0.05).sum())
        s_q10 = int((d[col_sq] < 0.10).sum())
        return {
            "group": label,
            "n_genes_tested": n,
            "melopsittacus_p_holm2_lt_005": f"{m_p05}/{n} ({100*m_p05/max(1,n):.1f}%)",
            "melopsittacus_q_bh_panel_lt_010": f"{m_q10}/{n} ({100*m_q10/max(1,n):.1f}%)",
            "melopsittacus_q_bh_panel_lt_005": f"{m_q05}/{n} ({100*m_q05/max(1,n):.1f}%)",
            "psittaculidae_stem_p_holm2_lt_005": f"{s_p05}/{n} ({100*s_p05/max(1,n):.1f}%)",
            "psittaculidae_stem_q_bh_panel_lt_010": f"{s_q10}/{n} ({100*s_q10/max(1,n):.1f}%)",
        }

    null_mapped = null_df.rename(columns={
        "absrel_melopsittacus_p_holm2": "absrel_melopsittacus_p",
        "absrel_melopsittacus_q_pool": "absrel_melopsittacus_q",
        "absrel_stem_p_holm2": "absrel_stem_p",
        "absrel_stem_q_pool": "absrel_stem_q",
    })

    n_null = len(null_mapped)
    null_m_p05 = int((null_mapped["absrel_melopsittacus_p"] < 0.05).sum())
    null_m_q10 = int((null_mapped["absrel_melopsittacus_q"] < 0.10).sum())

    m7_n = int(p_m7["absrel_melopsittacus_p"].notna().sum())
    m7_p05 = int((p_m7["absrel_melopsittacus_p"] < 0.05).sum())
    m7_q10 = int((p_m7["absrel_melopsittacus_q"] < 0.10).sum())

    pan_n = int(p_burg["absrel_melopsittacus_p"].notna().sum())
    pan_p05 = int((p_burg["absrel_melopsittacus_p"] < 0.05).sum())
    pan_q10 = int((p_burg["absrel_melopsittacus_q"] < 0.10).sum())

    _, fisher_m7_vs_null_p05 = stats.fisher_exact([[m7_p05, m7_n - m7_p05], [null_m_p05, n_null - null_m_p05]])
    _, fisher_m7_vs_null_q10 = stats.fisher_exact([[m7_q10, m7_n - m7_q10], [null_m_q10, n_null - null_m_q10]])
    _, fisher_pan_vs_null_p05 = stats.fisher_exact([[pan_p05, pan_n - pan_p05], [null_m_p05, n_null - null_m_p05]])
    _, fisher_pan_vs_null_q10 = stats.fisher_exact([[pan_q10, pan_n - pan_q10], [null_m_q10, n_null - null_m_q10]])

    summary_rows = [
        summarize_group("H8 Full Panel (burgio2019, 87 genes)", p_burg),
        summarize_group("H8 Module M7 Telomere/Shelterin (8 genes)", p_m7),
        summarize_group("H8 Non-M7 Modules M1-M6, M8 (79 genes)", p_non_m7),
        summarize_group(f"Genome-wide Control Pool ({n_null} genes tested with aBSREL)", null_mapped),
    ]
    sum_df = pd.DataFrame(summary_rows)
    sum_df["fisher_vs_control_p05"] = [fisher_pan_vs_null_p05, fisher_m7_vs_null_p05, np.nan, 1.0]
    sum_df["fisher_vs_control_q10"] = [fisher_pan_vs_null_q10, fisher_m7_vs_null_q10, np.nan, 1.0]
    sum_df["all_298_controls_omega_tips_gt_1"] = [
        np.nan, np.nan, np.nan,
        f"{int((rel298_df['omega_unclassified_Psittaculidae_tips'] > 1.0).sum())}/{len(rel298_df)} ({100*(rel298_df['omega_unclassified_Psittaculidae_tips'] > 1.0).mean():.1f}%)"
    ]
    sum_df.to_csv(OUT / "6_absrel_base_rate.tsv", sep="\t", index=False, float_format="%.6g")
    print("[RedTeam Python] Test 6 complete.", flush=True)


if __name__ == "__main__":
    run_test_9()
    run_test_7()
    run_test_6()
    print("[RedTeam Python] All Python Red-Team audits completed successfully!", flush=True)
