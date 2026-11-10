#!/usr/bin/env python3
"""Builds a self-contained interactive HTML raw-data dashboard at docs/parrots_raw_data_dashboard.html
containing all raw data across Phases 0-3, 126 assemblies / 71 validated species, Phase 2 H8 58 vs 71sp,
all 7,756 RERconverge genome-wide genes (with Tier A zero-branch vs Tier B informative filtering),
3D AlphaFold summaries, and all 9 Red-Team stress tests.
Zero external CDN dependencies (pure inline HTML + CSS + SVG + Vanilla JS).
"""
import json
import re
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def df_to_records(df):
    recs = []
    for _, row in df.iterrows():
        d = {}
        for k, v in row.items():
            if pd.isna(v):
                d[k] = None
            elif isinstance(v, (np.floating, float)):
                d[k] = round(float(v), 5)
            elif isinstance(v, (np.integer, int)):
                d[k] = int(v)
            elif isinstance(v, (np.bool_, bool)):
                d[k] = bool(v)
            else:
                d[k] = str(v)
        recs.append(d)
    return recs


def main():
    # 1. Species & Audit Raw Data
    ph = pd.read_csv(ROOT / "results/phase0/phenotypes.tsv", sep="\t")
    mat = pd.read_csv(ROOT / "results/phase0b/species_data_matrix_with_residuals.tsv", sep="\t")
    sg = pd.read_csv(ROOT / "results/phase1_audit/species_genomes.tsv", sep="\t")
    t71 = pd.read_csv(ROOT / "results/phase2_genes_71sp/rer_traits_71sp.tsv", sep="\t")
    aud = pd.read_csv(ROOT / "results/phase1_audit/audit_table.tsv", sep="\t")

    sp_df = ph[[
        "species", "life_expectancy_years", "max_longevity_years", "body_mass_g",
        "brain_mass_g", "clutch_size", "incubation_days"
    ]].merge(
        mat[["species", "family", "genus", "rel_log_life_expectancy", "rel_log_max_longevity_ssadj"]],
        on="species", how="left"
    ).merge(
        sg[["species", "group", "psittaculidae", "status", "accession", "source", "contig_n50", "busco_C", "identity"]],
        on="species", how="left"
    ).merge(
        t71[["species", "rel_log_brain", "rel_brain_given_le", "cohort"]],
        on="species", how="left"
    ).merge(
        aud[["species", "mito_overall", "mito_foreign_species", "mash_resid_vs_tree"]].drop_duplicates("species"),
        on="species", how="left"
    )
    sp_df["cohort"] = sp_df["cohort"].fillna(
        sp_df["status"].map(lambda s: "excluded" if pd.notna(s) and str(s).startswith("excluded") else "no_genome")
    )

    # 2. Phase 2 H8 & 71sp & Telomere RedTeam
    h8 = pd.read_csv(ROOT / "results/phase2_genes/h8_relax_by_gene.tsv", sep="\t")
    h8_b = h8[h8.tree == "burgio2019"].copy()
    comp71_p = ROOT / "results/phase2_genes_71sp/comparison_58_vs_71_confirmatory.tsv"
    if comp71_p.exists():
        comp71 = pd.read_csv(comp71_p, sep="\t")
    else:
        comp71 = pd.DataFrame(columns=["gene"])
    integ71 = pd.read_csv(ROOT / "results/phase2_genes_71sp/h8_integrity_71sp.tsv", sep="\t")
    s71 = pd.read_csv(ROOT / "results/phase2_genes_71sp/alignment_summary.tsv", sep="\t")
    tel7 = pd.read_csv(ROOT / "results/redteam_audit/7_telomere_alignment_and_rbh.tsv", sep="\t")

    p2_df = h8_b[[
        "gene", "module", "k", "relax_p", "relax_q", "absrel_stem_p", "absrel_melopsittacus_p", "absrel_melopsittacus_q"
    ]].merge(
        integ71[["gene", "disrupted_psittaculidae_71sp", "complete_psittaculidae_71sp", "disrupted_rest_71sp", "complete_rest_71sp"]],
        on="gene", how="left"
    ).merge(
        s71[["gene", "n_species", "aln_codons"]].rename(columns={"n_species": "n_sp_71", "aln_codons": "codons_71"}),
        on="gene", how="left"
    )
    if len(comp71) > 0 and "k_71sp" in comp71.columns:
        keep_cols = [c for c in ["gene", "k_71sp", "relax_p_71sp", "absrel_melo_p_71sp"] if c in comp71.columns]
        p2_df = p2_df.merge(comp71[keep_cols], on="gene", how="left")

    # 3. Phase 3 All 7,756 Genes + branch invariance + technical covariates
    le = pd.read_csv(ROOT / "results/phase3_genome/rer_burgio2019_rel_log_life_expectancy.tsv", sep="\t")
    br = pd.read_csv(ROOT / "results/phase3_genome/rer_burgio2019_rel_log_brain.tsv", sep="\t")
    br_p = pd.read_csv(ROOT / "results/phase3_genome/rer_burgio2019_rel_brain_given_le.tsv", sep="\t")
    r5k = pd.read_csv(ROOT / "results/redteam_audit/5_residual_brain_key_genes.tsv", sep="\t")

    tlen, n_nz, zfrac, n_br = {}, {}, {}, {}
    for line in open(ROOT / "results/phase3_genome/genetrees_burgio2019.txt"):
        g, nwk = line.rstrip("\n").split("\t", 1)
        els = [float(x) for x in re.findall(r":([0-9eE.+-]+)", nwk)]
        tlen[g] = sum(els)
        n_br[g] = len(els)
        n_nz[g] = sum(1 for e in els if e > 1e-6)
        zfrac[g] = 1.0 - (n_nz[g] / len(els)) if els else 1.0

    g_df = le[["gene", "Rho", "N", "P", "permP"]].rename(columns={"Rho": "Rho_LE", "P": "P_LE", "permP": "permP_LE"}).merge(
        br[["gene", "Rho", "P", "permP"]].rename(columns={"Rho": "Rho_Br", "P": "P_Br", "permP": "permP_Br"}),
        on="gene"
    ).merge(
        br_p[["gene", "Rho", "P"]].rename(columns={"Rho": "Rho_BrPGLS", "P": "P_BrPGLS"}),
        on="gene"
    )
    g_df["tree_len"] = g_df["gene"].map(tlen)
    g_df["n_nz"] = g_df["gene"].map(n_nz)
    g_df["zero_pct"] = (g_df["gene"].map(zfrac) * 100.0).round(1)
    g_df["tier"] = np.where(
        (g_df["zero_pct"] <= 80.0) & (g_df["n_nz"] >= 20), "B_informative", "A_ultraconserved"
    )

    # 4. Red-Team & 3D & Inventory payloads
    r1 = pd.read_csv(ROOT / "results/redteam_audit/1_reproducibility.tsv", sep="\t")
    r2 = pd.read_csv(ROOT / "results/redteam_audit/2_coupling_null_and_controls.tsv", sep="\t")
    r3c = pd.read_csv(ROOT / "results/redteam_audit/3_coupling_after_confounder_removal.tsv", sep="\t")
    r4 = pd.read_csv(ROOT / "results/redteam_audit/4_leave_one_clade_out.tsv", sep="\t")
    r5s = pd.read_csv(ROOT / "results/redteam_audit/5_residual_brain_summary.tsv", sep="\t")
    r6b = pd.read_csv(ROOT / "results/redteam_audit/6_absrel_base_rate.tsv", sep="\t")
    r8cv = pd.read_csv(ROOT / "results/redteam_audit/8_holdout_cv_key_genes.tsv", sep="\t")
    r9 = pd.read_csv(ROOT / "results/redteam_audit/9_statistical_honesty_inventory.tsv", sep="\t")
    s3d = pd.read_csv(ROOT / "results/phase3_genome/structural_3d_summary.tsv", sep="\t")
    tel_win = (ROOT / "results/redteam_audit/7_telomere_substitution_windows.txt").read_text()

    payload = {
        "species": df_to_records(sp_df),
        "panel": df_to_records(p2_df),
        "telomere_rt": df_to_records(tel7),
        "genome": df_to_records(g_df),
        "rt_repro": df_to_records(r1),
        "rt_null2": df_to_records(r2),
        "rt_conf": df_to_records(r3c),
        "rt_loco": df_to_records(r4),
        "rt_orth": df_to_records(r5s),
        "rt_orth_genes": df_to_records(r5k),
        "rt_base6": df_to_records(r6b),
        "rt_cv": df_to_records(r8cv),
        "rt_inv": df_to_records(r9),
        "struct_3d": df_to_records(s3d),
        "tel_windows_txt": tel_win,
    }

    html = HTML_TEMPLATE.replace("__DATA_JSON__", json.dumps(payload, separators=(",", ":")))
    out_path = ROOT / "docs/parrots_raw_data_dashboard.html"
    out_path.write_text(html)
    print(f"Wrote {out_path} ({out_path.stat().st_size / 1e6:.2f} MB)")


HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Parrots Comparative Genomics — Raw Data & Red-Team Audit Dashboard</title>
<style>
:root {
  --bg: #0b101b; --card: #131b2c; --card2: #182238; --border: #25334e;
  --text: #e6edf7; --muted: #93a4c3; --accent: #38bdf8; --red: #f87171;
  --green: #34d399; --amber: #fbbf24; --purple: #a78bfa;
}
* { box-sizing: border-box; }
body { margin: 0; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; background: var(--bg); color: var(--text); font-size: 13px; }
header { background: linear-gradient(120deg, #0f172a, #172554); border-bottom: 1px solid var(--border); padding: 14px 22px; display: flex; flex-wrap: wrap; justify-content: space-between; align-items: center; gap: 12px; }
header h1 { margin: 0; font-size: 18px; font-weight: 700; letter-spacing: -0.02em; }
header p { margin: 3px 0 0; color: var(--muted); font-size: 12px; }
.tabs { display: flex; gap: 6px; padding: 10px 22px 0; background: #0e1525; border-bottom: 1px solid var(--border); flex-wrap: wrap; }
.tab { padding: 8px 14px; border: 1px solid transparent; border-bottom: none; border-radius: 8px 8px 0 0; cursor: pointer; background: transparent; color: var(--muted); font-weight: 600; font-size: 12.5px; }
.tab.active { background: var(--card); color: var(--accent); border-color: var(--border); }
main { padding: 16px 22px 40px; max-width: 1520px; margin: 0 auto; }
.panel { display: none; }
.panel.active { display: block; }
.kpi-row { display: grid; grid-template-columns: repeat(auto-fit, minmax(185px, 1fr)); gap: 10px; margin-bottom: 14px; }
.kpi { background: var(--card); border: 1px solid var(--border); border-radius: 10px; padding: 11px 14px; }
.kpi .lbl { color: var(--muted); font-size: 11px; text-transform: uppercase; letter-spacing: .04em; }
.kpi .val { font-size: 20px; font-weight: 700; margin-top: 3px; color: var(--text); }
.kpi .sub { color: var(--muted); font-size: 11px; margin-top: 2px; }
.grid2 { display: grid; grid-template-columns: 1.15fr 1fr; gap: 14px; margin-bottom: 14px; }
@media (max-width: 1050px) { .grid2 { grid-template-columns: 1fr; } }
.box { background: var(--card); border: 1px solid var(--border); border-radius: 10px; padding: 13px 15px; margin-bottom: 14px; }
.box h3 { margin: 0 0 10px; font-size: 14px; color: #f1f5f9; display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 8px; }
.controls { display: flex; flex-wrap: wrap; gap: 8px; align-items: center; margin-bottom: 10px; }
input, select, button { background: var(--card2); color: var(--text); border: 1px solid var(--border); border-radius: 6px; padding: 5px 9px; font-size: 12px; }
button { cursor: pointer; font-weight: 600; }
button:hover { border-color: var(--accent); color: var(--accent); }
.tbl-wrap { max-height: 440px; overflow: auto; border: 1px solid var(--border); border-radius: 8px; }
table { width: 100%; border-collapse: collapse; font-size: 12px; }
th { position: sticky; top: 0; background: #1a243b; color: #cbd5e1; text-align: left; padding: 7px 8px; cursor: pointer; border-bottom: 1px solid var(--border); white-space: nowrap; }
th:hover { color: var(--accent); }
td { padding: 5px 8px; border-bottom: 1px solid #1e293b; white-space: nowrap; }
tr:hover td { background: rgba(56,189,248,0.07); }
.badge { display: inline-block; padding: 1px 6px; border-radius: 999px; font-size: 10.5px; font-weight: 600; }
.b-red { background: rgba(248,113,113,0.16); color: #fca5a5; }
.b-green { background: rgba(52,211,153,0.16); color: #6ee7b7; }
.b-blue { background: rgba(56,189,248,0.16); color: #7dd3fc; }
.b-amber { background: rgba(251,191,36,0.16); color: #fcd34d; }
svg.plot { width: 100%; height: 340px; background: #0c1322; border-radius: 8px; border: 1px solid var(--border); display: block; }
.tip { position: fixed; pointer-events: none; background: rgba(15,23,42,0.96); border: 1px solid var(--accent); border-radius: 8px; padding: 8px 10px; font-size: 11.5px; color: var(--text); max-width: 300px; z-index: 999; display: none; box-shadow: 0 8px 20px rgba(0,0,0,0.45); }
pre { background: #0b111e; border: 1px solid var(--border); border-radius: 8px; padding: 10px; font-size: 11px; overflow: auto; max-height: 360px; color: #cbd5e1; line-height: 1.4; }
.pager { display: flex; justify-content: space-between; align-items: center; margin-top: 8px; color: var(--muted); font-size: 11.5px; }
</style>
</head>
<body>
<header>
  <div>
    <h1>Parrots Comparative Genomics — Raw Data &amp; Red-Team Audit Explorer</h1>
    <p>Repository: <code>redkk123/Parrots</code> (<code>claude/new-session-lwh7xh</code>) &middot; Offline standalone dashboard with full raw matrices</p>
  </div>
  <div id="top-badges">
    <span class="badge b-blue">126 Assemblies / 71 Validated Spp</span>
    <span class="badge b-green">101 H8 Panel Genes (58 vs 71sp)</span>
    <span class="badge b-amber">7,756 RERconverge Orthologs</span>
    <span class="badge b-red">9 Red-Team Stress Audits</span>
  </div>
</header>

<div class="tabs">
  <button class="tab active" onclick="openTab('tab1', this)">1. Species, Phenotypes &amp; BUSCO Audit</button>
  <button class="tab" onclick="openTab('tab2', this)">2. Phase 2 H8 Panel (58 vs 71sp) &amp; Telomeres</button>
  <button class="tab" onclick="openTab('tab3', this)">3. Phase 3 Genome-Wide (All 7,756 Genes)</button>
  <button class="tab" onclick="openTab('tab4', this)">4. Red-Team Audit &amp; 3D Structures</button>
</div>

<main>
  <!-- TAB 1: SPECIES & AUDIT -->
  <div id="tab1" class="panel active">
    <div class="kpi-row">
      <div class="kpi"><div class="lbl">Phenotype Species</div><div class="val" id="k-sp-total">-</div><div class="sub">Phase 0/0b master matrix</div></div>
      <div class="kpi"><div class="lbl">Validated Genomes</div><div class="val">71</div><div class="sub">58 discovery + 13 holdout</div></div>
      <div class="kpi"><div class="lbl">Iridian Mito Contamination</div><div class="val">61 / 96</div><div class="sub">4 severe nuclear swaps/chimeras</div></div>
      <div class="kpi"><div class="lbl">Quality vs Longevity Bias</div><div class="val">NULL</div><div class="sub">BUSCO p=0.111; N50 p=0.249</div></div>
    </div>
    <div class="grid2">
      <div class="box">
        <h3>Interactive Species Scatter Plot
          <span>
            X: <select id="sp-x" onchange="renderSpeciesPlot()">
              <option value="log_mass">log10 Body Mass (g)</option>
              <option value="rel_log_life_expectancy">Rel. Longevity Residual</option>
              <option value="contig_n50_kb">log10 Contig N50 (kb)</option>
              <option value="busco_C">BUSCO Complete (%)</option>
            </select>
            Y: <select id="sp-y" onchange="renderSpeciesPlot()">
              <option value="log_le">log10 Life Expectancy (yr)</option>
              <option value="rel_log_brain" selected>Rel. Brain Residual</option>
              <option value="log_brain">log10 Brain Mass (g)</option>
              <option value="busco_C">BUSCO Complete (%)</option>
            </select>
          </span>
        </h3>
        <svg id="svg-sp" class="plot" viewBox="0 0 640 330"></svg>
      </div>
      <div class="box">
        <h3>Filters &amp; Summary
          <button onclick="exportCSV(filteredSpecies(), 'parrots_species_raw.csv')">Export Filtered CSV</button>
        </h3>
        <div class="controls">
          <input id="sp-q" placeholder="Search species / genus..." oninput="renderSpeciesTable()">
          <select id="sp-cohort" onchange="renderSpeciesTable(); renderSpeciesPlot()">
            <option value="all">All Cohorts</option>
            <option value="discovery_58">Discovery (58 sp)</option>
            <option value="holdout_13">Holdout (13 sp)</option>
            <option value="excluded">Excluded Assemblies</option>
          </select>
          <select id="sp-clade" onchange="renderSpeciesTable(); renderSpeciesPlot()">
            <option value="all">All Clades</option>
            <option value="Psittaculidae">Psittaculidae</option>
            <option value="other">Other Parrots</option>
          </select>
        </div>
        <div style="color:var(--muted);line-height:1.55;font-size:12px;">
          <b>Color Legend:</b> <span class="badge b-red">Psittaculidae</span> &middot; <span class="badge b-blue">Other Parrots (Included)</span> &middot; <span class="badge b-amber">Excluded / No Genome</span><br><br>
          <b>Key Takeaways from Raw Data:</b><br>
          &bull; <i>Melopsittacus undulatus</i> exhibits 20.8% (AZA) / 25.9% (UC Davis age-std) neoplasia vs 1.5&ndash;5.4% in lories.<br>
          &bull; Dropping the 13 Psittaculidae species with brain data reduces raw tip correlation <i>r</i>(Rel. LE, Rel. Brain) from <b>0.478</b> to <b>0.195</b>.<br>
          &bull; All 71 included species have <code>BUSCO_C &ge; 72.96%</code>.
        </div>
      </div>
    </div>
    <div class="box">
      <h3>Raw Species, Traits &amp; Genomic Audit Table (<span id="sp-count">-</span> rows)</h3>
      <div class="tbl-wrap"><table id="tbl-sp"><thead></thead><tbody></tbody></table></div>
    </div>
  </div>

  <!-- TAB 2: PHASE 2 H8 PANEL -->
  <div id="tab2" class="panel">
    <div class="kpi-row">
      <div class="kpi"><div class="lbl">H8 Panel Genes</div><div class="val">101</div><div class="sub">87 alignable (58sp) / 97 (71sp)</div></div>
      <div class="kpi"><div class="lbl">Confirmed Pseudogenes</div><div class="val">0 / 329</div><div class="sub">100% ref/GC/het artifacts</div></div>
      <div class="kpi"><div class="lbl">M7 Telomere aBSREL (Orig.)</div><div class="val">5 / 8 (62.5%)</div><div class="sub">vs 0/19 FDR in control pool</div></div>
      <div class="kpi"><div class="lbl">M7 After Multi-Hit D+T</div><div class="val">0 / 5</div><div class="sub">RTEL1/TERF2/TERT/TERF1 p&ge;0.97</div></div>
    </div>
    <div class="box">
      <h3>Telomeric (M7) Red-Team Sensitivity Audit &amp; Dual-Reference Orthology (7_telomere_alignment_and_rbh.tsv)</h3>
      <div class="tbl-wrap" style="max-height:240px;"><table id="tbl-tel"><thead></thead><tbody></tbody></table></div>
    </div>
    <div class="box">
      <h3>H8 Panel Gene Explorer (58 vs 71 Species RELAX, aBSREL &amp; Integrity)
        <span>
          <input id="p2-q" placeholder="Filter gene..." oninput="renderPanelTable()">
          <select id="p2-mod" onchange="renderPanelTable()">
            <option value="all">All Modules (M1-M8)</option>
            <option value="M1">M1 (p53/Cell Cycle)</option><option value="M2">M2 (PI3K/mTOR)</option>
            <option value="M3">M3 (Apoptosis)</option><option value="M4">M4 (DNA Repair)</option>
            <option value="M5">M5 (NRF2)</option><option value="M6">M6 (Antioxidant)</option>
            <option value="M7">M7 (Telomere)</option><option value="M8">M8 (Proteostasis)</option>
          </select>
          <button onclick="exportCSV(filteredPanel(), 'parrots_h8_panel_58_vs_71sp.csv')">Export CSV</button>
        </span>
      </h3>
      <div class="tbl-wrap"><table id="tbl-p2"><thead></thead><tbody></tbody></table></div>
    </div>
  </div>

  <!-- TAB 3: PHASE 3 ALL 7,756 GENES -->
  <div id="tab3" class="panel">
    <div class="kpi-row">
      <div class="kpi"><div class="lbl">Tested Orthologs</div><div class="val">7,756</div><div class="sub">58 sp, 1,000 BM perms each</div></div>
      <div class="kpi"><div class="lbl">Raw Coupling r(LE, Brain)</div><div class="val">0.6663</div><div class="sub">r_permQ&lt;0.10 = 0 genes</div></div>
      <div class="kpi"><div class="lbl">Tier A (&gt;80% Zero Edges)</div><div class="val">117 / 131</div><div class="sub">89.3% of shared hits!</div></div>
      <div class="kpi"><div class="lbl">Tier B (&ge;20 Nonzero Edges)</div><div class="val">14 genes</div><div class="sub">TKT, PHOSPHO1, SLC4A7, HUS1...</div></div>
    </div>
    <div class="grid2">
      <div class="box">
        <h3>Genome-Wide Scatter: Rho(Longevity) vs Rho(Brain)
          <span>
            Stratum: <select id="g-tier" onchange="renderGenomePlot(); gPage=0; renderGenomeTable()">
              <option value="all">All 7,756 Genes</option>
              <option value="A_ultraconserved">Tier A: Ultra-conserved (&gt;80% zero branches)</option>
              <option value="B_informative">Tier B: Informative (&ge;20 nonzero &amp; &le;80% zero)</option>
              <option value="shared_sig">Shared Significant (P&lt;0.05 &amp; permP&lt;0.05)</option>
            </select>
          </span>
        </h3>
        <svg id="svg-g" class="plot" viewBox="0 0 640 330"></svg>
      </div>
      <div class="box">
        <h3>True Partial Orthogonalization on Flagship Genes (Test 5 &amp; CV)
          <button onclick="exportCSV(filteredGenome(), 'parrots_rer7756_filtered.csv')">Export Filtered 7,756 CSV</button>
        </h3>
        <div class="controls">
          <input id="g-q" placeholder="Search gene symbol (PTEN, TKT...)" oninput="gPage=0; renderGenomeTable()">
          <label><input type="checkbox" id="g-sigonly" onchange="gPage=0; renderGenomeTable()"> Shared P&amp;permP &lt; 0.05 only</label>
        </div>
        <div class="tbl-wrap" style="max-height:290px;"><table id="tbl-orth-g"><thead></thead><tbody></tbody></table></div>
      </div>
    </div>
    <div class="box">
      <h3>All 7,756 Raw Genes (<span id="g-count">-</span> matching)</h3>
      <div class="tbl-wrap"><table id="tbl-g"><thead></thead><tbody></tbody></table></div>
      <div class="pager"><span id="g-page-info">-</span><span><button onclick="gPage=Math.max(0,gPage-1);renderGenomeTable()">&laquo; Prev 100</button> <button onclick="gPage++;renderGenomeTable()">Next 100 &raquo;</button></span></div>
    </div>
  </div>

  <!-- TAB 4: RED TEAM & 3D -->
  <div id="tab4" class="panel">
    <div class="grid2">
      <div class="box">
        <h3>Leave-One-Clade-Out &amp; Null Coupling Comparison</h3>
        <div class="tbl-wrap" style="max-height:280px;"><table id="tbl-loco"><thead></thead><tbody></tbody></table></div>
      </div>
      <div class="box">
        <h3>True Brain-Beyond-Longevity Orthogonalization Summary</h3>
        <div class="tbl-wrap" style="max-height:280px;"><table id="tbl-orth-s"><thead></thead><tbody></tbody></table></div>
      </div>
    </div>
    <div class="box">
      <h3>AlphaFold 3D Structural Domain Summary (Step 2)</h3>
      <div class="tbl-wrap" style="max-height:220px;"><table id="tbl-3d"><thead></thead><tbody></tbody></table></div>
    </div>
    <div class="box">
      <h3>Complete Statistical Honesty Inventory (78,802 Tests &amp; Null Findings)</h3>
      <div class="tbl-wrap" style="max-height:260px;"><table id="tbl-inv"><thead></thead><tbody></tbody></table></div>
    </div>
    <div class="box">
      <h3>Telomeric Alignment Windows Around <i>M. undulatus</i> Substitutions (7_telomere_substitution_windows.txt)</h3>
      <pre id="pre-win"></pre>
    </div>
  </div>
</main>

<div id="tip" class="tip"></div>

<script>
const D = __DATA_JSON__;
const tip = document.getElementById('tip');
let gPage = 0;

function openTab(id, btn) {
  document.querySelectorAll('.panel').forEach(p => p.classList.remove('active'));
  document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
  document.getElementById(id).classList.add('active');
  btn.classList.add('active');
}

function fmt(v, d=3) {
  if (v === null || v === undefined || v === '') return '—';
  if (typeof v === 'number') return (Math.abs(v) < 1e-3 && v !== 0) ? v.toExponential(2) : v.toFixed(d);
  return String(v);
}

function fillTable(id, rows, cols, sortState) {
  const tbl = document.getElementById(id);
  const thead = tbl.querySelector('thead');
  const tbody = tbl.querySelector('tbody');
  thead.innerHTML = '<tr>' + cols.map(c => `<th data-k="${c.k}">${c.lbl}</th>`).join('') + '</tr>';
  thead.querySelectorAll('th').forEach(th => {
    th.onclick = () => {
      const k = th.dataset.k;
      const asc = (sortState.k === k) ? !sortState.asc : true;
      sortState.k = k; sortState.asc = asc;
      rows.sort((a,b) => {
        const va = a[k], vb = b[k];
        if (va == null) return 1; if (vb == null) return -1;
        return (typeof va === 'number' && typeof vb === 'number') ? (asc ? va-vb : vb-va) : (asc ? String(va).localeCompare(String(vb)) : String(vb).localeCompare(String(va)));
      });
      fillTable(id, rows, cols, sortState);
    };
  });
  tbody.innerHTML = rows.map(r => '<tr>' + cols.map(c => `<td>${c.fn ? c.fn(r[c.k], r) : fmt(r[c.k], c.d ?? 3)}</td>`).join('') + '</tr>').join('');
}

function showTip(evt, html) {
  tip.innerHTML = html; tip.style.display = 'block';
  tip.style.left = Math.min(window.innerWidth - 320, evt.clientX + 12) + 'px';
  tip.style.top = Math.min(window.innerHeight - 160, evt.clientY + 12) + 'px';
}
function hideTip() { tip.style.display = 'none'; }

function exportCSV(rows, fname) {
  if (!rows.length) return;
  const keys = Object.keys(rows[0]);
  const csv = [keys.join(','), ...rows.map(r => keys.map(k => `"${String(r[k]??'').replace(/"/g,'""')}"`).join(','))].join('\n');
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([csv], {type:'text/csv'}));
  a.download = fname; a.click();
}

// TAB 1
function getVal(r, k) {
  if (k === 'log_mass') return r.body_mass_g ? Math.log10(r.body_mass_g) : null;
  if (k === 'log_le') return r.life_expectancy_years ? Math.log10(r.life_expectancy_years) : null;
  if (k === 'log_brain') return r.brain_mass_g ? Math.log10(r.brain_mass_g) : null;
  if (k === 'contig_n50_kb') return r.contig_n50 ? Math.log10(r.contig_n50 / 1000) : null;
  return r[k];
}
function filteredSpecies() {
  const q = document.getElementById('sp-q').value.toLowerCase();
  const coh = document.getElementById('sp-cohort').value;
  const cl = document.getElementById('sp-clade').value;
  return D.species.filter(r => {
    if (q && !(r.species||'').toLowerCase().includes(q) && !(r.family||'').toLowerCase().includes(q)) return false;
    if (coh !== 'all' && r.cohort !== coh) return false;
    if (cl === 'Psittaculidae' && !r.psittaculidae) return false;
    if (cl === 'other' && r.psittaculidae) return false;
    return true;
  });
}
function renderSpeciesPlot() {
  const xk = document.getElementById('sp-x').value, yk = document.getElementById('sp-y').value;
  const rows = filteredSpecies().map(r => ({r, x:getVal(r,xk), y:getVal(r,yk)})).filter(p => p.x!=null && p.y!=null && isFinite(p.x) && isFinite(p.y));
  const svg = document.getElementById('svg-sp');
  if (!rows.length) { svg.innerHTML = ''; return; }
  const xs = rows.map(p=>p.x), ys = rows.map(p=>p.y);
  const xmin = Math.min(...xs), xmax = Math.max(...xs), ymin = Math.min(...ys), ymax = Math.max(...ys);
  const pad = {l:48, r:16, t:16, b:34}, W=640, H=330;
  const sx = x => pad.l + (x - xmin) / (xmax - xmin || 1) * (W - pad.l - pad.r);
  const sy = y => H - pad.b - (y - ymin) / (ymax - ymin || 1) * (H - pad.t - pad.b);
  let h = `<line x1="${pad.l}" y1="${H-pad.b}" x2="${W-pad.r}" y2="${H-pad.b}" stroke="#334155"/><line x1="${pad.l}" y1="${pad.t}" x2="${pad.l}" y2="${H-pad.b}" stroke="#334155"/>`;
  h += `<text x="${W/2}" y="${H-8}" fill="#93a4c3" font-size="11" text-anchor="middle">${xk} (n=${rows.length})</text>`;
  h += `<text x="${pad.l}" y="${12}" fill="#93a4c3" font-size="10">${ymin.toFixed(2)} .. ${ymax.toFixed(2)} (${yk})</text>`;
  rows.forEach((p, idx) => {
    const col = p.r.psittaculidae ? '#f87171' : (p.r.cohort==='discovery_58'||p.r.cohort==='holdout_13' ? '#38bdf8' : '#fbbf24');
    const rad = p.r.species === 'Melopsittacus undulatus' ? 6.5 : 3.8;
    h += `<circle cx="${sx(p.x).toFixed(1)}" cy="${sy(p.y).toFixed(1)}" r="${rad}" fill="${col}" fill-opacity="0.82" stroke="#0b101b" stroke-width="0.8" data-i="${idx}"/>`;
  });
  svg.innerHTML = h;
  svg.querySelectorAll('circle').forEach(c => {
    const p = rows[+c.dataset.i].r;
    c.onmousemove = e => showTip(e, `<b><i>${p.species}</i></b> (${p.cohort})<br>Mass: ${fmt(p.body_mass_g,1)}g | LE: ${fmt(p.life_expectancy_years,1)}yr | Brain: ${fmt(p.brain_mass_g,2)}g<br>BUSCO: ${fmt(p.busco_C,1)}% | N50: ${p.contig_n50?fmt(p.contig_n50/1e3,1)+'kb':'—'}`);
    c.onmouseleave = hideTip;
  });
}
const spCols = [
  {k:'species', lbl:'Species', fn:v=>`<i>${v}</i>`},
  {k:'cohort', lbl:'Cohort'},
  {k:'psittaculidae', lbl:'Psitt.', fn:v=>v===true?'Yes':(v===false?'No':'—')},
  {k:'life_expectancy_years', lbl:'LE (yr)', d:1},
  {k:'body_mass_g', lbl:'Mass (g)', d:1},
  {k:'brain_mass_g', lbl:'Brain (g)', d:2},
  {k:'rel_log_life_expectancy', lbl:'Rel. LE', d:3},
  {k:'rel_log_brain', lbl:'Rel. Brain', d:3},
  {k:'accession', lbl:'Accession'},
  {k:'busco_C', lbl:'BUSCO%', d:1},
  {k:'contig_n50', lbl:'N50 (bp)', d:0},
  {k:'mito_overall', lbl:'Mito Audit'},
];
let spSort = {k:'species', asc:true};
function renderSpeciesTable() {
  const rows = filteredSpecies();
  document.getElementById('k-sp-total').textContent = D.species.length;
  document.getElementById('sp-count').textContent = rows.length;
  fillTable('tbl-sp', rows, spCols, spSort);
}

// TAB 2
function filteredPanel() {
  const q = document.getElementById('p2-q').value.toLowerCase();
  const m = document.getElementById('p2-mod').value;
  return D.panel.filter(r => (!q || r.gene.toLowerCase().includes(q)) && (m==='all' || r.module===m));
}
let p2Sort = {k:'absrel_melopsittacus_p', asc:true};
function renderPanelTable() {
  const cols = [
    {k:'gene', lbl:'Gene', fn:v=>`<b>${v}</b>`},
    {k:'module', lbl:'Mod.'},
    {k:'k', lbl:'k (58sp)', d:2},
    {k:'relax_q', lbl:'RELAX q (58)', d:4},
    {k:'absrel_melopsittacus_p', lbl:'aBSREL p (58)', d:4},
    {k:'k_71sp', lbl:'k (71sp)', d:2},
    {k:'relax_p_71sp', lbl:'RELAX p (71)', d:4},
    {k:'absrel_melo_p_71sp', lbl:'aBSREL p (71)', d:4},
    {k:'n_sp_71', lbl:'Spp (71)', d:0},
    {k:'disrupted_psittaculidae_71sp', lbl:'Disr. Psl', d:0},
    {k:'complete_psittaculidae_71sp', lbl:'Comp. Psl', d:0},
  ];
  fillTable('tbl-p2', filteredPanel(), cols, p2Sort);
}

// TAB 3: ALL 7,756 GENES
function filteredGenome() {
  const q = document.getElementById('g-q').value.toLowerCase();
  const tier = document.getElementById('g-tier').value;
  const sigOnly = document.getElementById('g-sigonly').checked;
  return D.genome.filter(r => {
    if (q && !r.gene.toLowerCase().includes(q)) return false;
    if (tier === 'A_ultraconserved' && r.tier !== 'A_ultraconserved') return false;
    if (tier === 'B_informative' && r.tier !== 'B_informative') return false;
    const isSig = r.P_LE<0.05 && r.permP_LE<0.05 && r.P_Br<0.05 && r.permP_Br<0.05 && r.Rho_LE*r.Rho_Br>0;
    if ((tier === 'shared_sig' || sigOnly) && !isSig) return false;
    return true;
  });
}
function renderGenomePlot() {
  const rows = filteredGenome();
  const svg = document.getElementById('svg-g');
  const W=640, H=330, pad={l:46,r:16,t:16,b:34};
  const sx = x => pad.l + (x + 0.65) / 1.3 * (W - pad.l - pad.r);
  const sy = y => H - pad.b - (y + 0.65) / 1.3 * (H - pad.t - pad.b);
  let h = `<line x1="${sx(0)}" y1="${pad.t}" x2="${sx(0)}" y2="${H-pad.b}" stroke="#334155" stroke-dasharray="3 3"/><line x1="${pad.l}" y1="${sy(0)}" x2="${W-pad.r}" y2="${sy(0)}" stroke="#334155" stroke-dasharray="3 3"/>`;
  // downsample to max 2500 points for crisp SVG responsiveness, always keeping all highlighted/sig hits
  const step = Math.max(1, Math.floor(rows.length / 2200));
  const drawn = [];
  rows.forEach((r, idx) => {
    const isKey = ['PTEN','MAD2L2','HSPA8','CTNNB1','TKT','PHOSPHO1','HUS1','SLC4A7'].includes(r.gene) || (r.P_LE<0.01 && r.P_Br<0.01);
    if (idx % step === 0 || isKey) drawn.push(r);
  });
  drawn.forEach((r, idx) => {
    const col = r.tier === 'B_informative' ? '#34d399' : '#f87171';
    const isTop = ['PTEN','MAD2L2','HSPA8','CTNNB1','TKT','HUS1'].includes(r.gene);
    h += `<circle cx="${sx(r.Rho_LE).toFixed(1)}" cy="${sy(r.Rho_Br).toFixed(1)}" r="${isTop?5.5:2.3}" fill="${isTop?'#fbbf24':col}" fill-opacity="${isTop?1:0.45}" data-i="${idx}"/>`;
    if (isTop) h += `<text x="${(sx(r.Rho_LE)+7).toFixed(1)}" y="${(sy(r.Rho_Br)+3).toFixed(1)}" fill="#fde68a" font-size="10" font-weight="700">${r.gene}</text>`;
  });
  h += `<text x="${W/2}" y="${H-8}" fill="#93a4c3" font-size="11" text-anchor="middle">Rho(Life Expectancy) [Green = Tier B Informative, Red = Tier A >80% Zero Branches]</text>`;
  svg.innerHTML = h;
  svg.querySelectorAll('circle').forEach(c => {
    const g = drawn[+c.dataset.i];
    c.onmousemove = e => showTip(e, `<b>${g.gene}</b> (${g.tier})<br>Rho_LE=${fmt(g.Rho_LE)} (permP=${fmt(g.permP_LE,4)})<br>Rho_Br=${fmt(g.Rho_Br)} (permP=${fmt(g.permP_Br,4)})<br>Zero branches: ${g.zero_pct}% (${g.n_nz} nonzero)`);
    c.onmouseleave = hideTip;
  });
}
let gSort = {k:'P_LE', asc:true};
function renderGenomeTable() {
  const rows = filteredGenome();
  document.getElementById('g-count').textContent = rows.length;
  const maxP = Math.max(0, Math.ceil(rows.length / 100) - 1);
  if (gPage > maxP) gPage = maxP;
  const slice = rows.slice(gPage * 100, gPage * 100 + 100);
  document.getElementById('g-page-info').textContent = `Showing ${gPage*100+1}-${Math.min(rows.length,(gPage+1)*100)} of ${rows.length}`;
  const cols = [
    {k:'gene', lbl:'Gene', fn:v=>`<b>${v}</b>`},
    {k:'tier', lbl:'Tier', fn:v=>`<span class="badge ${v==='B_informative'?'b-green':'b-red'}">${v==='B_informative'?'Tier B':'Tier A'}</span>`},
    {k:'zero_pct', lbl:'%Zero Br.', d:1},
    {k:'n_nz', lbl:'Nonzero', d:0},
    {k:'Rho_LE', lbl:'Rho_LE', d:3},
    {k:'P_LE', lbl:'P_LE', d:4},
    {k:'permP_LE', lbl:'permP_LE', d:4},
    {k:'Rho_Br', lbl:'Rho_Br', d:3},
    {k:'P_Br', lbl:'P_Br', d:4},
    {k:'permP_Br', lbl:'permP_Br', d:4},
    {k:'Rho_BrPGLS', lbl:'Rho_Br|LE(old)', d:3},
    {k:'tree_len', lbl:'Tree Len', d:4},
  ];
  fillTable('tbl-g', slice, cols, gSort);
}

// INIT ALL TABLES
renderSpeciesTable(); renderSpeciesPlot();
renderPanelTable();
renderGenomeTable(); renderGenomePlot();

let telSort = {k:'absrel_p_original_mafft', asc:true};
fillTable('tbl-tel', D.telomere_rt, [
  {k:'gene', lbl:'Gene', fn:v=>`<b>${v}</b>`},
  {k:'n_species', lbl:'Spp', d:0},
  {k:'melo_aa_subs_total', lbl:'Melo Subs', d:0},
  {k:'melo_aa_subs_clustered_or_gappy', lbl:'Noisy/Gap', d:0},
  {k:'absrel_p_original_mafft', lbl:'Orig p', d:4},
  {k:'absrel_p_linsi_realigned', lbl:'L-INS-i p', d:4},
  {k:'absrel_p_window_cleaned', lbl:'Cleaned p', d:4},
  {k:'absrel_p_MultiHit_DoubleTriple', lbl:'MultiHit p', d:4},
  {k:'dual_ref_kakapo_concordance', lbl:'Kakapo Orth.'},
], telSort);

let ogSort = {k:'P_Brain_branch_orth', asc:true};
fillTable('tbl-orth-g', D.rt_orth_genes, [
  {k:'gene', lbl:'Gene', fn:v=>`<b>${v}</b>`},
  {k:'Rho_LE_58sp', lbl:'Rho_LE', d:3},
  {k:'Rho_Brain_raw', lbl:'Rho_Br', d:3},
  {k:'Rho_Brain_OLS_tip_orth', lbl:'OLS_Br', d:3},
  {k:'permP_Brain_OLS_tip_orth', lbl:'OLS permP', d:3},
  {k:'Rho_Brain_branch_orth', lbl:'BranchOrth Br', d:3},
  {k:'P_Brain_branch_orth', lbl:'BranchOrth P', d:4},
], ogSort);

let lSort = {k:'genome_Rho_cor', asc:false};
fillTable('tbl-loco', D.rt_loco, [
  {k:'subset', lbl:'Clade Subset'},
  {k:'n_species', lbl:'N', d:0},
  {k:'tip_cor_pearson', lbl:'r_tip', d:3},
  {k:'branch_path_cor', lbl:'r_path', d:3},
  {k:'genome_Rho_cor', lbl:'r(Rho)', d:3},
  {k:'PTEN_Rho_LE', lbl:'PTEN LE', d:3},
  {k:'MAD2L2_Rho_LE', lbl:'MAD2L2 LE', d:3},
], lSort);

let oSort = {k:'genome_Rho_cor_with_LE_58sp', asc:false};
fillTable('tbl-orth-s', D.rt_orth, [
  {k:'method', lbl:'Orthogonalization Method'},
  {k:'tip_cor_with_LE', lbl:'r_tip(LE)', d:3},
  {k:'branch_path_cor_with_LE', lbl:'r_path(LE)', d:3},
  {k:'genome_Rho_cor_with_LE_58sp', lbl:'r(Rho)_58', d:3},
  {k:'genome_Rho_cor_with_LE_50sp', lbl:'r(Rho)_50', d:3},
], oSort);

let s3Sort = {k:'gene', asc:true};
fillTable('tbl-3d', D.struct_3d, [
  {k:'gene', lbl:'Protein', fn:v=>`<b>${v}</b>`},
  {k:'uniprot_id', lbl:'UniProt'},
  {k:'length_aa', lbl:'Len', d:0},
  {k:'global_plddt', lbl:'pLDDT', d:2},
  {k:'domain_boundaries', lbl:'Domains'},
  {k:'structural_finding', lbl:'3D Finding'},
], s3Sort);

let iSort = {k:'n_tests_run', asc:false};
fillTable('tbl-inv', D.rt_inv, [
  {k:'phase', lbl:'Phase'},
  {k:'n_tests_run', lbl:'Tests', d:0},
  {k:'significant_findings', lbl:'Positive Signals'},
  {k:'null_or_negative_findings', lbl:'Null / Negative Findings'},
], iSort);

document.getElementById('pre-win').textContent = D.tel_windows_txt;
</script>
</body>
</html>"""


if __name__ == "__main__":
    main()
