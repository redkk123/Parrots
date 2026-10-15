"""Cancer axis: budgerigar vs Loriini vs other parrots in zoo necropsies (Kapsetaki et al. 2024).

Source: Zenodo 10.5281/zenodo.11311870 (Supplementary_data_v2.csv, CC BY 4.0), 5,729 adult necropsies of
108 bird species from 25 AZA institutions; species with >= 20 necropsies. Counts are reconstructed as
prevalence x necropsies (exact to < 1e-6). Groups follow the project's taxonomy (genera), not the dataset's
inconsistent family column.

Order fixed in the protocol: Melopsittacus vs Loriini first, then the rest of Psittaculidae, with Nymphicus
(also domesticated) as the husbandry control.
  1. pooled counts: Clopper-Pearson 95 % CI and Fisher exact test vs the budgerigar;
  2. between-species heterogeneity: beta-binomial fitted to the comparator species; probability of a species
     with the budgerigar's necropsy count showing >= its number of cases (one-sided);
  3. domestication: domesticated vs non-domesticated parrots, excluding the budgerigar.
Writes results/phase1_cancer/necropsy_parrots.tsv and results/phase1_cancer/necropsy_contrasts.tsv.
No ages per species are available here, so age at death cannot be modelled (see report).
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, special, stats

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data/external/kapsetaki2024/Supplementary_data_v2.csv"
OUT = ROOT / "results/phase1_cancer"
LORIINI = {"Trichoglossus", "Lorius", "Eos", "Vini", "Chalcopsitta", "Pseudeos", "Charmosyna", "Glossopsitta",
           "Psitteuteles", "Neopsittacus", "Oreopsittacus", "Phigys"}
PSITTACULIDAE_OTHER = {"Agapornis", "Neophema", "Loriculus", "Psittacula", "Platycercus", "Eclectus", "Lathamus",
                       "Cyanoramphus", "Polytelis", "Alisterus", "Aprosmictus", "Tanygnathus", "Psephotus"}


def group(genus):
    if genus == "Melopsittacus":
        return "Melopsittacus"
    if genus in LORIINI:
        return "Loriini"
    if genus in PSITTACULIDAE_OTHER:
        return "Psittaculidae_other"
    return "non_Psittaculidae"


def cp(k, n):
    lo = stats.beta.ppf(0.025, k, n - k + 1) if k > 0 else 0.0
    hi = stats.beta.ppf(0.975, k + 1, n - k) if k < n else 1.0
    return lo, hi


def fit_betabinom(k, n):
    """ML beta-binomial (mean mu, overdispersion rho = 1/(a+b+1)) for species-level counts."""
    k, n = np.asarray(k), np.asarray(n)

    def nll(t):
        mu, rho = 1 / (1 + np.exp(-t[0])), 1 / (1 + np.exp(-t[1]))
        s = 1 / rho - 1
        return -stats.betabinom.logpmf(k, n, mu * s, (1 - mu) * s).sum()

    r = optimize.minimize(nll, x0=[np.log(k.sum() / (n.sum() - k.sum()) + 1e-9), -3], method="Nelder-Mead",
                          options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 5000})
    mu, rho = 1 / (1 + np.exp(-r.x[0])), 1 / (1 + np.exp(-r.x[1]))
    return mu, rho


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(SRC)
    p = d[d.Orders == "Psittaciformes"].copy()
    p["group"] = p.Genus.map(group)
    p["domesticated"] = p["Domesticated.or.Semidomesticated"].eq("yes")
    rows = []
    for c in ("Neoplasia", "Malignancy"):
        p[f"n_{c.lower()}"] = (p[f"{c}Prevalence"] * p.Necropsies).round().astype(int)
    sp = p[["Species", "Genus", "group", "domesticated", "Necropsies", "n_neoplasia", "n_malignancy",
            "NeoplasiaPrevalence", "MalignancyPrevalence", "Maximum_lifespan"]].sort_values(["group", "Species"])
    sp.to_csv(OUT / "necropsy_parrots.tsv", sep="\t", index=False, float_format="%.4g")

    bud = p[p.Genus == "Melopsittacus"].iloc[0]
    comparisons = {
        "Loriini": p.group == "Loriini",
        "Psittaculidae_other (no Loriini, no budgerigar)": p.group == "Psittaculidae_other",
        "non_Psittaculidae": p.group == "non_Psittaculidae",
        "all other parrots": p.Genus != "Melopsittacus",
        "Nymphicus (domesticated control)": p.Genus == "Nymphicus",
        "non-domesticated parrots": (~p.domesticated) & (p.Genus != "Melopsittacus"),
    }
    for c in ("neoplasia", "malignancy"):
        kb, nb = int(bud[f"n_{c}"]), int(bud.Necropsies)
        for name, mask in comparisons.items():
            g = p[mask]
            k, n = int(g[f"n_{c}"].sum()), int(g.Necropsies.sum())
            odds, pf = stats.fisher_exact([[kb, nb - kb], [k, n - k]], alternative="two-sided")
            r = {"outcome": c, "comparator": name, "n_species": len(g), "necropsies": n, "cases": k,
                 "prevalence": k / n, "ci_low": cp(k, n)[0], "ci_high": cp(k, n)[1],
                 "budgerigar_cases": kb, "budgerigar_necropsies": nb, "budgerigar_prevalence": kb / nb,
                 "fisher_odds_ratio": odds, "fisher_p": pf}
            if len(g) >= 3:
                mu, rho = fit_betabinom(g[f"n_{c}"], g.Necropsies)
                s = 1 / rho - 1
                xs = np.arange(kb, nb + 1)  # upper tail summed on the log scale (sf underflows for Loriini)
                lp = float(special.logsumexp(stats.betabinom.logpmf(xs, nb, mu * s, (1 - mu) * s)))
                r.update(betabinom_mu=mu, betabinom_rho=rho, log10_p_budgerigar_from_group=lp / np.log(10))
            rows.append(r)
        # domestication among the other parrots
        g = p[p.Genus != "Melopsittacus"]
        dom, non = g[g.domesticated], g[~g.domesticated]
        a = [int(dom[f"n_{c}"].sum()), int(dom.Necropsies.sum())]
        b = [int(non[f"n_{c}"].sum()), int(non.Necropsies.sum())]
        odds, pf = stats.fisher_exact([[a[0], a[1] - a[0]], [b[0], b[1] - b[0]]])
        rows.append({"outcome": c, "comparator": "domesticated vs non-domesticated (no budgerigar)",
                     "n_species": len(dom), "necropsies": a[1], "cases": a[0], "prevalence": a[0] / a[1],
                     "ci_low": cp(*a)[0], "ci_high": cp(*a)[1], "budgerigar_cases": b[0],
                     "budgerigar_necropsies": b[1], "budgerigar_prevalence": b[0] / b[1],
                     "fisher_odds_ratio": odds, "fisher_p": pf,
                     "note": "here the 'budgerigar_*' columns hold the non-domesticated group"})
    r = pd.DataFrame(rows)
    r.to_csv(OUT / "necropsy_contrasts.tsv", sep="\t", index=False, float_format="%.4g")
    print(sp.to_string(index=False))
    print(r.drop(columns=[c for c in ("note",) if c in r]).to_string(index=False, float_format=lambda v: f"{v:.3g}"))


if __name__ == "__main__":
    main()
