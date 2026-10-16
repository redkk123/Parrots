"""Cancer axis: tumour prevalence in psittacine necropsies, Drury R. Reavill Pathology Database (UC Davis).

Public database (https://drury-reavill-database.vercel.app/, 2026); only aggregate COUNTS are queried
(/api/records with pageSize=1 returns the total for a filter), no records are downloaded. Requests are cached in
data/external/drury_reavill/counts_cache.json and spaced by 0.5 s.

Necropsies: charge_type "Post Mortem ..." or "CalAvian Necropsy". Tumour: full-text match (whole words, OR) on the
diagnosis text for named neoplasms; generic words ("neoplasia", "tumor") are excluded because they also appear in
negative statements. Malignant = carcinomas, sarcomas, lymphoma/leukaemia, melanoma, seminoma, mesothelioma,
nephroblastoma; benign = adenomas, papillomas, lipomas, fibromas, haemangiomas etc. Groups pool the database's
"breed" values by common name within category "Avian, Psittacine".
Age bins (years): 0-2, 2-5, 5-10, 10-20, >=20. Age-standardised prevalence uses the pooled psittacine necropsy age
distribution as the standard (direct standardisation).
Writes results/phase1_cancer/drury_reavill_necropsy.tsv
Cite: Drury R. Reavill Pathology Database at UC Davis. Davis (CA): University of California, Davis, School of
Veterinary Medicine; 2026.
"""
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DIR = ROOT / "data/external/drury_reavill"
CACHE = DIR / "counts_cache.json"
OUT = ROOT / "results/phase1_cancer/drury_reavill_necropsy.tsv"
API = "https://drury-reavill-database.vercel.app/api/records?"
PM = ["Post Mortem, Jar", "Post mortem, Medium", "Post Mortem, Small", "Post Mortem Sm Jar", "Post Mortem, small",
      "Post Mortem, Large", "CalAvian Necropsy", "Post Mortem Jar-Lab", "Post Mortem Jar-Very",
      "Post Mortem Sm Jarsk", "Post Mortem, Small1)"]
MALIGNANT = ["carcinoma", "carcinomas", "adenocarcinoma", "adenocarcinomas", "cystadenocarcinoma",
             "cholangiocarcinoma", "carcinosarcoma", "sarcoma", "sarcomas", "fibrosarcoma", "liposarcoma",
             "hemangiosarcoma", "leiomyosarcoma", "osteosarcoma", "rhabdomyosarcoma", "chondrosarcoma",
             "lymphosarcoma", "lymphoma", "lymphomas", "leukemia", "leukaemia", "melanoma", "seminoma",
             "mesothelioma", "nephroblastoma", "dysgerminoma", "thymoma", "myeloproliferative"]
BENIGN = ["adenoma", "adenomas", "cystadenoma", "cholangioma", "papilloma", "papillomas", "lipoma", "lipomas",
          "myelolipoma", "hemangioma", "fibroma", "leiomyoma", "osteoma", "chondroma", "hepatoma"]
GROUPS = {  # common-name prefixes of the database's "breed" values
    "Budgerigar": ["Budgerigar"], "Cockatiel": ["Cockatiel"], "Lovebird": ["Lovebird"],
    "African grey": ["Parrot, African Grey", "Parrot, Congo African Grey", "Parrot, Timneh Grey"],
    "Amazon": ["Amazon"], "Macaw": ["Macaw"], "Cockatoo": ["Cockatoo"], "Eclectus": ["Parrot, Eclectus"],
    "Conure": ["Conure"], "Lory/lorikeet": ["Lory", "Lorikeet"],
    "Ringneck/Psittacula": ["Parakeet, Ringneck", "Parakeet, Moustache"],
}
AGE_BINS = [(0, 2), (2, 5), (5, 10), (10, 20), (20, 200)]


def load_cache():
    return json.loads(CACHE.read_text()) if CACHE.exists() else {}


def count(cache, **params):
    p = {"dataset": "main", "category": "Avian, Psittacine", "charge_type": "||".join(PM), **params,
         "page": "1", "pageSize": "1"}
    url = API + urllib.parse.urlencode(sorted(p.items()))
    if url not in cache:
        req = urllib.request.Request(url, headers={"User-Agent": "parrot-longevity-study (research; counts only)"})
        cache[url] = json.load(urllib.request.urlopen(req, timeout=60))["total"]
        CACHE.write_text(json.dumps(cache))
        time.sleep(0.5)
    return cache[url]


def breeds_for(prefixes):
    meta = json.loads((DIR / "meta.json").read_text())["facets"]["breed"]
    return [b["value"] for b in meta if any(b["value"] == p or b["value"].startswith(p + ",")
                                            or b["value"].startswith(p + " ") for p in prefixes)]


def cp(k, n):
    if n == 0:
        return np.nan, np.nan
    return (stats.beta.ppf(0.025, k, n - k + 1) if k else 0.0, stats.beta.ppf(0.975, k + 1, n - k) if k < n else 1.0)


def main():
    cache = load_cache()
    q_any = " OR ".join(MALIGNANT + BENIGN)
    q_mal = " OR ".join(MALIGNANT)
    rows = []
    groups = {g: breeds_for(p) for g, p in GROUPS.items()}
    groups["All psittacines"] = []  # no breed filter
    for g, br in groups.items():
        f = {"breed": "||".join(br)} if br else {}
        r = {"group": g, "n_breed_values": len(br)}
        r["necropsies"] = count(cache, **f)
        r["any_tumour"] = count(cache, q=q_any, cols="diagnosis", **f)
        r["malignant"] = count(cache, q=q_mal, cols="diagnosis", **f)
        r["classified_tumour_primary"] = count(cache, specific_lesions="Tumor", **f)
        for lo, hi in AGE_BINS:
            a = {"age_min": str(lo), "age_max": str(hi - 0.001)}
            r[f"n_{lo}_{hi}"] = count(cache, **f, **a)
            r[f"t_{lo}_{hi}"] = count(cache, q=q_any, cols="diagnosis", **f, **a)
        rows.append(r)
    t = pd.DataFrame(rows)
    bins = [f"{lo}_{hi}" for lo, hi in AGE_BINS]
    std = t.set_index("group").loc["All psittacines", [f"n_{b}" for b in bins]].to_numpy(float)
    w = std / std.sum()
    for c in ("any_tumour", "malignant"):
        t[f"pct_{c}"] = 100 * t[c] / t.necropsies
        ci = [cp(k, n) for k, n in zip(t[c], t.necropsies)]
        t[f"ci_low_{c}"] = [100 * x[0] for x in ci]
        t[f"ci_high_{c}"] = [100 * x[1] for x in ci]
    t["n_with_age"] = t[[f"n_{b}" for b in bins]].sum(axis=1)

    def age_std(r):
        ok = np.array([r[f"n_{b}"] >= 5 for b in bins])
        rate = np.array([r[f"t_{b}"] / r[f"n_{b}"] if r[f"n_{b}"] else np.nan for b in bins])
        return 100 * np.sum(w[ok] * rate[ok]) / w[ok].sum() if ok.any() else np.nan

    t["age_std_pct_any_tumour"] = [age_std(r) for _, r in t.iterrows()]
    t["age_bins_used"] = [sum(r[f"n_{b}"] >= 5 for b in bins) for _, r in t.iterrows()]
    t["median_age_bin"] = [bins[int(np.searchsorted(np.cumsum([r[f"n_{b}"] for b in bins]),
                                                      sum(r[f"n_{b}"] for b in bins) / 2))]
                           for _, r in t.iterrows()]
    bud = t.set_index("group").loc["Budgerigar"]
    t["fisher_p_vs_budgerigar"] = [stats.fisher_exact([[bud.any_tumour, bud.necropsies - bud.any_tumour],
                                                       [k, n - k]])[1] for k, n in zip(t.any_tumour, t.necropsies)]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    t.to_csv(OUT, sep="\t", index=False, float_format="%.4g")
    show = ["group", "necropsies", "any_tumour", "pct_any_tumour", "ci_low_any_tumour", "ci_high_any_tumour",
            "pct_malignant", "age_std_pct_any_tumour", "median_age_bin", "fisher_p_vs_budgerigar"]
    print(t[show].sort_values("pct_any_tumour", ascending=False).to_string(index=False,
                                                                          float_format=lambda v: f"{v:.3g}"))
    print("requests cached:", len(cache))


if __name__ == "__main__":
    main()
