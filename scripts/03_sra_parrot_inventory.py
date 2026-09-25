"""Inventory of public parrot sequencing data from the SRA metadata snapshot on AWS Open Data.

NCBI Datasets/E-utilities are not reachable from this environment, but NCBI publishes the
full SRA run metadata as Parquet in s3://sra-pub-metadata-us-east-1/sra/metadata/.
Only the `organism` column is scanned for every row group; the remaining columns are
read only for row groups that contain Psittaciformes records.

Outputs (data/processed/sra/):
  - sra_parrot_runs.tsv: one row per run
  - sra_parrot_species.tsv: per-species summary of whole-genome (WGS) data
"""
import re
import subprocess
import sys
import urllib.parse
from concurrent.futures import ThreadPoolExecutor
from datetime import date
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

sys.path.insert(0, str(Path(__file__).parent / "lib"))
from http_parquet import HTTPRangeFile  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data/processed/sra"
OUT.mkdir(parents=True, exist_ok=True)
BUCKET = "https://sra-pub-metadata-us-east-1.s3.amazonaws.com"
GENOME_GB = 1.2  # typical parrot genome size (Gb), used only for a rough coverage estimate
COLS = ["acc", "assay_type", "center_name", "instrument", "librarylayout", "libraryselection", "librarysource",
        "platform", "biosample", "organism", "sra_study", "bioproject", "releasedate", "mbases", "avgspotlen"]

# Genera recognised in recent taxonomies but absent from the Burgio et al. 2019 table.
EXTRA_GENERA = {"Palaeornis", "Alexandrinus", "Himalayapsitta", "Belocercus", "Nicopsitta", "Saudareos",
                "Glossoptilus", "Charminetta", "Hypocharmosyna", "Charmosynopsis", "Synorhacma", "Parvipsitta",
                "Thectocercus", "Pyrilia", "Lophochroa", "Zanda", "Psittacara", "Eupsittula", "Mascarinus",
                "Conuropsis", "Ognorhynchus", "Leptosittaca", "Bolbopsittacus", "Eunymphicus"}


def list_keys():
    keys, token = [], None
    while True:
        q = {"list-type": "2", "prefix": "sra/metadata/", "max-keys": "1000"}
        if token:
            q["continuation-token"] = token
        xml = subprocess.run(["curl", "-s", f"{BUCKET}/?{urllib.parse.urlencode(q)}"],
                             capture_output=True, text=True, check=True).stdout
        keys += re.findall(r"<Key>([^<]+)</Key>", xml)
        m = re.search(r"<NextContinuationToken>([^<]+)</NextContinuationToken>", xml)
        if not m:
            return [k for k in keys if not k.endswith("/")]
        token = m.group(1)


def parrot_genera():
    tax = pd.read_csv(ROOT / "data/external/parrottraits/data-raw/parrot-traits-database.csv", header=1,
                      encoding="utf-8-sig")
    return set(tax["genus"].dropna().str.strip()) | EXTRA_GENERA


def scan(key, genera_re):
    f = HTTPRangeFile(f"{BUCKET}/{key}")
    pf = pq.ParquetFile(f)
    hits = []
    for rg in range(pf.num_row_groups):
        org = pf.read_row_group(rg, columns=["organism"]).column("organism")
        mask = pc.match_substring_regex(org, genera_re)
        mask = pc.fill_null(mask, False)
        if pc.any(mask).as_py():
            t = pf.read_row_group(rg, columns=COLS).filter(mask)
            hits.append(t)
    return key, (pa.concat_tables(hits) if hits else None), f.bytes_fetched


def main():
    keys = list_keys()
    genera = parrot_genera()
    genera_re = r"^(" + "|".join(sorted(genera)) + r") [a-z]"
    print(f"{len(keys)} parquet files; {len(genera)} parrot genera", flush=True)
    tables, fetched = [], 0
    with ThreadPoolExecutor(6) as ex:
        for key, t, nbytes in ex.map(lambda k: scan(k, genera_re), keys):
            fetched += nbytes
            n = 0 if t is None else t.num_rows
            print(f"  {key.split('/')[-1][:40]}  hits={n}  MB={nbytes / 1e6:.1f}", flush=True)
            if t is not None:
                tables.append(t)
    runs = pa.concat_tables(tables).to_pandas()
    runs = runs.drop_duplicates("acc").sort_values(["organism", "releasedate", "acc"])
    runs.to_csv(OUT / "sra_parrot_runs.tsv", sep="\t", index=False)
    print(f"total fetched: {fetched / 1e6:.0f} MB; parrot runs: {len(runs)}", flush=True)

    wgs = runs[(runs["assay_type"] == "WGS") & (runs["librarysource"] == "GENOMIC")].copy()
    wgs["long_read"] = wgs["platform"].isin(["PACBIO_SMRT", "OXFORD_NANOPORE"])
    g = wgs.groupby("organism")
    summ = pd.DataFrame({
        "n_wgs_runs": g.size(),
        "n_biosamples": g["biosample"].nunique(),
        "n_bioprojects": g["bioproject"].nunique(),
        "total_gbases": g["mbases"].sum() / 1000,
        "short_read_gbases": wgs[~wgs["long_read"]].groupby("organism")["mbases"].sum() / 1000,
        "long_read_gbases": wgs[wgs["long_read"]].groupby("organism")["mbases"].sum() / 1000,
        "max_gbases_single_biosample": wgs.groupby(["organism", "biosample"])["mbases"].sum()
                                          .groupby(level=0).max() / 1000,
        "platforms": g["platform"].agg(lambda s: ";".join(sorted(set(s)))),
        "bioprojects": g["bioproject"].agg(lambda s: ";".join(sorted(set(s.dropna())))),
        "centers": g["center_name"].agg(lambda s: ";".join(sorted(set(s.dropna()))[:5])),
        "first_release": g["releasedate"].min(),
    }).fillna({"short_read_gbases": 0, "long_read_gbases": 0})
    summ["est_cov_best_individual_x"] = (summ["max_gbases_single_biosample"] / GENOME_GB).round(1)
    summ["snapshot"] = keys[0].split("/")[-1][:8]
    summ["retrieved"] = date.today().isoformat()
    summ.reset_index().rename(columns={"organism": "sra_organism"}) \
        .to_csv(OUT / "sra_parrot_species.tsv", sep="\t", index=False, float_format="%.2f")
    print(f"species with WGS: {len(summ)}; with >=20x in one individual: "
          f"{(summ['est_cov_best_individual_x'] >= 20).sum()}")


if __name__ == "__main__":
    main()
