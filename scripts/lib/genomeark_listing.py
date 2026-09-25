"""List Psittaciformes assemblies on GenomeArk (VGP) and fetch their gfastats summaries.

Usage: python3 scripts/lib/genomeark_listing.py <out_dir>

Parrot genera are taken from the Burgio et al. 2019 trait database, so the filter
does not depend on a hand-written list.
"""
import csv
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

BUCKET = "https://genomeark.s3.amazonaws.com"
TRAITS = Path("data/external/parrottraits/data-raw/parrot-traits-database.csv")


def s3_list(prefix, delimiter="/"):
    """Return (common_prefixes, keys) for an S3 prefix, following pagination."""
    prefixes, keys, token = [], [], None
    while True:
        q = {"list-type": "2", "prefix": prefix, "max-keys": "1000"}
        if delimiter:
            q["delimiter"] = delimiter
        if token:
            q["continuation-token"] = token
        with urllib.request.urlopen(f"{BUCKET}/?{urllib.parse.urlencode(q)}", timeout=60) as r:
            xml = r.read().decode()
        prefixes += re.findall(r"<Prefix>([^<]+)</Prefix>", xml)[1:] if delimiter else []
        keys += re.findall(r"<Key>([^<]+)</Key>", xml)
        m = re.search(r"<NextContinuationToken>([^<]+)</NextContinuationToken>", xml)
        if not m:
            return [p for p in prefixes if p != prefix], keys
        token = m.group(1)


def parrot_genera():
    with open(TRAITS, encoding="utf-8-sig") as f:
        rows = list(csv.reader(f))
    header = rows[1]
    gi = header.index("genus")
    return {r[gi].strip() for r in rows[2:] if len(r) > gi and r[gi].strip()}


def last_block(text):
    """gfastats files may hold several reports; keep the last one."""
    blocks, cur = [], {}
    for line in text.splitlines():
        if "\t" not in line:
            continue
        k, v = line.split("\t", 1)
        if k == "# scaffolds" and cur:
            blocks.append(cur)
            cur = {}
        cur[k] = v.replace(",", "")
    if cur:
        blocks.append(cur)
    return blocks[-1] if blocks else {}, len(blocks)


def main(out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    genera = parrot_genera()
    species, _ = s3_list("species/")
    parrots = sorted(p.split("/")[1] for p in species if p.split("/")[1].split("_")[0] in genera)
    rows = []
    for sp in parrots:
        tolids, _ = s3_list(f"species/{sp}/")
        for tp in tolids:
            tolid = tp.rstrip("/").split("/")[-1]
            folders, _ = s3_list(tp)
            _, keys = s3_list(f"{tp}assembly_curated/", delimiter=None)
            fastas = sorted((k for k in keys if re.search(r"\.(fa|fasta)\.gz$", k) and "/intermediates/" not in k
                            and not re.search(r"\.(alt|dup|hap2|pat|h2)\.", k) and ".MT." not in k),
                            key=lambda k: (".cur." in k, k))  # prefer curated over decontaminated
            stats = sorted(k for k in keys if "gfastats" in k and not re.search(r"\.(alt|hap2)\.", k))
            row = {
                "species": sp.replace("_", " "), "tolid": tolid,
                "assembly_folders": ";".join(f.rstrip("/").split("/")[-1] for f in folders),
                "latest_curated_fasta": fastas[-1].split("/")[-1] if fastas else "",
                "gfastats_file": stats[-1].split("/")[-1] if stats else "",
            }
            if stats:
                raw = urllib.request.urlopen(f"{BUCKET}/{urllib.parse.quote(stats[-1])}", timeout=120).read()
                if stats[-1].endswith(".gz"):
                    import gzip
                    raw = gzip.decompress(raw)
                (out / row["gfastats_file"].removesuffix(".gz")).write_bytes(raw)
                s, nblocks = last_block(raw.decode())
                row.update({
                    "gfastats_blocks": nblocks,
                    "total_scaffold_length_bp": s.get("Total scaffold length", ""),
                    "n_scaffolds": s.get("# scaffolds", ""),
                    "scaffold_n50_bp": s.get("Scaffold N50", ""),
                    "n_contigs": s.get("# contigs", ""),
                    "contig_n50_bp": s.get("Contig N50", ""),
                    "gc_percent": s.get("GC content %", ""),
                })
            rows.append(row)
            print(f"genomeark {sp} {tolid} {row['latest_curated_fasta'] or '-'}")
    cols = ["species", "tolid", "assembly_folders", "latest_curated_fasta", "gfastats_file", "gfastats_blocks",
            "total_scaffold_length_bp", "n_scaffolds", "scaffold_n50_bp", "n_contigs", "contig_n50_bp", "gc_percent"]
    with open(out / "genomeark_parrots.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main(sys.argv[1])
