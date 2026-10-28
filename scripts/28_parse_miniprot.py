"""Per-gene calls from a miniprot GFF of the Melopsittacus canonical proteome (data/genome_genes, which contains the
101 H8 panel proteins of data/pilot_genes with the same headers) on one genome.

Usage: 28_parse_miniprot.py <accession> <species> <genome.fna> <miniprot.gff>
For each query protein, the best alignment (Rank=1, highest score) gives the CDS, spliced from its CDS features and
translated (table 1). Call, as fixed in the plan (2026-10-16: a gene enters a species only with a complete ORF):
  complete   >= 90 % of the query aligned, CDS length a multiple of 3, no frameshift, no internal stop codon
  disrupted  frameshift or internal stop (candidate pseudogene: not called a loss without read-level checks, §4)
  partial    < 90 % of the query aligned, no disruption (assembly fragment; missing data)
  paralog    best alignment has < 60 % identity, or its locus is the best locus of another query protein with
             higher identity (e.g. BCL2L1 landing on BCL2 when its own locus is missing): treated as missing
  absent     no alignment
Coordinates are reported on the original sequences (window offsets from seqkit sliding are removed).
Writes, for the panel, results/phase2_genes/hits/<acc>.tsv (one row per gene) and results/phase2_genes/cds/<acc>.fna
(complete only); for the whole proteome, results/phase3_genome/calls/<acc>.tsv.gz and data/genome_cds/<acc>.fna.gz.
"""
import gzip
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase2_genes"
GENOME_OUT = ROOT / "results/phase3_genome/calls"
GENOME_CDS = ROOT / "data/genome_cds"
QUERY = ROOT / "data/genome_genes/canonical_proteins.faa"
PANEL = ROOT / "data/pilot_genes/canonical_proteins.faa"
MIN_QCOV = 0.90
MIN_IDENTITY = 0.60  # orthologs across parrots: >= 0.75 in the Loriculus test; paralog hits 0.36-0.52
CODE = {a + b + c: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"[16 * i + 4 * j + k]
        for i, a in enumerate("TCAG") for j, b in enumerate("TCAG") for k, c in enumerate("TCAG")}
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def fasta(path, keep=None):
    name, buf = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None and (keep is None or name in keep):
                    yield name, "".join(buf)
                name, buf = line[1:].split()[0], []
            elif keep is None or name in keep:
                buf.append(line.strip())
    if name is not None and (keep is None or name in keep):
        yield name, "".join(buf)


def translate(s):
    return "".join(CODE.get(s[i:i + 3].upper(), "X") for i in range(0, len(s) - 2, 3))


def genomic(ctg, start, end):
    """Window name from seqkit sliding ('chr_sliding:4500001-9500000') -> original sequence coordinates."""
    if "_sliding:" in ctg:
        base, win = ctg.rsplit("_sliding:", 1)
        off = int(win.split("-")[0]) - 1
        return base, start + off, end + off
    return ctg, start, end


def parse_gff(path):
    aln = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] == "mRNA":
            q, ts, te = attr["Target"].split()
            aln[attr["ID"]] = {"query": q, "ctg": f[0], "strand": f[6], "start": int(f[3]), "end": int(f[4]),
                               "score": float(f[5]), "rank": int(attr.get("Rank", 1)), "tstart": int(ts),
                               "tend": int(te), "identity": float(attr.get("Identity", 0)),
                               "frameshift": int(attr.get("Frameshift", 0)), "stop": int(attr.get("StopCodon", 0)),
                               "cds": []}
        elif f[2] == "CDS":
            aln[attr["Parent"]]["cds"].append((int(f[3]), int(f[4])))
    return aln


def main():
    acc, species, genome, gff = sys.argv[1:5]
    qlen = {n.split("|")[0]: len(s) for n, s in fasta(QUERY)}
    by_q = {}
    for a in parse_gff(gff).values():
        by_q.setdefault(a["query"].split("|")[0], []).append(a)
    best = {g: max(v, key=lambda a: (-a["rank"], a["score"])) for g, v in by_q.items()}
    seqs = dict(fasta(genome, keep={a["ctg"] for a in best.values()}))
    loci = {g: genomic(a["ctg"], a["start"], a["end"]) for g, a in best.items()}
    taken = set()  # queries whose best locus overlaps the best locus of another query with higher identity
    active, last_ctg = [], None
    for g, (c, x, y) in sorted(loci.items(), key=lambda kv: kv[1]):
        if c != last_ctg:
            active, last_ctg = [], c
        active = [h for h in active if loci[h][2] >= x]
        for h in active:
            if best[h]["identity"] > best[g]["identity"]:
                taken.add(g)
            elif best[g]["identity"] > best[h]["identity"]:
                taken.add(h)
        active.append(g)
    rows, out = [], []
    for g in qlen:
        a = best.get(g)
        r = {"gene": g, "n_alignments": len(by_q.get(g, []))}
        if a is None:
            rows.append({**r, "call": "absent"})
            continue
        s = seqs[a["ctg"]]
        cds = "".join(s[x - 1:y] for x, y in sorted(a["cds"]))
        if a["strand"] == "-":
            cds = cds.translate(COMP)[::-1]
        prot = translate(cds)
        if prot.endswith("*"):
            prot, cds = prot[:-1], cds[:-3]
        internal = prot.count("*")
        qcov = (a["tend"] - a["tstart"] + 1) / qlen[g]
        disrupted = a["frameshift"] > 0 or internal > 0 or len(cds) % 3 != 0
        call = "disrupted" if disrupted else ("complete" if qcov >= MIN_QCOV else "partial")
        if a["identity"] < MIN_IDENTITY or g in taken:
            call = "paralog"
        c, x, y = loci[g]
        rows.append({**r, "call": call, "query_cov": round(qcov, 3), "identity": a["identity"],
                     "frameshifts": a["frameshift"], "internal_stops": internal, "cds_len": len(cds),
                     "starts_met": prot[:1] == "M", "location": f"{c}:{x}-{y}{a['strand']}",
                     "n_exons": len(a["cds"])})
        if call == "complete":
            out.append((g, f">{species.replace(' ', '_')}|{g}|{acc}\n{cds}\n"))
    panel = {n.split("|")[0] for n, _ in fasta(PANEL)}
    cols = ["gene", "call", "query_cov", "identity", "frameshifts", "internal_stops", "cds_len", "starts_met",
            "n_exons", "n_alignments", "location"]

    def table(rs):
        return "\t".join(cols) + "\n" + "".join("\t".join(str(r.get(c, "")) for c in cols) + "\n" for r in rs)

    for d in (OUT / "hits", OUT / "cds", GENOME_OUT, GENOME_CDS):
        d.mkdir(parents=True, exist_ok=True)
    (OUT / "cds" / f"{acc}.fna").write_text("".join(x for g, x in out if g in panel))
    (OUT / "hits" / f"{acc}.tsv").write_text(table(r for r in rows if r["gene"] in panel))
    with gzip.open(GENOME_CDS / f"{acc}.fna.gz", "wt") as f:
        f.write("".join(x for _, x in out))
    with gzip.open(GENOME_OUT / f"{acc}.tsv.gz", "wt") as f:
        f.write(table(rows))
    for name, rs in (("panel", [r for r in rows if r["gene"] in panel]), ("genome", rows)):
        calls = [r["call"] for r in rs]
        print(acc, species, name, {c: calls.count(c) for c in ("complete", "partial", "disrupted", "paralog", "absent")})


if __name__ == "__main__":
    main()
