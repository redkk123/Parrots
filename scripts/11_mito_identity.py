"""Specimen identity and foreign-mitochondrion check for one assembly.

Mitochondrial-like contigs are found with COX1/CYTB/ND2 queries from five distant parrots, then each candidate
contig is compared to independent GenBank parrot sequences (records by Hains/Pirro/Iridian excluded).

Usage: python scripts/11_mito_identity.py <accession> <genome.fna> "<claimed species>"
Writes results/phase1_audit/mito/<accession>.json

Per (contig, gene), alignments >= 400 bp:
  own      the claimed species is the top hit (or within 1 pp of it) at >= 98 % identity
  foreign  another species is the top hit at >= 98 % and the claimed species, which has a reference for
           that gene, is > 2 pp lower  -> a second mitogenome in the assembly (cross-sample contamination
           or wrong specimen)
  numt     top hit < 95 % (diverged nuclear copy)
Overall:
  confirmed                  own signal, no foreign signal
  confirmed_with_foreign_mito own signal plus foreign signal (identity OK; contamination flag)
  conflict                   foreign signal without own signal (possible wrong specimen)
  no_own_reference           close (>= 98 %) hits only to other species, claimed species lacks references
  no_reliable_mito           no candidate >= 95 % to any parrot (NUMTs only) -> identity must come from reads
  no_mito_found              no hits at all
"""
import json
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / ".bioenv/bin"
REF = ROOT / "data/external/mito_ref"
QUERY_SPECIES = ("Melopsittacus_undulatus", "Ara_macao", "Cacatua_sulphurea", "Nestor_notabilis", "Psittacus_erithacus")
MIN_ALN, TOL, N_CAND = 400, 1.0, 40
CLOSE, FOREIGN_GAP, NUMT_MAX = 98.0, 2.0, 95.0


def run(*a):
    return subprocess.run([str(x) for x in a], check=True, capture_output=True, text=True).stdout


def read_fasta(p):
    seqs, name = {}, None
    for line in open(p):
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def acc_id(s):
    return s.strip("|").split("|")[-1]  # 'gb|ACC|' (makeblastdb -parse_seqids) -> 'ACC'


def classify(d, claimed_u, has_ref):
    top_sp = max(d, key=d.get)
    top = d[top_sp]
    c = d.get(claimed_u)
    if top < NUMT_MAX:
        cls = "numt"
    elif c is not None and c >= top - TOL and c >= CLOSE:
        cls = "own"
    elif top >= CLOSE and has_ref and (c is None or c < top - FOREIGN_GAP):
        cls = "foreign"
    elif top >= CLOSE and not has_ref:
        cls = "close_other_no_own_ref"
    else:
        cls = "ambiguous"
    return {"class": cls, "top_species": top_sp, "top_pid": top, "claimed_pid": c}


def main(acc, genome, claimed):
    claimed_u = claimed.replace(" ", "_")
    out = ROOT / f"results/phase1_audit/mito/{acc}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(dir=ROOT / "data/external"))
    genes = read_fasta(REF / "indep_genes.fasta")
    own_ref = {k.split("|")[1] for k in genes if k.split("|")[2] == claimed_u}
    picked = {}
    for k in genes:
        _, g, sp, _ = k.split("|")
        if sp in QUERY_SPECIES and (sp, g) not in picked and len(genes[k]) >= 900:
            picked[(sp, g)] = k
    with open(tmp / "q.fa", "w") as f:
        for k in picked.values():
            f.write(f">{k}\n{genes[k]}\n")
    run(BIN / "makeblastdb", "-in", genome, "-dbtype", "nucl", "-parse_seqids", "-out", tmp / "g")
    hits = run(BIN / "blastn", "-task", "dc-megablast", "-query", tmp / "q.fa", "-db", tmp / "g", "-evalue", "1e-30",
               "-outfmt", "6 qseqid sseqid pident length slen bitscore", "-max_target_seqs", "100", "-num_threads", "8")
    score, slen = defaultdict(float), {}
    for line in hits.splitlines():
        q, s, pid, ln, sl, bs = line.split("\t")
        s = acc_id(s)
        score[s] += float(bs)
        slen[s] = int(sl)
    res = {"accession": acc, "claimed": claimed, "own_reference_genes": sorted(own_ref), "contigs": {}}
    if score:
        cand = sorted(score, key=score.get, reverse=True)[:N_CAND]
        (tmp / "ids").write_text("\n".join(cand) + "\n")
        (tmp / "m.fa").write_text(run(BIN / "blastdbcmd", "-db", tmp / "g", "-entry_batch", tmp / "ids"))
        gc = {k: round(100 * sum(v.upper().count(x) for x in "GC") / max(1, len(v)), 1)
              for k, v in read_fasta(tmp / "m.fa").items()}
        hits = run(BIN / "blastn", "-task", "dc-megablast", "-query", tmp / "m.fa", "-db", REF / "db/indep",
                   "-evalue", "1e-50", "-outfmt", "6 qseqid sseqid pident length", "-max_target_seqs", "10000",
                   "-num_threads", "8")
        by = defaultdict(lambda: defaultdict(dict))  # contig -> gene -> species -> best pident
        for line in hits.splitlines():
            q, sid, pid, ln = line.split("\t")
            if int(ln) < MIN_ALN:
                continue
            _, g, sp, _ = sid.split("|")
            d = by[acc_id(q)][g]
            d[sp] = max(d.get(sp, 0.0), float(pid))
        for c, gs in by.items():
            res["contigs"][c] = {"len": slen.get(c), "gc": gc.get(c),
                                 "genes": {g: classify(d, claimed_u, g in own_ref) for g, d in gs.items()}}
    calls = [(c, g, v) for c, r in res["contigs"].items() for g, v in r["genes"].items()]
    own = [x for x in calls if x[2]["class"] == "own"]
    foreign = [x for x in calls if x[2]["class"] == "foreign"]
    res["own_genes"] = sorted({g for _, g, _ in own})
    res["foreign_species"] = sorted({v["top_species"] for _, _, v in foreign})
    res["foreign_contigs"] = sorted({c for c, _, _ in foreign})
    res["best_own_pid"] = max((v["claimed_pid"] for _, _, v in own), default=None)
    if not calls:
        res["overall"] = "no_mito_found"
    elif own:
        res["overall"] = "confirmed_with_foreign_mito" if foreign else "confirmed"
    elif foreign:
        res["overall"] = "conflict"
    elif any(v["class"] == "close_other_no_own_ref" for _, _, v in calls):
        res["overall"] = "no_own_reference"
    elif all(v["class"] == "numt" for _, _, v in calls):
        res["overall"] = "no_reliable_mito"
    else:
        res["overall"] = "ambiguous"
    json.dump(res, open(out, "w"), indent=1)
    run("rm", "-rf", tmp)
    return res


if __name__ == "__main__":
    r = main(sys.argv[1], sys.argv[2], sys.argv[3])
    print(r["accession"], r["claimed"], r["overall"], "own:", r["own_genes"], r["best_own_pid"],
          "foreign:", r["foreign_species"], r["foreign_contigs"])
