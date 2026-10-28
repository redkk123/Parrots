"""Canonical proteome of Melopsittacus undulatus (RefSeq GCF_012275295.1): one protein per gene, the genome-wide
universe for H6 and the H8 null (analysis plan, 2026-10-28).

Same isoform rule as the H8 panel (scripts/25_prepare_reference_geneset.py: longest valid CDS, then curated over
predicted, lower isoform index, accession). The 101 panel proteins are kept exactly as in data/pilot_genes (same
headers), and their RefSeq genes are not added twice. Genes without a valid CDS are skipped.
Second reference for gene-loss calls (plan §10, Nível 1, item 1: "contra Melopsittacus e uma espécie longeva com
montagem cromossômica"): the panel genes of the kakapo (Strigops habroptila, RefSeq GCF_004027225.2), same isoform
rule, appended as kakapo_<SYMBOL> (94 of 101 found; the rest have other names in that annotation).
Writes data/genome_genes/{canonical_proteins.faa, canonical_cds.fna, canonical_genes.tsv, kakapo_panel.tsv}.
"""
import importlib.util
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REF = ROOT / "data/external/ref_melopsittacus/pkg/ncbi_dataset/data/GCF_012275295.1"
KAKAPO = ROOT / "data/external/ref_strigops/pkg/ncbi_dataset/data/GCF_004027225.2"
PANEL = ROOT / "data/pilot_genes"
OUT = ROOT / "data/genome_genes"
spec = importlib.util.spec_from_file_location("refset", ROOT / "scripts/25_prepare_reference_geneset.py")
refset = importlib.util.module_from_spec(spec)
spec.loader.exec_module(refset)


def load(ref):
    cds_by_protein = {}
    for desc, seq in refset.parse_fasta(ref / "cds_from_genomic.fna").values():
        m = re.search(r"\[protein_id=([^\]]+)\]", desc)
        if m:
            cds_by_protein[m.group(1)] = (desc, seq)
    return (cds_by_protein, refset.parse_fasta(ref / "protein.faa")) + refset.parse_refseq_gff(ref / "genomic.gff")


def kakapo_panel():
    if not KAKAPO.exists():
        raise SystemExit(f"missing {KAKAPO}: datasets download genome accession GCF_004027225.2 --include protein,cds,gff3")
    cds_by_protein, proteins, mrna_info, gene_to_transcripts = load(KAKAPO)
    out = []
    for d in refset.GENE_DEFINITIONS:
        try:
            g = refset.select_canonical_isoforms([d], mrna_info, gene_to_transcripts, cds_by_protein, proteins)[0]
        except ValueError:
            continue
        g["prot"], g["cds"] = proteins[g["protein_accession"]][1], cds_by_protein[g["protein_accession"]][1]
        out.append(g)
    return out


def main():
    cds_by_protein, proteins, mrna_info, gene_to_transcripts = load(REF)
    panel = pd.read_csv(PANEL / "canonical_genes.tsv", sep="\t")
    panel_ids = set(panel.gene_id.astype(str))
    names = {i["gene_name"].upper(): i["gene_name"] for i in mrna_info.values() if i["gene_name"]}
    rows, skipped = [], 0
    for key in sorted(gene_to_transcripts):
        if mrna_info[gene_to_transcripts[key][0]]["gene_id"] in panel_ids:
            continue
        try:
            g = refset.select_canonical_isoforms([(names[key], "genome", names[key])], mrna_info,
                                                 gene_to_transcripts, cds_by_protein, proteins)[0]
        except ValueError:
            skipped += 1
            continue
        rows.append(g)
    OUT.mkdir(parents=True, exist_ok=True)
    prot = (PANEL / "canonical_proteins.faa").read_text()
    cds = (PANEL / "canonical_cds.fna").read_text()
    for g in rows:
        h = f">{g['gene_symbol']}|{g['gene_id']}|{g['transcript_accession']}|{g['protein_accession']}\n"
        prot += h + proteins[g["protein_accession"]][1] + "\n"
        cds += h + cds_by_protein[g["protein_accession"]][1] + "\n"
    kk = kakapo_panel()
    for g in kk:
        h = f">kakapo_{g['gene_symbol']}|{g['gene_id']}|{g['transcript_accession']}|{g['protein_accession']}\n"
        prot += h + g["prot"] + "\n"
        cds += h + g["cds"] + "\n"
    pd.DataFrame(kk)[["gene_symbol", "module", "gene_id", "transcript_accession", "protein_accession", "cds_length",
                      "protein_length"]].to_csv(OUT / "kakapo_panel.tsv", sep="\t", index=False)
    (OUT / "canonical_proteins.faa").write_text(prot)
    (OUT / "canonical_cds.fna").write_text(cds)
    cols = ["gene_symbol", "gene_id", "transcript_accession", "protein_accession", "cds_length", "protein_length",
            "exon_count"]
    t = pd.concat([panel.assign(in_panel=True), pd.DataFrame(rows)[cols].assign(module="", in_panel=False)])
    t.to_csv(OUT / "canonical_genes.tsv", sep="\t", index=False)
    print(f"{len(t)} genes ({len(panel)} panel + {len(rows)} others); {skipped} genes without a valid CDS; "
          f"{t.gene_symbol.str.startswith('LOC').sum()} LOC symbols; kakapo panel {len(kk)}")


if __name__ == "__main__":
    main()
