"""intogen_literature_weight_audit.py -- manuscript1.1 (IntOGen-primary draft) ONLY.
Literature-weight audit of unrecovered IntOGen driver genes, mirroring
cgc_literature_weight_audit.py exactly except the final gene-cancer pair set is
IntOGen's 1,740-pair conservative mapping instead of CGC's 781 Tier 1 annotations.
The engine-query logic (BM25 search, extended entity-table extraction, supporting-
paper/corpus-paper/relevance extraction for every gene that surfaces at all) is
identical and untouched -- only which gene-cancer pairs get written to the output
differs, since the engine itself has no concept of CGC or IntOGen.

Requires the local OncoDigger corpora (not in this repository).

Run: python -X utf8 scripts/analysis/intogen_literature_weight_audit.py
"""
from __future__ import annotations

import argparse
import csv
import sys
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ONCODIGGER = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET"))
OUT_CSV = ROOT / "data" / "derived" / "intogen_rescoring" / "intogen_literature_weight_audit.csv"
INTOGEN_DRIVERS = ROOT / "data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_query_planes import (  # noqa: E402
    CANCER_MAP, GENE_EXCLUDE, GENE_REMAP, OPTIONS, longest_match_lookups, open_corpus,
    suppress_ambiguous_symbols,
)

ENTITY_LIMIT = 3000
SEARCH_POOL = 1250
ANALYSIS_PAPERS = 1000

# Identical conservative mapping used throughout the manuscript1.1 IntOGen conversion
# (intogen_full_rescoring.py, intogen_ablation_panel.py, intogen_discordant_cross_check.py).
INTOGEN_CANCER_MAP = {
    "Breast Cancer": {"BRCA"}, "Lung Cancer": {"LUAD", "LUSC", "NSCLC", "SCLC"},
    "Colorectal Cancer": {"COAD", "READ", "COADREAD"}, "Prostate Cancer": {"PRAD"},
    "Melanoma": {"MEL", "SKCM", "UM"}, "Bladder Cancer": {"BLCA", "UTUC"},
    "Kidney Cancer": {"CCRCC", "CHRCC", "PRCC", "RCC", "WT"},
    "Pancreatic Cancer": {"PAAD", "PANET"}, "Liver Cancer": {"HCC", "CHOL", "LIHB"},
    "Stomach Cancer": {"STAD"}, "Esophageal Cancer": {"ESCA", "ESCC"},
    "Ovarian Cancer": {"OVT"}, "Endometrial Cancer": {"UCEC", "UCS"},
    "Cervical Cancer": {"CESC", "CEAD"}, "Thyroid Cancer": {"WDTC"},
    "Brain Cancer": {"GB", "GBM", "HGGNOS", "LGGNOS", "PAST", "MBL", "EPM", "ATRT"},
    "Leukemia": {"ALL", "AML", "CLLSLL", "CML", "MDS"},
    "Lymphoma": {"BL", "DLBCLNOS", "NHL", "MLYM"}, "Myeloma": {"PCM"},
}


def build_intogen_pairs():
    """cancer -> set of genes, mirroring build_cgc_pairs()'s shape exactly."""
    pairs = {}
    with open(INTOGEN_DRIVERS, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            code = row["CANCER_TYPE"]
            for cancer, codes in INTOGEN_CANCER_MAP.items():
                if code in codes:
                    pairs.setdefault(cancer, set()).add(row["SYMBOL"])
    return pairs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--oncodigger-root", type=Path, default=DEFAULT_ONCODIGGER)
    ap.add_argument("--cancers", nargs="+")
    args = ap.parse_args()

    sys.path.insert(0, str(args.oncodigger_root / "src"))
    from oncodigger.search import discovery
    from oncodigger.search.discovery import associated_entity_tables
    from oncodigger.search.lexical import bm25_search, normalize_query
    from oncodigger.search.lexicons import load_lexicons

    lexicons = load_lexicons(args.oncodigger_root)
    genes = suppress_ambiguous_symbols(lexicons["Genes"])
    genes = longest_match_lookups(genes, discovery)
    lookups = {"Genes": genes}

    intogen_pairs = build_intogen_pairs()
    db_base = args.oncodigger_root / "data" / "pubmed_processed"

    rows = []
    for key, label in CANCER_MAP:
        if args.cancers and label not in args.cancers:
            continue
        db = next((p for p in (db_base / key / "lexical_index.db", db_base / f"{key}_all" / "lexical_index.db")
                   if p.exists()), None)
        if db is None:
            print(f"missing index for {label}, skipping", flush=True)
            continue
        con = open_corpus(db, discovery)
        tokens = normalize_query("mutation mutations")
        hits = bm25_search(con, tokens, SEARCH_POOL, OPTIONS)[:ANALYSIS_PAPERS]
        tables = associated_entity_tables(con, hits, "mutation mutations", (2000, 2026), tokens, lookups,
                                           limit=ENTITY_LIMIT, discovery_mode="Primary Evidence", corpus_query="")
        df = tables.get("Genes")
        engine_rank = {}
        support = {}
        corpus_n = {}
        relevance = {}
        if df is not None and not df.empty:
            seen = set()
            r = 0
            for _, row in df.iterrows():
                gene = GENE_REMAP.get(row["Gene"], row["Gene"])
                if gene in GENE_EXCLUDE or gene in seen:
                    continue
                try:
                    n = con.execute("SELECT COUNT(DISTINCT pmid) FROM postings WHERE token = ?",
                                     (gene.lower(),)).fetchone()[0]
                except Exception:
                    n = 0
                if n == 0:
                    continue
                seen.add(gene)
                r += 1
                engine_rank[gene] = r
                support[gene] = row.get("Supporting papers", "")
                corpus_n[gene] = row.get("Corpus papers", "")
                relevance[gene] = row.get("Relevance", "")

        intogen_genes = intogen_pairs.get(label, set())
        for gene in sorted(intogen_genes):
            rows.append({
                "cancer": label, "gene": gene,
                "engine_rank": engine_rank.get(gene, ""),
                "supporting_papers": support.get(gene, ""),
                "corpus_papers": corpus_n.get(gene, ""),
                "relevance_score": relevance.get(gene, ""),
                "surfaced_at_all": gene in engine_rank,
            })
        con.close()
        n_surfaced = sum(1 for g in intogen_genes if g in engine_rank)
        print(f"{label:20} IntOGen genes={len(intogen_genes):3} surfaced={n_surfaced:3} "
              f"max_rank_seen={len(engine_rank)}", flush=True)

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["cancer", "gene", "engine_rank", "supporting_papers",
                                            "corpus_papers", "relevance_score", "surfaced_at_all"])
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {OUT_CSV} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
