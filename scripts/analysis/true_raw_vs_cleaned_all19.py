"""
true_raw_vs_cleaned_all19.py -- the correct vocabulary-QC counterfactual: Cleaned (this
paper's actual pipeline, confirmed identical to the live site's "Vetted" for Leukemia and
Melanoma) vs. TRUE Raw (gene_quality.py's genuinely unfiltered term list -- no
useful_gene_term() filter, no exclusions of any kind), not the partially-filtered
"less-cleaned" condition used in the earlier 2x2 (which had useful_gene_term() baked into
load_lexicons() regardless of the suppress_ambiguous_symbols step). "No QC" has no other
real alternative besides true Raw, so this is the comparison that actually answers "what
does vocabulary QC contribute."

Same real BM25 retrieval pool (confirmed vocabulary-independent, Phase 0) and real
rank-weighted+enrichment formula throughout -- only the vocabulary varies.

Run: python -X utf8 scripts/analysis/true_raw_vs_cleaned_all19.py
"""
from __future__ import annotations

import csv
import sys
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ONCODIGGER = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET"))
sys.path.insert(0, os.environ.get("GENE_QUALITY_TEST_ENGINE", "GENE_QUALITY_TEST_ENGINE_NOT_SET"))

import sqlite3

from oncodigger_search import gene_quality as gq
from oncodigger_search import discovery
from oncodigger_search.discovery import associated_entity_tables
from oncodigger_search.lexical import bm25_search, normalize_query
from oncodigger_search.lexicons import load_lexicons

QUERY = "mutation, mutations"
SEARCH_POOL, ANALYSIS_PAPERS = 1250, 1000
OPTIONS = {"year_min": 2000, "year_max": 2026, "title_match_multiplier": 1.2,
           "mesh_match_multiplier": 1.15, "recency_boost": True, "discovery_mode": "Primary Evidence"}
GENE_EXCLUDE = {"MICE", "DCR", "PC", "C2", "FBN1", "SGCG", "WDHD1"}
SUPPRESS_SYMBOLS = {"GC", "HR", "HCCS", "BPIFA4P", "MB", "ELN", "CP", "UBC", "CAMP", "EFS", "IVD", "SCT", "SON"}
GENE_REMAP = {"P53": "TP53", "HER2": "ERBB2"}
OUT_DIR = ROOT / "data" / "derived" / "intogen_rescoring"

CANCER_MAP = [
    ("breast_cancer_inclusive", "Breast Cancer"), ("lung_cancer_inclusive", "Lung Cancer"),
    ("colorectal_cancer_inclusive", "Colorectal Cancer"), ("prostate_cancer_inclusive", "Prostate Cancer"),
    ("melanoma_cancer_inclusive", "Melanoma"), ("bladder_cancer_inclusive", "Bladder Cancer"),
    ("kidney_cancer_inclusive", "Kidney Cancer"), ("pancreas_cancer_inclusive", "Pancreatic Cancer"),
    ("liver_cancer_inclusive", "Liver Cancer"), ("stomach_cancer_inclusive", "Stomach Cancer"),
    ("esophageal_cancer_inclusive", "Esophageal Cancer"), ("ovarian_cancer_inclusive", "Ovarian Cancer"),
    ("endometrial_cancer_inclusive", "Endometrial Cancer"), ("cervical_cancer_inclusive", "Cervical Cancer"),
    ("thyroid_cancer_inclusive", "Thyroid Cancer"), ("brain_cancer_inclusive", "Brain Cancer"),
    ("leukemia_cancer_inclusive", "Leukemia"), ("lymphoma_cancer_inclusive", "Lymphoma"),
    ("myeloma_cancer_inclusive", "Myeloma"),
]

INTOGEN_CANCER_MAP = {
    "Breast Cancer": {"BRCA"}, "Lung Cancer": {"LUAD", "LUSC", "NSCLC", "SCLC"},
    "Colorectal Cancer": {"COAD", "READ", "COADREAD"}, "Prostate Cancer": {"PRAD"},
    "Melanoma": {"MEL", "SKCM", "UM"}, "Bladder Cancer": {"BLCA", "UTUC"},
    "Kidney Cancer": {"CCRCC", "CHRCC", "PRCC", "RCC", "WT"}, "Pancreatic Cancer": {"PAAD", "PANET"},
    "Liver Cancer": {"HCC", "CHOL", "LIHB"}, "Stomach Cancer": {"STAD"},
    "Esophageal Cancer": {"ESCA", "ESCC"}, "Ovarian Cancer": {"OVT"},
    "Endometrial Cancer": {"UCEC", "UCS"}, "Cervical Cancer": {"CESC", "CEAD"},
    "Thyroid Cancer": {"WDTC"}, "Brain Cancer": {"GB", "GBM", "HGGNOS", "LGGNOS", "PAST", "MBL", "EPM", "ATRT"},
    "Leukemia": {"ALL", "AML", "CLLSLL", "MDS"}, "Lymphoma": {"BL", "DLBCLNOS", "NHL", "MLYM"},
    "Myeloma": {"PCM"},
}


def suppress_ambiguous_symbols(genes_lexicon: dict) -> dict:
    single = {k: v for k, v in genes_lexicon["single"].items() if k.upper() not in SUPPRESS_SYMBOLS}
    return {"single": single, "phrases": genes_lexicon["phrases"]}


def rank_genes(con, vocab_lexicon: dict) -> list[str]:
    discovery._sqlite_doc_counts.clear()
    discovery._sqlite_caches.clear()
    tokens = normalize_query(QUERY)
    hits = bm25_search(con, tokens, SEARCH_POOL, OPTIONS)[:ANALYSIS_PAPERS]
    tables = associated_entity_tables(con, hits, QUERY, (2000, 2026), tokens, {"Genes": vocab_lexicon},
                                      limit=100, discovery_mode="Primary Evidence", corpus_query="")
    genes_df = tables.get("Genes")
    ranked, seen = [], set()
    if genes_df is None or genes_df.empty:
        return ranked
    for _, row in genes_df.iterrows():
        gene = GENE_REMAP.get(row["Gene"], row["Gene"])
        if gene in GENE_EXCLUDE or gene in seen:
            continue
        seen.add(gene)
        ranked.append(gene)
        if len(ranked) >= 10:
            break
    return ranked


def load_intogen():
    genes_any, exact_pairs = set(), set()
    path = ROOT / "data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            sym = row["SYMBOL"].strip().upper()
            genes_any.add(sym)
            code = row["CANCER_TYPE"].strip()
            for cancer, codes in INTOGEN_CANCER_MAP.items():
                if code in codes:
                    exact_pairs.add((sym, cancer))
    return genes_any, exact_pairs


def main() -> None:
    print("Building vocabularies...")
    lex_manuscript_raw = load_lexicons(DEFAULT_ONCODIGGER)["Genes"]
    lex_cleaned = suppress_ambiguous_symbols(lex_manuscript_raw)
    raw_gq, vetted_gq, audit = gq.gene_lookups(DEFAULT_ONCODIGGER)
    print(f"gene_quality.py vetted audit: {len(audit)} actions, "
          f"{sum(1 for a in audit if a['Action']=='Removed alias')} alias collisions, "
          f"{sum(1 for a in audit if a['Action']=='Removed term')} term removals")

    genes_any, exact_pairs = load_intogen()

    rows = []
    cleaned_any = cleaned_exact = raw_any = raw_exact = 0
    for corpus_key, cancer in CANCER_MAP:
        db = DEFAULT_ONCODIGGER / "data" / "pubmed_processed" / f"{corpus_key}_all" / "lexical_index.db"
        con = sqlite3.connect(str(db))
        top_cleaned = rank_genes(con, lex_cleaned)
        top_raw = rank_genes(con, raw_gq)
        c_any = sum(1 for g in top_cleaned if g in genes_any)
        c_exact = sum(1 for g in top_cleaned if (g, cancer) in exact_pairs)
        r_any = sum(1 for g in top_raw if g in genes_any)
        r_exact = sum(1 for g in top_raw if (g, cancer) in exact_pairs)
        cleaned_any += c_any; cleaned_exact += c_exact
        raw_any += r_any; raw_exact += r_exact
        print(f"{cancer:20} cleaned={c_any}/10 any, {c_exact}/10 exact | "
              f"true_raw={r_any}/10 any, {r_exact}/10 exact")
        rows.append({"cancer": cancer, "cleaned_top10": "|".join(top_cleaned), "cleaned_any": c_any,
                     "cleaned_exact": c_exact, "true_raw_top10": "|".join(top_raw), "raw_any": r_any,
                     "raw_exact": r_exact})

    print()
    print(f"POOLED: cleaned {cleaned_any}/190 = {100*cleaned_any/190:.1f}% any-cancer, "
          f"{cleaned_exact}/190 = {100*cleaned_exact/190:.1f}% exact")
    print(f"POOLED: true_raw {raw_any}/190 = {100*raw_any/190:.1f}% any-cancer, "
          f"{raw_exact}/190 = {100*raw_exact/190:.1f}% exact")
    print(f"TRUE vocabulary-QC effect: {100*(cleaned_any-raw_any)/190:+.1f}pp any-cancer, "
          f"{100*(cleaned_exact-raw_exact)/190:+.1f}pp exact concordance")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUT_DIR / "true_raw_vs_cleaned_all19.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {OUT_DIR / 'true_raw_vs_cleaned_all19.csv'}")


if __name__ == "__main__":
    main()
