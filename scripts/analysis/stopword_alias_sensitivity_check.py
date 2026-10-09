"""
stopword_alias_sensitivity_check.py -- isolates how much of the "vocabulary QC" effect
(cleaned 93.2%/75.8% vs true-raw 44.2%/27.4%, BM25 pool, all 19 cancers) is
attributable ONLY to removing the 13 single-term aliases that gene_quality.py's
own audit log flags as "Short or common-word alias filter" removals -- i.e.
raw gene-symbol aliases that are also ordinary English function words:
  of, in, for, was, all, not, as, an, cml, phl
mapped to: BRIP1, CD44, SPI1, WWOX, WAS, BCR(x3), NR4A2, UBE3A, CCDC102B,
PAX6, DIAPH3.

This is a direct test of the external review's hole #2: "If five or six
terrible aliases account for most of the interaction, your finding is mainly
about catastrophic entity-recognition failures, not necessarily a broad
property of corpus retrieval and vocabulary cleaning." Diagnostic only, not
part of the frozen 2x2 -- logged here, to be folded into the frozen spec's
amendment log once reviewed.

Three vocabulary conditions, same real BM25 pool + real formula throughout:
  A. True raw (gene_quality.py raw, zero filtering) -- known: 44.2%/27.4%
  B. Raw MINUS ONLY the 13 stopword-alias term entries (minimal fix)
  C. Full cleaned (this paper's production vocabulary) -- known: 93.2%/75.8%

Run: python -X utf8 scripts/analysis/stopword_alias_sensitivity_check.py
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
GENE_REMAP = {"P53": "TP53", "HER2": "ERBB2"}
OUT_DIR = ROOT / "data" / "derived" / "intogen_rescoring"

STOPWORD_ALIASES = {"of", "in", "for", "was", "all", "not", "as", "an", "cml", "phl"}

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


def strip_stopword_aliases(lexicon: dict) -> dict:
    single = {k: v for k, v in lexicon["single"].items() if k not in STOPWORD_ALIASES}
    return {"single": single, "phrases": lexicon["phrases"]}


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
    raw_gq, vetted_gq, audit = gq.gene_lookups(DEFAULT_ONCODIGGER)
    minimal_fix = strip_stopword_aliases(raw_gq)
    removed = len(raw_gq["single"]) - len(minimal_fix["single"])
    print(f"Minimal fix removed {removed} single-token alias entries from {len(raw_gq['single'])} raw entries "
          f"(expected 10-13 incl. duplicate-keyed stopwords)")

    genes_any, exact_pairs = load_intogen()

    rows = []
    minfix_any = minfix_exact = 0
    for corpus_key, cancer in CANCER_MAP:
        db = DEFAULT_ONCODIGGER / "data" / "pubmed_processed" / f"{corpus_key}_all" / "lexical_index.db"
        con = sqlite3.connect(str(db))
        top_minfix = rank_genes(con, minimal_fix)
        m_any = sum(1 for g in top_minfix if g in genes_any)
        m_exact = sum(1 for g in top_minfix if (g, cancer) in exact_pairs)
        minfix_any += m_any; minfix_exact += m_exact
        print(f"{cancer:20} minimal_fix={m_any}/10 any, {m_exact}/10 exact  top10={top_minfix}")
        rows.append({"cancer": cancer, "minimal_fix_top10": "|".join(top_minfix),
                     "minfix_any": m_any, "minfix_exact": m_exact})

    print()
    print(f"POOLED: minimal_fix (raw minus 13 stopword aliases only) "
          f"{minfix_any}/190 = {100*minfix_any/190:.1f}% any-cancer, "
          f"{minfix_exact}/190 = {100*minfix_exact/190:.1f}% exact")
    print("Reference: true-raw = 84/190 = 44.2% any / 52/190 = 27.4% exact")
    print("Reference: full-cleaned = 177/190 = 93.2% any / 144/190 = 75.8% exact")
    print(f"Minimal-fix recovers {100*(minfix_any-84)/(177-84):.1f}% of the any-cancer gap, "
          f"{100*(minfix_exact-52)/(144-52):.1f}% of the exact-concordance gap, "
          f"using only 13 of the full cleaning policy's removals.")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUT_DIR / "stopword_alias_sensitivity_check.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {OUT_DIR / 'stopword_alias_sensitivity_check.csv'}")


if __name__ == "__main__":
    main()
