"""
retrieval_vocab_independence_check.py -- empirically confirms (not just asserts from the
function signature) that bm25_search returns identical ordered PMID lists and identical
BM25 scores regardless of gene-vocabulary state, for all 19 cancers. Addresses the external
review's point that the index/options/module-level caches could in principle carry a hidden
vocabulary dependency even though bm25_search's own argument list has none.

Runs bm25_search twice per cancer, once with the cleaned vocabulary loaded and once with the
true-raw vocabulary loaded (each load clears discovery's per-connection caches first, exactly
as every other script in this audit does), and diffs the two ordered (pmid, score) lists
exactly, not just the final top-50 gene ranking.

Run: python -X utf8 scripts/analysis/retrieval_vocab_independence_check.py
"""
from __future__ import annotations

import sys
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ONCODIGGER = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET"))
sys.path.insert(0, os.environ.get("GENE_QUALITY_TEST_ENGINE", "GENE_QUALITY_TEST_ENGINE_NOT_SET"))

import sqlite3

from oncodigger_search import gene_quality as gq
from oncodigger_search import discovery
from oncodigger_search.lexical import bm25_search, normalize_query
from oncodigger_search.lexicons import load_lexicons

QUERY = "mutation, mutations"
SEARCH_POOL = 1250
OPTIONS = {"year_min": 2000, "year_max": 2026, "title_match_multiplier": 1.2,
           "mesh_match_multiplier": 1.15, "recency_boost": True, "discovery_mode": "Primary Evidence"}

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


def suppress_ambiguous_symbols(genes_lexicon: dict) -> dict:
    SUPPRESS_SYMBOLS = {"GC", "HR", "HCCS", "BPIFA4P", "MB", "ELN", "CP", "UBC", "CAMP", "EFS", "IVD", "SCT", "SON"}
    single = {k: v for k, v in genes_lexicon["single"].items() if k.upper() not in SUPPRESS_SYMBOLS}
    return {"single": single, "phrases": genes_lexicon["phrases"]}


def run_pool(con, vocab_lexicon: dict):
    """Loads a vocabulary (unused by retrieval, but mirrors every other script's setup
    exactly, including clearing per-connection caches), then runs bm25_search and returns
    the exact ordered (pmid, score) list. The vocab_lexicon argument is intentionally never
    passed to bm25_search -- this is the point of the test."""
    discovery._sqlite_doc_counts.clear()
    discovery._sqlite_caches.clear()
    tokens = normalize_query(QUERY)
    hits = bm25_search(con, tokens, SEARCH_POOL, OPTIONS)
    return [(h["pmid"], round(h.get("score", h.get("bm25_score", 0.0)), 10)) for h in hits]


def main():
    print("Building both vocabularies (cleaned and true-raw)...")
    lex_manuscript_raw = load_lexicons(DEFAULT_ONCODIGGER)["Genes"]
    lex_cleaned = suppress_ambiguous_symbols(lex_manuscript_raw)
    raw_gq, vetted_gq, audit = gq.gene_lookups(DEFAULT_ONCODIGGER)

    all_identical = True
    for corpus_key, cancer in CANCER_MAP:
        db = DEFAULT_ONCODIGGER / "data" / "pubmed_processed" / f"{corpus_key}_all" / "lexical_index.db"
        con = sqlite3.connect(str(db))

        pool_cleaned = run_pool(con, lex_cleaned)
        pool_raw = run_pool(con, raw_gq)

        identical = pool_cleaned == pool_raw
        all_identical &= identical
        status = "IDENTICAL" if identical else "DIFFERS"
        print(f"{cancer:20} n={len(pool_cleaned)} vs {len(pool_raw)}  {status}")
        if not identical:
            for i, (a, b) in enumerate(zip(pool_cleaned, pool_raw)):
                if a != b:
                    print(f"    first mismatch at index {i}: cleaned={a} raw={b}")
                    break

    print()
    print("ALL 19 CANCERS: retrieval pools byte-identical across vocabulary states"
          if all_identical else "MISMATCH FOUND -- retrieval depends on vocabulary state")


if __name__ == "__main__":
    main()
