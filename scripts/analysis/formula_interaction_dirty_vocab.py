"""
Recovered 2026-10-08 via a full repo-wide CSV inventory: this script's output,
formula_interaction_dirty_vocab.csv, backs Results 5.8's "the combined formula (44.2%)"
figure (confirmed by independently recomputing real_any's pooled sum from the committed CSV:
84/190 = 44.2%, exact match) but had not been committed to this repository until now.

Known limitation: this script imports `oncodigger_search.gene_quality`, the live
oncodigger.com application's own gene-quality policy module (Methods 3.3's "OncoDigger's
gene_quality.py policy version 2026-10-05-v2"), located via the GENE_QUALITY_TEST_ENGINE
environment variable. That module is not part of this repository, the manuscript's vendored
engine copy, or the public OncoDigger-app GitHub clone as currently checked out -- it could
not be located to independently re-run this script end-to-end as part of this recovery. Only
the already-committed output was verified, by cross-checking its number against the
manuscript's own cited figure.

formula_interaction_dirty_vocab.py -- tests whether the scoring formula's
contribution (real rank-weighted+enrichment formula vs. flat mention-count)
depends on vocabulary state, the open question flagged during this audit
("how are you measuring the 3rd variable, gene enrichment?"). The formula's
effect was previously measured ONLY under cleaned vocabulary (+1.6pp any-
cancer, matched_frequency_baseline.py: 93.2% real vs 91.6% flat). This script
computes the SAME comparison (real formula vs flat count) under TRUE-RAW
vocabulary instead, in a single pass so flat-count and real-formula rankings
come from the IDENTICAL BM25 pool and vocabulary construction (avoids any
cross-script pairing risk between true_raw_flat_count_2x2.py and
true_raw_vs_cleaned_all19.py, which used separately-constructed pools even
though both targeted the same query/options).

Real formula: associated_entity_tables's own ranking (rank-weighted mention
position x log(1+enrichment) x log(1+supporting papers)).
Flat count: same evidence/support-count data, re-ranked by raw mention count
only, independent (alphabetical) tie-breaking.

Both rankings computed from one BM25 pool fetch + one entity-matching pass
per cancer -- not two separate runs -- so this is NOT vulnerable to subtle
pool differences between scripts. Diagnostic only.

Run: python -X utf8 scripts/analysis/formula_interaction_dirty_vocab.py
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

QUERY = "mutation, mutations"
SEARCH_POOL, ANALYSIS_PAPERS = 1250, 1000
OPTIONS = {"year_min": 2000, "year_max": 2026, "title_match_multiplier": 1.2,
           "mesh_match_multiplier": 1.15, "recency_boost": True, "discovery_mode": "Primary Evidence"}
GENE_EXCLUDE = {"MICE", "DCR", "PC", "C2", "FBN1", "SGCG", "WDHD1"}
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


def rank_real_and_flat(con, vocab_lexicon: dict):
    """Single BM25 pool fetch + single matching pass; returns (real_top10, flat_top10)."""
    discovery._sqlite_doc_counts.clear()
    discovery._sqlite_caches.clear()
    tokens = normalize_query(QUERY)
    hits = bm25_search(con, tokens, SEARCH_POOL, OPTIONS)[:ANALYSIS_PAPERS]
    tables = associated_entity_tables(con, hits, QUERY, (2000, 2026), tokens, {"Genes": vocab_lexicon},
                                      limit=200, discovery_mode="Primary Evidence", corpus_query="")
    genes_df = tables.get("Genes")
    if genes_df is None or genes_df.empty:
        return [], []

    # Real formula: table's own order (already ranked by relevance).
    real_ranked, seen = [], set()
    for _, row in genes_df.iterrows():
        gene = GENE_REMAP.get(row["Gene"], row["Gene"])
        if gene in GENE_EXCLUDE or gene in seen:
            continue
        seen.add(gene)
        real_ranked.append(gene)
        if len(real_ranked) >= 10:
            break

    # Flat count: same rows, re-sorted by mention/support count, alphabetical tie-break.
    count_col = None
    for candidate in ("Supporting papers", "Supporting Papers", "Mentions", "Count", "Papers"):
        if candidate in genes_df.columns:
            count_col = candidate
            break
    if count_col is None:
        raise RuntimeError(f"Could not find a supporting-papers count column. Columns: {list(genes_df.columns)}")

    flat_rows = []
    seen2 = set()
    for _, row in genes_df.iterrows():
        gene = GENE_REMAP.get(row["Gene"], row["Gene"])
        if gene in GENE_EXCLUDE or gene in seen2:
            continue
        seen2.add(gene)
        count_val = row[count_col] if count_col else 0
        flat_rows.append((gene, count_val))
    flat_rows.sort(key=lambda gc: (-gc[1], gc[0]))
    flat_ranked = [g for g, _ in flat_rows[:10]]

    return real_ranked, flat_ranked, count_col


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
    print("Building true-raw vocabulary via gene_quality.py gene_lookups()...")
    raw_gq, vetted_gq, audit = gq.gene_lookups(DEFAULT_ONCODIGGER)
    genes_any, exact_pairs = load_intogen()

    rows = []
    real_any = real_exact = flat_any = flat_exact = 0
    count_col_seen = None
    for corpus_key, cancer in CANCER_MAP:
        db = DEFAULT_ONCODIGGER / "data" / "pubmed_processed" / f"{corpus_key}_all" / "lexical_index.db"
        con = sqlite3.connect(str(db))
        real_top10, flat_top10, count_col = rank_real_and_flat(con, raw_gq)
        count_col_seen = count_col
        r_any = sum(1 for g in real_top10 if g in genes_any)
        r_exact = sum(1 for g in real_top10 if (g, cancer) in exact_pairs)
        f_any = sum(1 for g in flat_top10 if g in genes_any)
        f_exact = sum(1 for g in flat_top10 if (g, cancer) in exact_pairs)
        real_any += r_any; real_exact += r_exact
        flat_any += f_any; flat_exact += f_exact
        print(f"{cancer:20} real={r_any}/10 any {r_exact}/10 exact  flat={f_any}/10 any {f_exact}/10 exact  "
              f"(count_col={count_col})")
        print(f"  real_top10: {real_top10}")
        print(f"  flat_top10: {flat_top10}")
        rows.append({"cancer": cancer, "real_top10": "|".join(real_top10), "real_any": r_any, "real_exact": r_exact,
                     "flat_top10": "|".join(flat_top10), "flat_any": f_any, "flat_exact": f_exact})

    print()
    print(f"count column used for flat ranking: {count_col_seen}")
    print(f"POOLED (true-raw vocab, BM25 pool, n=190):")
    print(f"  REAL formula: {real_any}/190 = {100*real_any/190:.1f}% any-cancer, "
          f"{real_exact}/190 = {100*real_exact/190:.1f}% exact")
    print(f"  FLAT count:   {flat_any}/190 = {100*flat_any/190:.1f}% any-cancer, "
          f"{flat_exact}/190 = {100*flat_exact/190:.1f}% exact")
    print(f"  Formula effect under DIRTY (true-raw) vocab: "
          f"{100*(real_any-flat_any)/190:+.1f}pp any-cancer, {100*(real_exact-flat_exact)/190:+.1f}pp exact")
    print()
    print("Reference, formula effect under CLEANED vocab (matched_frequency_baseline.py, already established):")
    print("  REAL (production): 93.2% any / 75.8% exact")
    print("  FLAT (matched baseline): 91.6% any / 76.3% exact")
    print("  Formula effect under CLEAN vocab: +1.6pp any-cancer, -0.5pp exact")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUT_DIR / "formula_interaction_dirty_vocab.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {OUT_DIR / 'formula_interaction_dirty_vocab.csv'}")


if __name__ == "__main__":
    main()
