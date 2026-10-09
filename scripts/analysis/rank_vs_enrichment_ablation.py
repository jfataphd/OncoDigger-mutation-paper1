"""
rank_vs_enrichment_ablation.py -- isolates rank-weighting vs. enrichment within
the real scoring formula (ChatGPT review's requested cheap ablation), under
true-raw vocabulary + the same frozen BM25 pool used throughout this audit.
Diagnostic only.

Real formula: S_g = W_g * log(1+E_g) * log(1+n_g), where
  W_g = sum over supporting docs of 1/sqrt(rank_in_pool + 1)   (rank-weighted mention score)
  E_g = (n_g/relevant_papers) / (background_g/corpus_papers)  (enrichment)
  n_g = distinct supporting-paper count (same evidence set used throughout)
Reproduced directly from discovery.py (lines ~500-553) and association.py's
enrichment_score, not re-derived -- same formula, same BM25 pool, same
matching loop, just computed once per cancer and re-ranked four ways:

  1. Flat count:        rank by n_g alone                          (reference: 28/190, 10/190)
  2. Rank-weighting only: rank by W_g * log(1+n_g), E_g forced = 1  (drops enrichment)
  3. Enrichment only:     rank by log(1+E_g) * log(1+n_g), W_g = 1  (drops rank-weighting)
  4. Full formula:        rank by W_g * log(1+E_g) * log(1+n_g)     (reference: 84/190, 52/190)

Min-support threshold (n_g >= 2 when pool >= 10 docs) applied identically to
all four, matching production. Alphabetical tie-breaking throughout.

Run: python -X utf8 scripts/analysis/rank_vs_enrichment_ablation.py
"""
from __future__ import annotations

import csv
import math
import sys
import os
from pathlib import Path
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ONCODIGGER = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET"))
sys.path.insert(0, os.environ.get("GENE_QUALITY_TEST_ENGINE", "GENE_QUALITY_TEST_ENGINE_NOT_SET"))

import sqlite3

from oncodigger_search import gene_quality as gq
from oncodigger_search import discovery
from oncodigger_search.association import entity_background_papers
from oncodigger_search.discovery import _docs_for_results
from oncodigger_search.lexical import bm25_search, normalize_query

QUERY = "mutation, mutations"
SEARCH_POOL, ANALYSIS_PAPERS = 1250, 1000
OPTIONS = {"year_min": 2000, "year_max": 2026, "title_match_multiplier": 1.2,
           "mesh_match_multiplier": 1.15, "recency_boost": True, "discovery_mode": "Primary Evidence"}
GENE_EXCLUDE = {"MICE", "DCR", "PC", "C2", "FBN1", "SGCG", "WDHD1"}
GENE_REMAP = {"P53": "TP53", "HER2": "ERBB2"}
OUT_DIR = ROOT / "data" / "derived" / "intogen_rescoring"
MIN_SUPPORT = 2

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


def build_combined(genes_lexicon: dict):
    combined_single = defaultdict(list)
    combined_phrases = defaultdict(list)
    for token, entities in genes_lexicon.get("single", {}).items():
        combined_single[token].extend(entities)
    for anchor, phrase_entities in genes_lexicon.get("phrases", {}).items():
        combined_phrases[anchor].extend(phrase_entities)
    return combined_single, combined_phrases


def compute_raw_terms(con, vocab_lexicon: dict):
    """Returns {gene_id: (weighted_score, support_count, background_count)} for one cancer."""
    discovery._sqlite_doc_counts.clear()
    discovery._sqlite_caches.clear()
    tokens = normalize_query(QUERY)
    query_terms = {t.replace("-", " ") for t in tokens}
    hits = bm25_search(con, tokens, SEARCH_POOL, OPTIONS)[:ANALYSIS_PAPERS]
    docs_map = _docs_for_results(con, hits)
    pool_docs = []
    for h in hits:
        d = docs_map.get(h["pmid"])
        if d is not None:
            d = dict(d)
            d["pmid"] = h["pmid"]
            pool_docs.append(d)

    combined_single, combined_phrases = build_combined(vocab_lexicon)
    document_tokens = discovery.document_tokens
    full_document_text = discovery.full_document_text

    scores = defaultdict(float)
    evidence = defaultdict(set)
    matched_terms = defaultdict(set)

    for rank, doc in enumerate(pool_docs):
        text = full_document_text(doc)
        padded = f" {text} "
        token_set = set(document_tokens(doc))
        weight = 1.0 / math.sqrt(rank + 1)
        seen = set()
        for token in token_set:
            if token in query_terms:
                continue
            for entity in combined_single.get(token, []):
                gene = entity["id"]
                seen.add(gene)
                matched_terms[gene].add(token)
        for anchor in token_set:
            for phrase, entity in combined_phrases.get(anchor, []):
                if phrase in query_terms:
                    continue
                if f" {phrase} " in padded:
                    gene = entity["id"]
                    seen.add(gene)
                    matched_terms[gene].add(phrase)
        for gene in seen:
            scores[gene] += weight
            evidence[gene].add(doc["pmid"])

    corpus_doc_count = discovery._doc_count(con)
    results_out = {}
    for gene, w_score in scores.items():
        support_count = len(evidence[gene])
        if support_count < MIN_SUPPORT and len(pool_docs) >= 10:
            continue
        background_count = entity_background_papers(con, list(matched_terms[gene])) or support_count
        results_out[gene] = (w_score, support_count, background_count, len(pool_docs), corpus_doc_count)
    return results_out


def enrichment_value(support_count, relevant_papers, background_count, corpus_papers):
    corpus_papers = max(corpus_papers, 1)
    background_count = max(background_count, support_count, 1)
    relevant_fraction = support_count / relevant_papers
    background_fraction = background_count / corpus_papers
    return relevant_fraction / max(background_fraction, 1 / corpus_papers)


def rank_top10(gene_terms: dict, mode: str):
    scored = []
    for gene, (w_score, support_count, background_count, relevant_papers, corpus_papers) in gene_terms.items():
        if gene in GENE_EXCLUDE:
            continue
        gene_r = GENE_REMAP.get(gene, gene)
        e_val = enrichment_value(support_count, relevant_papers, background_count, corpus_papers)
        if mode == "flat":
            key = (support_count, 0.0)
        elif mode == "rank_only":
            key = (w_score * math.log1p(support_count), 0.0)
        elif mode == "enrichment_only":
            key = (math.log1p(e_val) * math.log1p(support_count), 0.0)
        elif mode == "full":
            key = (w_score * math.log1p(e_val) * math.log1p(support_count), 0.0)
        else:
            raise ValueError(mode)
        scored.append((gene_r, key[0]))
    # merge duplicate remapped genes (e.g. P53->TP53) by max score
    best = {}
    for gene_r, val in scored:
        if gene_r not in best or val > best[gene_r]:
            best[gene_r] = val
    scored2 = sorted(best.items(), key=lambda gv: (-gv[1], gv[0]))
    return [g for g, _ in scored2[:10]]


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


def main():
    print("Building true-raw vocabulary...")
    raw_gq, vetted_gq, audit = gq.gene_lookups(DEFAULT_ONCODIGGER)
    genes_any, exact_pairs = load_intogen()

    totals = {m: [0, 0] for m in ("flat", "rank_only", "enrichment_only", "full")}
    rows = []
    for corpus_key, cancer in CANCER_MAP:
        db = DEFAULT_ONCODIGGER / "data" / "pubmed_processed" / f"{corpus_key}_all" / "lexical_index.db"
        con = sqlite3.connect(str(db))
        terms = compute_raw_terms(con, raw_gq)
        row = {"cancer": cancer}
        line = f"{cancer:20}"
        for mode in ("flat", "rank_only", "enrichment_only", "full"):
            top10 = rank_top10(terms, mode)
            a = sum(1 for g in top10 if g in genes_any)
            e = sum(1 for g in top10 if (g, cancer) in exact_pairs)
            totals[mode][0] += a
            totals[mode][1] += e
            row[f"{mode}_any"] = a
            row[f"{mode}_exact"] = e
            row[f"{mode}_top10"] = "|".join(top10)
            line += f"  {mode}={a}/10,{e}/10"
        print(line)
        rows.append(row)

    print()
    print("POOLED (true-raw vocab, BM25 pool, n=190):")
    for mode in ("flat", "rank_only", "enrichment_only", "full"):
        a, e = totals[mode]
        print(f"  {mode:16} any={a}/190={100*a/190:.1f}%  exact={e}/190={100*e/190:.1f}%")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUT_DIR / "rank_vs_enrichment_ablation.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {OUT_DIR / 'rank_vs_enrichment_ablation.csv'}")


if __name__ == "__main__":
    main()
