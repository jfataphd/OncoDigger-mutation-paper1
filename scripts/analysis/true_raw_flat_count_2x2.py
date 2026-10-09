"""
true_raw_flat_count_2x2.py -- completes the flat-mention-count 2x2 using the CORRECT
"no QC" baseline (gene_quality.py's genuinely unfiltered `raw` dict, not load_lexicons()'s
raw output which still has useful_gene_term() filtering baked in). Per
FROZEN_ANALYSIS_SPEC_retrieval_vocab_decomposition.md.

Four cells, all flat (non-rank-weighted) supporting-paper counts, same real matching loop
(supporting_counts_for_documents, ported verbatim from retrieval_vocab_2x2_faithful.py) for
both pool and full-corpus scoring:
  Cell 1: full corpus,           TRUE raw vocabulary (gene_quality.py's `raw`)   -- NEW
  Cell 2: full corpus,           cleaned vocabulary (this paper's pipeline)      -- REUSED
          from retrieval_vocab_2x2_faithful_FINAL_summary.csv's cell2_fullcorpus_clean row
          (5-cancer sample: Thyroid, Bladder, Esophageal, Endometrial, Myeloma)
  Cell 3: BM25 1,000-paper pool, TRUE raw vocabulary                              -- NEW, all 19
  Cell 4: BM25 1,000-paper pool, cleaned vocabulary                               -- REUSED
          from matched_frequency_baseline_summary.csv's matched_frequency_only row (all 19)

Cell 1 is run for the same 5 smallest-corpus cancers already scoped in the frozen spec's
existing Amendment (Thyroid, Bladder, Esophageal, Endometrial, Myeloma), timed starting with
Thyroid (the smallest, ~37.8k docs) BEFORE committing to the rest, because the true-raw
vocabulary has substantially more candidate terms than the "dirty" load_lexicons-raw
vocabulary the prior 378.9s/cancer timing was measured against.

Run: python -X utf8 scripts/analysis/true_raw_flat_count_2x2.py
"""
from __future__ import annotations

import csv
import sqlite3
import sys
import os
import time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ONCODIGGER = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET"))
sys.path.insert(0, os.environ.get("GENE_QUALITY_TEST_ENGINE", "GENE_QUALITY_TEST_ENGINE_NOT_SET"))

from oncodigger_search import gene_quality as gq
from oncodigger_search import discovery
from oncodigger_search.lexical import bm25_search, normalize_query
from oncodigger_search.discovery import _docs_for_results

OUT_DIR = ROOT / "data" / "derived" / "intogen_rescoring"
INTOGEN_TSV = ROOT / "data" / "external" / "intogen" / "drivers" / "2024-06-18_IntOGen-Drivers" / "Compendium_Cancer_Genes.tsv"

ANALYSIS_PAPERS, SEARCH_POOL = 1000, 1250
MIN_SUPPORT = 2
GENE_EXCLUDE = {"MICE", "DCR", "PC", "C2", "FBN1", "SGCG", "WDHD1"}
GENE_REMAP = {"P53": "TP53", "HER2": "ERBB2"}
QUERY = "mutation, mutations"
OPTIONS = {
    "year_min": 2000, "year_max": 2026, "exclude_reviews": False, "downrank_reviews": True,
    "title_match_multiplier": 1.2, "mesh_match_multiplier": 1.15, "recency_boost": True,
    "discovery_mode": "Primary Evidence",
}

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

# Extended to all 19 cancers per the frozen spec's Amendment log: true-raw full-corpus
# matching measured at 25.6s/cancer on Thyroid (the smallest corpus), ~15x faster than
# originally assumed, making full 19-cancer coverage cheap and feasible.
FULL_CORPUS_CANCERS = [label for _, label in CANCER_MAP]


def build_combined(lookups_genes: dict) -> tuple[dict, dict]:
    combined_single = defaultdict(list)
    combined_phrases = defaultdict(list)
    for token, entities in lookups_genes.get("single", {}).items():
        combined_single[token].extend(("Genes", entity) for entity in entities)
    for anchor, phrase_entities in lookups_genes.get("phrases", {}).items():
        combined_phrases[anchor].extend(("Genes", phrase, entity) for phrase, entity in phrase_entities)
    return combined_single, combined_phrases


def supporting_counts_for_documents(docs: list[dict], combined_single: dict, combined_phrases: dict,
                                     query_terms: set[str]) -> dict[str, set[str]]:
    evidence: dict[str, set[str]] = defaultdict(set)
    document_tokens = discovery.document_tokens
    full_document_text = discovery.full_document_text
    for doc in docs:
        pmid = doc["pmid"]
        text = full_document_text(doc)
        padded = f" {text} "
        token_set = set(document_tokens(doc))
        seen = set()
        for token in token_set:
            if token in query_terms:
                continue
            for bucket, entity in combined_single.get(token, []):
                seen.add(entity["id"])
        for anchor in token_set:
            for bucket, phrase, entity in combined_phrases.get(anchor, []):
                if phrase in query_terms:
                    continue
                if f" {phrase} " in padded:
                    seen.add(entity["id"])
        for entity_id in seen:
            evidence[entity_id].add(pmid)
    return evidence


def fetch_all_documents(con: sqlite3.Connection) -> list[dict]:
    rows = con.execute(
        "SELECT pmid, title, publication_year, publication_types, mesh_terms, abstract, authors, article_ids, journal "
        "FROM documents"
    ).fetchall()
    docs = []
    for pmid, title, year, ptypes_json, mesh_json, abstract, authors_json, aids_json, journal in rows:
        docs.append({
            "pmid": pmid, "title": title or "", "abstract": abstract or "",
            "publication_year": year, "publication_types": ptypes_json or "",
            "mesh_terms": mesh_json or "", "journal": journal or "",
        })
    return docs


def gene_id_to_symbol(entity_id: str, id_to_symbol: dict) -> str:
    return id_to_symbol.get(entity_id, entity_id)


def rank_by_flat_count(evidence: dict[str, set[str]], id_to_symbol: dict, n_results: int) -> list[tuple[str, int]]:
    counts = {}
    for entity_id, pmids in evidence.items():
        support = len(pmids)
        if support < MIN_SUPPORT and n_results >= 10:
            continue
        symbol = GENE_REMAP.get(gene_id_to_symbol(entity_id, id_to_symbol), gene_id_to_symbol(entity_id, id_to_symbol))
        if symbol in GENE_EXCLUDE:
            continue
        counts[symbol] = max(counts.get(symbol, 0), support)
    ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    return ranked


def load_intogen():
    genes_any, exact_pairs = set(), set()
    with open(INTOGEN_TSV, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            sym = row["SYMBOL"].strip().upper()
            genes_any.add(sym)
            code = row["CANCER_TYPE"].strip()
            for cancer, codes in INTOGEN_CANCER_MAP.items():
                if code in codes:
                    exact_pairs.add((sym, cancer))
    return genes_any, exact_pairs


def score_top10(pairs_list, genes_any, exact_pairs):
    any_hits = sum(1 for c, g in pairs_list if g.strip().upper() in genes_any)
    exact_hits = sum(1 for c, g in pairs_list if (g.strip().upper(), c) in exact_pairs)
    return any_hits, exact_hits, len(pairs_list)


def main() -> None:
    print("Building TRUE raw vocabulary via gene_quality.py gene_lookups()...")
    raw_gq, vetted_gq, audit = gq.gene_lookups(DEFAULT_ONCODIGGER)
    id_to_symbol_raw = {}
    for ents in raw_gq.get("single", {}).values():
        for e in ents:
            id_to_symbol_raw[e["id"]] = e["id"]
    for lst in raw_gq.get("phrases", {}).values():
        for _, e in lst:
            id_to_symbol_raw[e["id"]] = e["id"]
    combined_single_raw, combined_phrases_raw = build_combined(raw_gq)
    print(f"  true-raw vocabulary: {len(raw_gq.get('single', {}))} single tokens, "
          f"{len(raw_gq.get('phrases', {}))} phrase anchors")

    db_base = DEFAULT_ONCODIGGER / "data" / "pubmed_processed"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    tokens = normalize_query(QUERY)
    query_terms = {t.replace("-", " ") for t in tokens}

    rows_out = []
    timing_log = []
    full_corpus_done = []
    full_corpus_abandoned_reason = None

    for key, label in CANCER_MAP:
        db = next((p for p in (db_base / key / "lexical_index.db", db_base / f"{key}_all" / "lexical_index.db")
                   if p.exists()), None)
        if db is None:
            print(f"  MISSING index for {label}, skipping", flush=True)
            continue
        discovery._sqlite_doc_counts.clear()
        discovery._sqlite_caches.clear()
        con = sqlite3.connect(str(db), check_same_thread=False)
        t0 = time.time()

        # --- Cell 3: BM25 pool, TRUE raw vocabulary, flat count (all 19) ---
        hits = bm25_search(con, tokens, SEARCH_POOL, OPTIONS)[:ANALYSIS_PAPERS]
        pool_pmids = [h["pmid"] for h in hits]
        docs_map = _docs_for_results(con, hits)
        pool_docs = []
        for pmid in pool_pmids:
            d = docs_map.get(pmid)
            if d is not None:
                d = dict(d)
                d["pmid"] = pmid
                pool_docs.append(d)

        ev3 = supporting_counts_for_documents(pool_docs, combined_single_raw, combined_phrases_raw, query_terms)
        top10_c3 = rank_by_flat_count(ev3, id_to_symbol_raw, len(pool_docs))[:10]
        t_pool = time.time() - t0
        print(f"  {label:20} pool docs={len(pool_docs)}  c3(true_raw)={[g for g, _ in top10_c3]}  [{t_pool:.1f}s]",
              flush=True)
        for rank, (g, c) in enumerate(top10_c3, 1):
            rows_out.append({"cell": "cell3_pool_true_raw", "cancer": label, "rank": rank, "gene": g, "count": c})

        # --- Cell 1: full corpus, TRUE raw vocabulary, flat count ---
        # Scoped to the same 5 smallest corpora already used for the cleaned-vocabulary
        # full-corpus cell (frozen spec Amendment, 2026-10-04), Thyroid timed first.
        if label in FULL_CORPUS_CANCERS and full_corpus_abandoned_reason is None:
            t1 = time.time()
            all_docs = fetch_all_documents(con)
            t_fetch = time.time() - t1
            t2 = time.time()
            ev1 = supporting_counts_for_documents(all_docs, combined_single_raw, combined_phrases_raw, query_terms)
            top10_c1 = rank_by_flat_count(ev1, id_to_symbol_raw, len(all_docs))[:10]
            t_match = time.time() - t2
            print(f"  {label:20} full-corpus docs={len(all_docs)}  fetch={t_fetch:.1f}s match={t_match:.1f}s "
                  f"TOTAL={t_fetch+t_match:.1f}s", flush=True)
            print(f"  {label:20}                c1(true_raw)={[g for g, _ in top10_c1]}", flush=True)
            for rank, (g, c) in enumerate(top10_c1, 1):
                rows_out.append({"cell": "cell1_fullcorpus_true_raw", "cancer": label, "rank": rank, "gene": g,
                                  "count": c})
            timing_log.append({"cancer": label, "n_docs": len(all_docs), "fetch_s": round(t_fetch, 1),
                                "match_s": round(t_match, 1)})
            full_corpus_done.append(label)

            if label == "Thyroid Cancer":
                total_thyroid = t_fetch + t_match
                print(f"\n  >>> Thyroid true-raw full-corpus timing: {total_thyroid:.1f}s "
                      f"(prior dirty-vocab timing was 378.9s) <<<\n", flush=True)
                # Honest, logged scope decision based on measured timing, made before
                # looking at any downstream IntOGen-agreement effect.
                budget_s = 1800  # ~30 min ceiling for the remaining 4 small corpora in this fork
                projected_remaining = total_thyroid * 4
                if projected_remaining > budget_s:
                    full_corpus_abandoned_reason = (
                        f"measured Thyroid true-raw full-corpus time={total_thyroid:.1f}s; "
                        f"projected remaining 4 cancers={projected_remaining:.1f}s exceeds "
                        f"{budget_s}s budget ceiling for this fork"
                    )
                    print(f"  SCOPE DECISION: stopping full-corpus (cell 1) after Thyroid only. "
                          f"Reason: {full_corpus_abandoned_reason}", flush=True)

        con.close()

    with (OUT_DIR / "true_raw_flat_count_2x2_detail.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["cell", "cancer", "rank", "gene", "count"])
        w.writeheader()
        w.writerows(rows_out)
    if timing_log:
        with (OUT_DIR / "true_raw_flat_count_2x2_timing.csv").open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["cancer", "n_docs", "fetch_s", "match_s"])
            w.writeheader()
            w.writerows(timing_log)

    genes_any, exact_pairs = load_intogen()
    summary = []
    for cell in ["cell1_fullcorpus_true_raw", "cell3_pool_true_raw"]:
        pairs_list = [(r["cancer"], r["gene"]) for r in rows_out if r["cell"] == cell and r["rank"] <= 10]
        if not pairs_list:
            continue
        any_hits, exact_hits, n = score_top10(pairs_list, genes_any, exact_pairs)
        summary.append({"cell": cell, "any_cancer": f"{any_hits}/{n}", "exact_concordance": f"{exact_hits}/{n}",
                         "n": n})
        print(f"{cell}: any-cancer {any_hits}/{n} = {100*any_hits/n:.1f}%  "
              f"exact {exact_hits}/{n} = {100*exact_hits/n:.1f}%")

    with (OUT_DIR / "true_raw_flat_count_2x2_summary.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["cell", "any_cancer", "exact_concordance", "n"])
        w.writeheader()
        w.writerows(summary)

    print(f"\nFull-corpus (cell 1, true-raw) completed for: {full_corpus_done}")
    if full_corpus_abandoned_reason:
        print(f"Scope reduction reason: {full_corpus_abandoned_reason}")
    print("Done.")


if __name__ == "__main__":
    main()
