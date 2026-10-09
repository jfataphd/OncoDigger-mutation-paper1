"""
retrieval_vocab_2x2_faithful.py -- manuscript1.1 (IntOGen-primary draft) ONLY. Production-
faithful replacement for retrieval_vocab_2x2.py, which used a hand-reimplemented retrieval
scorer that omitted the real "Primary Evidence" document multiplier. This script uses the
REAL vendored engine throughout (load_lexicons -> suppress_ambiguous_symbols ->
longest_match_lookups for vocabulary; bm25_search + associated_entity_tables's own entity-
matching loop for both retrieval-pool and full-corpus gene recognition), per
FROZEN_ANALYSIS_SPEC_retrieval_vocab_decomposition.md.

Four cells, flat (non-rank-weighted) supporting-paper counts, independent (alphabetical)
tie-breaking:
  Cell 1: full corpus,          less-cleaned vocabulary (raw load_lexicons() output)
  Cell 2: full corpus,          fully cleaned vocabulary (suppress + longest-match)
  Cell 3: BM25 1,000-paper pool, less-cleaned vocabulary
  Cell 4: BM25 1,000-paper pool, fully cleaned vocabulary (should closely match
          matched_frequency_baseline.py's existing 91.6%/76.3% result -- same real engine path)

For full-corpus cells, re-implements associated_entity_tables' core per-document matching
loop directly (same combined_single/combined_phrases lexicon structures, same evidence-set
construction), rather than calling associated_entity_tables with a synthetic 'results' list
of the whole corpus, because that function's SQL fetch is IN-clause-batched for pool-sized
(<=1250) requests, not full corpora (up to 245k documents for Breast Cancer). The matching
RULES are identical (same lexicon dicts built the same way); only the document-fetch and
looping is inlined for full-corpus scale. Supporting-paper count is a simple distinct-PMID
tally either way, independent of rank/order, so this does not change the matching semantics.

Run: python -X utf8 scripts/analysis/retrieval_vocab_2x2_faithful.py [--cancers "Lung Cancer" ...] [--skip-full-corpus]
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sqlite3
import sys
import time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ONCODIGGER = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET"))
OUT_DIR = ROOT / "data" / "derived" / "intogen_rescoring"
INTOGEN_TSV = ROOT / "data" / "external" / "intogen" / "drivers" / "2024-06-18_IntOGen-Drivers" / "Compendium_Cancer_Genes.tsv"

ANALYSIS_PAPERS, SEARCH_POOL = 1000, 1250
MIN_SUPPORT = 2
GENE_EXCLUDE = {"MICE", "DCR", "PC", "C2", "FBN1", "SGCG", "WDHD1"}
SUPPRESS_SYMBOLS = {"GC", "HR", "HCCS", "BPIFA4P", "MB", "ELN", "CP", "UBC", "CAMP", "EFS", "IVD", "SCT", "SON"}
GENE_REMAP = {"P53": "TP53", "HER2": "ERBB2"}
QUERY = "mutation mutations"
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


def suppress_ambiguous_symbols(genes_lexicon: dict) -> dict:
    single = {token: list(entities) for token, entities in genes_lexicon["single"].items()}
    for symbol in SUPPRESS_SYMBOLS:
        token = symbol.lower()
        entries = single.get(token)
        if not entries:
            continue
        kept = [e for e in entries if e["id"] != symbol]
        if kept:
            single[token] = kept
        elif token in single:
            del single[token]
    return {"single": single, "phrases": genes_lexicon["phrases"]}


def longest_match_lookups(genes_lexicon: dict, discovery) -> dict:
    import re
    names: dict[str, list] = {}
    for anchor, lst in genes_lexicon["phrases"].items():
        for term, entity in lst:
            names.setdefault(term, []).append(entity)
    ids = {t: {e["id"] for e in es} for t, es in names.items()}
    containers = set()
    for q in names:
        w = q.split()
        for i in range(len(w)):
            for j in range(i + 1, len(w) + 1):
                sub = " ".join(w[i:j])
                if (i, j) != (0, len(w)) and sub in ids and ids[sub] - ids[q]:
                    containers.add(q)
    pattern = re.compile(r"(?<![a-z0-9'])(?:" + "|".join(re.escape(q) for q in sorted(containers, key=len, reverse=True))
                         + r")(?![a-z0-9'])")
    phrases = {k: list(v) for k, v in genes_lexicon["phrases"].items()}
    for q in containers:
        for e in names[q]:
            phrases.setdefault(q.split()[0], []).append((q.replace(" ", "_"), e))
    original = discovery.full_document_text

    def masked_text(doc: dict) -> str:
        source = " ".join([doc.get("title", ""), doc.get("abstract", "")])
        if doc.get("_manuscript_masked_source") != source:
            doc["_oncodigger_full_text"] = pattern.sub(lambda m: m.group(0).replace(" ", "_"), original(doc))
            doc["_manuscript_masked_source"] = source
        return doc["_oncodigger_full_text"]

    discovery.full_document_text = masked_text
    return {"single": genes_lexicon["single"], "phrases": phrases}


def build_combined(lookups_genes: dict) -> tuple[dict, dict]:
    """Same structure associated_entity_tables builds internally from {'Genes': lookups_genes}."""
    combined_single = defaultdict(list)
    combined_phrases = defaultdict(list)
    for token, entities in lookups_genes.get("single", {}).items():
        combined_single[token].extend(("Genes", entity) for entity in entities)
    for anchor, phrase_entities in lookups_genes.get("phrases", {}).items():
        combined_phrases[anchor].extend(("Genes", phrase, entity) for phrase, entity in phrase_entities)
    return combined_single, combined_phrases


def supporting_counts_for_documents(discovery, docs: list[dict], combined_single: dict, combined_phrases: dict,
                                     query_terms: set[str]) -> dict[str, set[str]]:
    """Faithful reimplementation of associated_entity_tables' matching loop (discovery.py
    lines ~501-539), rank-weight omitted since only the flat evidence (supporting-paper) set
    is needed for this experiment. Returns {gene_id: {pmid, ...}}."""
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
    # Independent tie-break: alphabetical, not relevance-order.
    ranked = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    return ranked


def load_intogen():
    genes_any = set()
    pairs = set()
    with open(INTOGEN_TSV, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genes_any.add(row["SYMBOL"].strip().upper())
            code = row["CANCER_TYPE"].strip()
            for cancer, codes in INTOGEN_CANCER_MAP.items():
                if code in codes:
                    pairs.add((row["SYMBOL"].strip().upper(), cancer))
    return genes_any, pairs


def score_top10(pairs_list: list[tuple[str, str]], genes_any: set, exact_pairs: set) -> tuple[int, int, int]:
    any_hits = sum(1 for c, g in pairs_list if g.strip().upper() in genes_any)
    exact_hits = sum(1 for c, g in pairs_list if (g.strip().upper(), c) in exact_pairs)
    return any_hits, exact_hits, len(pairs_list)


def main() -> None:
    if os.environ.get("PYTHONHASHSEED") != "0":
        import subprocess
        env = dict(os.environ, PYTHONHASHSEED="0")
        sys.exit(subprocess.run([sys.executable, "-X", "utf8"] + sys.argv, env=env).returncode)

    ap = argparse.ArgumentParser()
    ap.add_argument("--oncodigger-root", type=Path, default=DEFAULT_ONCODIGGER)
    ap.add_argument("--cancers", nargs="+", default=None)
    ap.add_argument("--skip-full-corpus", action="store_true", help="only run cells 3/4 (BM25 pool)")
    args = ap.parse_args()

    sys.path.insert(0, str(args.oncodigger_root / "src"))
    from oncodigger.search import discovery
    from oncodigger.search.discovery import associated_entity_tables
    from oncodigger.search.lexical import bm25_search, normalize_query

    def fresh_lexicons():
        from oncodigger.search.lexicons import load_lexicons
        return load_lexicons(args.oncodigger_root)

    # Dirty = raw load_lexicons() output, no suppress/longest-match.
    lex_dirty_raw = fresh_lexicons()
    genes_dirty = lex_dirty_raw["Genes"]
    id_to_symbol_dirty = {}
    for ents in genes_dirty["single"].values():
        for e in ents:
            id_to_symbol_dirty[e["id"]] = e["id"]
    for lst in genes_dirty["phrases"].values():
        for _, e in lst:
            id_to_symbol_dirty[e["id"]] = e["id"]

    # Clean = suppress_ambiguous_symbols -> longest_match_lookups, same as production.
    lex_clean_raw = fresh_lexicons()
    genes_suppressed = suppress_ambiguous_symbols(lex_clean_raw["Genes"])
    genes_clean = longest_match_lookups(genes_suppressed, discovery)
    id_to_symbol_clean = dict(id_to_symbol_dirty)

    combined_single_dirty, combined_phrases_dirty = build_combined(genes_dirty)
    combined_single_clean, combined_phrases_clean = build_combined(genes_clean)

    db_base = args.oncodigger_root / "data" / "pubmed_processed"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    tokens = normalize_query(QUERY)
    query_terms = {t.replace("-", " ") for t in tokens}

    rows_out = []
    timing_log = []

    for key, label in CANCER_MAP:
        if args.cancers and label not in args.cancers:
            continue
        db = next((p for p in (db_base / key / "lexical_index.db", db_base / f"{key}_all" / "lexical_index.db")
                   if p.exists()), None)
        if db is None:
            print(f"  MISSING index for {label}, skipping", flush=True)
            continue
        discovery._sqlite_doc_counts.clear()
        discovery._sqlite_caches.clear()
        con = sqlite3.connect(str(db), check_same_thread=False)
        t0 = time.time()

        # --- Cells 3 & 4: BM25 pool (real retrieval incl. Primary Evidence multiplier) ---
        hits = bm25_search(con, tokens, SEARCH_POOL, OPTIONS)[:ANALYSIS_PAPERS]
        pool_pmids = [h["pmid"] for h in hits]
        from oncodigger.search.discovery import _docs_for_results
        pool_docs = list(_docs_for_results(con, hits).values())
        for d in pool_docs:
            d.setdefault("pmid", None)
        # _docs_for_results doesn't set 'pmid' key on the dict itself; re-key via pool order.
        pool_docs = []
        docs_map = _docs_for_results(con, hits)
        for pmid in pool_pmids:
            d = docs_map.get(pmid)
            if d is not None:
                d = dict(d)
                d["pmid"] = pmid
                pool_docs.append(d)

        ev3 = supporting_counts_for_documents(discovery, pool_docs, combined_single_dirty, combined_phrases_dirty, query_terms)
        ev4 = supporting_counts_for_documents(discovery, pool_docs, combined_single_clean, combined_phrases_clean, query_terms)
        top10_c3 = rank_by_flat_count(ev3, id_to_symbol_dirty, len(pool_docs))[:10]
        top10_c4 = rank_by_flat_count(ev4, id_to_symbol_clean, len(pool_docs))[:10]
        t_pool = time.time() - t0

        print(f"  {label:20} pool docs={len(pool_docs)}  c3(dirty)={[g for g,_ in top10_c3]}", flush=True)
        print(f"  {label:20}                c4(clean)={[g for g,_ in top10_c4]}  [{t_pool:.1f}s]", flush=True)

        for rank, (g, c) in enumerate(top10_c3, 1):
            rows_out.append({"cell": "cell3_pool_dirty", "cancer": label, "rank": rank, "gene": g, "count": c})
        for rank, (g, c) in enumerate(top10_c4, 1):
            rows_out.append({"cell": "cell4_pool_clean", "cancer": label, "rank": rank, "gene": g, "count": c})

        # --- Cells 1 & 2: full corpus ---
        if not args.skip_full_corpus:
            t1 = time.time()
            all_docs = fetch_all_documents(con)
            t_fetch = time.time() - t1
            t2 = time.time()
            ev1 = supporting_counts_for_documents(discovery, all_docs, combined_single_dirty, combined_phrases_dirty, query_terms)
            ev2 = supporting_counts_for_documents(discovery, all_docs, combined_single_clean, combined_phrases_clean, query_terms)
            top10_c1 = rank_by_flat_count(ev1, id_to_symbol_dirty, len(all_docs))[:10]
            top10_c2 = rank_by_flat_count(ev2, id_to_symbol_clean, len(all_docs))[:10]
            t_match = time.time() - t2
            print(f"  {label:20} full-corpus docs={len(all_docs)}  fetch={t_fetch:.1f}s match={t_match:.1f}s", flush=True)
            print(f"  {label:20}                c1(dirty)={[g for g,_ in top10_c1]}", flush=True)
            print(f"  {label:20}                c2(clean)={[g for g,_ in top10_c2]}", flush=True)
            for rank, (g, c) in enumerate(top10_c1, 1):
                rows_out.append({"cell": "cell1_fullcorpus_dirty", "cancer": label, "rank": rank, "gene": g, "count": c})
            for rank, (g, c) in enumerate(top10_c2, 1):
                rows_out.append({"cell": "cell2_fullcorpus_clean", "cancer": label, "rank": rank, "gene": g, "count": c})
            timing_log.append({"cancer": label, "n_docs": len(all_docs), "fetch_s": round(t_fetch, 1), "match_s": round(t_match, 1)})

        con.close()

    with (OUT_DIR / "retrieval_vocab_2x2_faithful_detail.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["cell", "cancer", "rank", "gene", "count"])
        w.writeheader()
        w.writerows(rows_out)
    if timing_log:
        with (OUT_DIR / "retrieval_vocab_2x2_faithful_timing.csv").open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["cancer", "n_docs", "fetch_s", "match_s"])
            w.writeheader()
            w.writerows(timing_log)

    genes_any, exact_pairs = load_intogen()
    summary = []
    for cell in ["cell1_fullcorpus_dirty", "cell2_fullcorpus_clean", "cell3_pool_dirty", "cell4_pool_clean"]:
        pairs_list = [(r["cancer"], r["gene"]) for r in rows_out if r["cell"] == cell and r["rank"] <= 10]
        if not pairs_list:
            continue
        any_hits, exact_hits, n = score_top10(pairs_list, genes_any, exact_pairs)
        summary.append({"cell": cell, "any_cancer": f"{any_hits}/{n}", "exact_concordance": f"{exact_hits}/{n}", "n": n})
        print(f"{cell}: any-cancer {any_hits}/{n} = {100*any_hits/n:.1f}%  exact {exact_hits}/{n} = {100*exact_hits/n:.1f}%")

    with (OUT_DIR / "retrieval_vocab_2x2_faithful_summary.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["cell", "any_cancer", "exact_concordance", "n"])
        w.writeheader()
        w.writerows(summary)
    print("Done.")


if __name__ == "__main__":
    main()
