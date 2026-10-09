"""
matched_frequency_baseline.py -- manuscript1.1 (IntOGen-primary draft) ONLY. Builds a
frequency-only ranking control that holds BM25 retrieval AND gene vocabulary IDENTICAL to
the real pipeline, changing ONLY the gene-ranking formula, to isolate what rank-weighting
and enrichment specifically contribute. This replaces the existing full-corpus frequency
baseline (data/derived/freq_baseline_top50_aliases.csv), whose 43.2% figure is confounded:
it counts gene mentions across each cancer's ENTIRE corpus (tens to hundreds of thousands
of abstracts), not the same 1,000-paper BM25-selected analysis pool the real pipeline uses,
so its gap from 93.2% mixes retrieval, weighting, enrichment and vocabulary QC together.

Method: run the exact same retrieval and gene-matching step the real pipeline uses
(run_query_planes.py's rank_genes(), via associated_entity_tables()) for each of the 19
cancers, under the mutation plane. associated_entity_tables() already computes, for every
gene passing its >=2-supporting-paper floor, both:
  - "Supporting papers": a flat count of distinct pool papers mentioning the gene (no
    rank-weighting, no enrichment) -- exactly the frequency-only quantity needed here.
  - "Relevance": the real pipeline's rank-weighted, enrichment-normalized score.
Both come from the identical retrieval pool and identical cleaned vocabulary (load_lexicons
-> suppress_ambiguous_symbols -> longest_match_lookups, GENE_EXCLUDE applied), so re-ranking
the SAME computed table by "Supporting papers" instead of "Relevance" isolates the ranking
formula's own contribution with nothing else varying.

Sanity check: re-deriving the real pipeline's own ranking from this same run (sorting by
Relevance, same as run_query_planes.py does) must reproduce data/canonical/
precision_at_k_detail.csv's top-10 exactly, confirming this script's retrieval and matching
are faithful to the real pipeline before trusting the frequency-only re-ranking.

Outputs (new files, does not touch or overwrite any existing canonical/derived file):
  data/derived/intogen_rescoring/matched_frequency_baseline_detail.csv
  data/derived/intogen_rescoring/matched_frequency_baseline_summary.csv
  data/derived/intogen_rescoring/matched_pipeline_sanity_check.csv

Run: python -X utf8 scripts/analysis/matched_frequency_baseline.py [--cancers "Lung Cancer" ...]
"""
from __future__ import annotations

import argparse
import csv
import os
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ONCODIGGER = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET"))
OUT_DIR = ROOT / "data" / "derived" / "intogen_rescoring"
INTOGEN_TSV = ROOT / "data" / "external" / "intogen" / "drivers" / "2024-06-18_IntOGen-Drivers" / "Compendium_Cancer_Genes.tsv"

TOP_K, LIMIT, ANALYSIS_PAPERS, SEARCH_POOL = 50, 100, 1000, 1250
MIN_SUPPORT = 2  # same >=2-supporting-paper floor used throughout this paper
GENE_EXCLUDE = {"MICE", "DCR", "PC", "C2", "FBN1", "SGCG", "WDHD1"}
SUPPRESS_SYMBOLS = {"GC", "HR", "HCCS", "BPIFA4P", "MB", "ELN", "CP", "UBC", "CAMP", "EFS", "IVD", "SCT", "SON"}
GENE_REMAP = {"P53": "TP53", "HER2": "ERBB2"}
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
OPTIONS = {
    "year_min": 2000, "year_max": 2026, "exclude_reviews": False, "downrank_reviews": True,
    "title_match_multiplier": 1.2, "mesh_match_multiplier": 1.15, "recency_boost": True,
    "discovery_mode": "Primary Evidence",
}
QUERY = "mutation mutations"

# IntOGen cancer-code mapping, identical to the one used throughout this paper
# (Supplementary Table S8 / intogen_full_rescoring.py).
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


def open_corpus(db_path, discovery) -> sqlite3.Connection:
    discovery._sqlite_doc_counts.clear()
    discovery._sqlite_caches.clear()
    return sqlite3.connect(str(db_path), check_same_thread=False)


def rank_both_ways(con, engine) -> tuple[list[str], list[str], list[dict]]:
    """Returns (relevance_ranked_top10, frequency_ranked_top10, full_detail_rows)."""
    bm25_search, normalize_query, associated_entity_tables, gene_lookups = engine
    tokens = normalize_query(QUERY)
    hits = bm25_search(con, tokens, SEARCH_POOL, OPTIONS)[:ANALYSIS_PAPERS]
    tables = associated_entity_tables(con, hits, QUERY, (2000, 2026), tokens, gene_lookups,
                                       limit=LIMIT, discovery_mode="Primary Evidence", corpus_query="")
    genes_df = tables.get("Genes")
    if genes_df is None or genes_df.empty:
        return [], [], []

    detail = []
    seen_genes = set()
    for _, row in genes_df.iterrows():
        gene = GENE_REMAP.get(row["Gene"], row["Gene"])
        if gene in GENE_EXCLUDE or gene in seen_genes:
            continue
        seen_genes.add(gene)
        detail.append({
            "gene": gene,
            "relevance": float(row["Relevance"]),
            "supporting_papers": int(row["Supporting papers"]),
        })

    # Relevance ranking: genes_df already arrives Relevance-sorted (associated_entity_tables'
    # own enriched_entities.sort by Relevance descending); dedup above preserves that order.
    relevance_ranked = [d["gene"] for d in detail][:TOP_K]

    # Frequency-only ranking: identical gene set and supporting_papers counts, re-sorted by
    # that flat count alone. Ties broken by relevance rank (stable sort preserves original
    # relative order for equal supporting_papers), a deterministic, documented tie rule.
    freq_ranked = [d["gene"] for d in sorted(detail, key=lambda d: d["supporting_papers"], reverse=True)][:TOP_K]

    return relevance_ranked[:10], freq_ranked[:10], detail


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


def main() -> None:
    if os.environ.get("PYTHONHASHSEED") != "0":
        import subprocess
        env = dict(os.environ, PYTHONHASHSEED="0")
        sys.exit(subprocess.run([sys.executable, "-X", "utf8"] + sys.argv, env=env).returncode)

    ap = argparse.ArgumentParser()
    ap.add_argument("--oncodigger-root", type=Path, default=DEFAULT_ONCODIGGER)
    ap.add_argument("--cancers", nargs="+", default=None)
    args = ap.parse_args()

    sys.path.insert(0, str(args.oncodigger_root / "src"))
    from oncodigger.search import discovery
    from oncodigger.search.discovery import associated_entity_tables
    from oncodigger.search.lexical import bm25_search, normalize_query
    from oncodigger.search.lexicons import load_lexicons

    lexicons = load_lexicons(args.oncodigger_root)
    genes_suppressed = suppress_ambiguous_symbols(lexicons["Genes"])
    genes = longest_match_lookups(genes_suppressed, discovery)
    engine = (bm25_search, normalize_query, associated_entity_tables, {"Genes": genes})

    db_base = args.oncodigger_root / "data" / "pubmed_processed"
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    sanity_rows = []
    detail_rows = []
    for key, label in CANCER_MAP:
        if args.cancers and label not in args.cancers:
            continue
        db = next((p for p in (db_base / key / "lexical_index.db", db_base / f"{key}_all" / "lexical_index.db")
                   if p.exists()), None)
        if db is None:
            print(f"  MISSING index for {label}, skipping", flush=True)
            continue
        con = open_corpus(db, discovery)
        rel_top10, freq_top10, detail = rank_both_ways(con, engine)
        con.close()

        for rank, gene in enumerate(rel_top10, 1):
            sanity_rows.append({"cancer": label, "rank": rank, "gene": gene})
        for rank, gene in enumerate(freq_top10, 1):
            detail_rows.append({"cancer": label, "rank": rank, "gene": gene})

        print(f"  {label:20} relevance-top10: {' '.join(rel_top10)}", flush=True)
        print(f"  {label:20} freq-only-top10: {' '.join(freq_top10)}", flush=True)

    with (OUT_DIR / "matched_pipeline_sanity_check.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["cancer", "rank", "gene"])
        w.writeheader()
        w.writerows(sanity_rows)

    with (OUT_DIR / "matched_frequency_baseline_detail.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["cancer", "rank", "gene"])
        w.writeheader()
        w.writerows(detail_rows)

    # Score both the sanity-check (relevance) ranking and the frequency-only ranking against
    # IntOGen, using the identical mapping and methodology as the rest of this paper.
    genes_any, exact_pairs = load_intogen()

    def score(rows):
        top10 = [(r["cancer"], r["gene"]) for r in rows if r["rank"] <= 10]
        any_hits = sum(1 for c, g in top10 if g.strip().upper() in genes_any)
        exact_hits = sum(1 for c, g in top10 if (g.strip().upper(), c) in exact_pairs)
        n = len(top10)
        return any_hits, exact_hits, n

    rel_any, rel_exact, rel_n = score(sanity_rows)
    freq_any, freq_exact, freq_n = score(detail_rows)

    with (OUT_DIR / "matched_frequency_baseline_summary.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["ranking", "precision_at_10_any", "exact_concordance", "n"])
        w.writerow(["real_pipeline_sanity_check", f"{rel_any}/{rel_n}", f"{rel_exact}/{rel_n}", rel_n])
        w.writerow(["matched_frequency_only", f"{freq_any}/{freq_n}", f"{freq_exact}/{freq_n}", freq_n])

    print()
    print(f"Sanity check (relevance ranking, should match canonical 177/190): {rel_any}/{rel_n} = {100*rel_any/rel_n:.1f}% any-cancer, {rel_exact}/{rel_n} = {100*rel_exact/rel_n:.1f}% exact concordance")
    print(f"Matched frequency-only baseline: {freq_any}/{freq_n} = {100*freq_any/freq_n:.1f}% any-cancer, {freq_exact}/{freq_n} = {100*freq_exact/freq_n:.1f}% exact concordance")


if __name__ == "__main__":
    main()
