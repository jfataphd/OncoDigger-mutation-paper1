"""
intogen_parameter_sensitivity.py -- manuscript1.1 (IntOGen-primary draft) ONLY. Mirrors
parameter_sensitivity.py exactly (same 24-config one-at-a-time sweep over k1, b,
title_match_multiplier, mesh_match_multiplier, review_multiplier, analysis-paper cutoff,
over all 19 cancers) except Precision@10 is graded against IntOGen's any-cancer driver
set instead of COSMIC CGC Tier 1. Does not modify or depend on parameter_sensitivity.py
or cosmic_census_2026-05-23.csv; only which gene set counts as a "hit" differs.

Run:  python -X utf8 scripts/analysis/intogen_parameter_sensitivity.py
Output: data/derived/intogen_rescoring/intogen_parameter_sensitivity.csv
        (one row per config: param, value, precision_at_10, n_hits, n_total)
"""
from __future__ import annotations

import csv
import os
import re
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ONCODIGGER = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET"))
OUT = ROOT / "data" / "derived" / "intogen_rescoring" / "intogen_parameter_sensitivity.csv"
INTOGEN_DRIVERS = ROOT / "data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"

GENE_EXCLUDE = {
    "GC", "HR", "MICE", "DCR", "BPIFA4P", "PC", "HCCS", "WDHD1", "SGCG", "C2", "FBN1",
    # 2026-09-30 false-alias sweep -- kept in sync with run_query_planes.py's GENE_EXCLUDE.
    "MB", "ELN", "CP", "UBC", "CAMP", "EFS", "IVD", "SCT", "SON",
}
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
QUERY = "mutation, mutations"
BASELINE_OPTIONS = {
    "year_min": 2000, "year_max": 2026, "exclude_reviews": False, "downrank_reviews": True,
    "title_match_multiplier": 1.2, "mesh_match_multiplier": 1.15, "recency_boost": True,
    "discovery_mode": "Primary Evidence",
}
BASELINE_K1, BASELINE_B = 1.5, 0.75
BASELINE_ANALYSIS_PAPERS, BASELINE_SEARCH_POOL = 1000, 1250

SWEEP = {
    "k1": [1.0, 1.2, BASELINE_K1, 1.8, 2.0],
    "b": [0.5, 0.65, BASELINE_B, 0.85, 1.0],
    "title_match_multiplier": [1.0, 1.1, 1.2, 1.3, 1.5],
    "mesh_match_multiplier": [1.0, 1.05, 1.15, 1.25, 1.35],
    "review_multiplier": [0.5, 0.65, 0.8, 1.0],
    "analysis_papers": [500, 750, 1000, 1500, 2000],
}


def load_intogen_any_genes() -> set:
    genes = set()
    with open(INTOGEN_DRIVERS, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genes.add(row["SYMBOL"].upper())
    return genes


def build_patched_lexical(oncodigger_root: Path, sweep_params: dict):
    """Exec a private copy of lexical.py with k1/b sourced from sweep_params (mutable dict)."""
    src_path = oncodigger_root / "src" / "oncodigger" / "search" / "lexical.py"
    source = src_path.read_text(encoding="utf-8")
    patched = source.replace("k1 = 1.5", "k1 = sweep_params['k1']")
    patched = patched.replace("b = 0.75", "b = sweep_params['b']")
    assert patched.count("sweep_params['k1']") == 2, "expected 2 k1 patch sites"
    assert patched.count("sweep_params['b']") == 2, "expected 2 b patch sites"
    patched = patched.replace(
        "from .discovery_modes import document_mode_multiplier",
        "from oncodigger.search.discovery_modes import document_mode_multiplier",
    )
    namespace = {"sweep_params": sweep_params, "__name__": "oncodigger.search._lexical_patched"}
    exec(compile(patched, str(src_path) + " (patched)", "exec"), namespace)
    return namespace["bm25_search"], namespace["normalize_query"]


def longest_match_lookups(genes_lexicon: dict, discovery) -> dict:
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


def open_corpus(db_path, discovery):
    import sqlite3
    discovery._sqlite_doc_counts.clear()
    discovery._sqlite_caches.clear()
    return sqlite3.connect(str(db_path), check_same_thread=False)


def rank_top10(con, bm25_search, normalize_query, associated_entity_tables, genes_lexicon,
                options: dict, analysis_papers: int, search_pool: int) -> list[str]:
    tokens = normalize_query(QUERY)
    hits = bm25_search(con, tokens, search_pool, options)[:analysis_papers]
    tables = associated_entity_tables(con, hits, QUERY, (2000, 2026), tokens, {"Genes": genes_lexicon},
                                      limit=100, discovery_mode="Primary Evidence", corpus_query="")
    genes_df = tables.get("Genes")
    ranked, seen = [], set()
    if genes_df is None or genes_df.empty:
        return ranked
    for _, row in genes_df.iterrows():
        gene = GENE_REMAP.get(row["Gene"], row["Gene"])
        if gene in GENE_EXCLUDE or gene in seen:
            continue
        try:
            n = con.execute("SELECT COUNT(DISTINCT pmid) FROM postings WHERE token = ?", (gene.lower(),)).fetchone()[0]
        except Exception:
            n = 0
        if n == 0:
            continue
        seen.add(gene)
        ranked.append(gene)
        if len(ranked) >= 10:
            break
    return ranked


def run_config(label: str, value, options: dict, sweep_params: dict, analysis_papers: int, search_pool: int,
                bm25_search, normalize_query, associated_entity_tables, genes_lexicon, tier1_genes: set,
                db_base: Path, discovery) -> dict:
    hits = 0
    total = 0
    per_cancer = []
    for key, cancer_label in CANCER_MAP:
        db = next((p for p in (db_base / key / "lexical_index.db", db_base / f"{key}_all" / "lexical_index.db")
                   if p.exists()), None)
        if db is None:
            sys.exit(f"missing index for {cancer_label}")
        con = open_corpus(db, discovery)
        top10 = rank_top10(con, bm25_search, normalize_query, associated_entity_tables, genes_lexicon,
                            options, analysis_papers, search_pool)
        con.close()
        n_hit = sum(1 for g in top10 if g.upper() in tier1_genes)
        hits += n_hit
        total += len(top10)
        per_cancer.append((cancer_label, n_hit, len(top10)))
    precision = hits / total if total else 0.0
    print(f"  {label}={value!r}: Precision@10 = {precision:.4f} ({hits}/{total})", flush=True)
    return {"param": label, "value": value, "hits": hits, "total": total, "precision_at_10": precision}


def main() -> None:
    if os.environ.get("PYTHONHASHSEED") != "0":
        import subprocess
        env = dict(os.environ, PYTHONHASHSEED="0")
        sys.exit(subprocess.run([sys.executable, "-X", "utf8"] + sys.argv, env=env).returncode)

    oncodigger_root = DEFAULT_ONCODIGGER
    sys.path.insert(0, str(oncodigger_root / "src"))
    from oncodigger.search import discovery
    from oncodigger.search.discovery import associated_entity_tables
    from oncodigger.search.lexicons import load_lexicons

    tier1_genes = load_intogen_any_genes()
    print(f"Loaded {len(tier1_genes)} IntOGen driver gene symbols (any cancer type)", flush=True)

    lexicons = load_lexicons(oncodigger_root)
    genes_lexicon = longest_match_lookups(lexicons["Genes"], discovery)
    db_base = oncodigger_root / "data" / "pubmed_processed"

    sweep_params = {"k1": BASELINE_K1, "b": BASELINE_B}
    bm25_search, normalize_query = build_patched_lexical(oncodigger_root, sweep_params)

    results = []
    t_start = time.time()

    # Baseline (shared reference point for every parameter's sweep)
    print("\n=== Baseline ===", flush=True)
    sweep_params["k1"], sweep_params["b"] = BASELINE_K1, BASELINE_B
    base_options = dict(BASELINE_OPTIONS)
    baseline_result = run_config("baseline", "current", base_options, sweep_params,
                                  BASELINE_ANALYSIS_PAPERS, BASELINE_SEARCH_POOL,
                                  bm25_search, normalize_query, associated_entity_tables,
                                  genes_lexicon, tier1_genes, db_base, discovery)
    results.append(baseline_result)

    for value in SWEEP["k1"]:
        if value == BASELINE_K1:
            continue
        sweep_params["k1"], sweep_params["b"] = value, BASELINE_B
        results.append(run_config("k1", value, dict(BASELINE_OPTIONS), sweep_params,
                                   BASELINE_ANALYSIS_PAPERS, BASELINE_SEARCH_POOL,
                                   bm25_search, normalize_query, associated_entity_tables,
                                   genes_lexicon, tier1_genes, db_base, discovery))
    sweep_params["k1"] = BASELINE_K1

    for value in SWEEP["b"]:
        if value == BASELINE_B:
            continue
        sweep_params["k1"], sweep_params["b"] = BASELINE_K1, value
        results.append(run_config("b", value, dict(BASELINE_OPTIONS), sweep_params,
                                   BASELINE_ANALYSIS_PAPERS, BASELINE_SEARCH_POOL,
                                   bm25_search, normalize_query, associated_entity_tables,
                                   genes_lexicon, tier1_genes, db_base, discovery))
    sweep_params["b"] = BASELINE_B

    for value in SWEEP["title_match_multiplier"]:
        if value == BASELINE_OPTIONS["title_match_multiplier"]:
            continue
        opts = dict(BASELINE_OPTIONS, title_match_multiplier=value)
        results.append(run_config("title_match_multiplier", value, opts, sweep_params,
                                   BASELINE_ANALYSIS_PAPERS, BASELINE_SEARCH_POOL,
                                   bm25_search, normalize_query, associated_entity_tables,
                                   genes_lexicon, tier1_genes, db_base, discovery))

    for value in SWEEP["mesh_match_multiplier"]:
        if value == BASELINE_OPTIONS["mesh_match_multiplier"]:
            continue
        opts = dict(BASELINE_OPTIONS, mesh_match_multiplier=value)
        results.append(run_config("mesh_match_multiplier", value, opts, sweep_params,
                                   BASELINE_ANALYSIS_PAPERS, BASELINE_SEARCH_POOL,
                                   bm25_search, normalize_query, associated_entity_tables,
                                   genes_lexicon, tier1_genes, db_base, discovery))

    for value in SWEEP["review_multiplier"]:
        if value == 0.65:
            continue
        opts = dict(BASELINE_OPTIONS, review_multiplier=value)
        results.append(run_config("review_multiplier", value, opts, sweep_params,
                                   BASELINE_ANALYSIS_PAPERS, BASELINE_SEARCH_POOL,
                                   bm25_search, normalize_query, associated_entity_tables,
                                   genes_lexicon, tier1_genes, db_base, discovery))

    for value in SWEEP["analysis_papers"]:
        if value == BASELINE_ANALYSIS_PAPERS:
            continue
        pool = int(value * 1.25)
        results.append(run_config("analysis_papers", value, dict(BASELINE_OPTIONS), sweep_params,
                                   value, pool,
                                   bm25_search, normalize_query, associated_entity_tables,
                                   genes_lexicon, tier1_genes, db_base, discovery))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["param", "value", "hits", "total", "precision_at_10"])
        w.writeheader()
        w.writerows(results)
    print(f"\nSaved {OUT} ({len(results)} configs) in {time.time()-t_start:.0f}s total", flush=True)


if __name__ == "__main__":
    main()
