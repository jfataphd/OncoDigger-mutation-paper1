"""
uncleaned_vocabulary_sensitivity.py -- does the HGNC vocabulary QC layer (alias-collision
removal, common-word filtering, longest-match nested-name masking, the 7 exclusions + 4
alias substitutions) materially change the two resulting top-10 gene lists, vetted vs. raw?

An "entirely uncleaned" sensitivity test already exists for the frequency-only baseline
control (20.5% uncleaned vs 52.6% vetted; Methods 2.2/2.3), but no equivalent test exists
for the main BM25 + enrichment pipeline that actually produces the paper's headline number.
This script builds one, reusing the same frozen engine, corpora, query ("mutation,
mutations"), and retrieval settings (k1=1.5, b=0.75, title x1.2, MeSH x1.15, Primary Evidence
mode) as the primary analysis -- the ONLY thing that changes is the gene vocabulary itself:

  Vetted (current manuscript pipeline): preferred_alias substitution (GC->VDBP etc.),
    GENE_EXCLUDE (7 symbols dropped), GENE_REMAP (P53->TP53, HER2->ERBB2), the generic
    length/common-word alias filter, and longest-match nested-gene-name masking.

  Raw/uncleaned: every gene indexed under its own literal HGNC approved symbol, with
    every alias from the HGNC alias column, no filtering, no exclusions, no
    substitutions, and no longest-match correction -- i.e. the same policy the public
    oncodigger.com site runs (per run_query_planes.py's own docstring note).

Retired 2026-10-08: this script originally also scored both top-10 lists against COSMIC
Cancer Gene Census Tier 1 (`hits_vetted`/`hits_raw` columns, `cosmic_census_2026-05-23.csv`,
not included in this repository). That COSMIC-based score is not used by manuscript1.1's
IntOGen-calibrated analysis (intogen_ablation_panel.py rescores the gene lists against
IntOGen instead) and has been removed; only the two gene lists (top10_vetted, top10_raw)
are still used.

Run:  python -X utf8 scripts/analysis/uncleaned_vocabulary_sensitivity.py
Output: data/derived/uncleaned_vocabulary_precision.csv (cancer, top10_vetted, top10_raw)
"""
from __future__ import annotations

import csv
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ONCODIGGER = Path(r"C:\Users\van0g\Documents\Codex\2026-05-03\can-you-get-access-to-my\OncoDigger")
OUT = ROOT / "data" / "derived" / "uncleaned_vocabulary_precision.csv"

GENE_EXCLUDE = {
    "GC", "HR", "MICE", "DCR", "BPIFA4P", "PC", "HCCS", "WDHD1", "SGCG", "C2", "FBN1",
    # 2026-09-30 false-alias sweep -- kept in sync with run_query_planes.py's GENE_EXCLUDE (vetted side only).
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
OPTIONS = {
    "year_min": 2000, "year_max": 2026, "exclude_reviews": False, "downrank_reviews": True,
    "title_match_multiplier": 1.2, "mesh_match_multiplier": 1.15, "recency_boost": True,
    "discovery_mode": "Primary Evidence",
}
ANALYSIS_PAPERS, SEARCH_POOL = 1000, 1250


def build_raw_lexicon(root: Path, add_entity, build_lookup, norm_text) -> dict:
    """Every gene under its own literal HGNC symbol, every alias, no filtering at all."""
    import pandas as pd
    frame = pd.read_csv(root / "Human Genes.csv").fillna("")
    genes = {}
    for row in frame.to_dict("records"):
        symbol = str(row.get("symbol2") or row.get("symbol")).strip().upper()
        if not symbol:
            continue
        name = str(row.get("Approved name", ""))
        aliases = [a.strip() for a in str(row.get("Alias symbols", "")).split(",") if a.strip()]
        terms = sorted({norm_text(t) for t in [symbol, name, *aliases] if norm_text(t)})
        if symbol in genes:
            terms = sorted(set(terms) | set(genes[symbol]["terms"]))
        add_entity(genes, symbol, symbol, "Gene", terms, source_id=symbol, description=name)
    return build_lookup(genes)


def open_corpus(db_path, discovery):
    import sqlite3
    discovery._sqlite_doc_counts.clear()
    discovery._sqlite_caches.clear()
    return sqlite3.connect(str(db_path), check_same_thread=False)


def rank_top10(con, bm25_search, normalize_query, associated_entity_tables, genes_lexicon,
                apply_exclude_remap: bool) -> list[str]:
    tokens = normalize_query(QUERY)
    hits = bm25_search(con, tokens, SEARCH_POOL, OPTIONS)[:ANALYSIS_PAPERS]
    tables = associated_entity_tables(con, hits, QUERY, (2000, 2026), tokens, {"Genes": genes_lexicon},
                                      limit=100, discovery_mode="Primary Evidence", corpus_query="")
    genes_df = tables.get("Genes")
    ranked, seen = [], set()
    if genes_df is None or genes_df.empty:
        return ranked
    for _, row in genes_df.iterrows():
        gene = row["Gene"]
        if apply_exclude_remap:
            gene = GENE_REMAP.get(gene, gene)
            if gene in GENE_EXCLUDE:
                continue
        if gene in seen:
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


def main() -> None:
    if os.environ.get("PYTHONHASHSEED") != "0":
        import subprocess
        env = dict(os.environ, PYTHONHASHSEED="0")
        sys.exit(subprocess.run([sys.executable, "-X", "utf8"] + sys.argv, env=env).returncode)

    oncodigger_root = DEFAULT_ONCODIGGER
    sys.path.insert(0, str(oncodigger_root / "src"))
    from oncodigger.search import discovery
    from oncodigger.search.discovery import associated_entity_tables
    from oncodigger.search.lexical import bm25_search, normalize_query
    from oncodigger.search.lexicons import load_lexicons, add_entity, build_lookup, norm_text

    vetted_lexicons = load_lexicons(oncodigger_root)
    genes_vetted = vetted_lexicons["Genes"]
    genes_raw = build_raw_lexicon(oncodigger_root, add_entity, build_lookup, norm_text)
    print(f"Vetted vocabulary: {len(genes_vetted['single']) + len(genes_vetted['phrases'])} term buckets", flush=True)
    print(f"Raw vocabulary:    {len(genes_raw['single']) + len(genes_raw['phrases'])} term buckets", flush=True)

    db_base = oncodigger_root / "data" / "pubmed_processed"
    rows = []
    t_start = time.time()

    for key, cancer_label in CANCER_MAP:
        db = next((p for p in (db_base / key / "lexical_index.db", db_base / f"{key}_all" / "lexical_index.db")
                   if p.exists()), None)
        if db is None:
            sys.exit(f"missing index for {cancer_label}")

        con = open_corpus(db, discovery)
        top10_vetted = rank_top10(con, bm25_search, normalize_query, associated_entity_tables,
                                    genes_vetted, apply_exclude_remap=True)
        con.close()

        con = open_corpus(db, discovery)
        top10_raw = rank_top10(con, bm25_search, normalize_query, associated_entity_tables,
                                 genes_raw, apply_exclude_remap=False)
        con.close()

        rows.append({
            "cancer": cancer_label,
            "top10_vetted": " ".join(top10_vetted),
            "top10_raw": " ".join(top10_raw),
        })
        print(f"  {cancer_label:20} vetted {top10_vetted}", flush=True)
        if top10_vetted != top10_raw:
            print(f"    raw:    {top10_raw}", flush=True)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["cancer", "top10_vetted", "top10_raw"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nSaved {OUT} in {time.time()-t_start:.0f}s total", flush=True)
    print("Scored against IntOGen by scripts/analysis/intogen_ablation_panel.py")


if __name__ == "__main__":
    main()
