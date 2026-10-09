"""
vocab_pool_contamination_trace.py -- Stage 3 of the vocabulary-collision audit. Full-
corpus contamination (Stage 1/2) may or may not matter for a gene's actual rank: BM25
retrieval for the "mutation, mutations" query already selects the 1,000-paper analysis
pool, which could incidentally concentrate genuine gene literature (since real gene papers
are more likely to also be mutation-focused) or could pull in exactly the same noise. This
reproduces the REAL retrieval pool (same bm25_search call, same options, as
run_query_planes.py) for each flagged gene's cancer, and reports: how many of its
supporting papers are actually inside that pool (not just the full corpus), and of those,
how many the frozen Stage-2 rule classifies as genuine.

Run: python -X utf8 scripts/analysis/vocab_pool_contamination_trace.py
"""
from __future__ import annotations

import json
import os
import re
import sqlite3
import sys
from pathlib import Path

if os.environ.get("PYTHONHASHSEED") != "0":
    import subprocess
    env = dict(os.environ, PYTHONHASHSEED="0")
    sys.exit(subprocess.run([sys.executable, "-X", "utf8"] + sys.argv, env=env).returncode)

ONCODIGGER_ROOT = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET"))
PUBMED_RAW = ONCODIGGER_ROOT / "data" / "pubmed_raw"
PUBMED_PROCESSED = ONCODIGGER_ROOT / "data" / "pubmed_processed"
sys.path.insert(0, str(ONCODIGGER_ROOT / "src"))

QUERY = "mutation, mutations"
SEARCH_POOL, ANALYSIS_PAPERS = 1250, 1000
OPTIONS = {
    "year_min": 2000, "year_max": 2026, "exclude_reviews": False, "downrank_reviews": True,
    "title_match_multiplier": 1.2, "mesh_match_multiplier": 1.15, "recency_boost": True,
    "discovery_mode": "Primary Evidence",
}

CODON_PATTERN = re.compile(r"(arg|lys|leu|val|ile|thr)\W{0,15}met\b", re.IGNORECASE)
GENE_RULES = {
    "MET": {
        "corpus": "lung_cancer_inclusive", "cancer": "Lung Cancer",
        "positive": ["c-met", "proto-oncogene", "protooncogene", "receptor tyrosine kinase",
                     "hepatocyte growth factor", "hgf", "amplification", "exon 14", "crizotinib",
                     "capmatinib", "tepotinib", "savolitinib", "kinase inhibitor", "juxtamembrane",
                     "met receptor", "met gene", "met oncogene", "met mutation", "met inhibitor",
                     "met amplif", "met kinase", "met signaling", "met signalling", "met pathway",
                     "met alterations", "copy number"],
        "negative": ["mesenchymal-to-epithelial transition", "mesenchymal to epithelial transition",
                     "mesenchymal-epithelial transition", "methioninase", "pseudomonas putida",
                     "met-pet", "methionine pet", "11c-methionine", "c-11 methionine", "mucoepidermoid",
                     "metronomic"],
        "codon": True,
    },
    "KIT": {
        "corpus": "melanoma_cancer_inclusive", "cancer": "Melanoma",
        "positive": ["c-kit", "kit mutation", "kit aberration", "kit alteration", "kit expression",
                     "kit inhibit", "kit locus", "kit gene", "kit protein", "kit amplif", "cd117",
                     "kit-positive", "kit positive", "imatinib", "sunitinib", "dasatinib", "nilotinib",
                     "stem cell factor"],
        "negative": ["elisa kit", "assay kit", "detection kit", "extraction kit", "diagnostic kit",
                     "test kit", "commercial kit", "starter kit", "kit according to the manufacturer",
                     "kit was used", "using a kit", "isolation kit", "purification kit", "sequencing kit",
                     "library kit", "cloning kit"],
        "codon": False,
    },
}


def classify(text: str, symbol: str, rule: dict) -> str:
    caps = re.compile(r"\b" + re.escape(symbol) + r"\b")
    if not caps.search(text):
        return "false_no_caps"
    low = text.lower()
    has_pos = any(t in low for t in rule["positive"])
    has_neg = any(t in low for t in rule["negative"])
    if rule["codon"] and CODON_PATTERN.search(text):
        has_neg = True
    if has_pos and not has_neg:
        return "genuine"
    if has_neg and not has_pos:
        return "false_negative_term"
    if has_pos and has_neg:
        return "ambiguous"
    return "unresolved"


def main() -> None:
    from oncodigger.search.lexical import bm25_search, normalize_query

    for symbol, rule in GENE_RULES.items():
        print(f"\n{'=' * 70}\n{symbol} ({rule['cancer']})\n{'=' * 70}")
        db_path = PUBMED_PROCESSED / f"{rule['corpus']}_all" / "lexical_index.db"
        con = sqlite3.connect(str(db_path), check_same_thread=False)

        tokens = normalize_query(QUERY)
        hits = bm25_search(con, tokens, SEARCH_POOL, OPTIONS)[:ANALYSIS_PAPERS]
        pool_pmids = {str(h["pmid"]) for h in hits}
        print(f"actual BM25 analysis-pool size: {len(pool_pmids)}")

        gene_pmids = {str(p[0]) for p in con.execute(
            "SELECT DISTINCT pmid FROM postings WHERE token = ?", (symbol.lower(),)).fetchall()}
        in_pool = gene_pmids & pool_pmids
        print(f"{symbol}'s full-corpus postings: {len(gene_pmids)}")
        print(f"{symbol}'s postings actually inside the real analysis pool: {len(in_pool)} "
              f"({100 * len(in_pool) / len(gene_pmids):.1f}% of full-corpus postings)")

        raw_path = PUBMED_RAW / rule["corpus"] / "abstracts.jsonl"
        raw = {}
        with open(raw_path, encoding="utf-8") as f:
            for line in f:
                rec = json.loads(line)
                raw[str(rec.get("pmid", ""))] = (rec.get("title") or "") + " " + (rec.get("abstract") or "")

        counts = {}
        for pmid in in_pool:
            c = classify(raw.get(pmid, ""), symbol, rule)
            counts[c] = counts.get(c, 0) + 1
        total = len(in_pool)
        print(f"\nClassification of the {total} in-pool postings:")
        for k, v in sorted(counts.items(), key=lambda kv: -kv[1]):
            print(f"  {k}: {v} ({100 * v / total:.1f}%)" if total else f"  {k}: {v}")

        genuine_in_pool = counts.get("genuine", 0)
        if total:
            reduction = 100 * (1 - genuine_in_pool / total)
            print(f"\n>>> {symbol}: supporting-paper count would change from {total} "
                  f"(current, uncorrected, everything in the pool) to {genuine_in_pool} "
                  f"(corrected, genuine-only) if this rule were applied -- a "
                  f"{reduction:.1f}% reduction in supporting evidence.")
        else:
            print(f"\n>>> {symbol}: no postings found in the analysis pool.")


if __name__ == "__main__":
    main()
