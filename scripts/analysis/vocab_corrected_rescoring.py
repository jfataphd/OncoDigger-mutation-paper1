"""
vocab_corrected_rescoring.py -- re-runs the REAL scoring engine (same setup as
run_query_planes.py: load_lexicons, suppress_ambiguous_symbols, longest_match_lookups,
identical OPTIONS/SEARCH_POOL/ANALYSIS_PAPERS/TOP_K) for Kidney Cancer, Lung Cancer and
Thyroid Cancer specifically, substituting the corrected corpus copies built by
vocab_corrected_rerun.py (MET/REST/FH false positives physically removed from the postings
table) in place of the originals. Every other cancer, the vocabulary, and the engine code
itself are identical to the original run -- only these three corpora's postings for
'met'/'rest'/'fh' differ.

Not an approximation: this calls the same bm25_search / associated_entity_tables /
load_lexicons functions run_query_planes.py uses, against real (corrected) SQLite data,
not a hand-recomputed relevance score.

Run: python -X utf8 scripts/analysis/vocab_corrected_rescoring.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

DEFAULT_ONCODIGGER = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET"))
CORRECTED_DIR = ROOT / "data" / "derived" / "intogen_rescoring" / "corrected_corpora"

# (corpus_key, label, corrected_db_path)
TARGETS = [
    ("kidney_cancer_inclusive", "Kidney Cancer", CORRECTED_DIR / "kidney_cancer_inclusive_all" / "lexical_index.db"),
    ("lung_cancer_inclusive", "Lung Cancer", CORRECTED_DIR / "lung_cancer_inclusive_all" / "lexical_index.db"),
    ("thyroid_cancer_inclusive", "Thyroid Cancer", CORRECTED_DIR / "thyroid_cancer_inclusive_all" / "lexical_index.db"),
]


def main() -> None:
    if os.environ.get("PYTHONHASHSEED") != "0":
        import subprocess
        env = dict(os.environ, PYTHONHASHSEED="0")
        sys.exit(subprocess.run([sys.executable, "-X", "utf8", __file__], env=env).returncode)

    from run_query_planes import (  # noqa: E402
        suppress_ambiguous_symbols, longest_match_lookups, open_corpus, rank_genes, PLANES,
    )

    sys.path.insert(0, str(DEFAULT_ONCODIGGER / "src"))
    from oncodigger.search import discovery
    from oncodigger.search.discovery import associated_entity_tables
    from oncodigger.search.lexical import bm25_search, normalize_query
    from oncodigger.search.lexicons import load_lexicons

    for db_path in (t[2] for t in TARGETS):
        assert db_path.exists(), f"missing corrected db: {db_path}"

    lexicons = load_lexicons(DEFAULT_ONCODIGGER)
    genes_suppressed = suppress_ambiguous_symbols(lexicons["Genes"])
    genes = longest_match_lookups(genes_suppressed, discovery)
    engine = (bm25_search, normalize_query, associated_entity_tables, {"Genes": genes})

    results = {}
    for corpus_key, label, db_path in TARGETS:
        con = open_corpus(db_path, discovery)
        ranked = rank_genes(con, PLANES["mutation"], engine)
        con.close()
        results[label] = ranked
        print(f"{label}: top-10 = {ranked[:10]}")
        print(f"{label}: full top-{len(ranked)} = {ranked}")

    out = ROOT / "data" / "derived" / "intogen_rescoring" / "vocab_corrected_rankings.csv"
    import csv
    with out.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["cancer", "rank", "gene"])
        for label, ranked in results.items():
            for i, g in enumerate(ranked, 1):
                w.writerow([label, i, g])
    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
