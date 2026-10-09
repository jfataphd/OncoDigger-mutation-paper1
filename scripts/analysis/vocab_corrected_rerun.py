"""
vocab_corrected_rerun.py -- applies the frozen, validated MET/REST/FH false-positive rule
(from vocab_collision_pipeline.py's NEGATIVE_PHRASES + capitalization logic, unchanged) to
physically correct Lung Cancer's, Thyroid Cancer's and Kidney Cancer's corpora, then re-runs
the REAL scoring engine (not an approximated relevance-score recompute) to get the true
corrected top-50 rankings for those three cancers.

Corrected 2026-10-07 after a reproducibility audit found this script's original version did
not do what its own docstring claimed. The original version deleted the false-positive PMIDs'
rows from the postings table only. The engine's gene-mention detection (associated_entity_tables
in discovery.py) derives mentions by re-tokenizing documents.title/documents.abstract directly
and never queries postings, so deleting postings rows left every false-positive mention fully
intact in scoring. Verified empirically: MET's "Supporting papers" count was byte-identical
(65) before and after the old postings-only "correction," and a specific confirmed
false-positive PMID remained in MET's supporting-paper set both before and after. Deleting
postings did have one real, unintended side effect: it shrank the enrichment denominator
(background rate), which inflated relevance scores for the "corrected" genes for a reason
unrelated to removing false mentions, and for REST specifically, whose false-positive rate for
this cancer happened to be 100% (136/136), it zeroed out its postings count entirely, which
triggered an unrelated downstream gate in rank_genes() that drops any gene with zero remaining
postings regardless of actual text-based support. That gate, not a genuine demonstration of
false-mention removal, is what caused REST's prior disappearance from the ranking.

This version instead redacts each false-positive PMID's gene-token occurrences directly in the
documents table's title and abstract text (case-insensitive, word-boundary, replaced with a
placeholder that cannot re-match any lexicon term), in addition to removing its postings row,
so the engine's real text-based entity detection genuinely stops finding that mention. Re-run
against the real engine, the same three qualitative conclusions reported in the manuscript
still held under this corrected methodology: MET retained its Lung Cancer top-10 position
(Supporting papers 65 -> 56, rank unchanged, still 10th), FH retained its Kidney Cancer rank-5
position (Supporting papers 88 -> 86, rank unchanged), and REST was genuinely removed from the
Thyroid Cancer ranking (0 remaining text-based mentions, confirmed not a gate artifact),
replaced at rank 50 by MAPK3, exactly as previously reported. The numbers changed; the
conclusions did not.

Kidney Cancer/FH was previously missing from this script's main() entirely; the corrected
database used to generate the manuscript's reported FH result was built by a one-off manual
command, not by this script, so a fresh checkout could not reproduce it. Added here.

Does not touch the original corpus databases under the Codex OncoDigger pubmed_processed
directory -- works on copies only, written to
data/derived/intogen_rescoring/corrected_corpora/.

Run: python -X utf8 scripts/analysis/vocab_corrected_rerun.py
"""
from __future__ import annotations

import json
import re
import shutil
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from vocab_collision_pipeline import NEGATIVE_PHRASES, CANCER_KEY, BASE_DB, BASE_RAW  # noqa: E402

OUT_DIR = ROOT / "data" / "derived" / "intogen_rescoring" / "corrected_corpora"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SENTENCE_START = re.compile(r"(^|[.!?]\s+)\s*$")
REDACTION_PLACEHOLDER = "ZZREDACTEDZZ"


def false_positive_pmids(gene: str, cancer: str) -> set[str]:
    """Exact same logic as vocab_collision_pipeline.capitalization_and_negative_filter,
    but returns the set of PMIDs that FAIL (false positives) instead of just a count."""
    corpus_key = CANCER_KEY[cancer]
    db = BASE_DB / f"{corpus_key}_all" / "lexical_index.db"
    con = sqlite3.connect(str(db))
    token = gene.lower()
    pmids = {str(p[0]) for p in con.execute(
        "SELECT DISTINCT pmid FROM postings WHERE token = ?", (token,)).fetchall()}
    raw = BASE_RAW / corpus_key / "abstracts.jsonl"
    word_boundary = chr(92) + "b"
    pattern = re.compile(word_boundary + re.escape(token) + word_boundary, re.IGNORECASE)
    neg_phrases = NEGATIVE_PHRASES.get(gene.upper(), [])

    false_pos = set()
    found = set()
    with open(raw, encoding="utf-8") as f:
        for line in f:
            rec = json.loads(line)
            pmid = str(rec.get("pmid", ""))
            if pmid not in pmids:
                continue
            found.add(pmid)
            text = (rec.get("title") or "") + " " + (rec.get("abstract") or "")
            is_caps = any(
                m.group(0) != token and not SENTENCE_START.search(text[:m.start()])
                for m in pattern.finditer(text)
            )
            if not is_caps:
                false_pos.add(pmid)
                continue
            low = text.lower()
            if neg_phrases and any(p in low for p in neg_phrases):
                false_pos.add(pmid)
    missing = pmids - found
    print(f"  [{gene}/{cancer}] postings pmids={len(pmids)} found_in_raw={len(found)} "
          f"missing_from_raw={len(missing)} false_positive={len(false_pos)} "
          f"genuine_kept={len(found) - len(false_pos)}")
    return false_pos


def make_corrected_db(cancer: str, gene: str, fp_pmids: set[str]) -> Path:
    """Build a corrected corpus copy that genuinely removes the false-positive mentions.

    Redacts the gene token from the false-positive PMIDs' title/abstract text in the copy's
    documents table (so the engine's real, text-based entity detection stops finding them),
    and removes their postings rows (so background/enrichment counts reflect the correction
    too). Genuine (non-false-positive) PMIDs are left completely untouched.
    """
    corpus_key = CANCER_KEY[cancer]
    src = BASE_DB / f"{corpus_key}_all" / "lexical_index.db"
    dst_dir = OUT_DIR / f"{corpus_key}_all"
    dst_dir.mkdir(parents=True, exist_ok=True)
    dst = dst_dir / "lexical_index.db"
    shutil.copy2(src, dst)
    con = sqlite3.connect(str(dst))
    token = gene.lower()
    word_boundary = chr(92) + "b"
    pattern = re.compile(word_boundary + re.escape(token) + word_boundary, re.IGNORECASE)

    before = con.execute("SELECT COUNT(*) FROM postings WHERE token = ?", (token,)).fetchone()[0]

    redacted = 0
    for pmid in fp_pmids:
        row = con.execute("SELECT title, abstract FROM documents WHERE pmid = ?", (pmid,)).fetchone()
        if row is None:
            continue
        title, abstract = row
        new_title = pattern.sub(REDACTION_PLACEHOLDER, title or "")
        new_abstract = pattern.sub(REDACTION_PLACEHOLDER, abstract or "")
        con.execute("UPDATE documents SET title = ?, abstract = ? WHERE pmid = ?",
                    (new_title, new_abstract, pmid))
        redacted += 1

    fp_list = list(fp_pmids)
    placeholders = ",".join(["?"] * len(fp_list))
    q = "DELETE FROM postings WHERE token = ? AND pmid IN (" + placeholders + ")"
    con.execute(q, [token, *fp_list])
    con.commit()
    after = con.execute("SELECT COUNT(*) FROM postings WHERE token = ?", (token,)).fetchone()[0]
    con.close()
    print(f"  [{cancer}] redacted text for {redacted}/{len(fp_pmids)} false-positive PMIDs; "
          f"postings for '{token}': {before} -> {after} (removed {before - after})")
    return dst


def main():
    print("Step 1: identify false-positive PMIDs using the frozen validated rule")
    met_fp = false_positive_pmids("MET", "Lung Cancer")
    rest_fp = false_positive_pmids("REST", "Thyroid Cancer")
    fh_fp = false_positive_pmids("FH", "Kidney Cancer")

    print()
    print("Step 2: build corrected corpus copies (originals untouched), text genuinely redacted")
    make_corrected_db("Lung Cancer", "MET", met_fp)
    make_corrected_db("Thyroid Cancer", "REST", rest_fp)
    make_corrected_db("Kidney Cancer", "FH", fh_fp)

    print()
    print("Done. Original corpora are untouched.")
    print("Corrected copies at", OUT_DIR)


if __name__ == "__main__":
    main()
