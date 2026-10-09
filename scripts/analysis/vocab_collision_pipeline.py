"""
vocab_collision_pipeline.py -- diagnostic only, not yet wired into the scoring pipeline.
Single, reusable, uniformly-applied vocabulary-collision check for any flagged gene symbol.
Replaces the two earlier batch scripts (vocab_collision_check_batch1.py,
vocab_collision_check_batch2.py), whose checks are reproduced here via one shared function
so the method is identical for every gene, not bespoke per-gene code.

For each gene, run() does, in this fixed order, every time:
  1. Look up the gene's real cancer(s)/rank(s) directly from
     data/canonical/precision_at_k_detail.csv. Never assumed or hardcoded -- this is the
     step that was skipped once for REST (checked against Kidney Cancer, when REST's only
     actual ranking is Thyroid Cancer rank 50) and must not recur.
  2. Pull the real postings count for the gene's best-ranked cancer: prefer
     data/derived/intogen_rescoring/intogen_literature_weight_audit.csv's supporting_papers
     column (the in-pool, score-relevant count) if that gene/cancer pair exists there, else
     fall back to full-corpus SQLite postings count, clearly labelled which source was used.
  3. Capitalization filter on the RAW, case-preserved title+abstract text: a match counts as
     "not all-lowercase" if ANY occurrence of the whole word (case-insensitive) is not purely
     lowercase. This accepts both human nomenclature (ALL-CAPS: FAS) and mouse nomenclature
     (Title case: Fas) -- the correction found necessary after the strict all-caps-only
     filter undercounted FAS (22.1% -> 98.5% once corrected) and KIT (79.0% -> 84.6%).
  4. Negative-phrase filter: of the capitalization-filtered set, discard any paper containing
     one of that gene's pre-derived alternate-meaning phrases (NEGATIVE_PHRASES below). The
     phrase LIST is gene-specific; the lookup/application CODE is identical for every gene.
  5. Report: gene, cancer(s)/rank(s), postings count + source, capitalized rate, final
     genuine count/rate.

Run: python -X utf8 scripts/analysis/vocab_collision_pipeline.py GENE1 GENE2 ...
     (no args = runs the full validated list from both prior batches, for reproducibility)
"""
from __future__ import annotations

import json
import re
import sqlite3
import sys
import os
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CANONICAL = ROOT / "data" / "canonical" / "precision_at_k_detail.csv"
AUDIT_CSV = ROOT / "data" / "derived" / "intogen_rescoring" / "intogen_literature_weight_audit.csv"

BASE_DB = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET")) / "data" / "pubmed_processed"
BASE_RAW = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET")) / "data" / "pubmed_raw"

# Verified against the actual pubmed_processed directory listing (19 folders, one per
# study cancer). "Pancreatic Cancer" maps to "pancreas_cancer_inclusive", not
# "pancreatic_cancer_inclusive" -- the one irregular case.
CANCER_KEY = {
    "Breast Cancer": "breast_cancer_inclusive", "Lung Cancer": "lung_cancer_inclusive",
    "Colorectal Cancer": "colorectal_cancer_inclusive", "Prostate Cancer": "prostate_cancer_inclusive",
    "Melanoma": "melanoma_cancer_inclusive", "Bladder Cancer": "bladder_cancer_inclusive",
    "Kidney Cancer": "kidney_cancer_inclusive", "Pancreatic Cancer": "pancreas_cancer_inclusive",
    "Liver Cancer": "liver_cancer_inclusive", "Stomach Cancer": "stomach_cancer_inclusive",
    "Esophageal Cancer": "esophageal_cancer_inclusive", "Ovarian Cancer": "ovarian_cancer_inclusive",
    "Endometrial Cancer": "endometrial_cancer_inclusive", "Cervical Cancer": "cervical_cancer_inclusive",
    "Thyroid Cancer": "thyroid_cancer_inclusive", "Brain Cancer": "brain_cancer_inclusive",
    "Leukemia": "leukemia_cancer_inclusive", "Lymphoma": "lymphoma_cancer_inclusive",
    "Myeloma": "myeloma_cancer_inclusive",
}

# One central, documented data structure. Each gene's negative-phrase list is derived from
# actually reading real matched-paper text (see the two batch scripts' docstrings/samples
# for the provenance of each entry), not guessed. Hard-negative = confirmed unrelated
# meaning. Genes with an empty list were checked and found to have no confirmed false
# meaning worth filtering (capitalization alone was sufficient).
NEGATIVE_PHRASES: dict[str, list[str]] = {
    "MET": ["mesenchymal-to-epithelial transition", "mesenchymal to epithelial transition",
            "mesenchymal-epithelial transition", "emt/met", "emt / met", "methioninase",
            "pseudomonas putida", "met-pet", "methionine pet", "11c-methionine",
            "c-11 methionine", "mucoepidermoid", "metronomic",
            # 2026-10-05 additions, found via manual adjudication of a 40-paper random sample
            # from the "capitalized, no recognized positive or negative signal" bucket
            # (Methods 3.3): confirmed false-positive patterns the original 11-phrase list
            # did not cover. "met+"/"met-"/"nmet" were tried and DROPPED after verification
            # showed they also match genuine usage at high volume ("MET-amplified",
            # "MET-positive", "MET-overexpressing", and the unrelated English word "unmet"),
            # so a single abbreviated-metastasis-status false-positive paper (PMID 38433721)
            # is left as a disclosed, unresolved limitation rather than patched with an
            # unsafe blanket rule.
            "metformin", "met/met", "met/val", "val/met", "met-5a", "met 5a"],
    "FAS": ["fetal alcohol syndrome", "functional assessment", "full analysis set"],
    # 2026-10-05 addition, not part of the original 23-gene IntOGen-driver audit (FH is not
    # an IntOGen-confirmed driver in Kidney Cancer, Results 5.3) but checked following an
    # external review prompt. Manual adjudication of 40 of FH's 297 Kidney Cancer postings
    # found 36/40 (90%) genuine fumarate hydratase gene references and 4/40 (10%) false,
    # all "favorable histology (FH)", a Wilms tumor clinical-staging abbreviation, not
    # "family history" as initially hypothesized.
    "FH": ["favorable histology", "unfavorable histology"],
    "KIT": ["assay kit", "test kit", "elisa kit", "detection kit", "extraction kit",
            "diagnostic kit", "starter kit", "commercial kit"],
    "MAF": [],
    "REST": ["nephrogenic rest", "rest 2008", "representational state transfer",
             "rest api", "restful"],
    "CIC": [],
    "MPL": [],
    "FAT1": [],
    "EGFR": ["estimated glomerular filtration rate", "glomerular filtration rate",
             "creatinine clearance", "renal function", "chronic kidney disease", "ckd stage"],
    "ESR1": ["endoplasmic reticulum", "unfolded protein response", "er stress"],
    "POLE": [],
    "RET": ["reverse electron transport", "resistance exercise training"],
    "APC": ["antigen-presenting cell", "antigen presenting cell", "activated protein c",
            "atrial premature contraction"],
    "FLT3": [],
    "GNAS": [],
    "MEN1": [],
    "NF1": [],
    "PMS2": ["pms2cl"],
    "RB1": ["right bundle branch block", "rbbb"],
    "TSC1": ["lymphangioleiomyomatosis"],
    "TSC2": [],
    "VHL": [],
}


def lookup_rankings(gene: str) -> list[tuple[str, int]]:
    """Step 1: real cancer(s)/rank(s) from the canonical top-50 rankings. Never assumed."""
    r = pd.read_csv(CANONICAL)
    sub = r[r["gene"].str.upper() == gene.upper()].sort_values("rank")
    return list(zip(sub["cancer"], sub["rank"]))


def lookup_postings(gene: str, cancer: str) -> tuple[int, str]:
    """Step 2: prefer in-pool supporting_papers from the audit CSV, else full-corpus count."""
    audit = pd.read_csv(AUDIT_CSV)
    row = audit[(audit["gene"] == gene) & (audit["cancer"] == cancer)]
    if len(row) and pd.notna(row.iloc[0]["supporting_papers"]):
        return int(row.iloc[0]["supporting_papers"]), "in-pool supporting_papers (audit CSV)"
    corpus_key = CANCER_KEY[cancer]
    db = BASE_DB / f"{corpus_key}_all" / "lexical_index.db"
    con = sqlite3.connect(str(db))
    n = con.execute("SELECT COUNT(DISTINCT pmid) FROM postings WHERE token = ?",
                     (gene.lower(),)).fetchone()[0]
    return n, "full-corpus postings (SQLite, audit CSV had no row for this pair)"


def capitalization_and_negative_filter(gene: str, cancer: str) -> dict:
    """Steps 3-4: corrected capitalization filter, then gene-specific negative-phrase filter,
    both applied against the real raw case-preserved text for the gene's full-corpus matches."""
    corpus_key = CANCER_KEY[cancer]
    db = BASE_DB / f"{corpus_key}_all" / "lexical_index.db"
    con = sqlite3.connect(str(db))
    token = gene.lower()
    pmids = {str(p[0]) for p in con.execute(
        "SELECT DISTINCT pmid FROM postings WHERE token = ?", (token,)).fetchall()}
    raw = BASE_RAW / corpus_key / "abstracts.jsonl"
    pattern = re.compile(r"\b" + re.escape(token) + r"\b", re.IGNORECASE)
    neg_phrases = NEGATIVE_PHRASES.get(gene.upper(), [])

    # A match preceded only by whitespace/start-of-string or by sentence-ending punctuation
    # ('. ', '? ', '! ') is merely capitalized because it starts a sentence -- indistinguishable
    # from ordinary-word use, NOT a genuine signal of gene nomenclature (Title case or ALL-CAPS
    # mid-sentence). Found via REST: 5 of 6 "capitalized" matches were just "Rest of patients...",
    # sentence-initial, not the gene. Excluding these is required for any gene whose lowercase
    # form is also a plain English word.
    sentence_start = re.compile(r"(^|[.!?]\s+)\s*$")

    total = not_lower = genuine = 0
    with open(raw, encoding="utf-8") as f:
        for line in f:
            rec = json.loads(line)
            if str(rec.get("pmid", "")) not in pmids:
                continue
            total += 1
            text = (rec.get("title") or "") + " " + (rec.get("abstract") or "")
            is_caps = any(
                m.group(0) != token and not sentence_start.search(text[:m.start()])
                for m in pattern.finditer(text)
            )
            if not is_caps:
                continue
            not_lower += 1
            low = text.lower()
            if neg_phrases and any(p in low for p in neg_phrases):
                continue
            genuine += 1
    return {"total": total, "not_lower": not_lower, "genuine": genuine}


def run(gene: str) -> None:
    rankings = lookup_rankings(gene)
    if not rankings:
        print(f"{gene}: not in top-50 anywhere, skipping (no corpus/cancer to check against).")
        return
    best_cancer, best_rank = rankings[0]
    postings, source = lookup_postings(gene, best_cancer)
    stats = capitalization_and_negative_filter(gene, best_cancer)

    print(f"=== {gene} ===")
    print(f"  cancer(s)/rank(s): {rankings}")
    print(f"  postings for {best_cancer}: {postings} [{source}]")
    print(f"  full-corpus matches checked: {stats['total']}")
    if stats["total"] == 0:
        print("  (no full-corpus matches found -- postings count may be from a different corpus slice)")
        print()
        return
    cg = 100 * stats["not_lower"] / stats["total"]
    gr = 100 * stats["genuine"] / stats["total"]
    print(f"  not-all-lowercase (C_g): {stats['not_lower']} ({cg:.1f}%)")
    print(f"  genuine after negative-phrase filter: {stats['genuine']} ({gr:.1f}% of total, "
          f"{100 * stats['genuine'] / stats['not_lower']:.1f}% of capitalized)" if stats['not_lower'] else "  (none capitalized)")
    print()


VALIDATED_RERUN = ["FAS", "KIT", "MAF", "CIC", "MPL", "FAT1", "EGFR", "ESR1", "POLE", "RET", "REST"]
NEW_GENES = ["APC", "FLT3", "GNAS", "MEN1", "NF1", "PMS2", "RB1", "TSC1", "TSC2", "VHL"]


def main() -> None:
    genes = sys.argv[1:] or (VALIDATED_RERUN + NEW_GENES)
    for g in genes:
        run(g)


if __name__ == "__main__":
    main()
