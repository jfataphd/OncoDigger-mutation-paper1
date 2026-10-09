"""
intogen_full_corpus_presence_check.py -- manuscript1.1 (IntOGen-primary draft) ONLY. Direct
test of whether the 986 "never-surfaced" IntOGen driver pairs (Results 5.5) are truly absent
from their cancer's corpus, or merely outside the 1,000-paper / 3,000-entity analysis window
the literature-weight audit used (intogen_literature_weight_audit.py: ANALYSIS_PAPERS=1000,
ENTITY_LIMIT=3000). This bypasses BM25 retrieval and entity-table extraction entirely and
queries each gene's raw postings count directly against the FULL corpus index (every paper in
that cancer's collection, not just the BM25-selected analysis pool), using the gene's single-
token symbol.

This is a presence/absence check only (does the term appear in ANY paper in the full corpus),
not a re-ranking: it does not reproduce OncoDigger's scoring, title/MeSH weighting or
longest-match correction, so a nonzero count here is necessary but not sufficient for the gene
to have ranked highly under the actual pipeline.

One false-alias collision was found by this check itself and is excluded here, mirroring the
existing vocabulary QC sweep (Methods 3.3): WAS (Brain Cancer) returned 68,416 postings, two
orders of magnitude above every other result, because "was" is a common English auxiliary verb,
not a genuine mention of the WAS (Wiskott-Aldrich syndrome) gene. No other result in this check
approaches that magnitude; MAX (Breast Cancer 670, Brain Cancer 427) is retained since it is
plausibly genuine (MAX is a real, if less prominent, MYC-dimerisation-partner proto-oncogene)
but is flagged as not independently vetted against live PubMed the way the primary top-50
rankings were.

Input: data/derived/intogen_rescoring/intogen_recall_missed_genes.csv (986 pairs)
Output: data/derived/intogen_rescoring/intogen_full_corpus_presence_check.csv (986 rows, includes WAS with a flag)
        data/derived/intogen_rescoring/intogen_full_corpus_presence_summary.csv

Run: python -X utf8 scripts/analysis/intogen_full_corpus_presence_check.py
"""
from __future__ import annotations

import sys
import os
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ONCODIGGER = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET"))
OUT = ROOT / "data" / "derived" / "intogen_rescoring" / "intogen_full_corpus_presence_check.csv"

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


def main():
    import sqlite3
    missed = pd.read_csv(ROOT / "data" / "derived" / "intogen_rescoring" / "intogen_recall_missed_genes.csv")
    db_base = DEFAULT_ONCODIGGER / "data" / "pubmed_processed"

    rows = []
    for cancer, group in missed.groupby("cancer"):
        key = CANCER_KEY[cancer]
        db = next((p for p in (db_base / key / "lexical_index.db", db_base / f"{key}_all" / "lexical_index.db")
                   if p.exists()), None)
        if db is None:
            print(f"missing index for {cancer}, skipping", flush=True)
            continue
        con = sqlite3.connect(str(db))
        n_total = 0
        n_zero = 0
        for gene in group["gene"]:
            try:
                n = con.execute("SELECT COUNT(DISTINCT pmid) FROM postings WHERE token = ?",
                                 (gene.lower(),)).fetchone()[0]
            except Exception:
                n = 0
            rows.append({"cancer": cancer, "gene": gene, "full_corpus_postings": n})
            n_total += 1
            if n == 0:
                n_zero += 1
        con.close()
        print(f"{cancer:20} {n_total:4} missed genes, {n_zero:4} truly zero full-corpus postings "
              f"({100*n_zero/n_total:.1f}%)", flush=True)

    df = pd.DataFrame(rows)
    df["false_alias_excluded"] = (df["cancer"] == "Brain Cancer") & (df["gene"] == "WAS")
    df.to_csv(OUT, index=False)
    print(f"\nWrote {OUT} ({len(df)} rows)")

    clean = df[~df["false_alias_excluded"]].copy()
    n_total_pairs = 1740  # full IntOGen pair set this 986-row missed list is drawn from
    n_never_surfaced = len(df)  # 986, unchanged regardless of the WAS exclusion (WAS stays in the "missed" set)
    zero = (clean["full_corpus_postings"] == 0).sum()
    # WAS is counted under "present" for totals (it has a nonzero raw count), but that count is
    # known to be dominated by the "was" word collision and is not used for any magnitude claim.
    present = (clean["full_corpus_postings"] > 0).sum() + 1
    print(f"Excluding WAS (false alias) from magnitude stats, counting it as present (nonzero, unreliable magnitude) "
          f"for totals: {len(clean)} of {len(df)} pairs magnitude-scored")
    print(f"  Truly zero full-corpus postings: {zero} ({100*zero/n_total_pairs:.1f}% of all 1,740 pairs; "
          f"{100*zero/n_never_surfaced:.1f}% of the 986 never-surfaced)")
    print(f"  Present but outside the 1,000-paper window: {present} ({100*present/n_total_pairs:.1f}% of all 1,740; "
          f"{100*present/n_never_surfaced:.1f}% of the 986)")
    nonzero = clean[clean["full_corpus_postings"] > 0]
    print(f"  median postings (nonzero) = {nonzero['full_corpus_postings'].median():.0f}")
    print(f"  max postings (nonzero) = {nonzero['full_corpus_postings'].max():.0f}")
    print(f"  >=10 postings: {(nonzero['full_corpus_postings']>=10).sum()}")
    print(f"  >=50 postings: {(nonzero['full_corpus_postings']>=50).sum()}")
    print(f"  >=100 postings: {(nonzero['full_corpus_postings']>=100).sum()}")

    summary = pd.DataFrame([
        {"category": "surfaced_in_window", "n": 754, "pct_of_1740": 100 * 754 / n_total_pairs},
        {"category": "present_but_outside_window", "n": int(present), "pct_of_1740": 100 * present / n_total_pairs},
        {"category": "truly_absent_from_corpus", "n": int(zero), "pct_of_1740": 100 * zero / n_total_pairs},
    ])
    summary.to_csv(OUT.parent / "intogen_full_corpus_presence_summary.csv", index=False)
    print(f"\n{summary.to_string(index=False)}")


if __name__ == "__main__":
    main()
