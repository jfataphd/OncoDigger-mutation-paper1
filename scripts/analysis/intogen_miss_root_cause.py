"""
intogen_miss_root_cause.py -- manuscript1.1 (IntOGen-primary draft) ONLY. Root-cause
classification of the pairs still never surfaced by OncoDigger's ranking after restricting
to IntOGen's higher-confidence driver calls (Supplementary Results S11: present in the full
corpus AND backed by >=2 of IntOGen's seven driver-detection methods).

For each such pair, classifies the explanation using only data already computed (full-corpus
mention count, Supplementary Results S10) plus a new, direct sample of the actual corpus text
(titles of up to 20 papers mentioning the gene in that cancer's corpus, queried directly
against the postings/documents tables, bypassing BM25 retrieval and ranking entirely), checked
for whether the sampled titles are mutation-framed (regex match on "mutat").

Categories (checked in this order):
  1. low_volume          -- fewer than 20 full-corpus mentions. Consistent with the broader
                             "present but outside window" finding (Supplementary Results S10):
                             real but sparse literature, unlikely to rank regardless of pool size.
  2. vocabulary_artefact  -- manually confirmed, by direct inspection of sampled titles, to be
                             dominated by a false-alias text collision (the gene symbol matches
                             an unrelated common word/phrase), the same class of bug already
                             documented for MB, ELN, CP, WAS etc. (Methods 3.3, Supplementary
                             Results S10). Confirmed cases are hardcoded below after manual
                             review; this script does not auto-detect them.
  3. non_mutation_framed  -- >=20 full-corpus mentions, but fewer than 20% of the sampled
                             titles contain "mutat": real, on-topic cancer literature not framed
                             around mutation status (expression, transcription-factor biology,
                             immune/HLA typing, germline risk polymorphisms, etc.), so it does
                             not compete well in a mutation-focused BM25 query regardless of
                             corpus size.
  4. unresolved           -- >=20 full-corpus mentions and >=20% of the sampled titles are
                             mutation-framed: a real, mutation-relevant literature presence that
                             still did not rank. Reported as a genuine open question.

Requires the local OncoDigger corpus databases (not in this repository).

Inputs: data/derived/intogen_rescoring/intogen_literature_weight_audit.csv,
        intogen_full_corpus_presence_check.csv, intogen_confidence_breakdown_detail.csv
Outputs: data/derived/intogen_rescoring/intogen_high_volume_mutation_titling.csv (title-sample detail)
         data/derived/intogen_rescoring/intogen_miss_root_cause_detail.csv (per-pair classification)
         data/derived/intogen_rescoring/intogen_miss_root_cause_summary.csv (category counts)

Run: python -X utf8 scripts/analysis/intogen_miss_root_cause.py
"""
from __future__ import annotations

import os
import re
import sqlite3
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_ONCODIGGER = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET"))
DER = ROOT / "data" / "derived" / "intogen_rescoring"

VOLUME_THRESHOLD = 20
MUTATION_TITLE_THRESHOLD = 20  # percent
SAMPLE_SIZE = 20

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

# Manually confirmed by direct inspection of sampled corpus titles (see module docstring);
# not auto-detected. MAX/Brain Cancer: 427 "postings" were overwhelmingly "maximum
# [dose/tolerated/...]" and the soybean species "Glycine max", not the MAX gene.
CONFIRMED_VOCAB_ARTEFACTS = {("Brain Cancer", "MAX")}


def title_sample_stats(db_base: Path, cancer: str, gene: str) -> tuple[int, float | None]:
    key = CANCER_KEY[cancer]
    db = next((p for p in (db_base / key / "lexical_index.db", db_base / f"{key}_all" / "lexical_index.db")
               if p.exists()), None)
    con = sqlite3.connect(str(db))
    token = gene.lower()
    pmids = [p[0] for p in con.execute("SELECT DISTINCT pmid FROM postings WHERE token = ?", (token,)).fetchall()]
    n = min(SAMPLE_SIZE, len(pmids))
    if n == 0:
        con.close()
        return 0, None
    placeholders = ",".join(["?"] * n)
    titles = [t[0] for t in con.execute(
        f"SELECT title FROM documents WHERE pmid IN ({placeholders})", pmids[:n]).fetchall()]
    con.close()
    mut_hits = sum(1 for t in titles if t and re.search(r"mutat", t))
    return n, 100 * mut_hits / n


def main():
    audit = pd.read_csv(DER / "intogen_literature_weight_audit.csv")
    presence = pd.read_csv(DER / "intogen_full_corpus_presence_check.csv")
    conf = pd.read_csv(DER / "intogen_confidence_breakdown_detail.csv")

    absent_pairs = set(zip(presence[presence["full_corpus_postings"] == 0]["cancer"],
                            presence[presence["full_corpus_postings"] == 0]["gene"]))
    audit["is_absent"] = audit.apply(lambda r: (r["cancer"], r["gene"]) in absent_pairs, axis=1)
    audit = audit.merge(conf[["cancer", "gene", "n_methods"]], on=["cancer", "gene"], how="left")
    audit = audit.merge(presence[["cancer", "gene", "full_corpus_postings"]], on=["cancer", "gene"], how="left")

    clean = audit[(~audit["is_absent"]) & (audit["n_methods"] >= 2)]
    missed = clean[clean["engine_rank"].isna()].copy()
    print(f"Still-missed pairs (detectable + >=2 methods): {len(missed)}")

    high_volume = missed[missed["full_corpus_postings"] >= VOLUME_THRESHOLD].copy()
    print(f"Of these, >= {VOLUME_THRESHOLD} full-corpus mentions: {len(high_volume)} -- sampling titles...")
    high_volume.sort_values("full_corpus_postings", ascending=False).to_csv(
        DER / "intogen_high_volume_misses.csv", index=False)

    db_base = DEFAULT_ONCODIGGER / "data" / "pubmed_processed"
    titling_rows = []
    for _, r in high_volume.iterrows():
        n_sampled, pct = title_sample_stats(db_base, r["cancer"], r["gene"])
        titling_rows.append({"cancer": r["cancer"], "gene": r["gene"], "n_sampled": n_sampled,
                              "pct_mutation_titled": pct})
    titling = pd.DataFrame(titling_rows)
    titling = titling.merge(high_volume[["cancer", "gene", "full_corpus_postings"]], on=["cancer", "gene"])
    titling = titling.sort_values("full_corpus_postings", ascending=False)
    titling.to_csv(DER / "intogen_high_volume_mutation_titling.csv", index=False)

    missed = missed.merge(titling[["cancer", "gene", "pct_mutation_titled"]], on=["cancer", "gene"], how="left")

    def classify(r):
        if (r["cancer"], r["gene"]) in CONFIRMED_VOCAB_ARTEFACTS:
            return "vocabulary_artefact"
        if r["full_corpus_postings"] < VOLUME_THRESHOLD:
            return "low_volume"
        if pd.notna(r["pct_mutation_titled"]) and r["pct_mutation_titled"] >= MUTATION_TITLE_THRESHOLD:
            return "unresolved"
        return "non_mutation_framed"

    missed["category"] = missed.apply(classify, axis=1)
    missed.to_csv(DER / "intogen_miss_root_cause_detail.csv", index=False)

    order = ["low_volume", "non_mutation_framed", "vocabulary_artefact", "unresolved"]
    summ = missed["category"].value_counts().reindex(order, fill_value=0).rename("n").reset_index()
    summ.columns = ["category", "n"]
    summ["pct_of_total"] = 100 * summ["n"] / len(missed)
    summ.to_csv(DER / "intogen_miss_root_cause_summary.csv", index=False)

    print()
    print(summ.to_string(index=False))
    print()
    print("Unresolved pairs:")
    print(missed[missed["category"] == "unresolved"][["cancer", "gene", "full_corpus_postings", "pct_mutation_titled"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
