"""
cancer_specific_drivers_pubmed_validation.py -- empirical PubMed cross-validation for the
91 corpus-separation-recovered genes that are confirmed IntOGen drivers (Results 5.x:
"median supporting PubMed count was 368 papers, IQR 179-799").

For each (cancer, gene) pair in corpus_specific_genes_full_list.csv with is_intogen=True,
queries NCBI esearch for:

  ("<cancer MeSH term>"[mh]) AND (<gene>[tiab]) AND (2000/01/01:2026/12/31[pdat])

A high paper count confirms the gene has substantial dedicated literature in that cancer,
not a sporadic mention -- a publication-volume cross-validation independent of IntOGen's own
driver-detection methodology.

Reconstructed 2026-10-08: this script's output had been committed without a generating
script anywhere in this repo (adapted here from a COSMIC-era predecessor,
validate_cancer_specific_drivers_pubmed.py, found in the author's archived early manuscript
drafts, updated to source its (cancer, gene) pairs from the current IntOGen-based
corpus_specific_genes_full_list.csv instead of a COSMIC Tier 1 list). PubMed's index grows
continuously, so re-querying today will not reproduce the exact historical counts byte-for-
byte, only approximately -- the same caveat already documented in PROVENANCE.md for
DISEASES/PubTator. Re-running produced counts of the same order of magnitude and the same
median/IQR as the committed file (see run log).

Input:  data/derived/intogen_rescoring/corpus_specific_genes_full_list.csv
Output: data/derived/intogen_rescoring/cancer_specific_drivers_pubmed_validation_intogen.csv

Run: python -X utf8 scripts/analysis/cancer_specific_drivers_pubmed_validation.py
"""
from __future__ import annotations

import csv
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[2]
DER = ROOT / "data" / "derived" / "intogen_rescoring"
SPEC_FILE = DER / "corpus_specific_genes_full_list.csv"
OUT_FILE = DER / "cancer_specific_drivers_pubmed_validation_intogen.csv"

CANCER_TO_MESH = {
    "Bladder Cancer": "Urinary Bladder Neoplasms",
    "Brain Cancer": "Brain Neoplasms",
    "Breast Cancer": "Breast Neoplasms",
    "Cervical Cancer": "Uterine Cervical Neoplasms",
    "Colorectal Cancer": "Colorectal Neoplasms",
    "Endometrial Cancer": "Endometrial Neoplasms",
    "Esophageal Cancer": "Esophageal Neoplasms",
    "Kidney Cancer": "Kidney Neoplasms",
    "Leukemia": "Leukemia",
    "Liver Cancer": "Liver Neoplasms",
    "Lung Cancer": "Lung Neoplasms",
    "Lymphoma": "Lymphoma",
    "Melanoma": "Melanoma",
    "Myeloma": "Multiple Myeloma",
    "Ovarian Cancer": "Ovarian Neoplasms",
    "Pancreatic Cancer": "Pancreatic Neoplasms",
    "Prostate Cancer": "Prostatic Neoplasms",
    "Stomach Cancer": "Stomach Neoplasms",
    "Thyroid Cancer": "Thyroid Neoplasms",
}

ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
SLEEP_SEC = 0.4  # ~3 req/sec; NCBI guideline without an API key


def esearch_count(mesh: str, gene: str, retries: int = 3) -> int:
    term = f'("{mesh}"[mh]) AND ({gene}[tiab]) AND (2000/01/01:2026/12/31[pdat])'
    for attempt in range(retries):
        try:
            r = requests.get(
                ESEARCH_URL,
                params={"db": "pubmed", "term": term, "retmode": "json", "retmax": 0},
                timeout=30,
            )
            r.raise_for_status()
            return int(r.json()["esearchresult"]["count"])
        except Exception as e:
            print(f"  Retry {attempt+1}/{retries} for {gene}/{mesh}: {e}")
            time.sleep(2)
    return -1


def main() -> None:
    pairs = []
    with open(SPEC_FILE, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["is_intogen"] == "True":
                pairs.append((row["cancer"], row["gene"]))
    print(f"Cancer-specific IntOGen-driver gene-cancer pairs: {len(pairs)}")

    rows = []
    for i, (cancer, gene) in enumerate(pairs):
        mesh = CANCER_TO_MESH.get(cancer)
        if mesh is None:
            print(f"  [{i+1}/{len(pairs)}] {cancer}: NO MESH MAPPING")
            continue
        count = esearch_count(mesh, gene)
        rows.append({"cancer": cancer, "gene": gene, "pubmed_count": count})
        print(f"  [{i+1}/{len(pairs)}] {cancer:22s} {gene:8s} -> {count:6d} papers", flush=True)
        time.sleep(SLEEP_SEC)

    with open(OUT_FILE, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["cancer", "gene", "pubmed_count"])
        w.writeheader()
        w.writerows(rows)
    counts = sorted(r["pubmed_count"] for r in rows)
    n = len(counts)
    median = counts[n // 2] if n % 2 else (counts[n // 2 - 1] + counts[n // 2]) / 2
    print(f"\nWrote {OUT_FILE}: {n} pairs, median={median}")


if __name__ == "__main__":
    main()
