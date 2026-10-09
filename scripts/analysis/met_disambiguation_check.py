"""
met_disambiguation_check.py -- diagnostic only, not yet wired into the scoring pipeline.
MET is a confirmed high-confidence IntOGen driver (q as low as 5.7e-6 in LUAD, 2.97e-7 in
NSCLC, up to 32% cohort mutation rate), so this checks whether targeted secondary filters
can separate genuine MET gene mentions from confirmed false-positive sources, rather than
excluding the symbol outright the way the 7 GENE_EXCLUDE symbols are handled.

Confirmed false-positive sources (by manual inspection of raw, case-preserved text):
  1. The ordinary English word "met" (lowercase -- separated by requiring standalone caps).
  2. "Mesenchymal-to-epithelial transition" (a cancer-biology process acronym).
  3. "Methioninase" (a *Pseudomonas putida* bacterial gene used in gene therapy).
  4. "MET-PET" (a radiolabelled-methionine PET imaging tracer).
  5. "Mucoepidermoid tumor" (abbreviated MET in at least one sampled paper).
  6. "Metronomic [therapy]" (a dosing-regimen term abbreviated MET in some papers).
  7. The amino acid methionine in codon-substitution notation (e.g. "ATG(MET)", "Arg->Met").

Deliberately NOT used: a "co-occurs with other known gene symbols" heuristic (e.g. "EGFR,
MET, HER2, and AKT"). An earlier version of this script included that as a positive signal,
but it is a generic cross-gene rule, not a MET-specific confirmed pattern, and no other gene
in this paper's vocabulary QC gets that treatment -- applying it only to MET would be an
inconsistent, bespoke exception. Every term below is instead a documented pattern specific
to MET's own confirmed false-positive sources, the same evidentiary standard used for every
other symbol fix in this paper (MB, WAS, MAX, GC/HR/HCCS, etc.).

Run: python -X utf8 scripts/analysis/met_disambiguation_check.py
"""
from __future__ import annotations

import os

import json
import re
import sqlite3
from pathlib import Path

DB = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET")) / "data" / "pubmed_processed" / "lung_cancer_inclusive_all" / "lexical_index.db"
RAW = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET")) / "data" / "pubmed_raw" / "lung_cancer_inclusive" / "abstracts.jsonl"

CAPS_MET = re.compile(r"\bMET\b")

POSITIVE_TERMS = [
    "c-met", "proto-oncogene", "protooncogene", "receptor tyrosine kinase",
    "hepatocyte growth factor", "hgf", "amplification", "exon 14", "crizotinib",
    "capmatinib", "tepotinib", "savolitinib", "kinase inhibitor", "juxtamembrane",
    "met receptor", "met gene", "met oncogene", "met mutation", "met inhibitor",
    "met amplif", "met kinase", "met signaling", "met signalling", "met pathway",
    "met alterations", "copy number",
]
NEGATIVE_TERMS = [
    "mesenchymal-to-epithelial transition", "mesenchymal to epithelial transition",
    "mesenchymal-epithelial transition", "methioninase", "pseudomonas putida",
    "met-pet", "methionine pet", "11c-methionine", "c-11 methionine", "mucoepidermoid",
    "metronomic",
]
# Arg->Met / AGG(ARG)...ATG(MET) style codon-substitution notation for the amino acid.
CODON_PATTERN = re.compile(r"(arg|lys|leu|val|ile|thr)\W{0,15}met\b", re.IGNORECASE)


def classify(text: str) -> str:
    low = text.lower()
    has_pos = any(t in low for t in POSITIVE_TERMS)
    has_neg = any(t in low for t in NEGATIVE_TERMS) or bool(CODON_PATTERN.search(text))

    if has_neg and not has_pos:
        return "false_positive"
    if has_pos and not has_neg:
        return "genuine"
    if has_pos and has_neg:
        return "ambiguous_both_signals"
    return "ambiguous_neither_signal"


def main() -> None:
    con = sqlite3.connect(str(DB))
    met_pmids = {str(p[0]) for p in con.execute("SELECT DISTINCT pmid FROM postings WHERE token = 'met'").fetchall()}
    print(f"met-tagged pmids (lowercase token match): {len(met_pmids)}")

    results = {"genuine": [], "false_positive": [],
               "ambiguous_both_signals": [], "ambiguous_neither_signal": [], "no_caps_met": []}
    with open(RAW, encoding="utf-8") as f:
        for line in f:
            rec = json.loads(line)
            pmid = str(rec.get("pmid", ""))
            if pmid not in met_pmids:
                continue
            text = (rec.get("title") or "") + " " + (rec.get("abstract") or "")
            if not CAPS_MET.search(text):
                results["no_caps_met"].append(pmid)
                continue
            results[classify(text)].append(pmid)

    total = sum(len(v) for v in results.values())
    print(f"\ntotal checked: {total}")
    for bucket, pmids in results.items():
        print(f"  {bucket}: {len(pmids)} ({100 * len(pmids) / total:.1f}%)")

    confident_genuine = len(results["genuine"])
    print(f"\nconfident genuine count: {confident_genuine} ({100 * confident_genuine / total:.1f}%)")
    print(f"remaining unresolved (ambiguous + no-caps excluded as non-genuine): "
          f"{len(results['ambiguous_both_signals']) + len(results['ambiguous_neither_signal'])} ambiguous")


if __name__ == "__main__":
    main()
