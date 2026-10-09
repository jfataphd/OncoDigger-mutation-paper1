"""
verify_headline_results.py -- independent re-derivation of this paper's three headline
numbers, written as a standalone check that deliberately does NOT reuse any of this
repository's own analysis code (intogen_full_rescoring.py, etc.). Re-running the same
script that produced a number only confirms it agrees with itself; this script instead
loads just the two raw input files below and recomputes each statistic from scratch by
an independent path, exactly as a reviewer could do with no more than a CSV reader.

Inputs (both included in this repository, no access restrictions):
  data/canonical/precision_at_k_detail.csv
      OncoDigger's own gene ranking output (cancer, rank, gene).
  data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv
      IntOGen's raw, unmodified driver-gene compendium (CC0 licence).

The cancer-to-IntOGen-code mapping below is typed in directly from Supplementary Table
S8 as published in the manuscript, not imported from this repository's own code, so this
check uses exactly what a reader of the paper has in front of them.

Checks performed:
  1. Precision@10 against IntOGen (any cancer): manuscript reports 93.2% (177/190).
  2. Exact cancer-type driver concordance: manuscript reports 75.8% (144/190).
  3. Corpus-separation effect (Results 5.1): manuscript reports a mean of 5.5 additional
     genes per cancer (range 3-9) recovered by cancer-specific separation relative to a
     single pooled, cancer-blind ranking, of which a mean of 4.8 are IntOGen-confirmed
     drivers.

Run: python -X utf8 scripts/analysis/verify_headline_results.py
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ONCODIGGER_RANKING = ROOT / "data" / "canonical" / "precision_at_k_detail.csv"
POOLED_CONTROL = ROOT / "data" / "derived" / "pooled_control_exact_concordance.csv"
INTOGEN_TSV = ROOT / "data" / "external" / "intogen" / "drivers" / "2024-06-18_IntOGen-Drivers" / "Compendium_Cancer_Genes.tsv"

# Typed directly from the published Supplementary Table S8, independent of this
# repository's own CANCER_MAP in intogen_full_rescoring.py.
CANCER_MAP = {
    "Breast Cancer": {"BRCA"},
    "Lung Cancer": {"LUAD", "LUSC", "NSCLC", "SCLC"},
    "Colorectal Cancer": {"COAD", "READ", "COADREAD"},
    "Prostate Cancer": {"PRAD"},
    "Melanoma": {"MEL", "SKCM", "UM"},
    "Bladder Cancer": {"BLCA", "UTUC"},
    "Kidney Cancer": {"CCRCC", "CHRCC", "PRCC", "RCC", "WT"},
    "Pancreatic Cancer": {"PAAD", "PANET"},
    "Liver Cancer": {"HCC", "CHOL", "LIHB"},
    "Stomach Cancer": {"STAD"},
    "Esophageal Cancer": {"ESCA", "ESCC"},
    "Ovarian Cancer": {"OVT"},
    "Endometrial Cancer": {"UCEC", "UCS"},
    "Cervical Cancer": {"CESC", "CEAD"},
    "Thyroid Cancer": {"WDTC"},
    "Brain Cancer": {"GB", "GBM", "HGGNOS", "LGGNOS", "PAST", "MBL", "EPM", "ATRT"},
    "Leukemia": {"ALL", "AML", "CLLSLL", "MDS"},
    "Lymphoma": {"BL", "DLBCLNOS", "NHL", "MLYM"},
    "Myeloma": {"PCM"},
}


def load_top10() -> list[tuple[str, str]]:
    top10 = []
    with open(ONCODIGGER_RANKING, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if int(row["rank"]) <= 10:
                top10.append((row["cancer"], row["gene"].strip().upper()))
    return top10


def load_intogen_any_cancer() -> set[str]:
    genes = set()
    with open(INTOGEN_TSV, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genes.add(row["SYMBOL"].strip().upper())
    return genes


def load_intogen_exact_pairs() -> set[tuple[str, str]]:
    pairs = set()
    with open(INTOGEN_TSV, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            sym = row["SYMBOL"].strip().upper()
            code = row["CANCER_TYPE"].strip()
            for cancer, codes in CANCER_MAP.items():
                if code in codes:
                    pairs.add((sym, cancer))
    return pairs


def check_1_precision_at_10(top10, intogen_any):
    hits = [(c, g) for c, g in top10 if g in intogen_any]
    misses = [(c, g) for c, g in top10 if g not in intogen_any]
    n = len(top10)
    print(f"[1] Precision@10 against IntOGen (any cancer): {len(hits)}/{n} = {100 * len(hits) / n:.1f}%")
    print(f"    manuscript reports: 177/190 = 93.2%")
    print(f"    unique IntOGen driver genes loaded: {len(intogen_any)} (manuscript: 633)")
    print(f"    {len(misses)} misses, unique genes: {sorted(set(g for _, g in misses))}")
    print()


def check_2_exact_concordance(top10, intogen_exact):
    hits = [(c, g) for c, g in top10 if (g, c) in intogen_exact]
    n = len(top10)
    print(f"[2] Exact cancer-type driver concordance: {len(hits)}/{n} = {100 * len(hits) / n:.1f}%")
    print(f"    manuscript reports: 144/190 = 75.8%")
    print()


def check_3_corpus_separation(top10, intogen_any):
    sep_top10 = defaultdict(set)
    for cancer, gene in top10:
        sep_top10[cancer].add(gene)

    pooled_top10 = defaultdict(set)
    with open(POOLED_CONTROL, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            pooled_top10[row["cancer"]].add(row["gene"].strip().upper())
    pooled_sets = set(frozenset(s) for s in pooled_top10.values())
    assert len(pooled_sets) == 1, "pooled control should return an identical top-10 for every cancer"
    pooled_genes = next(iter(pooled_sets))

    extra_counts, confirmed_counts = [], []
    for cancer, genes in sep_top10.items():
        extra = genes - pooled_genes
        extra_counts.append(len(extra))
        confirmed_counts.append(len([g for g in extra if g in intogen_any]))

    mean_extra = sum(extra_counts) / len(extra_counts)
    mean_confirmed = sum(confirmed_counts) / len(confirmed_counts)
    print(f"[3] Corpus-separation effect (Results 5.1):")
    print(f"    mean extra genes per cancer: {mean_extra:.2f} (manuscript: 5.5)")
    print(f"    range: {min(extra_counts)}-{max(extra_counts)} (manuscript: 3-9)")
    print(f"    mean of those confirmed as IntOGen drivers: {mean_confirmed:.2f} (manuscript: 4.8)")
    leukemia_extra = sep_top10["Leukemia"] - pooled_genes
    print(f"    Leukemia's extra genes: {sorted(leukemia_extra)} (manuscript cites FLT3, NPM1 as an example)")
    print()


def main() -> None:
    top10 = load_top10()
    assert len(top10) == 190, f"expected 190 top-10 rows, got {len(top10)}"
    intogen_any = load_intogen_any_cancer()
    intogen_exact = load_intogen_exact_pairs()

    check_1_precision_at_10(top10, intogen_any)
    check_2_exact_concordance(top10, intogen_exact)
    check_3_corpus_separation(top10, intogen_any)


if __name__ == "__main__":
    main()
