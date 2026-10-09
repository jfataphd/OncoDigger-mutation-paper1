"""
intogen_confidence_breakdown.py -- manuscript1.1 (IntOGen-primary draft) ONLY. Checks whether
OncoDigger's 986 "never-surfaced" IntOGen driver pairs (Results 5.5, Supplementary Results S6)
are concentrated in IntOGen's own lower-confidence driver calls, using signals already present
in IntOGen's raw compendium: how many of its seven independent statistical methods agreed on
each call (METHODS column), the combined q-value (QVALUE_COMBINATION), and the fraction of the
cohort carrying a mutation in that gene (%_SAMPLES_COHORT). None of these three fields were used
anywhere else in this paper's scoring or mapping -- they are IntOGen's own internal confidence
signals, read here only to characterise the surfaced/never-surfaced split, not to re-score or
filter anything.

For each of the 1,740 (gene, study-cancer) driver pairs (Methods 2.3), the matching IntOGen
cohort row(s) are aggregated by taking the strongest evidence across cohorts: max method count,
min q-value, max %-of-cohort-mutated. Pairs are then split by whether OncoDigger's ranking ever
surfaced them (data/derived/intogen_rescoring/intogen_literature_weight_audit.csv), and the two
groups are compared with Mann-Whitney U tests (none of these three quantities is normally
distributed).

Inputs: data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv,
        data/derived/intogen_rescoring/intogen_literature_weight_audit.csv
Outputs: data/derived/intogen_rescoring/intogen_confidence_breakdown_detail.csv (1,740 rows)
         data/derived/intogen_rescoring/intogen_confidence_breakdown_summary.csv

Run: python -X utf8 scripts/analysis/intogen_confidence_breakdown.py
"""
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
INTOGEN_DRIVERS = ROOT / "data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"
OUT_DIR = ROOT / "data/derived/intogen_rescoring"

INTOGEN_CANCER_MAP = {
    "Breast Cancer": {"BRCA"}, "Lung Cancer": {"LUAD", "LUSC", "NSCLC", "SCLC"},
    "Colorectal Cancer": {"COAD", "READ", "COADREAD"}, "Prostate Cancer": {"PRAD"},
    "Melanoma": {"MEL", "SKCM", "UM"}, "Bladder Cancer": {"BLCA", "UTUC"},
    "Kidney Cancer": {"CCRCC", "CHRCC", "PRCC", "RCC", "WT"},
    "Pancreatic Cancer": {"PAAD", "PANET"}, "Liver Cancer": {"HCC", "CHOL", "LIHB"},
    "Stomach Cancer": {"STAD"}, "Esophageal Cancer": {"ESCA", "ESCC"},
    "Ovarian Cancer": {"OVT"}, "Endometrial Cancer": {"UCEC", "UCS"},
    "Cervical Cancer": {"CESC", "CEAD"}, "Thyroid Cancer": {"WDTC"},
    "Brain Cancer": {"GB", "GBM", "HGGNOS", "LGGNOS", "PAST", "MBL", "EPM", "ATRT"},
    "Leukemia": {"ALL", "AML", "CLLSLL", "CML", "MDS"},
    "Lymphoma": {"BL", "DLBCLNOS", "NHL", "MLYM"}, "Myeloma": {"PCM"},
}


def main():
    rows = []
    with open(INTOGEN_DRIVERS, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rows.append(row)
    df = pd.DataFrame(rows)
    df["n_methods"] = df["METHODS"].apply(lambda s: 1 if s == "combination" else len(s.split(",")))
    df["qval"] = pd.to_numeric(df["QVALUE_COMBINATION"], errors="coerce")
    df["pct_samples"] = pd.to_numeric(df["%_SAMPLES_COHORT"], errors="coerce")

    code_to_cancer = {code: cancer for cancer, codes in INTOGEN_CANCER_MAP.items() for code in codes}
    df["study_cancer"] = df["CANCER_TYPE"].map(code_to_cancer)
    mapped = df.dropna(subset=["study_cancer"])

    agg = (mapped.groupby(["SYMBOL", "study_cancer"])
           .agg(n_methods=("n_methods", "max"), qval=("qval", "min"), pct_samples=("pct_samples", "max"))
           .reset_index())
    agg.columns = ["gene", "cancer", "n_methods", "qval", "pct_samples"]

    audit = pd.read_csv(OUT_DIR / "intogen_literature_weight_audit.csv")[["cancer", "gene", "engine_rank"]]
    audit["surfaced"] = audit["engine_rank"].notna()

    merged = agg.merge(audit[["cancer", "gene", "surfaced"]], on=["cancer", "gene"], how="inner")
    assert len(merged) == 1740, f"expected 1740 pairs, got {len(merged)}"
    merged["group"] = merged["surfaced"].map({True: "surfaced", False: "never_surfaced"})

    merged.to_csv(OUT_DIR / "intogen_confidence_breakdown_detail.csv", index=False)

    surf = merged[merged["surfaced"]]
    miss = merged[~merged["surfaced"]]

    summary_rows = []
    for col in ["n_methods", "qval", "pct_samples"]:
        u, p = stats.mannwhitneyu(surf[col].dropna(), miss[col].dropna(), alternative="two-sided")
        summary_rows.append({
            "metric": col,
            "surfaced_median": surf[col].median(), "surfaced_mean": surf[col].mean(), "surfaced_n": surf[col].notna().sum(),
            "never_surfaced_median": miss[col].median(), "never_surfaced_mean": miss[col].mean(), "never_surfaced_n": miss[col].notna().sum(),
            "mannwhitney_u": u, "p_value": p,
        })
    summary_rows.append({
        "metric": "share_single_method",
        "surfaced_median": None, "surfaced_mean": (surf["n_methods"] == 1).mean(), "surfaced_n": len(surf),
        "never_surfaced_median": None, "never_surfaced_mean": (miss["n_methods"] == 1).mean(), "never_surfaced_n": len(miss),
        "mannwhitney_u": None, "p_value": None,
    })
    summ = pd.DataFrame(summary_rows)
    summ.to_csv(OUT_DIR / "intogen_confidence_breakdown_summary.csv", index=False)

    print(summ.to_string(index=False))
    print()
    print("n_methods distribution, surfaced:")
    print((surf["n_methods"].value_counts(normalize=True).sort_index() * 100).round(1))
    print("n_methods distribution, never surfaced:")
    print((miss["n_methods"].value_counts(normalize=True).sort_index() * 100).round(1))


if __name__ == "__main__":
    main()
