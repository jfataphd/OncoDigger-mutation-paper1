"""
intogen_ablation_panel.py -- manuscript1.1 (IntOGen-primary draft) ONLY. Scores every
control variant already computed for the CGC-based manuscript1.md ablation story
(full BM25+enrichment pipeline, frequency-only baselines, vocabulary-QC on/off)
against IntOGen instead, using the identical gene-level rankings already on disk
(no corpus, retrieval, or scoring-formula changes -- only the reference standard
used to grade the output changes). Does not modify or depend on any CGC-based
script. All outputs are new files under data/derived/intogen_rescoring/.

Inputs (already produced for the CGC-based analysis; reused unchanged):
  data/canonical/precision_at_k_detail.csv          -- full pipeline, OncoDigger
  data/derived/freq_baseline_top50_aliases.csv      -- frequency-only, vetted + aliases
  data/derived/freq_baseline_top50_primary.csv      -- frequency-only, primary symbols only
  data/derived/uncleaned_vocabulary_precision.csv   -- full pipeline, vetted vs. raw vocab

Run: python -X utf8 scripts/analysis/intogen_ablation_panel.py
"""
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
INTOGEN_DRIVERS = ROOT / "data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"
OUT_DIR = ROOT / "data/derived/intogen_rescoring"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CANCER_MAP = {
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


def load_intogen():
    genes_any, pairs = set(), set()
    with open(INTOGEN_DRIVERS, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genes_any.add(row["SYMBOL"])
            for cancer, codes in CANCER_MAP.items():
                if row["CANCER_TYPE"] in codes:
                    pairs.add((row["SYMBOL"], cancer))
    return genes_any, pairs


def score(df: pd.DataFrame, genes_any: set, pairs: set) -> pd.DataFrame:
    df = df.copy()
    df["intogen_any"] = df["gene"].isin(genes_any)
    df["intogen_exact"] = df.apply(lambda r: (r["gene"], r["cancer"]) in pairs, axis=1)
    return df


def summarize(name: str, df: pd.DataFrame) -> dict:
    n = len(df)
    return {
        "variant": name, "n": n,
        "any_hits": int(df["intogen_any"].sum()), "any_pct": 100 * df["intogen_any"].mean(),
        "exact_hits": int(df["intogen_exact"].sum()), "exact_pct": 100 * df["intogen_exact"].mean(),
    }


def main():
    genes_any, pairs = load_intogen()
    rows = []
    detail_frames = {}

    full = pd.read_csv(ROOT / "data/canonical/precision_at_k_detail.csv")
    full = full[full["rank"] <= 10][["cancer", "gene"]]
    full_scored = score(full, genes_any, pairs)
    detail_frames["full_pipeline"] = full_scored
    rows.append(summarize("full_pipeline (BM25+enrichment, vetted vocab)", full_scored))

    freq_al = pd.read_csv(ROOT / "data/derived/freq_baseline_top50_aliases.csv")
    freq_al = freq_al[freq_al["rank"] <= 10][["cancer", "gene"]]
    freq_al_scored = score(freq_al, genes_any, pairs)
    detail_frames["frequency_only_vetted_aliases"] = freq_al_scored
    rows.append(summarize("frequency_only (vetted vocab + aliases, no BM25/enrichment)", freq_al_scored))

    freq_pr = pd.read_csv(ROOT / "data/derived/freq_baseline_top50_primary.csv")
    freq_pr = freq_pr[freq_pr["rank"] <= 10][["cancer", "gene"]]
    freq_pr_scored = score(freq_pr, genes_any, pairs)
    detail_frames["frequency_only_primary_symbols"] = freq_pr_scored
    rows.append(summarize("frequency_only (primary HGNC symbols only, strictest)", freq_pr_scored))

    unc = pd.read_csv(ROOT / "data/derived/uncleaned_vocabulary_precision.csv")
    vetted_rows = [{"cancer": r["cancer"], "gene": g} for _, r in unc.iterrows() for g in r["top10_vetted"].split()]
    raw_rows = [{"cancer": r["cancer"], "gene": g} for _, r in unc.iterrows() for g in r["top10_raw"].split()]
    vetted_scored = score(pd.DataFrame(vetted_rows), genes_any, pairs)
    raw_scored = score(pd.DataFrame(raw_rows), genes_any, pairs)
    detail_frames["bm25_vetted_no_longest_match_fix"] = vetted_scored
    detail_frames["bm25_raw_uncleaned_vocab"] = raw_scored
    rows.append(summarize("bm25_pipeline (vetted vocab, no longest-match fix)", vetted_scored))
    rows.append(summarize("bm25_pipeline (raw, uncleaned vocab)", raw_scored))

    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_DIR / "ablation_panel_summary.csv", index=False)
    for name, df in detail_frames.items():
        df.to_csv(OUT_DIR / f"ablation_{name}_detail.csv", index=False)

    print(summary.to_string(index=False))

    # Paired tests: full pipeline vs. each control, on exact concordance (per-cancer)
    def per_cancer_exact(df):
        return df.groupby("cancer")["intogen_exact"].mean()

    full_pc = per_cancer_exact(full_scored)
    print("\nPaired comparisons vs. full pipeline (exact concordance, per-cancer):")
    comparisons = {
        "frequency_only_vetted_aliases": freq_al_scored,
        "bm25_raw_uncleaned_vocab": raw_scored,
    }
    for name, df in comparisons.items():
        other_pc = per_cancer_exact(df)
        cancers = sorted(full_pc.index)
        a, b = full_pc.loc[cancers].values, other_pc.loc[cancers].values
        tt, pt = stats.ttest_rel(a, b)
        print(f"  full vs {name}: mean diff = {100*(a-b).mean():.1f} pp, paired t-test p={pt:.2e}")


if __name__ == "__main__":
    main()
