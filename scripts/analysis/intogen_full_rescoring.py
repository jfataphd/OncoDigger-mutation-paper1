"""
intogen_full_rescoring.py -- NEW, EXPLORATORY ANALYSIS. Does not modify or depend on
any manuscript-facing CGC-based script or output. All outputs are written to new
files under data/derived/intogen_rescoring/ so nothing already used by
OncoDigger_manuscript1.md is touched.

Purpose: score all four systems already benchmarked in the manuscript (OncoDigger,
CancerMine, ChatGPT, and the cancer-blind pooled control) against IntOGen
(Martinez-Jimenez et al. 2020, release 2024.09.20) instead of COSMIC CGC, using
the exact same case-bootstrap statistical methodology already used throughout this
project (resample the 19 cancers with replacement, 2000 replicates, seed 42), so
the two reference standards (CGC and IntOGen) can be compared on genuinely equal
methodological footing before any decision is made about using this in the
manuscript.

For each system, computes:
  1. IntOGen-Precision@10 equivalent: is the gene an IntOGen driver in ANY of its
     86 cancer-type codes (analogous to CGC Tier 1 "anywhere" Precision@10).
  2. IntOGen exact concordance: is the gene an IntOGen driver specifically within
     a cohort mapped to the SAME study cancer (analogous to CGC exact tumour-type
     concordance).
Both are reported pooled (190 pairs), per-cancer, and with a case-bootstrap 95% CI.
Paired comparisons (OncoDigger vs each other system) use paired t-test, Wilcoxon
signed-rank, and a paired case-bootstrap CI of the mean difference, mirroring
head_to_head_bootstrap.py / cancer_aware_concordance_bootstrap.py exactly.

Run: python -X utf8 scripts/analysis/intogen_full_rescoring.py
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
INTOGEN_DRIVERS = ROOT / "data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"
OUT_DIR = ROOT / "data/derived/intogen_rescoring"
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_BOOT = 2000
SEED = 42

# Identical conservative mapping used in intogen_discordant_cross_check.py, kept
# in sync deliberately so the two analyses agree; duplicated here rather than
# imported so this script has no dependency on the earlier one.
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
    "Leukemia": {"ALL", "AML", "CLLSLL", "CML", "MDS"},
    "Lymphoma": {"BL", "DLBCLNOS", "NHL", "MLYM"},
    "Myeloma": {"PCM"},
}

SYSTEMS = {
    "oncodigger": ROOT / "data/canonical/precision_at_k_detail.csv",
    "cancermine": ROOT / "data/derived/cancermine_top10.csv",
    "chatgpt": ROOT / "data/derived/chatgpt_top10.csv",
    "pooled_control": ROOT / "data/derived/pooled_control_exact_concordance.csv",
}


def load_intogen():
    genes_any = set()
    pairs = set()
    with open(INTOGEN_DRIVERS, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genes_any.add(row["SYMBOL"])
            code = row["CANCER_TYPE"]
            for cancer, codes in CANCER_MAP.items():
                if code in codes:
                    pairs.add((row["SYMBOL"], cancer))
    return genes_any, pairs


def load_top10(system: str, path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if system == "oncodigger":
        df = df[df["rank"] <= 10][["cancer", "gene"]].reset_index(drop=True)
    else:
        df = df[["cancer", "gene"]].reset_index(drop=True)
    assert len(df) == 190, f"{system}: expected 190 rows, got {len(df)}"
    return df


def score(df: pd.DataFrame, genes_any: set, pairs: set) -> pd.DataFrame:
    df = df.copy()
    df["intogen_any"] = df["gene"].isin(genes_any)
    df["intogen_exact"] = df.apply(lambda r: (r["gene"], r["cancer"]) in pairs, axis=1)
    return df


def case_bootstrap_ci(per_cancer_rate: np.ndarray, n_boot=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    n = len(per_cancer_rate)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        boots[i] = per_cancer_rate[idx].mean()
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return per_cancer_rate.mean(), lo, hi


def paired_bootstrap_ci(a: np.ndarray, b: np.ndarray, n_boot=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    diffs = a - b
    n = len(diffs)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        boots[i] = diffs[idx].mean()
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return diffs.mean(), lo, hi, (boots <= 0).mean()


def per_cancer_summary(scored: pd.DataFrame) -> pd.DataFrame:
    g = scored.groupby("cancer").agg(
        n=("gene", "size"),
        p10_intogen_any=("intogen_any", "mean"),
        p10_intogen_exact=("intogen_exact", "mean"),
    ).reset_index()
    return g


def main():
    print("Loading IntOGen compendium...")
    genes_any, pairs = load_intogen()
    print(f"  {len(genes_any)} distinct driver genes (any cancer type)")
    print(f"  {len(pairs)} distinct (gene, study-cancer) driver pairs after conservative mapping\n")

    summaries = {}
    scored_all = {}
    for system, path in SYSTEMS.items():
        df = load_top10(system, path)
        scored = score(df, genes_any, pairs)
        scored.to_csv(OUT_DIR / f"{system}_intogen_detail.csv", index=False)
        summ = per_cancer_summary(scored)
        summ.to_csv(OUT_DIR / f"{system}_intogen_summary.csv", index=False)
        summaries[system] = summ
        scored_all[system] = scored

        pooled_any = scored["intogen_any"].mean()
        pooled_exact = scored["intogen_exact"].mean()
        mean_any, lo_any, hi_any = case_bootstrap_ci(summ["p10_intogen_any"].values)
        mean_exact, lo_exact, hi_exact = case_bootstrap_ci(summ["p10_intogen_exact"].values)
        print(f"=== {system} ===")
        print(f"  IntOGen-Precision@10 (any cancer):  {scored['intogen_any'].sum()}/190 = {100*pooled_any:.1f}%  "
              f"(95% CI {100*lo_any:.1f}-{100*hi_any:.1f}%)")
        print(f"  IntOGen exact concordance:          {scored['intogen_exact'].sum()}/190 = {100*pooled_exact:.1f}%  "
              f"(95% CI {100*lo_exact:.1f}-{100*hi_exact:.1f}%)")
        print()

    # Paired OncoDigger-vs-X comparisons, both metrics
    print("=" * 70)
    print("Paired comparisons vs OncoDigger (case-bootstrap over 19 cancers)")
    print("=" * 70)
    od = summaries["oncodigger"].set_index("cancer")
    rows = []
    for metric, col in [("IntOGen-Precision@10 (any)", "p10_intogen_any"),
                         ("IntOGen exact concordance", "p10_intogen_exact")]:
        for other in ["cancermine", "chatgpt", "pooled_control"]:
            o = summaries[other].set_index("cancer")
            cancers = sorted(od.index)
            a = od.loc[cancers, col].values
            b = o.loc[cancers, col].values
            wins, ties, losses = int((a > b).sum()), int((a == b).sum()), int((a < b).sum())
            tt, pt = stats.ttest_rel(a, b)
            try:
                ww, pw = stats.wilcoxon(a, b)
            except ValueError:
                pw = float("nan")
            mean_diff, lo, hi, frac_le0 = paired_bootstrap_ci(a, b)
            print(f"\n[{metric}] OncoDigger vs {other}:")
            print(f"  wins {wins}, ties {ties}, losses {losses} (of {len(cancers)} cancers)")
            print(f"  paired t-test p={pt:.4g}, Wilcoxon p={pw:.4g}")
            print(f"  mean diff = {mean_diff*100:.1f} pp, 95% bootstrap CI = [{lo*100:.1f}, {hi*100:.1f}] pp")
            rows.append({"metric": metric, "comparator": other, "wins": wins, "ties": ties, "losses": losses,
                         "p_ttest": pt, "p_wilcoxon": pw, "mean_diff_pp": mean_diff * 100,
                         "ci_low_pp": lo * 100, "ci_high_pp": hi * 100, "frac_resamples_le0": frac_le0})
    pd.DataFrame(rows).to_csv(OUT_DIR / "paired_comparisons_vs_oncodigger.csv", index=False)

    # Pooled-level summary table across all systems/metrics for quick reference
    summary_rows = []
    for system in SYSTEMS:
        scored = scored_all[system]
        summ = summaries[system]
        for metric, col_pair, col_summ in [("intogen_any", "intogen_any", "p10_intogen_any"),
                                            ("intogen_exact", "intogen_exact", "p10_intogen_exact")]:
            mean_v, lo, hi = case_bootstrap_ci(summ[col_summ].values)
            summary_rows.append({
                "system": system, "metric": metric,
                "n_hits": int(scored[col_pair].sum()), "n_total": len(scored),
                "pooled_pct": 100 * scored[col_pair].mean(),
                "ci_low_pct": 100 * lo, "ci_high_pct": 100 * hi,
            })
    pd.DataFrame(summary_rows).to_csv(OUT_DIR / "all_systems_summary.csv", index=False)
    print(f"\nAll outputs written to {OUT_DIR}")


if __name__ == "__main__":
    main()
