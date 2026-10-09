"""
intogen_figure3_data.py -- manuscript1.1 (IntOGen-primary draft) ONLY. Produces the three
IntOGen-scored data files Figure 3 needs, mirroring the CGC-based files used by
build_figure3_validation.py exactly in shape, scored against IntOGen's any-cancer driver
set instead of COSMIC CGC Tier 1. None of the underlying rankings are re-run -- same
OncoDigger full ranking (data/derived/query_planes/mutation.csv), same frequency-only
baseline ranking (data/derived/freq_baseline_top50_aliases.csv), same pooled-corpus
control ranking (data/derived/pooled_control_mutation.csv) -- only which gene set counts
as a hit changes.

Outputs (data/derived/intogen_rescoring/):
  freq_baseline_per_cancer_detail_intogen.csv   mirrors freq_baseline_per_cancer_detail.csv
  freq_baseline_precision_at_k_aliases_intogen.csv  mirrors freq_baseline_precision_at_k_aliases.csv
  corpus_specificity_clean_intogen.csv          mirrors corpus_specificity_clean.csv

Run: python -X utf8 scripts/analysis/intogen_figure3_data.py
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
DER = ROOT / "data" / "derived"
OUT_DIR = DER / "intogen_rescoring"
OUT_DIR.mkdir(parents=True, exist_ok=True)
INTOGEN_DRIVERS = ROOT / "data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"
KS = [5, 10, 15, 20, 30, 40, 50]
N_BOOT, SEED = 2000, 42


def load_genes_any() -> set:
    genes = set()
    with open(INTOGEN_DRIVERS, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genes.add(row["SYMBOL"].upper())
    return genes


def cum_precision(rk: pd.DataFrame, genes_any: set) -> pd.DataFrame:
    """rk: columns cancer, rank, gene (sorted by cancer, rank). Returns cancer, rank, precision_at_k."""
    rk = rk.sort_values(["cancer", "rank"]).copy()
    rk["hit"] = rk["gene"].str.upper().isin(genes_any)
    rk["cum"] = rk.groupby("cancer")["hit"].cumsum()
    rk["precision_at_k"] = rk["cum"] / rk["rank"]
    return rk


def case_bootstrap_ci(x: np.ndarray, n_boot=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    n = len(x)
    boots = np.array([x[rng.integers(0, n, n)].mean() for _ in range(n_boot)])
    return np.percentile(boots, [2.5, 97.5])


def main():
    genes_any = load_genes_any()

    # Reads the same frozen canonical ranking every other figure/table in this manuscript
    # uses (data/canonical/precision_at_k_detail.csv), not a live rerun of the engine --
    # a live rerun (query_planes/mutation.csv) reproduces this exactly through k=40 but can
    # drift by ~0.1pp at k=50 as the live corpus grows past the frozen snapshot date.
    od_rk = pd.read_csv(ROOT / "data" / "canonical" / "precision_at_k_detail.csv")[["cancer", "rank", "gene"]]
    fq_rk = pd.read_csv(DER / "freq_baseline_top50_aliases.csv")[["cancer", "rank", "gene"]]

    od = cum_precision(od_rk, genes_any)
    fq = cum_precision(fq_rk, genes_any)

    # per-cancer wide detail
    rows = []
    for c in sorted(od["cancer"].unique()):
        row = {"cancer": c}
        for k in KS:
            od_p = od[(od["cancer"] == c) & (od["rank"] == k)]["precision_at_k"].iloc[0]
            fq_p = fq[(fq["cancer"] == c) & (fq["rank"] == k)]["precision_at_k"].iloc[0]
            row[f"OD@{k}"] = od_p
            row[f"Freq@{k}"] = fq_p
            row[f"Δ@{k}"] = 100 * (od_p - fq_p)
        rows.append(row)
    detail = pd.DataFrame(rows)
    detail.to_csv(OUT_DIR / "freq_baseline_per_cancer_detail_intogen.csv", index=False)

    # pooled summary with bootstrap CI + paired tests
    summ_rows = []
    for k in KS:
        a = detail[f"OD@{k}"].values
        b = detail[f"Freq@{k}"].values
        delta = a - b
        lo, hi = case_bootstrap_ci(delta)
        tt, pt = stats.ttest_rel(a, b)
        try:
            ww, pw = stats.wilcoxon(a, b)
        except ValueError:
            pw = float("nan")
        summ_rows.append({
            "k": k, "OD_mean": a.mean(), "Freq_mean": b.mean(), "delta_pp": 100 * delta.mean(),
            "ci_low_pp": 100 * lo, "ci_high_pp": 100 * hi, "p_t": pt, "p_wilcoxon": pw,
            "od_wins": int((a > b).sum()), "ties": int((a == b).sum()), "freq_wins": int((a < b).sum()),
        })
    pd.DataFrame(summ_rows).to_csv(OUT_DIR / "freq_baseline_precision_at_k_aliases_intogen.csv", index=False)

    # corpus specificity: OD top-10 genes per cancer absent from the pooled-corpus top-10
    pooled = pd.read_csv(DER / "pooled_control_mutation.csv").sort_values("rank")
    global_top10 = set(pooled.loc[pooled["rank"] <= 10, "gene"])
    top10 = {c: list(g.sort_values("rank")["gene"]) for c, g in od_rk[od_rk["rank"] <= 10].groupby("cancer")}

    spec_rows = []
    for c, genes in top10.items():
        s = set(genes)
        od_p10 = sum(1 for g in s if g.upper() in genes_any) / 10
        specific = sorted(s - global_top10)
        spec_rows.append({
            "cancer": c, "OD_p10": od_p10, "overlap_with_global": len(s & global_top10),
            "cancer_specific_in_top10": len(specific),
            "cancer_specific_intogen": sum(1 for g in specific if g.upper() in genes_any),
            "cancer_specific_genes": "; ".join(specific),
        })
    pd.DataFrame(spec_rows).sort_values("cancer").to_csv(OUT_DIR / "corpus_specificity_clean_intogen.csv", index=False)

    print("Wrote freq_baseline_per_cancer_detail_intogen.csv, freq_baseline_precision_at_k_aliases_intogen.csv, "
          "corpus_specificity_clean_intogen.csv")
    for k in KS:
        r = [r for r in summ_rows if r["k"] == k][0]
        print(f"  k={k}: OD={r['OD_mean']:.3f} Freq={r['Freq_mean']:.3f} delta={r['delta_pp']:.1f}pp "
              f"p_t={r['p_t']:.2e}")


if __name__ == "__main__":
    main()
