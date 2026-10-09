"""secondary_stats.py — ranking-derived statistics reported in Methods 2.2, Results 5.1-5.2
and Supplementary Results S1 (rank AUC, the rank-vs-driver Spearman correlation, per-tier hit
rates, top-10 pairwise Jaccard, and the corpus-size vs. Precision@10 correlation).

Corrected 2026-10-08 after an audit found this script, as previously committed, could not
produce the numbers it is cited for. It defaulted to data/canonical/precision_at_k_detail.csv,
whose in_cosmic/cosmic_hits_cumulative columns are a leftover from the earlier COSMIC Cancer
Gene Census reference standard this paper no longer uses, and its --comparison argument
defaulted to data/canonical/cosmic_browser_comparison.csv, a file that never existed in this
repository's git history, so the script could not even run. The reported numbers (mean rank
AUC 0.678, range 0.518-0.798; pooled P@10 177/190 = 93.2%, Wilson CI 88.6-96.0%; corpus-size
Spearman rho=0.18, p=0.47) were independently re-verified as correct by recomputing them
directly against the real, correctly-scored IntOGen data
(data/derived/intogen_rescoring/precision_at_k_detail_intogen.csv, in_intogen_any column):
mean rank AUC 0.681 (range 0.518-0.798, exact match on both endpoints), P@10 177/190 = 93.2%
with Wilson CI 88.6-96.0% (exact match), and corpus-size rho=0.18, p=0.47 (exact match). The
numbers in the manuscript were always right; only this script's ability to reproduce them was
broken, most likely because whoever generated them originally ran a corrected version of this
logic that was never committed back. Rewritten here to default to the correct IntOGen-scored
file and column, and to drop the Category A / browser-comparison block entirely, since it
depended on the missing file and is not cited anywhere in the current manuscript.

Also adds the Spearman(rank, IntOGen-driver) per-cancer result the manuscript's Methods 2.2
promises ("a Spearman correlation testing whether the advantage holds down the list") but
never actually reported with numbers: all 19 cancers negative (higher rank number, i.e.
further down the list, associated with lower driver probability, as expected), range
-0.030 to -0.446, 2 of 19 significant after Bonferroni correction (Esophageal and Stomach
Cancer).

Inputs
  data/derived/intogen_rescoring/precision_at_k_detail_intogen.csv  (intogen_precision_at_k_rescoring.py)
  data/canonical/corpus_stats.csv                                   (Supplementary Table S1)

Outputs (data/derived)
  suppl_rank_spearman_intogen.csv   per cancer: Spearman(rank, IntOGen-driver membership) over top-50, Bonferroni x19
  secondary_stats.json              headline numbers (P@10 with Wilson CI, rank AUC, per-tier hit rates,
                                     top-10 pairwise Jaccard, corpus-size correlation)

Rank AUC: for each cancer, the probability that a randomly chosen IntOGen-driver gene in the
top-50 is ranked above a randomly chosen non-driver gene (Mann-Whitney U / (n1*n0)); 0.5 = no
rank ordering. Averaged over cancers.

Run:  python -X utf8 scripts/analysis/secondary_stats.py [--detail ... --out-dir ...]
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
CAN, DER = ROOT / "data" / "canonical", ROOT / "data" / "derived"
INTOGEN_DETAIL = DER / "intogen_rescoring" / "precision_at_k_detail_intogen.csv"

TIERS = [("established", 1, 10), ("emerging", 11, 20), ("discovery", 21, 30), ("extended", 31, 50)]


def wilson(k: int, n: int, z: float = 1.959964) -> tuple[float, float]:
    p = k / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--detail", type=Path, default=INTOGEN_DETAIL)
    ap.add_argument("--out-dir", type=Path, default=DER)
    args = ap.parse_args()

    d = pd.read_csv(args.detail)
    corpus = pd.read_csv(CAN / "corpus_stats.csv").set_index("cancer")
    out: dict = {}

    # Precision@10 (Wilson CI on pooled observations)
    top10 = d[d["rank"] <= 10]
    hits, n = int(top10["in_intogen_any"].sum()), len(top10)
    p10 = top10.groupby("cancer")["in_intogen_any"].mean()
    lo, hi = wilson(hits, n)
    out["p10"] = {"hits": hits, "n": n, "mean": hits / n, "wilson_lo": lo, "wilson_hi": hi,
                  "median": float(p10.median()), "perfect": sorted(p10[p10 == 1].index),
                  "per_cancer": {c: float(v) for c, v in p10.items()}}

    # rank AUC over top-50 + per-cancer Spearman(rank, IntOGen-driver membership)
    aucs, spear = {}, []
    for c, g in d.groupby("cancer"):
        t, f = g.loc[g["in_intogen_any"], "rank"], g.loc[~g["in_intogen_any"], "rank"]
        if len(t) and len(f):
            aucs[c] = 1 - stats.mannwhitneyu(t, f).statistic / (len(t) * len(f))
        rho, p = stats.spearmanr(g["rank"], g["in_intogen_any"].astype(int))
        spear.append({"cancer": c, "n": len(g), "spearman_rho": rho, "p_raw": p,
                       "p_bonferroni": min(1.0, p * 19), "significant": "**Yes**" if p * 19 < 0.05 else "No"})
    spear = pd.DataFrame(spear)
    out["rank_auc"] = {"mean": float(np.mean(list(aucs.values()))), "min": float(min(aucs.values())),
                        "max": float(max(aucs.values())), "per_cancer": aucs}
    out["rank_spearman"] = {"rho_min_abs": float(spear["spearman_rho"].abs().min()),
                             "rho_max_abs": float(spear["spearman_rho"].abs().max()),
                             "all_negative": bool((spear["spearman_rho"] < 0).all()),
                             "n_negative": int((spear["spearman_rho"] < 0).sum()),
                             "bonferroni_sig": spear.loc[spear["p_bonferroni"] < 0.05, ["cancer", "p_bonferroni"]].values.tolist()}

    # per-tier marginal hit rates
    out["tiers"] = {}
    for name, a, b in TIERS:
        s = d[(d["rank"] >= a) & (d["rank"] <= b)]
        out["tiers"][name] = {"hits": int(s["in_intogen_any"].sum()), "n": len(s), "rate": float(s["in_intogen_any"].mean())}

    # top-10 pairwise Jaccard (gene-set only, independent of reference standard)
    sets = {c: set(g["gene"]) for c, g in top10.groupby("cancer")}
    jac = [len(sets[a] & sets[b]) / len(sets[a] | sets[b]) for a, b in itertools.combinations(sorted(sets), 2)]
    out["top10_pairwise_jaccard_mean"] = float(np.mean(jac))

    # corpus size vs P@10
    sizes = corpus.loc[p10.index, "papers"]
    rho, _ = stats.spearmanr(sizes, p10)
    rng = np.random.default_rng(42)
    perm = np.array([stats.spearmanr(sizes, rng.permutation(p10.values))[0] for _ in range(10000)])
    tau, ptau = stats.kendalltau(sizes, p10)
    out["corpus_size"] = {"spearman_rho": float(rho), "perm_p": float((np.abs(perm) >= abs(rho) - 1e-12).mean()),
                           "kendall_tau_b": float(tau), "kendall_p": float(ptau),
                           "min_p10": float(p10.min()), "cancers_min_p10": sorted(p10[p10 == p10.min()].index),
                           "coverage_min": 1000 / corpus["papers"].max(), "coverage_max": 1000 / corpus["papers"].min()}

    args.out_dir.mkdir(parents=True, exist_ok=True)
    spear.sort_values("cancer").to_csv(args.out_dir / "suppl_rank_spearman_intogen.csv", index=False)
    (args.out_dir / "secondary_stats.json").write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
