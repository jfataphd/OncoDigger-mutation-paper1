"""
intogen_precision_at_k_rescoring.py -- manuscript1.1 (IntOGen-primary draft) ONLY.
Rescores the full rank-1-50 OncoDigger gene ranking (data/canonical/precision_at_k_detail.csv,
950 rows = 19 cancers x 50 ranks) against IntOGen's any-cancer driver set instead of COSMIC
CGC Tier 1. The underlying gene/rank/cancer ranking is untouched (same engine output) --
only which gene set counts as a hit, and the resulting cumulative precision_at_k, changes.
Does not modify or depend on precision_at_k_detail.csv or any CGC-based script.

Run: python -X utf8 scripts/analysis/intogen_precision_at_k_rescoring.py
Output: data/derived/intogen_rescoring/precision_at_k_detail_intogen.csv (same shape as the
        CGC original, with in_cosmic/cosmic_hits_cumulative renamed to
        in_intogen_any/intogen_hits_cumulative)
"""
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
INTOGEN_DRIVERS = ROOT / "data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"
OUT = ROOT / "data/derived/intogen_rescoring/precision_at_k_detail_intogen.csv"


def load_genes_any() -> set:
    genes = set()
    with open(INTOGEN_DRIVERS, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genes.add(row["SYMBOL"].upper())
    return genes


def main():
    genes_any = load_genes_any()
    d = pd.read_csv(ROOT / "data/canonical/precision_at_k_detail.csv")
    d = d.drop(columns=["in_cosmic", "cosmic_hits_cumulative", "precision_at_k"])
    d["in_intogen_any"] = d["gene"].str.upper().isin(genes_any)
    d = d.sort_values(["cancer", "rank"])
    d["intogen_hits_cumulative"] = d.groupby("cancer")["in_intogen_any"].cumsum()
    d["precision_at_k"] = d["intogen_hits_cumulative"] / d["rank"]

    OUT.parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(OUT, index=False)
    print(f"Wrote {OUT} ({len(d)} rows)")

    for k in (5, 10, 15, 20, 30, 40, 50):
        sub = d[d["rank"] == k]
        print(f"  k={k}: mean precision = {sub['precision_at_k'].mean():.4f}")


if __name__ == "__main__":
    main()
