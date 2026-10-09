"""
intogen_diseases_pubtator_rescoring.py -- manuscript1.1 (IntOGen-primary draft) ONLY.
Rescores the already-computed DISEASES 2.0 and PubTator 3.0 subtype-aggregated top-10
rankings (diseases_top10.csv, pubtator_v2_top10.csv -- gene lists and ranks are untouched,
no retrieval re-run) against IntOGen's any-cancer driver set instead of COSMIC CGC Tier 1.
Does not modify or depend on diseases_benchmark.py or pubtator_benchmark_v2_subtype_aggregated.py.

Run: python -X utf8 scripts/analysis/intogen_diseases_pubtator_rescoring.py
"""
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
INTOGEN_DRIVERS = ROOT / "data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"
OUT_DIR = ROOT / "data/derived/intogen_rescoring"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def load_genes_any() -> set:
    genes = set()
    with open(INTOGEN_DRIVERS, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genes.add(row["SYMBOL"].upper())
    return genes


def rescore(name: str, path: Path) -> None:
    genes_any = load_genes_any()
    df = pd.read_csv(path)
    df["intogen_any"] = df["gene"].str.upper().isin(genes_any)
    df.to_csv(OUT_DIR / f"{name}_intogen_detail.csv", index=False)
    per_cancer = df.groupby("cancer")["intogen_any"].mean()
    hits, n = int(df["intogen_any"].sum()), len(df)
    print(f"{name}: {hits}/{n} = {hits/n:.1%} against IntOGen (any cancer)")
    per_cancer.to_csv(OUT_DIR / f"{name}_intogen_summary.csv")


def main():
    rescore("diseases", ROOT / "data/derived/diseases_top10.csv")
    rescore("pubtator_v2", ROOT / "data/derived/pubtator_v2_top10.csv")


if __name__ == "__main__":
    main()
