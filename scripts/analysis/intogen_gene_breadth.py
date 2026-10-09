"""
intogen_gene_breadth.py -- for every gene in IntOGen's compendium, the number of distinct
IntOGen CANCER_TYPE codes it is called a driver in (pan-cancer breadth), supporting the
Results 5.x point that the pooled control's inflated score traces to famous pan-cancer genes
such as TP53 and KRAS.

Reconstructed 2026-10-08: this script's output had been committed without a generating
script anywhere in this repo. Trivially reconstructed as a direct count from the IntOGen
compendium file already used throughout this repo; re-running this reconstruction reproduced
the committed intogen_gene_breadth.csv exactly before this script was added.

Input:  data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv
Output: data/derived/intogen_rescoring/intogen_gene_breadth.csv

Run: python -X utf8 scripts/analysis/intogen_gene_breadth.py
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
INTOGEN_DRIVERS = ROOT / "data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"
OUT = ROOT / "data" / "derived" / "intogen_rescoring" / "intogen_gene_breadth.csv"


def main() -> None:
    breadth: dict[str, set[str]] = defaultdict(set)
    with open(INTOGEN_DRIVERS, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            breadth[row["SYMBOL"].strip().upper()].add(row["CANCER_TYPE"].strip())

    rows = sorted(
        ({"SYMBOL": sym, "n_cancer_types_intogen": len(types)} for sym, types in breadth.items()),
        key=lambda r: r["SYMBOL"],
    )
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["SYMBOL", "n_cancer_types_intogen"])
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {OUT}: {len(rows)} genes")


if __name__ == "__main__":
    main()
