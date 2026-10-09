"""
corpus_specific_genes_full_list.py -- explodes corpus_specificity_clean_intogen.csv's
"cancer_specific_genes" column (semicolon-joined) into one (cancer, gene, is_intogen) row per
gene, backing Results 5.x's "104 separation-recovered genes... 91 IntOGen pairs" claim.

Reconstructed 2026-10-08: this script's output had been committed without a generating script
anywhere in this repo. Recognized it as a row-per-gene explosion of
corpus_specificity_clean_intogen.csv's cancer_specific_genes column (intogen_figure3_data.py),
already present and verified, with is_intogen re-derived from IntOGen any-cancer driver
membership (an initial attempt using exact cancer-type-matched membership did not reproduce
the committed file; any-cancer membership did, exactly). Re-running this reconstruction
reproduced the committed corpus_specific_genes_full_list.csv exactly before this script was
added.

Input:  data/derived/intogen_rescoring/corpus_specificity_clean_intogen.csv
Output: data/derived/intogen_rescoring/corpus_specific_genes_full_list.csv

Run: python -X utf8 scripts/analysis/corpus_specific_genes_full_list.py
"""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DER = ROOT / "data" / "derived" / "intogen_rescoring"


INTOGEN_DRIVERS = ROOT / "data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"


def load_genes_any() -> set[str]:
    genes = set()
    with open(INTOGEN_DRIVERS, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genes.add(row["SYMBOL"].strip().upper())
    return genes


def main() -> None:
    genes_any = load_genes_any()
    rows = []
    with open(DER / "corpus_specificity_clean_intogen.csv", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            cancer = row["cancer"]
            genes = [g.strip() for g in row["cancer_specific_genes"].split(";") if g.strip()]
            for gene in genes:
                rows.append({"cancer": cancer, "gene": gene, "is_intogen": gene.upper() in genes_any})

    out = DER / "corpus_specific_genes_full_list.csv"
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["cancer", "gene", "is_intogen"])
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {out}: {len(rows)} genes, {sum(r['is_intogen'] for r in rows)} IntOGen drivers")


if __name__ == "__main__":
    main()
