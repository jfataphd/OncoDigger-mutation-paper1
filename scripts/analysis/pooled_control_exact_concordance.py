"""
pooled_control_exact_concordance.py — builds the (pooled-control gene x cancer) pair table:
the pooled-control's own top-10 genes, each applied against all 19 cancers in turn, for
intogen_full_rescoring.py to score against IntOGen's exact tumour-type-matched driver pairs
(Results 5.1's "pooled_control" exact-concordance number).

A reviewer could reasonably ask whether the pooled top 10's genes are also exact matches for
every cancer they're scored against, since the same 10 genes are applied identically to all
19 cancers -- this builds all 190 (pooled gene, cancer) pairs for that check, rather than only
the subset that happens to also appear in each cancer's own top 10.

Retired 2026-10-08: this script originally also scored each pair against COSMIC Cancer Gene
Census Tier 1 (`tier1`/`exact_match` columns, via cancer_aware_precision.py's retired
load_tumour_types/cancer_matches). That COSMIC-based score is not used by manuscript1.1's
IntOGen-calibrated analysis and has been removed; only the (cancer, gene) pairs themselves
are still used, scored against IntOGen by intogen_full_rescoring.py.

Run: python -X utf8 scripts/analysis/pooled_control_exact_concordance.py
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

CANCERS = [
    "Bladder Cancer", "Brain Cancer", "Breast Cancer", "Cervical Cancer", "Colorectal Cancer",
    "Endometrial Cancer", "Esophageal Cancer", "Kidney Cancer", "Leukemia", "Liver Cancer",
    "Lung Cancer", "Lymphoma", "Melanoma", "Myeloma", "Ovarian Cancer", "Pancreatic Cancer",
    "Prostate Cancer", "Stomach Cancer", "Thyroid Cancer",
]


def main() -> None:
    pooled = pd.read_csv(ROOT / "data" / "derived" / "pooled_control_mutation.csv")
    pooled_top10 = pooled[pooled["rank"] <= 10]["gene"].tolist()

    rows = [{"cancer": c, "gene": g} for c in CANCERS for g in pooled_top10]

    df = pd.DataFrame(rows)
    print(f"Pooled top-10 genes: {pooled_top10}")
    print(f"Total (gene, cancer) pairs: {len(df)}")
    out = ROOT / "data" / "derived" / "pooled_control_exact_concordance.csv"
    df.to_csv(out, index=False)
    print(f"Saved {out.relative_to(ROOT)}")
    print("Scored against IntOGen by scripts/analysis/intogen_full_rescoring.py")


if __name__ == "__main__":
    main()
