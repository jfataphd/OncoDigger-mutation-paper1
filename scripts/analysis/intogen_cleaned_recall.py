"""
intogen_cleaned_recall.py -- manuscript1.1 (IntOGen-primary draft) ONLY. Recomputes OncoDigger's
recall of IntOGen's driver catalogue (Results 5.5, Supplementary Results S6) after removing two
classes of pairs that do not reflect a ranking failure:

  (1) Pairs with zero mentions anywhere in the full corpus (Supplementary Results S10): these are
      structurally undetectable regardless of ranking quality, not a recall failure.
  (2) Pairs backed by only one of IntOGen's seven driver-detection methods (Supplementary Results
      S9): these are overwhelmingly CGC-rescue cases (a CGC gene can enter IntOGen's compendium via
      a lenient CGC-restricted threshold with zero independently significant methods, verified
      against the public `intogen-plus` v2024 pipeline source: `combination/intogen_combination
      /stouffer_script.py`'s `select_significant_bidders(q<=0.1)` and `core/intogen_core
      /postprocess/drivers/vetting.py`'s `get_drivers()`), not independently corroborated IntOGen
      evidence. A non-CGC gene needs agreement from >=2 such methods to enter the compendium at
      all, so requiring >=2 here selects IntOGen's own higher-confidence tier rather than an
      arbitrary cut.

Both quantities (full-corpus postings, method count) were already computed for Supplementary
Results S9/S10 and are reused unchanged here; this script only recombines them against the
literature-weight audit's per-pair rank data to report recall at standard depths under four
nested filters.

Inputs: data/derived/intogen_rescoring/intogen_literature_weight_audit.csv,
        intogen_full_corpus_presence_check.csv, intogen_confidence_breakdown_detail.csv
Output: data/derived/intogen_rescoring/intogen_cleaned_recall_summary.csv

Run: python -X utf8 scripts/analysis/intogen_cleaned_recall.py
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DER = ROOT / "data" / "derived" / "intogen_rescoring"
DEPTHS = [10, 20, 50]


def recall_row(label: str, df: pd.DataFrame) -> dict:
    n = len(df)
    row = {"filter": label, "n_pairs": n,
           "any_rank_pct": 100 * df["engine_rank"].notna().sum() / n}
    for d in DEPTHS:
        row[f"rank{d}_pct"] = 100 * (df["engine_rank"] <= d).sum() / n
    return row


def main():
    audit = pd.read_csv(DER / "intogen_literature_weight_audit.csv")
    presence = pd.read_csv(DER / "intogen_full_corpus_presence_check.csv")
    conf = pd.read_csv(DER / "intogen_confidence_breakdown_detail.csv")

    absent = presence[presence["full_corpus_postings"] == 0]
    absent_pairs = set(zip(absent["cancer"], absent["gene"]))
    audit["is_absent"] = audit.apply(lambda r: (r["cancer"], r["gene"]) in absent_pairs, axis=1)
    audit = audit.merge(conf[["cancer", "gene", "n_methods"]], on=["cancer", "gene"], how="left")

    rows = [
        recall_row("No filter (all 1,740 pairs)", audit),
        recall_row("Remove zero-mention pairs", audit[~audit["is_absent"]]),
        recall_row("Remove zero-mention + require >=2 methods",
                   audit[(~audit["is_absent"]) & (audit["n_methods"] >= 2)]),
    ]
    summary = pd.DataFrame(rows)
    summary.to_csv(DER / "intogen_cleaned_recall_summary.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
