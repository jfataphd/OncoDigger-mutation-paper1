"""
interaction_all19.py -- the retrieval x vocabulary 2x2 factorial interaction (Results 5.8,
Table 2, Abstract: "+45.3 percentage points"), for all 19 cancers.

Reconstructed 2026-10-08: this script's output (interaction_all19_final.csv) had been
committed without the script that produced it. Reconstructed by recognizing the file is a
straightforward extension, to all 19 cancers, of interaction_bootstrap_5cancer.py's own
5-cancer diagnostic logic and its two already-committed, already-scripted input files below;
re-running this reconstruction reproduced the committed interaction_all19_final.csv exactly,
all 152 (19 cancers x 8 cell/metric values) data points, before this script was added.

Inputs:
  data/derived/intogen_rescoring/true_raw_flat_count_2x2_detail.csv       (true_raw_flat_count_2x2.py)
  data/derived/intogen_rescoring/retrieval_vocab_2x2_faithful_detail.csv  (retrieval_vocab_2x2_faithful.py)

The four cells of the factorial:
  cell1 (full corpus, true-raw vocabulary)   <- true_raw_flat_count_2x2_detail.csv, cell1_fullcorpus_true_raw
  cell2 (full corpus, cleaned vocabulary)    <- retrieval_vocab_2x2_faithful_detail.csv, cell2_fullcorpus_clean
  cell3 (BM25 pool, true-raw vocabulary)     <- true_raw_flat_count_2x2_detail.csv, cell3_pool_true_raw
  cell4 (BM25 pool, cleaned vocabulary)      <- retrieval_vocab_2x2_faithful_detail.csv, cell4_pool_clean

Interaction = cell4 - cell2 - cell3 + cell1 (observed combined effect minus the two main
effects' additive prediction), per cancer and pooled, against IntOGen any-cancer and exact
cancer-type-matched driver membership.

Output: data/derived/intogen_rescoring/interaction_all19_final.csv

Run: python -X utf8 scripts/analysis/interaction_all19.py
"""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DER = ROOT / "data" / "derived" / "intogen_rescoring"
INTOGEN_DRIVERS = ROOT / "data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"

INTOGEN_CANCER_MAP = {
    "Breast Cancer": {"BRCA"}, "Lung Cancer": {"LUAD", "LUSC", "NSCLC", "SCLC"},
    "Colorectal Cancer": {"COAD", "READ", "COADREAD"}, "Prostate Cancer": {"PRAD"},
    "Melanoma": {"MEL", "SKCM", "UM"}, "Bladder Cancer": {"BLCA", "UTUC"},
    "Kidney Cancer": {"CCRCC", "CHRCC", "PRCC", "RCC", "WT"}, "Pancreatic Cancer": {"PAAD", "PANET"},
    "Liver Cancer": {"HCC", "CHOL", "LIHB"}, "Stomach Cancer": {"STAD"},
    "Esophageal Cancer": {"ESCA", "ESCC"}, "Ovarian Cancer": {"OVT"},
    "Endometrial Cancer": {"UCEC", "UCS"}, "Cervical Cancer": {"CESC", "CEAD"},
    "Thyroid Cancer": {"WDTC"}, "Brain Cancer": {"GB", "GBM", "HGGNOS", "LGGNOS", "PAST", "MBL", "EPM", "ATRT"},
    "Leukemia": {"ALL", "AML", "CLLSLL", "CML", "MDS"}, "Lymphoma": {"BL", "DLBCLNOS", "NHL", "MLYM"},
    "Myeloma": {"PCM"},
}
CANCERS = list(INTOGEN_CANCER_MAP.keys())


def load_intogen():
    genes_any, exact_pairs = set(), set()
    with open(INTOGEN_DRIVERS, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            sym = row["SYMBOL"].strip().upper()
            genes_any.add(sym)
            code = row["CANCER_TYPE"].strip()
            for cancer, codes in INTOGEN_CANCER_MAP.items():
                if code in codes:
                    exact_pairs.add((sym, cancer))
    return genes_any, exact_pairs


def top10_from_detail(path: Path, cell_name: str, cancers: list[str]) -> dict[str, list[str]]:
    rows = list(csv.DictReader(open(path, encoding="utf-8")))
    out = {c: [] for c in cancers}
    for r in rows:
        if r["cell"] != cell_name or r["cancer"] not in out:
            continue
        if int(r["rank"]) <= 10:
            out[r["cancer"]].append((int(r["rank"]), r["gene"]))
    for c in out:
        out[c] = [g for _, g in sorted(out[c])]
    return out


def score(top10_by_cancer: dict, genes_any: set, exact_pairs: set) -> tuple[dict, dict]:
    any_hits, exact_hits = {}, {}
    for cancer, genes in top10_by_cancer.items():
        any_hits[cancer] = sum(1 for g in genes if g.upper() in genes_any)
        exact_hits[cancer] = sum(1 for g in genes if (g.upper(), cancer) in exact_pairs)
    return any_hits, exact_hits


def main() -> None:
    genes_any, exact_pairs = load_intogen()

    true_raw_detail = DER / "true_raw_flat_count_2x2_detail.csv"
    cleaned_detail = DER / "retrieval_vocab_2x2_faithful_detail.csv"

    cell1_top10 = top10_from_detail(true_raw_detail, "cell1_fullcorpus_true_raw", CANCERS)
    cell3_top10 = top10_from_detail(true_raw_detail, "cell3_pool_true_raw", CANCERS)
    cell2_top10 = top10_from_detail(cleaned_detail, "cell2_fullcorpus_clean", CANCERS)
    cell4_top10 = top10_from_detail(cleaned_detail, "cell4_pool_clean", CANCERS)

    cell1_any, cell1_exact = score(cell1_top10, genes_any, exact_pairs)
    cell2_any, cell2_exact = score(cell2_top10, genes_any, exact_pairs)
    cell3_any, cell3_exact = score(cell3_top10, genes_any, exact_pairs)
    cell4_any, cell4_exact = score(cell4_top10, genes_any, exact_pairs)

    rows = []
    for c in CANCERS:
        any_interaction = cell4_any[c] - cell2_any[c] - cell3_any[c] + cell1_any[c]
        exact_interaction = cell4_exact[c] - cell2_exact[c] - cell3_exact[c] + cell1_exact[c]
        rows.append({
            "cancer": c,
            "cell1_any": cell1_any[c], "cell2_any": cell2_any[c], "cell3_any": cell3_any[c], "cell4_any": cell4_any[c],
            "any_interaction": any_interaction,
            "cell1_exact": cell1_exact[c], "cell2_exact": cell2_exact[c], "cell3_exact": cell3_exact[c], "cell4_exact": cell4_exact[c],
            "exact_interaction": exact_interaction,
        })

    out = DER / "interaction_all19_final.csv"
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["cancer", "cell1_any", "cell2_any", "cell3_any", "cell4_any", "any_interaction",
                                           "cell1_exact", "cell2_exact", "cell3_exact", "cell4_exact", "exact_interaction"])
        w.writeheader()
        w.writerows(rows)

    any_total = sum(cell4_any.values()) - sum(cell2_any.values()) - sum(cell3_any.values()) + sum(cell1_any.values())
    exact_total = sum(cell4_exact.values()) - sum(cell2_exact.values()) - sum(cell3_exact.values()) + sum(cell1_exact.values())
    print(f"Wrote {out}")
    print(f"cell1={sum(cell1_any.values())} cell2={sum(cell2_any.values())} "
          f"cell3={sum(cell3_any.values())} cell4={sum(cell4_any.values())} (any, /190 each)")
    print(f"Any-cancer interaction: {100*any_total/190:+.1f}pp")
    print(f"Exact-concordance interaction: {100*exact_total/190:+.1f}pp")
    print("Positive in all 19 cancers:" if all(cell4_any[c]-cell2_any[c]-cell3_any[c]+cell1_any[c] > 0 for c in CANCERS) else "NOT positive in all 19:",
          [c for c in CANCERS if cell4_any[c]-cell2_any[c]-cell3_any[c]+cell1_any[c] <= 0] or "none")


if __name__ == "__main__":
    main()
