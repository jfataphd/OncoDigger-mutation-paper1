"""
intogen_figure4_headtohead_data.py -- manuscript1.1 (IntOGen-primary draft) ONLY. Produces
the IntOGen-scored data Figure 4 (head-to-head benchmark) needs: rescores the original
(single-parent-MeSH, pre-correction) PubTator ranking against IntOGen for the Panel C
before/after comparison, and computes paired t-test p-values (OncoDigger vs CancerMine,
DISEASES, PubTator) for Panel A, reusing the already-computed any-cancer per-cancer summary
files from intogen_full_rescoring.py and intogen_diseases_pubtator_rescoring.py.

Output: data/derived/intogen_rescoring/pubtator_original_intogen_detail.csv
        data/derived/intogen_rescoring/figure4_headtohead_pvalues.csv
"""
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
DER = ROOT / "data" / "derived"
OUT_DIR = DER / "intogen_rescoring"
INTOGEN_DRIVERS = ROOT / "data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"


def load_genes_any() -> set:
    genes = set()
    with open(INTOGEN_DRIVERS, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genes.add(row["SYMBOL"].upper())
    return genes


def main():
    genes_any = load_genes_any()

    # PubTator original (single-parent-MeSH, pre-correction) rescored against IntOGen
    df = pd.read_csv(DER / "pubtator_top10.csv")
    df["intogen_any"] = df["gene"].str.upper().isin(genes_any)
    df.to_csv(OUT_DIR / "pubtator_original_intogen_detail.csv", index=False)
    hits, n = int(df["intogen_any"].sum()), len(df)
    print(f"PubTator original (pre-correction) vs IntOGen: {hits}/{n} = {hits/n:.1%}")

    # Paired t-tests, OncoDigger vs each tool, any-cancer Precision@10 per cancer
    od = pd.read_csv(OUT_DIR / "oncodigger_intogen_summary.csv").set_index("cancer")["p10_intogen_any"]
    rows = []
    for name, path in [("cancermine", OUT_DIR / "cancermine_intogen_summary.csv"),
                        ("diseases", OUT_DIR / "diseases_intogen_summary.csv"),
                        ("pubtator_v2", OUT_DIR / "pubtator_v2_intogen_summary.csv")]:
        s = pd.read_csv(path).set_index("cancer")
        col = "p10_intogen_any" if "p10_intogen_any" in s.columns else "intogen_any"
        s = s[col]
        cancers = sorted(set(od.index) & set(s.index))
        a, b = od.loc[cancers].values, s.loc[cancers].values
        tt, pt = stats.ttest_rel(a, b)
        rows.append({"tool": name, "od_mean": a.mean(), "tool_mean": b.mean(), "p_ttest": pt, "n": len(cancers)})
        print(f"  OncoDigger vs {name}: OD={a.mean():.1%} tool={b.mean():.1%} p={pt:.3g}")

    pd.DataFrame(rows).to_csv(OUT_DIR / "figure4_headtohead_pvalues.csv", index=False)


if __name__ == "__main__":
    main()
