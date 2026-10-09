"""ChatGPT (free/anonymous tier) top-10 gene ranking per cancer, on the same 19 cancer
types used throughout this paper (see docs/head_to_head_plan.md).

Unlike CancerMine, DISEASES, and PubTator, ChatGPT is not a versioned,
bulk-downloadable text-mining dataset -- it is a single, non-reproducible
generative response to one prompt, sent once to the logged-out/anonymous
ChatGPT interface on 2026-09-27 (raw response, prompt, and access conditions
recorded in full in
data/external/chatgpt/chatgpt_free_response_2026-09-27.md). The gene lists
below are transcribed verbatim from that single response and are hardcoded
here (there is no bulk file to parse and no guarantee that re-running the
same prompt today would reproduce them) so that the *analysis* of that fixed
snapshot is at least itself reproducible, even though the underlying
generative query is not.

ChatGPT explicitly grounded its answer in TCGA/cBioPortal somatic mutation
*frequency* (via live web search), not literature-mined driver prominence --
a materially different ground truth than the other three tools, which is
reported as a methods caveat, not elided.

Retired 2026-10-08: this script originally also scored its own top-10 against COSMIC
Cancer Gene Census Tier 1 (`cosmic_tier1` column, `cosmic_census_2026-05-23.csv`, not
included in this repository) and ran a paired significance test against OncoDigger's
own COSMIC-era (not IntOGen) per-cancer Precision@10. Neither is used by manuscript1.1's
IntOGen-calibrated analysis and both have been removed; only the ranking itself (gene,
rank, cancer) is still used, scored against IntOGen by intogen_full_rescoring.py.
"""

import csv
from collections import defaultdict
from itertools import combinations
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
OUT_RANKINGS = REPO / "data" / "derived" / "chatgpt_top10.csv"

# Transcribed verbatim from data/external/chatgpt/chatgpt_free_response_2026-09-27.md
CHATGPT_TOP10 = {
    "Breast Cancer": ["PIK3CA", "TP53", "TTN", "CDH1", "MUC16", "GATA3", "KMT2C", "MAP3K1", "PIK3R1", "AKT1"],
    "Lung Cancer": ["TP53", "TTN", "MUC16", "CSMD3", "KRAS", "LRP1B", "RYR2", "KEAP1", "STK11", "EGFR"],
    "Colorectal Cancer": ["APC", "TP53", "TTN", "KRAS", "PIK3CA", "MUC16", "SYNE1", "FAT4", "RYR2", "OBSCN"],
    "Prostate Cancer": ["TP53", "SPOP", "TTN", "FOXA1", "KMT2D", "KMT2C", "APC", "ATM", "CDK12", "PIK3CA"],
    "Melanoma": ["BRAF", "NRAS", "NF1", "TP53", "TTN", "MUC16", "LRP1B", "CDKN2A", "FAT4", "RYR2"],
    "Bladder Cancer": ["TP53", "TTN", "KDM6A", "ARID1A", "KMT2D", "PIK3CA", "KMT2C", "MUC16", "RB1", "FGFR3"],
    "Kidney Cancer": ["VHL", "PBRM1", "TTN", "SETD2", "BAP1", "KDM5C", "MUC16", "MTOR", "KMT2C", "TP53"],
    "Pancreatic Cancer": ["KRAS", "TP53", "SMAD4", "CDKN2A", "TTN", "MUC16", "RNF43", "GNAS", "KDM6A", "ARID1A"],
    "Liver Cancer": ["TP53", "CTNNB1", "TTN", "MUC16", "ALB", "AXIN1", "ARID1A", "APOB", "FAT4", "KEAP1"],
    "Stomach Cancer": ["TP53", "TTN", "MUC16", "ARID1A", "CDH1", "LRP1B", "PIK3CA", "SYNE1", "FAT4", "KMT2C"],
    "Esophageal Cancer": ["TP53", "TTN", "MUC16", "FAT4", "CSMD3", "SYNE1", "LRP1B", "PIK3CA", "CDKN2A", "KMT2D"],
    "Ovarian Cancer": ["TP53", "TTN", "MUC16", "CSMD3", "BRCA1", "BRCA2", "FAT3", "NF1", "LRP1B", "RYR2"],
    "Endometrial Cancer": ["PTEN", "PIK3CA", "TP53", "ARID1A", "TTN", "PIK3R1", "KMT2D", "CTNNB1", "RPL22", "KRAS"],
    "Cervical Cancer": ["PIK3CA", "KMT2D", "EP300", "FBXW7", "MAPK1", "PTEN", "TP53", "HLA-B", "CASP8", "FAT1"],
    "Thyroid Cancer": ["BRAF", "NRAS", "HRAS", "TTN", "MUC16", "RYR2", "DNAH5", "FAT1", "DICER1", "LRP1B"],
    "Brain Cancer": ["TP53", "IDH1", "PTEN", "ATRX", "EGFR", "NF1", "TTN", "PIK3CA", "CIC", "PIK3R1"],
    "Leukemia": ["NPM1", "DNMT3A", "FLT3", "IDH2", "IDH1", "TET2", "RUNX1", "TP53", "CEBPA", "WT1"],
    "Lymphoma": ["KMT2D", "TP53", "CREBBP", "GNA13", "EZH2", "MYD88", "TNFRSF14", "B2M", "CD79B", "CARD11"],
    "Myeloma": ["KRAS", "NRAS", "TP53", "DIS3", "FAM46C", "BRAF", "CYLD", "TRAF3", "ATM", "RB1"],
}


def main():
    OUT_RANKINGS.parent.mkdir(parents=True, exist_ok=True)
    rank_records = []
    for cancer, genes in CHATGPT_TOP10.items():
        for rank, gene in enumerate(genes, start=1):
            rank_records.append({"cancer": cancer, "rank": rank, "gene": gene})

    with open(OUT_RANKINGS, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["cancer", "rank", "gene"])
        w.writeheader()
        w.writerows(rank_records)

    # Breadth / genericity: how often does each gene appear across the 19 top-10s
    gene_counts = defaultdict(int)
    top10_sets = {c: set(g.upper() for g in genes) for c, genes in CHATGPT_TOP10.items()}
    for genes in top10_sets.values():
        for g in genes:
            gene_counts[g] += 1
    most_common = sorted(gene_counts.items(), key=lambda kv: -kv[1])[:8]

    jaccards = []
    for a, b in combinations(top10_sets, 2):
        inter = len(top10_sets[a] & top10_sets[b])
        union = len(top10_sets[a] | top10_sets[b])
        jaccards.append(inter / union)
    mean_jaccard = sum(jaccards) / len(jaccards)

    print(f"Wrote {OUT_RANKINGS}")
    print("Scored against IntOGen by scripts/analysis/intogen_full_rescoring.py")
    print(f"Mean pairwise Jaccard among ChatGPT top-10s: {mean_jaccard:.3f}")
    print(f"Most broadly recurring genes: {most_common}")


if __name__ == "__main__":
    main()
