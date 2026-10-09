"""DISEASES 2.0 top-10 gene ranking per cancer, matched to the same 19 cancer
types used throughout this paper (Methods 2.11, 2.12; Results 3.3, 3.4; see
docs/head_to_head_plan.md).

DISEASES (Pletscher-Frankild et al., 2015; Grissa et al., 2022) provides a
z-score per (gene, DOID disease term) pair from PubMed/PMC text mining. The
"filtered" bulk file keeps only each gene's own top-scoring diseases (mostly
broad ontology-parent terms), which is unusable for ranking genes *within* a
given cancer, so the full bulk file was streamed and grep-filtered to cancer-
related rows only (`data/external/diseases/textmining_full_cancer_subset.tsv`;
not committed — see docs/head_to_head_plan.md for how to regenerate it) before
this script aggregates it into per-cancer top-10 gene lists exactly as
`cancermine_benchmark.py` does for CancerMine.

Column layout of the (headerless) source rows: ENSP id, gene symbol, DOID,
disease name, z-score, confidence (z-score/2, rounded), URL.

Retired 2026-10-08: this script originally also scored its own top-10 against COSMIC
Cancer Gene Census Tier 1 (`cosmic_tier1` column, `cosmic_census_2026-05-23.csv`, not
included in this repository). That COSMIC-based score is not used by manuscript1.1's
IntOGen-calibrated analysis and has been removed; only the ranking itself (gene, rank,
cancer, max_z_score) is still used, scored against IntOGen by
intogen_diseases_pubtator_rescoring.py.
"""

import csv
import re
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DISEASES_SUBSET_TSV = (
    REPO / "data" / "external" / "diseases" / "textmining_full_cancer_subset.tsv"
)
OUT_RANKINGS = REPO / "data" / "derived" / "diseases_top10.csv"

# Same organ/site substrings as scripts/analysis/cancermine_benchmark.py's
# CANCER_MAP, with one fix: "bladder" is given a negative lookbehind so it
# does not also match "gallbladder" (found by audit, see docs/head_to_head_plan.md).
CANCER_MAP = {
    "Breast Cancer": ["breast"],
    "Lung Cancer": ["lung"],
    "Colorectal Cancer": ["colorectal", "colon", "rectal", "rectum"],
    "Prostate Cancer": ["prostate"],
    "Melanoma": ["melanoma"],
    "Bladder Cancer": [r"(?<!gall)bladder"],
    "Kidney Cancer": ["kidney", "renal"],
    "Pancreatic Cancer": ["pancrea"],
    "Liver Cancer": ["liver", "hepato"],
    "Stomach Cancer": ["stomach", "gastric"],
    "Esophageal Cancer": ["esophag", "gastroesophageal"],
    "Ovarian Cancer": ["ovar"],
    "Endometrial Cancer": ["endometri", "uterine corpus"],
    "Cervical Cancer": ["cervi"],
    "Thyroid Cancer": ["thyroid"],
    "Brain Cancer": ["brain", "glioma", "glioblastoma", "astrocytoma"],
    "Leukemia": ["leukemia"],
    "Lymphoma": ["lymphoma"],
    "Myeloma": ["myeloma"],
}

# Explicit reassignment for terms confirmed by DOID definition lookup (OLS,
# docs/head_to_head_plan.md) to be a subtype of one specific cancer, not the
# one a naive substring match would suggest:
#   DOID:5828 "endometrioid ovary carcinoma" -- an ovary adenocarcinoma
#   DOID:6212 "ovarian endometrial cancer"    -- an ovary epithelial cancer
#   DOID:5830 "cervical endometrioid adenocarcinoma" -- a cervical adenocarcinoma
EXPLICIT_TERM_CANCERS = {
    "Endometrioid ovary carcinoma": {"Ovarian Cancer"},
    "Ovarian endometrial cancer": {"Ovarian Cancer"},
    "Cervical endometrioid adenocarcinoma": {"Cervical Cancer"},
}

# "Hereditary breast ovarian cancer syndrome" is deliberately left assigned to
# both Breast and Ovarian Cancer: this is real dual-organ hereditary biology,
# consistent with how the source manuscript treats BRCA1 (Category B in both
# Breast and Prostate/Pancreatic Cancer, Category A in Ovarian Cancer;
# Section 3.3), not a mapping bug.
DUAL_ORGAN_TERMS = {"Hereditary breast ovarian cancer syndrome"}


def classify_term(term: str) -> set[str]:
    """Return the set of cancer categories a DISEASES disease name belongs to,
    applying disease-type overrides before falling back to organ substring
    matching. See docs/head_to_head_plan.md for the audit that derived these
    rules from the 37 cross-matching terms found in this dataset."""
    if term in EXPLICIT_TERM_CANCERS:
        return EXPLICIT_TERM_CANCERS[term]

    lower = term.lower()
    has_lymphoma = "lymphoma" in lower
    has_leukemia = "leukemia" in lower
    has_melanoma = "melanoma" in lower

    if term in DUAL_ORGAN_TERMS:
        pass  # fall through to organ matching below (deliberately dual)
    elif has_lymphoma and has_leukemia:
        return {"Leukemia", "Lymphoma"}
    elif has_lymphoma:
        # A lymphoma named after the organ it arose in/spread to (e.g. "Liver
        # lymphoma", "Gastric lymphoma") is a lymphoma, not a primary organ
        # carcinoma.
        return {"Lymphoma"}
    elif has_melanoma:
        # Likewise "Cervix melanoma", "Gallbladder melanoma", etc. are
        # melanoma, not the named organ's primary carcinoma.
        return {"Melanoma"}

    matched = set()
    for cancer, terms in CANCER_MAP.items():
        pattern = re.compile("|".join(terms), re.IGNORECASE)
        if pattern.search(term):
            matched.add(cancer)
    return matched


def load_rows(path: Path):
    rows = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            rows.append(
                {
                    "ensp": parts[0],
                    "gene": parts[1],
                    "doid": parts[2],
                    "disease_name": parts[3],
                    "z_score": parts[4],
                }
            )
    return rows


def build_term_to_cancers(rows) -> dict[str, set[str]]:
    all_terms = {r["disease_name"] for r in rows}
    return {term: classify_term(term) for term in all_terms}


def aggregate_top10(rows, cancer, term_to_cancers):
    # Max z-score across matched disease terms, not sum: z-scores are already
    # a normalized statistic per (gene, term) pair, and cancers here often
    # have several near-synonymous DOID subtypes (e.g. "Breast cancer",
    # "Sporadic breast cancer", "Bilateral breast cancer"); summing would let
    # a gene's score scale with how many near-duplicate subtypes happen to
    # exist for that cancer rather than with its evidence strength, unlike
    # CancerMine's citation_count, which is an additive count of papers where
    # summing subtypes is the correct aggregation (docs/head_to_head_plan.md).
    best: dict[str, float] = defaultdict(float)
    for row in rows:
        if cancer not in term_to_cancers.get(row["disease_name"], set()):
            continue
        gene = row["gene"].strip().upper()
        try:
            z = float(row["z_score"])
        except ValueError:
            continue
        if z > best[gene]:
            best[gene] = z
    ranked = sorted(best.items(), key=lambda kv: (-kv[1], kv[0]))
    return ranked[:10]


def main():
    if not DISEASES_SUBSET_TSV.exists():
        raise SystemExit(
            f"Missing {DISEASES_SUBSET_TSV}. Regenerate it with:\n"
            '  curl -sL https://download.jensenlab.org/human_disease_textmining_full.tsv | '
            'grep -iE "cancer|leukemia|lymphoma|melanoma|myeloma|carcinoma|glioma|astrocytoma|'
            f'hepatocellular|adenocarcinoma" > {DISEASES_SUBSET_TSV}'
        )
    rows = load_rows(DISEASES_SUBSET_TSV)
    term_to_cancers = build_term_to_cancers(rows)

    OUT_RANKINGS.parent.mkdir(parents=True, exist_ok=True)
    rank_records = []
    for cancer in CANCER_MAP:
        top10 = aggregate_top10(rows, cancer, term_to_cancers)
        for rank, (gene, score) in enumerate(top10, start=1):
            rank_records.append({"cancer": cancer, "rank": rank, "gene": gene, "max_z_score": score})

    with open(OUT_RANKINGS, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["cancer", "rank", "gene", "max_z_score"])
        writer.writeheader()
        writer.writerows(rank_records)

    print(f"Wrote {OUT_RANKINGS}")
    print("Scored against IntOGen by scripts/analysis/intogen_diseases_pubtator_rescoring.py")


if __name__ == "__main__":
    main()
