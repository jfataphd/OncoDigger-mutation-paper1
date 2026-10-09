"""CancerMine top-10 gene ranking per cancer, matched to the same 19 cancer
types used throughout this paper (see docs/scope.md, docs/head_to_head_plan.md).

CancerMine (Lever et al. 2019a) classifies gene-cancer pairs as Driver,
Oncogene or Tumor_Suppressor from PubMed/PMC text mining, with a citation
count per (cancer, gene, role) triple aggregated across all papers reporting
that triple. It does not group cancer subtypes the way OncoDigger's 19
"Common cancers" corpora do, so subtypes are aggregated per cancer using the
substring mapping below. Within a cancer, a gene's aggregate score is the sum
of citation_count across every matching disease term and every role (Driver +
Oncogene + Tumor_Suppressor), and the top-10 genes by that score are the
ranking scored against IntOGen by intogen_full_rescoring.py for every
CancerMine number this manuscript reports.

Retired 2026-10-08: this script originally also scored its own top-10 against
COSMIC Cancer Gene Census Tier 1 (`cosmic_tier1` column, `cosmic_census_2026-05-23.csv`,
not included in this repository). That COSMIC-based score is not used by
manuscript1.1's IntOGen-calibrated analysis and has been removed; only the
ranking itself (gene, rank, cancer, citation_count_sum) is still used.

Source data: cancermine_collated.tsv (Zenodo record 7689627), fetched to
data/external/cancermine/cancermine_collated.tsv.
"""

import csv
import re
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CANCERMINE_TSV = REPO / "data" / "external" / "cancermine" / "cancermine_collated.tsv"
OUT_RANKINGS = REPO / "data" / "derived" / "cancermine_top10.csv"

# Substring match against CancerMine's `cancer_normalized` field (lowercased).
# Broad ("include all subtype") aggregation, matching the COSMIC Browser
# tissue-matching approach in Section 2.4 of the source manuscript.
CANCER_MAP = {
    "Breast Cancer": ["breast"],
    "Lung Cancer": ["lung"],
    "Colorectal Cancer": ["colorectal", "colon", "rectal", "rectum"],
    "Prostate Cancer": ["prostate"],
    "Melanoma": ["melanoma"],
    "Bladder Cancer": ["bladder"],
    "Kidney Cancer": ["kidney", "renal"],
    "Pancreatic Cancer": ["pancrea"],
    "Liver Cancer": ["liver", "hepato"],
    "Stomach Cancer": ["stomach", "gastric"],
    "Esophageal Cancer": ["esophag", "gastroesophageal"],
    "Ovarian Cancer": ["ovar"],
    "Endometrial Cancer": ["endometri", "uterine corpus"],
    "Cervical Cancer": ["cervi"],
    "Thyroid Cancer": ["thyroid"],
    "Brain Cancer": ["brain", "glioma", "glioblastoma", "astrocytoma", "ependymoma",
                      "medulloblastoma", "central nervous system"],
    "Leukemia": ["leukemia"],
    "Lymphoma": ["lymphoma"],
    "Myeloma": ["myeloma"],
}

# Substring-matching false positives found by manual audit (docs/head_to_head_plan.md):
# these disease terms are misassigned by a naive substring match and are excluded from
# the cancer whose substring they spuriously trigger. Terms combining two disease
# categories that are genuinely both (e.g. "adult T-cell leukemia/lymphoma") are left
# double-counted deliberately, since that reflects real diagnostic overlap, not a bug.
TERM_EXCLUSIONS = {
    "Endometrial Cancer": ["endometrioid ovary carcinoma"],  # an ovarian-cancer subtype
    "Liver Cancer": ["liver lymphoma"],  # a lymphoma, not a liver carcinoma
    "Ovarian Cancer": ["ovarian lymphoma"],  # a lymphoma, not an ovarian carcinoma
    # "central nervous system lymphoma" contains "central nervous system" (now matched
    # for Brain Cancer, added to recover ependymoma/medulloblastoma evidence found missing
    # by the red-team audit's cross-resource crosswalk check) but is itself a lymphoma,
    # not a glial/embryonal brain tumor; it is already correctly counted toward Lymphoma
    # via that category's own "lymphoma" substring, so it is excluded here only.
    "Brain Cancer": ["central nervous system lymphoma"],
}

# CancerMine entity-normalization artifact confirmed by manual audit of
# cancermine_sentences.tsv (docs/head_to_head_plan.md): every one of the 1,329
# supporting sentences for gene_entrez_id 5729 (PTGDR) actually mentions the literal
# string "AS1" as part of an unrelated lncRNA symbol (FEZF1-AS1, VIM-AS1, FOXD2-AS1,
# LMCD1-AS1, ...), not PTGDR. This is excluded as a confirmed false positive, not a
# judgment call.
EXCLUDED_GENE_ENTREZ_IDS = {"5729"}  # PTGDR


def load_cancermine(path: Path):
    rows = []
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows.append(row)
    return rows


def aggregate_top10(rows, cancer, terms):
    pattern = re.compile("|".join(re.escape(t) for t in terms), re.IGNORECASE)
    excluded_terms = set(TERM_EXCLUSIONS.get(cancer, []))
    scores: dict[str, float] = defaultdict(float)
    for row in rows:
        if row["cancer_normalized"] in excluded_terms:
            continue
        if row["gene_entrez_id"] in EXCLUDED_GENE_ENTREZ_IDS:
            continue
        if pattern.search(row["cancer_normalized"]):
            gene = row["gene_normalized"].strip().upper()
            try:
                scores[gene] += float(row["citation_count"])
            except ValueError:
                continue
    ranked = sorted(scores.items(), key=lambda kv: (-kv[1], kv[0]))
    return ranked[:10]


def main():
    if not CANCERMINE_TSV.exists():
        raise SystemExit(
            f"Missing {CANCERMINE_TSV}. Fetch it first:\n"
            "  curl -sL https://zenodo.org/records/7689627/files/cancermine_collated.tsv "
            f"-o {CANCERMINE_TSV}"
        )
    rows = load_cancermine(CANCERMINE_TSV)

    OUT_RANKINGS.parent.mkdir(parents=True, exist_ok=True)
    rank_records = []
    for cancer, terms in CANCER_MAP.items():
        top10 = aggregate_top10(rows, cancer, terms)
        for rank, (gene, score) in enumerate(top10, start=1):
            rank_records.append(
                {"cancer": cancer, "rank": rank, "gene": gene, "citation_count_sum": score}
            )

    with open(OUT_RANKINGS, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["cancer", "rank", "gene", "citation_count_sum"])
        writer.writeheader()
        writer.writerows(rank_records)

    print(f"Wrote {OUT_RANKINGS}")
    print("Scored against IntOGen by scripts/analysis/intogen_full_rescoring.py")


if __name__ == "__main__":
    main()
