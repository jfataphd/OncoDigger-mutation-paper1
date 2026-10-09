"""Cancer-to-subtype vocabulary mapping (CANCER_TERMS).

Retired 2026-10-08: this module originally also computed a COSMIC Cancer Gene Census
Tier 1-based cancer-specific Precision@10 metric (requiring a top-10 gene's own CGC
`Tumour Types(Somatic)`/`(Germline)` field to name the specific cancer being scored, not
just membership in Tier 1 anywhere). That metric is not used by manuscript1.1's
IntOGen-calibrated analysis and has been removed. The `cosmic_census_2026-05-23.csv`
file it depended on is not included in this repository.

What remains, and is still used: CANCER_TERMS, the curated substring mapping from each
of this study's 19 cancers to the vocabulary terms found in CGC's own tumour-type field
during the original audit. Supplementary Methods S4 (PubTator disease mapping) reuses
this same curated subtype terminology, converted to full disease-name phrases for
PubTator's entity-autocomplete API (`pubtator_benchmark_v2_subtype_aggregated.py`).
"""

import re

# Curated substring mapping, built from manual inspection of the actual CGC
# Tumour Types(Somatic)/(Germline) vocabulary for the 71 genes appearing in
# some cancer's OncoDigger top-10 (docs/head_to_head_plan.md audit note).
# Matching is case-insensitive substring search against each gene's combined
# tumour-type string.
CANCER_TERMS = {
    "Breast Cancer": ["breast"],
    "Lung Cancer": ["nsclc", "lung", "sclc"],
    "Colorectal Cancer": ["colorectal", "crc", "colon", "large intestine"],
    "Prostate Cancer": ["prostate"],
    "Melanoma": ["melanoma", "melanocytic", "spitzoid"],
    "Bladder Cancer": ["bladder", "urothelial"],
    "Kidney Cancer": [r"\brcc\b", "renal", "wilms", "kidney"],
    "Pancreatic Cancer": ["pancrea"],
    "Liver Cancer": ["hepat", "cholangiocarcinoma", "liver"],
    "Stomach Cancer": ["gastric", "stomach"],
    "Esophageal Cancer": ["esophag", "oesophag"],
    "Ovarian Cancer": ["ovar"],
    "Endometrial Cancer": ["endometri", "uterine"],
    "Cervical Cancer": ["cervi"],
    "Thyroid Cancer": ["thyroid"],
    "Brain Cancer": ["glioma", "glioblastoma", r"\bcns\b", "medulloblastoma",
                      "astrocytoma", "central nervous system", "dipg", r"\bgbm\b",
                      "ependymoma"],
    "Leukemia": [r"\baml\b", r"\ball\b", r"\bcml\b", r"\bcll\b", r"\bapl\b", "mds",
                  "myelodysplastic syndrome", "myelodysplastic syndromes",
                  r"\bcmml\b", r"\bjmml\b", r"\bacml\b", "leukaemia", "leukemia",
                  r"\bt-cll\b", r"\bt-pll\b", r"\bcnl\b"],
    "Lymphoma": ["lymphoma", r"\bnhl\b", r"\bdlbcl\b", r"\bmalt\b", r"\bwm\b",
                  "burkitt", "hodgkin", r"\balcl\b", r"\bmcl\b", r"\bpmbl\b", r"\bsmzl\b"],
    "Myeloma": [r"\bmm\b", "myeloma"],
}

# Two compound CGC tumour-type strings are genuinely ambiguous across this
# study's simplified 19-category scheme and were reviewed manually rather
# than left to regex collision (found while extending this mapping to the
# full 592-gene census for the retired CGC-Recall analysis):
#   "primary central nervous system lymphoma" -- clinically and pathologically
#   a lymphoma (managed by haematology-oncology), not a glioma/astrocytoma-type
#   primary brain tumour, despite its anatomical location; mapped to Lymphoma
#   only, excluded from Brain Cancer.
#   "primary central nervous system melanocytic neoplasms" -- a rare
#   leptomeningeal melanocytic tumour, distinct from both cutaneous/uveal
#   melanoma and from glioma-type brain cancer; excluded from both categories
#   as not confidently resolvable to either.
AMBIGUOUS_TERM_OVERRIDES = {
    "primary central nervous system lymphoma": {"Lymphoma"},
    "primary central nervous system melanocytic neoplasms": set(),
}


def matching_cancers(tumour_type_string: str) -> set:
    """All of the 19 study cancers this tumour-type string matches. The two
    ambiguous compound terms in AMBIGUOUS_TERM_OVERRIDES are masked out before
    normal regex matching (so their own substrings, e.g. "melanocytic" or
    "central nervous system", cannot trigger a default match) and their
    manually reviewed category overrides are unioned back in afterward --
    this must not affect matching of any *other* tumour-type term present for
    the same gene (e.g. GNA11 and GNAQ separately carry "uveal melanoma",
    which must still match Melanoma even though their CNS-melanocytic term is
    excluded)."""
    masked = tumour_type_string
    override_cancers = set()
    for phrase, override in AMBIGUOUS_TERM_OVERRIDES.items():
        if phrase in masked:
            masked = masked.replace(phrase, "")
            override_cancers |= override
    matched = {cancer for cancer, terms in CANCER_TERMS.items()
               if any(re.search(t, masked) for t in terms)}
    return matched | override_cancers


def cancer_matches(cancer: str, tumour_type_string: str) -> bool:
    return cancer in matching_cancers(tumour_type_string)
