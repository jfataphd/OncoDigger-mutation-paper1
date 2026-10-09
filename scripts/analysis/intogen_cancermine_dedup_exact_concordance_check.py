"""
intogen_cancermine_dedup_exact_concordance_check.py -- manuscript1.1 (IntOGen-primary
draft) ONLY. Mirrors cancermine_dedup_exact_concordance_check.py exactly except both the
any-cancer hit and the exact tumour-type concordance are graded against IntOGen's
conservative cancer mapping instead of COSMIC CGC. Follow-up to item 7's double-counting
check (intogen_cancermine_double_counting_check.py): confirms pooled gene-level
Precision@10 is identical whether CancerMine's top-10 is built by summed citation count
or by deduplicated unique-PMID count, and additionally checks whether the stricter exact
IntOGen concordance metric is also unaffected by the same choice, even though top-10
membership does change in 4/19 cancers under deduplication.

Builds the deduplicated top-10 (identical logic to intogen_cancermine_double_counting_check.py),
scores it against IntOGen's (gene, study-cancer) pairs using the same conservative mapping
as every other IntOGen script in this draft, and reports the resulting exact-concordance
count and per-cancer breakdown.

Run: python scripts/analysis/intogen_cancermine_dedup_exact_concordance_check.py
"""
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cancermine_benchmark import CANCER_MAP, TERM_EXCLUSIONS, EXCLUDED_GENE_ENTREZ_IDS

ROOT = Path(__file__).resolve().parents[2]
SENTENCES_TSV = ROOT / "data" / "external" / "cancermine" / "cancermine_sentences.tsv"
INTOGEN_DRIVERS = ROOT / "data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"

INTOGEN_CANCER_MAP = {
    "Breast Cancer": {"BRCA"}, "Lung Cancer": {"LUAD", "LUSC", "NSCLC", "SCLC"},
    "Colorectal Cancer": {"COAD", "READ", "COADREAD"}, "Prostate Cancer": {"PRAD"},
    "Melanoma": {"MEL", "SKCM", "UM"}, "Bladder Cancer": {"BLCA", "UTUC"},
    "Kidney Cancer": {"CCRCC", "CHRCC", "PRCC", "RCC", "WT"},
    "Pancreatic Cancer": {"PAAD", "PANET"}, "Liver Cancer": {"HCC", "CHOL", "LIHB"},
    "Stomach Cancer": {"STAD"}, "Esophageal Cancer": {"ESCA", "ESCC"},
    "Ovarian Cancer": {"OVT"}, "Endometrial Cancer": {"UCEC", "UCS"},
    "Cervical Cancer": {"CESC", "CEAD"}, "Thyroid Cancer": {"WDTC"},
    "Brain Cancer": {"GB", "GBM", "HGGNOS", "LGGNOS", "PAST", "MBL", "EPM", "ATRT"},
    "Leukemia": {"ALL", "AML", "CLLSLL", "CML", "MDS"},
    "Lymphoma": {"BL", "DLBCLNOS", "NHL", "MLYM"}, "Myeloma": {"PCM"},
}


def load_intogen():
    genes_any, pairs = set(), set()
    with open(INTOGEN_DRIVERS, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genes_any.add(row["SYMBOL"].upper())
            for cancer, codes in INTOGEN_CANCER_MAP.items():
                if row["CANCER_TYPE"] in codes:
                    pairs.add((row["SYMBOL"].upper(), cancer))
    return genes_any, pairs


def main() -> None:
    gene_pmids: dict[tuple[str, str], set] = {}
    with SENTENCES_TSV.open(encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["gene_entrez_id"] in EXCLUDED_GENE_ENTREZ_IDS:
                continue
            term = row["cancer_normalized"]
            for cancer, terms in CANCER_MAP.items():
                if term in set(TERM_EXCLUSIONS.get(cancer, [])):
                    continue
                if re.search("|".join(re.escape(t) for t in terms), term, re.IGNORECASE):
                    gene = row["gene_normalized"].strip().upper()
                    gene_pmids.setdefault((cancer, gene), set()).add(row["pmid"])

    genes_any, pairs = load_intogen()
    by_cancer = defaultdict(list)
    for (cancer, gene), pmids in gene_pmids.items():
        by_cancer[cancer].append((gene, len(pmids)))

    exact_total = 0
    any_total = 0
    rows = []
    for cancer in CANCER_MAP:
        ranked = sorted(by_cancer[cancer], key=lambda x: (-x[1], x[0]))[:10]
        for gene, _n in ranked:
            is_any = gene in genes_any
            exact = (gene, cancer) in pairs
            if is_any:
                any_total += 1
            if exact:
                exact_total += 1
            rows.append((cancer, gene, is_any, exact))

    print(f"Deduplicated-ranking IntOGen-any Precision@10: {any_total}/190 = {any_total/190:.1%}")
    print(f"Deduplicated-ranking exact IntOGen concordance: {exact_total}/190 = {exact_total/190:.1%}")

    out = ROOT / "data" / "derived" / "intogen_rescoring" / "cancermine_dedup_exact_concordance_detail.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["cancer", "gene", "intogen_any", "exact_concordance"])
        w.writerows(rows)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
