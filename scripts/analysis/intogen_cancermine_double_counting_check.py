"""
intogen_cancermine_double_counting_check.py -- manuscript1.1 (IntOGen-primary draft) ONLY.
Mirrors cancermine_double_counting_check.py exactly except hits are graded against
IntOGen's any-cancer driver set instead of COSMIC CGC Tier 1. Addresses red-team audit
item 7: does
cancermine_benchmark.py's aggregation (summing citation_count across matched disease
subtypes AND across the three roles Driver/Oncogene/Tumor_Suppressor) double-count
the same evidence, if the same PMID/sentence supports a gene under more than one role
or subtype for the same cancer?

This does not change CancerMine's own methodology or benchmark it any differently from
how its own interface presents results (citation-count summation is independently
confirmed elsewhere in this project to match CancerMine's own default gene sort) --
it only checks whether that choice happens to double-count evidence, using the
sentence-level file (`cancermine_sentences.tsv`, Zenodo 7689627) to build an
alternative, fully-deduplicated ranking by distinct supporting PMID per (cancer, gene),
and compares the two.

Source data: data/external/cancermine/cancermine_sentences.tsv (83MB, gitignored --
regenerate with:
  curl -s -L -o data/external/cancermine/cancermine_sentences.tsv \
    "https://zenodo.org/api/records/7689627/files/cancermine_sentences.tsv/content"

Run:  python scripts/analysis/cancermine_double_counting_check.py
Output: prints a per-cancer comparison and the pooled Precision@10 under both methods.
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
OLD_TOP10_CSV = ROOT / "data" / "derived" / "cancermine_top10.csv"


def load_intogen_drivers() -> set:
    genes = set()
    with open(INTOGEN_DRIVERS, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genes.add(row["SYMBOL"].upper())
    return genes


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

    intogen_drivers = load_intogen_drivers()
    by_cancer = defaultdict(list)
    for (cancer, gene), pmids in gene_pmids.items():
        by_cancer[cancer].append((gene, len(pmids)))

    old_top10 = defaultdict(list)
    with OLD_TOP10_CSV.open(encoding="utf-8") as f:
        for r in csv.DictReader(f):
            old_top10[r["cancer"]].append(r["gene"])

    old_total, new_total, diff_cancers = 0, 0, 0
    for cancer in CANCER_MAP:
        ranked = sorted(by_cancer[cancer], key=lambda x: (-x[1], x[0]))[:10]
        new_genes = [g for g, _ in ranked]
        old_genes = old_top10[cancer]
        old_hits = sum(1 for g in old_genes if g in intogen_drivers)
        new_hits = sum(1 for g in new_genes if g in intogen_drivers)
        old_total += old_hits
        new_total += new_hits
        if set(new_genes) != set(old_genes):
            diff_cancers += 1
            print(f"{cancer}: citation-sum={old_genes} [{old_hits}/10]")
            print(f"{cancer}: unique-PMID ={new_genes} [{new_hits}/10]")

    print(f"\nCancers with any top-10 membership difference: {diff_cancers} / 19")
    print(f"Pooled Precision@10, citation-count summation: {old_total}/190 = {old_total/190:.1%}")
    print(f"Pooled Precision@10, unique-PMID deduplication: {new_total}/190 = {new_total/190:.1%}")


if __name__ == "__main__":
    main()
