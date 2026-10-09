"""
compute_frequency_baseline_with_aliases.py
==========================================
Frequency-only baseline on OncoDigger's own gene vocabulary (Methods 2.9; Figure 3A,B).

Vocabulary: the single-token gene terms of the OncoDigger gene lexicon (approved symbols of any length and the
aliases that survive the engine's filters), exactly as used by the ranking engine; a token that maps to several
genes credits each of them, as in the engine. The same symbol exclusions (run_query_planes.GENE_EXCLUDE) are
applied. Multi-word full-name phrases are not used (postings are token-level).

Per-cancer per-gene count = |distinct PMIDs containing any of the gene's tokens| (union, not sum), over the whole
cancer corpus; no BM25, no query, no enrichment. Genes ranked by this count; ties broken alphabetically.

Outputs (data/derived):
  - freq_baseline_top50_aliases.csv             frequency-only top-50 per cancer with paper counts

Retired 2026-10-08: this script originally also scored its own top-50 against COSMIC Cancer
Gene Census Tier 1 (`tier1` column, `cosmic_census_2026-05-23.csv`, not included in this
repository) at several k values, comparing against OncoDigger's own legacy, COSMIC-scored
`precision_at_k` column (since removed from `precision_at_k_detail.csv`). Neither is used by
manuscript1.1's IntOGen-calibrated analysis and both have been removed; only the frequency-only
ranking itself (gene, rank, cancer, papers) is still used, scored against IntOGen by
intogen_ablation_panel.py.
"""
from __future__ import annotations

import sqlite3
import sys
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd

HGNC_FILE = Path(__file__).resolve().parents[2] / "data" / "reference" / "Human_Genes.csv"
CORPUS_BASE = Path(r"C:\Users\van0g\Documents\Codex\2026-05-03\can-you-get-access-to-my\OncoDigger\data\pubmed_processed")
OUT_DIR = Path(__file__).resolve().parents[2] / "data" / "derived"

CANCER_TO_CORPUS = {
    "Breast Cancer": "breast_cancer_inclusive_all",
    "Lung Cancer": "lung_cancer_inclusive_all",
    "Colorectal Cancer": "colorectal_cancer_inclusive_all",
    "Prostate Cancer": "prostate_cancer_inclusive_all",
    "Melanoma": "melanoma_cancer_inclusive_all",
    "Bladder Cancer": "bladder_cancer_inclusive_all",
    "Kidney Cancer": "kidney_cancer_inclusive_all",
    "Pancreatic Cancer": "pancreas_cancer_inclusive_all",
    "Liver Cancer": "liver_cancer_inclusive_all",
    "Stomach Cancer": "stomach_cancer_inclusive_all",
    "Esophageal Cancer": "esophageal_cancer_inclusive_all",
    "Ovarian Cancer": "ovarian_cancer_inclusive_all",
    "Endometrial Cancer": "endometrial_cancer_inclusive_all",
    "Cervical Cancer": "cervical_cancer_inclusive_all",
    "Thyroid Cancer": "thyroid_cancer_inclusive_all",
    "Brain Cancer": "brain_cancer_inclusive_all",
    "Leukemia": "leukemia_cancer_inclusive_all",
    "Lymphoma": "lymphoma_cancer_inclusive_all",
    "Myeloma": "myeloma_cancer_inclusive_all",
}

SOFT_EXCLUDE = {"gc","hr","mice","dcr","bpifa4p","pc","hccs","wdhd1","sgcg","c2","fbn1"}
PRIMARY_ENGLISH_COLLISIONS = {"impact", "nodal"}
ENGLISH_COLLISIONS = {
    "was","set","are","for","the","with","have","after","from","more",
    "test","case","type","form","side","rate","size","step","rest","spin",
    "fish","chip","find","stop","simple","base","cell","dna","rna",
    "max","min","best","junk","many","much","ever","even","fine",
    "rank","shape","real","true","main","rich","pick","mass",
    "pace","axis","core","hand","head","back","hat",
} | PRIMARY_ENGLISH_COLLISIONS


def main():
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import run_query_planes as rq
    sys.path.insert(0, str(rq.DEFAULT_ONCODIGGER / "src"))
    from oncodigger.search.lexicons import load_lexicons
    genes_lexicon = load_lexicons(rq.DEFAULT_ONCODIGGER)["Genes"]
    # Apply the same term-level ambiguous-symbol suppression as the primary ranking pipeline
    # (2026-09-30), so this "vetted" baseline uses the identical vocabulary policy the Methods
    # text claims ("same alias exclusions and substitutions used throughout this paper's primary
    # ranking"), rather than relying on the now-narrowed GENE_EXCLUDE (7 true exclusions only) to
    # also suppress the 9 bare-symbol collisions it no longer covers.
    single = rq.suppress_ambiguous_symbols(genes_lexicon)["single"]
    token_to_genes: dict[str, set[str]] = defaultdict(set)
    for token, entities in single.items():
        for e in entities:
            g = rq.GENE_REMAP.get(e["id"], e["id"])
            if g not in rq.GENE_EXCLUDE:
                token_to_genes[token].add(g)
    print(f"  Engine gene vocabulary: {len(token_to_genes):,} single tokens, "
          f"{len(set().union(*token_to_genes.values())):,} genes")

    # Step 4: per-cancer top-50 by |union(distinct PMIDs)| across gene tokens
    print()
    print("Querying each cancer corpus...")
    freq_rankings: dict[str, list[str]] = {}
    freq_counts: dict[str, list] = {}
    tokens_list = list(token_to_genes.keys())
    chunk_size = 800

    for cancer, corpus_name in CANCER_TO_CORPUS.items():
        db = CORPUS_BASE / corpus_name / "lexical_index.db"
        t0 = time.time()
        con = sqlite3.connect(str(db))
        cur = con.cursor()
        gene_pmids: dict[str, set] = defaultdict(set)
        for i in range(0, len(tokens_list), chunk_size):
            chunk = tokens_list[i:i + chunk_size]
            placeholders = ",".join("?" * len(chunk))
            cur.execute(
                f"SELECT token, pmid FROM postings WHERE token IN ({placeholders})",
                chunk,
            )
            for token, pmid in cur.fetchall():
                for g in token_to_genes[token]:
                    gene_pmids[g].add(pmid)
        con.close()
        ranked = sorted(((g, len(p)) for g, p in gene_pmids.items()), key=lambda x: (-x[1], x[0]))
        freq_rankings[cancer] = [g for g, _ in ranked[:50]]
        freq_counts[cancer] = ranked[:50]
        print(f"  {cancer:22s}  unique genes: {len(gene_pmids):5d}  [{time.time()-t0:.1f}s]")

    pd.DataFrame([{"cancer": c, "rank": i, "gene": g, "papers": n}
                  for c, lst in freq_counts.items() for i, (g, n) in enumerate(lst, 1)]
                 ).to_csv(OUT_DIR / "freq_baseline_top50_aliases.csv", index=False)

    print()
    print("=== Frequency-only top-10 lists (alias-inclusive) ===")
    for cancer in CANCER_TO_CORPUS:
        print(f"  {cancer:22s} {freq_rankings[cancer][:10]}")


if __name__ == "__main__":
    main()
