"""
compute_frequency_baseline_primary_only.py
===========================================
QC-vocabulary frequency-only baseline (reviewer-requested robustness check on Methods 2.9;
see docs/reviewer_revision_plan.md item 7). Restricts the vocabulary used by
compute_frequency_baseline_with_aliases.py to PRIMARY HGNC symbols only -- no aliases, no
approved (full) names -- to test whether the gap between OncoDigger and a frequency-only
baseline persists once alias-collision noise (Section 2.2) is structurally impossible
rather than merely down-weighted by the engine's alias filters.

Vocabulary: for each gene entity in the OncoDigger gene lexicon, only the single token that
equals the entity's own primary display symbol (the HGNC-approved symbol, or its
disambiguating preferred alias for the small number of entities indexed that way; see
oncodigger.search.lexicons.load_lexicons) is kept. All alias/full-name tokens are dropped.
The same symbol exclusions (run_query_planes.GENE_EXCLUDE) are applied.

Per-cancer per-gene count = |distinct PMIDs containing the gene's primary-symbol token|,
over the whole cancer corpus; no BM25, no query, no enrichment. Genes ranked by this count;
ties broken alphabetically.

Outputs (data/derived):
  - freq_baseline_top50_primary.csv             frequency-only top-50 per cancer with paper counts

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


def main():
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import run_query_planes as rq
    sys.path.insert(0, str(rq.DEFAULT_ONCODIGGER / "src"))
    from oncodigger.search.lexicons import load_lexicons, norm_text

    genes_lexicon = load_lexicons(rq.DEFAULT_ONCODIGGER)["Genes"]
    # Term-level ambiguous-symbol suppression (2026-09-30), same as the primary pipeline: for the
    # 9 symbols whose own primary display label IS the ambiguous bare symbol (MB, ELN, CP, UBC,
    # CAMP, EFS, IVD, SCT, SON), this removes their only single-token entry entirely, which is the
    # semantically correct outcome for this specific "primary-symbol-only, no aliases" variant --
    # a gene whose sole candidate token is inherently collision-prone cannot be meaningfully
    # represented once aliases are stripped, the same reason these symbols cannot rely on the
    # alias-substitution mechanism (GC->VDBP etc.) that resolves the other four. Genes with a
    # preferred-alias label (GC/HR/HCCS/BPIFA4P, already indexed under VDBP/KDM3D/CCHL/LATH) are
    # unaffected, since suppression only ever targets a bare, still-ambiguous primary symbol.
    single = rq.suppress_ambiguous_symbols(genes_lexicon)["single"]

    # Primary-symbol-only vocabulary: keep a (token -> entity) pair only when the token IS
    # that entity's own primary display symbol, i.e. drop every alias/full-name token.
    token_to_genes: dict[str, set[str]] = defaultdict(set)
    dropped_alias_tokens = 0
    for token, entities in single.items():
        for e in entities:
            if token != norm_text(e["label"]):
                dropped_alias_tokens += 1
                continue
            g = rq.GENE_REMAP.get(e["id"], e["id"])
            if g not in rq.GENE_EXCLUDE:
                token_to_genes[token].add(g)
    n_genes = len(set().union(*token_to_genes.values())) if token_to_genes else 0
    print(f"  Primary-symbol-only vocabulary: {len(token_to_genes):,} single tokens, "
          f"{n_genes:,} genes (dropped {dropped_alias_tokens:,} alias/full-name token-entity pairs)")

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
                 ).to_csv(OUT_DIR / "freq_baseline_top50_primary.csv", index=False)

    print()
    print("=== Frequency-only top-10 lists (primary-symbol-only) ===")
    for cancer in CANCER_TO_CORPUS:
        print(f"  {cancer:22s} {freq_rankings[cancer][:10]}")


if __name__ == "__main__":
    main()
