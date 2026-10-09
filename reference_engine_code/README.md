# Vendored engine source (reproducibility pin)

Addresses red-team audit item 15: the manuscript's Methods described the "Primary Evidence"
weighting function and the BM25 retrieval implementation as living in a separate application
repository (`github.com/jfataphd/OncoDigger`), which continues to be actively developed. Citing
that repo without an exact commit pin means a future reader cloning it would get a different,
evolved version of this code, not the one that actually produced this paper's numbers.

These four files are an exact, byte-for-byte copy of the engine source this paper's entire
analysis pipeline (`scripts/analysis/run_query_planes.py` and everything downstream of it)
actually imports and calls, at the exact commit the manuscript's own frozen release depends on:

- **Source repository:** `github.com/jfataphd/OncoDigger`
- **Exact commit:** `4e39fc0681024dc4cef5b28978386f40878bcb4e`
- **Commit date:** 2026-05-31
- **Commit message:** "Add alias substitution for ambiguous gene symbols and gene validation UI"
- **Verified:** working tree at this commit has zero uncommitted local modifications (`git status`
  clean) as of the date these files were copied into this repository, 2026-09-30.

Files (original path `src/oncodigger/search/` in the source repository):

| File | Role in this paper's pipeline |
|---|---|
| `lexical.py` | BM25 retrieval (k1=1.5, b=0.75) and metadata multipliers (title/MeSH/recency/review) |
| `discovery_modes.py` | The "Primary Evidence" content-type weighting function (Methods 2.1) |
| `discovery.py` | Entity aggregation and gene enrichment scoring |
| `lexicons.py` | HGNC gene vocabulary loading, including the `useful_gene_term` common-word filter |

These are reference copies for reproducibility only; they are not re-imported by anything in
this repository (the analysis scripts still import live from `--oncodigger-root`, defaulting to
the frozen local working copy at the commit above). If that working copy is ever lost or
modified, these four files are the fallback source of truth for exactly what code produced every
number in this paper.

Also included here, extracted from the same exact commit for the same reason: the gene,
drug, phytochemical and compound vocabulary files `lexicons.py` loads (`Human Genes.csv`,
`kegg_drug_list_lowercase.csv`, `phytochemicals.csv`, `kegg_compounds_lowercase.csv`). Only
the gene vocabulary is used by the mutation-focused ranking this paper reports, but
`load_lexicons()` loads all four unconditionally, so all four are included together.

Re-running `run_query_planes.py` end to end additionally requires the 19 cancer-specific
SQLite corpora it builds its ranking from (Methods 2.1) — not included here, since they total
over 1.6 million PubMed abstracts; they are reproducible from the freely available NCBI
E-utilities API (exact queries: Supplementary Table S2). Every result built on top of the
frozen ranking this script already produced (`data/canonical/precision_at_k_detail.csv`) —
the IntOGen scoring, all comparator benchmarks, all ablations and controls — is fully
reproducible directly from the files in this repository, independent of the raw corpora.
