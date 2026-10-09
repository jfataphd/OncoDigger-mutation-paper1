# OncoDigger: IntOGen-calibrated validation (companion repository)

This repository contains the code and data needed to reproduce every
numerical result, table and figure in the manuscript:

> OncoDigger Reveals Distinct Cancer-Specific Mutation Landscapes from the
> Published Literature

validated against IntOGen (Martínez-Jiménez et al., 2020) as the primary
cancer-type driver reference standard.

## About the author

Jimmie Fata is a Professor of Biology at the College of Staten Island, CUNY.
His background is in molecular biology, mammary gland biology and cancer
research, with training at the University of Guelph, the University of
Toronto and Lawrence Berkeley National Laboratory. His current interests
include cancer bioinformatics, literature mining, AI-assisted biology and
educational technology. He created OncoDigger to explore how computational
tools can support biological research.

## Contents

- `manuscript/` — the manuscript source (`OncoDigger_manuscript1.1.md`) and the
  rendered PDF.
- `scripts/analysis/` — every analysis script cited in the manuscript,
  including the OncoDigger ranking-generation script (`run_query_planes.py`).
- `scripts/analysis/pubtator_v1_single_parent_reconstruction.py` — a live-API
  reconstruction of the original, pre-correction single-parent-MeSH PubTator
  methodology, used to confirm the original "before" benchmark percentage
  could not be independently re-verified (Supplementary Table S3).
- `data/external/intogen/` — the IntOGen driver-gene compendium (release
  2024.09.20, CC0 licence), included with no access restrictions.
- `data/external/chatgpt/` — the saved ChatGPT free-tier response used as a
  practical baseline.
- `data/external/PROVENANCE.md` — exact release identity, filenames and
  checksums for every external comparator dataset (CancerMine, DISEASES,
  PubTator).
- `data/derived/`, `data/canonical/`, `data/derived/intogen_rescoring/` — every
  derived result table the manuscript's text, figures and supplementary
  material cite, as actually produced by the scripts above.
- `figures/png/` — the rendered PNG for every main and supplementary figure.
- `reference_engine_code/` — an exact, byte-for-byte vendored copy of the
  OncoDigger ranking engine source this paper's pipeline depends on, pinned to
  commit `4e39fc0681024dc4cef5b28978386f40878bcb4e` (2026-05-31) of
  `github.com/jfataphd/OncoDigger`, together with the gene/drug/phytochemical/
  compound vocabulary files that engine loads, extracted from that exact same
  commit. See `reference_engine_code/README.md`.

## What this does and does not include

The 19 cancer-specific PubMed corpora themselves (over 1.6 million abstracts)
are not bundled here — they are reproducible from the freely available NCBI
E-utilities API with no subscription or institutional access required (exact
queries: Supplementary Table S2 in the manuscript). `run_query_planes.py`
documents exactly how the primary ranking (`data/canonical/
precision_at_k_detail.csv`) was generated from those corpora plus the vendored
engine and vocabulary files above; everything built on top of that frozen
ranking (the IntOGen scoring, CancerMine/DISEASES/PubTator comparisons, all
ablations and controls) is fully reproducible directly from the files in this
repository, independent of the raw corpora.

Two large, continuously-updated external bulk files (DISEASES' full
text-mining release and PubTator's full relation-extraction release) are not
committed for size reasons; `data/external/PROVENANCE.md` records their exact
checksums and access dates, and the relevant scripts document how to
regenerate the filtered subsets actually used.

`scripts/analysis/vocab_corrected_rerun.py` physically corrects three cancers'
working-copy corpora (Lung, Thyroid, Kidney) to test a vocabulary-collision
rule found during the audit, then `vocab_corrected_rescoring.py` re-runs the
real scoring engine on the copies and scores the result against IntOGen. The
corrected copies themselves (multi-gigabyte) are not committed, only the
resulting small output table (`vocab_corrected_rankings.csv`); the copies are
fully reproducible from the original corpora plus the documented, frozen
correction rule in `vocab_collision_pipeline.py`.

The figure-generation and PDF-building scripts (PowerPoint/ReportLab
presentation tooling, not scientific computation) and the editable PowerPoint
figure sources are kept in the authors' development repository rather than
here; every number, table and statistic they draw from is in this repository
and independently reproducible from the scripts and data above. The rendered
figures themselves are included as PNGs.

## Reproducing the results

Install dependencies (see the installation note in `requirements.txt` for
why this is two steps, not a plain `pip install -r requirements.txt`):

```
pip install numpy==2.4.6 pandas==2.2.3 requests==2.34.2
pip install --no-deps scipy==1.14.1
```

Then, confirmed via a clean-checkout test (2026-10-07) to reproduce the
manuscript's reported numbers exactly:

```
python scripts/analysis/intogen_full_rescoring.py
python scripts/analysis/intogen_cleaned_recall.py
```

`scripts/analysis/intogen_miss_root_cause.py` reproduces the same pipeline
logic and classification rules, but its title-sampling step queries the raw
per-cancer PubMed corpus databases directly (not bundled here, see "What
this does and does not include" above), so it will not run standalone from
only the files in this repository; its own docstring documents this.

Each script's own docstring documents its exact inputs and outputs.

## License and citation

Code and data in this repository are released under the MIT License (see
`LICENSE`). See `CITATION.cff` for how to cite this work; a Zenodo DOI will
be added once the accompanying GitHub release is cut.
