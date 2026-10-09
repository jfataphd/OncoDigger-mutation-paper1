<!--
BUILD NOTES FOR CLAUDE CODE (NOT PART OF THE MANUSCRIPT)

Target journal: Bioinformatics (Oxford University Press)
Article type: Original Paper
Target category: Data and text mining
Journal guidance checked 2026-09-27:
- Original Papers: up to 7 journal pages, approximately 5,000 words excluding figures.
- Structured abstract: Motivation; Results; Availability and Implementation; Contact; Supplementary Information; recommended maximum 150 words.
- Recommended manuscript sequence: Title page; Structured Abstract; Introduction; System and methods; Algorithm; Implementation; Discussion; References.
- Initial submission is format-free, but the text below is organized to be close to Bioinformatics conventions.
- Keep main figures embedded at the indicated markers when generating the review PDF.
- Supplementary material is included after the references so no study data or interpretations are lost while the main paper remains concise.
- Do not invent or alter numerical values, citations, datasets, URLs, software versions, or conclusions.
- Do not silently replace the frozen manuscript pipeline with output from the live OncoDigger website.
-->

# OncoDigger Reveals Distinct Cancer-Specific Mutation Landscapes from the Published Literature

**Jimmie E. Fata** (ORCID 0000-0001-5797-6619)

Department of Biology, College of Staten Island, City University of New York, Staten Island, New York, USA

**Corresponding author.** jimmie.fata@csi.cuny.edu

**Keywords.** literature mining, cancer genes, mutations, BM25, IntOGen

---

## Abstract

**Motivation.** Global biomedical literature can make broadly studied cancer genes dominate disease-specific searches. OncoDigger instead ranks genes within separate cancer-specific PubMed collections.

**Results.** Across 19 cancers, a mutation-focused query yielded 93.2% Precision@10 against IntOGen (177/190 gene-cancer pairs), a cancer-type-resolved driver catalogue from real tumour-sequencing cohorts (Martínez-Jiménez et al., 2020). Document retrieval and gene-vocabulary quality control interact substantially rather than independently (interaction +45.3 percentage points beyond their additive effects, positive in all 19 cancers), against an unfiltered-vocabulary, full-corpus baseline of 10.5% Precision@10. The cancer-aware strategy reached 75.8% exact cancer-type concordance, versus 58.4% for a pooled, cancer-blind index, higher than CancerMine (60.0%). OncoDigger combines cancer-specific literature organization with strong driver-gene concordance.

**Availability and Implementation.** Code, data and the frozen analysis release are at the companion repository, github.com/jfataphd/OncoDigger-mutation-paper1 (archived at Zenodo, DOI: 10.5281/zenodo.23269250). The web application is at oncodigger.com.

**Contact.** jimmie.fata@csi.cuny.edu

**Supplementary Information.** Supplementary data are available online.

---

## 1 Introduction

OncoDigger (oncodigger.com) is a cancer literature-mining platform that searches a chosen cancer's PubMed abstracts and ranks results specific to that cancer, instead of one pooled list where broadly-studied genes dominate regardless of which cancer is searched. It maintains 19 separate cancer-specific collections (1,640,505 abstracts, 2000 through September 2026, Methods 2.1), ranking documents within each cancer-specific collection using BM25 before gene enrichment and ranking. Cancer-specific separation recovers genes absent from a pooled top-10, while the complete method achieves 93.2% Precision@10 against IntOGen, an external, cohort-computed driver gene reference.

This paper evaluates **mutation-focused cancer-gene prioritization** only. Other query types, such as expression, pathways or drug response, need their own reference standards.

The same `mutation mutations` query was applied to all 19 cancers, scored against IntOGen (633 driver genes, release 2024.09.20, Martínez-Jiménez et al., 2020), used as an **external, cohort-computed calibration reference**, not ground truth. IntOGen's gene membership, cancer-type assignments and significance calls were never used to tune OncoDigger's ranking algorithm, scoring parameters or vocabulary rules to improve agreement with it. One later, disclosed audit (Methods 3.3) uses IntOGen membership only to select genes for additional vocabulary-collision scrutiny, a search strategy, not a tuning decision. Each correction was determined solely by independently reading real PubMed text, not by whether it would improve agreement with IntOGen. A secondary `driver mutation driver mutations` query served as a wording sensitivity check.

OncoDigger was compared with CancerMine (Lever et al., 2019a), DISEASES 2.0 (Pletscher-Frankild et al., 2015, and Grissa et al., 2022) and PubTator 3.0 (Wei et al., 2024). DisGeNET and CIViCmine were excluded as scope-mismatched (Piñero et al., 2017, and Lever et al., 2019b). These three draw on one combined, cancer-general literature base and differ from OncoDigger in gene recognition, source data and scoring, so a score difference reflects the whole system, not corpus separation specifically. A separate matched pooled-vs-separated comparison characterizes the cancer-aware versus cancer-blind strategy instead (Discussion). A single ChatGPT free-tier response was evaluated separately as a practical baseline, not a peer system.

The study asks three questions. Does separating cancer literature produce different rankings? Does the scoring method recover established cancer genes better than simple counting? How does the complete system compare with published literature-mining resources (Table 1)?

The cancer-aware strategy recovered a mean of 5.5 additional top-10 genes per cancer relative to the pooled, cancer-blind control (Results 5.1). The complete pipeline achieved 93.2% Precision@10 compared with 43.2% for the publication-counting control, a 50.0-percentage-point difference reflecting changes in retrieval, vocabulary handling and scoring together (Results 5.2). The matched factorial analysis instead identified a substantial interaction between document selection and vocabulary cleaning, with only a 1.6-percentage-point improvement in any-cancer Precision@10 from the full scoring formula after vocabulary cleaning (Results 5.8).

<!-- INSERT TABLE 1 HERE -->

---

## 2 System and methods

### 2.1 Cancer-specific literature collections

The study included 19 cancer types (Breast, Lung, Colorectal, Prostate, Melanoma, Bladder, Kidney, Pancreatic, Liver, Stomach, Esophageal, Ovarian, Endometrial, Cervical, Thyroid, Brain Cancer, Leukemia, Lymphoma, Myeloma). Each collection contained English-language PubMed abstracts published 2000 through September 2026, selected via a single cancer-specific MeSH term through the NCBI PubMed E-utilities API. Collection sizes ranged from 22,079 to 250,061 abstracts (1,640,505 papers across all 19 before cross-collection deduplication, Supplementary Table S1, with queries in Supplementary Table S2).

Document scoring included title and MeSH matches, publication recency and OncoDigger's "Primary Evidence" discovery mode, which up-weights original-research/clinical-evidence terminology over reviews, retained but down-weighted by 0.65, not excluded.

The same `mutation mutations` query, ranked by BM25 (Methods 3.1), selected the top 1,000 highest-ranked abstracts per cancer for gene scoring. Pool-size sensitivity testing across 250-2,000 papers confirmed this is a stable operating point (Supplementary Methods S1).

### 2.2 IntOGen calibration and internal controls

For each cancer, Precision@10 is the percentage of the top 10 ranked genes confirmed as drivers anywhere in IntOGen's 633-gene compendium. Statistical confidence uses case-bootstrap (19 cancer types resampled with replacement, 2,000 times), the primary CI throughout since the cancer type, not the individual gene, is the real experimental unit. A Wilson CI treating all 190 slots as independent is also reported for comparison but understates cross-cancer variation (Supplementary Results S1). As orientation only, random gene selection would hit a true IntOGen driver about 3.3% of the time (633/19,250 protein-coding genes, Seal et al., 2026), well below expected literature-mining performance.

Beyond the top 10, ranking quality across each cancer's top 50 was summarised by rank AUC (probability a random driver gene ranks above a random non-driver gene, 0.5=chance, 1.0=perfect, averaged across cancers) and a Spearman correlation testing whether the advantage holds down the list, not just at the top (`scripts/analysis/secondary_stats.py`).

Two controls test which part of the pipeline drives results. **Pooled-literature control** merges all 19 collections into one combined pool (1,514,148 unique papers) and runs the same procedure once. Comparing this single cancer-blind ranking with the separate cancer-specific analyses evaluates the combined effect of cancer-aware corpus organization and retrieval, rather than isolating physical index separation alone.

**Full pipeline versus publication-counting control.** The 19 collections were kept separate, but BM25, rank weighting and enrichment were replaced with a raw mention count. This changes retrieval, gene-recognition and ranking together, so it is a complete-pipeline-versus-counting comparison, not an isolated formula test (Results 5.8 gives the properly isolated decomposition). A stricter variant counted only HGNC primary symbols, no aliases (`scripts/analysis/intogen_ablation_panel.py`).

### 2.3 Exact cancer-type driver concordance

Beyond the any-cancer benchmark (Methods 2.2), a stricter question asks whether external evidence supports a gene as a driver specifically within the cancer it was ranked in, answered against IntOGen (release 2024.09.20, cancer-type-resolved, computed from real tumour-sequencing cohorts via several combined statistical methods, Martínez-Jiménez et al., 2020). IntOGen's cancer-type codes were conservatively mapped to the 19 study categories (ambiguous or duplicate codes excluded, full mapping in Supplementary Methods S5). This **exact cancer-type driver concordance** is not cancer-specific biological truth. IntOGen detects somatic point-mutation signal only, so germline- or fusion-driven genes can be real drivers without reaching significance there (Results 5.3). The same rule was applied unchanged to CancerMine, the pooled control, and ChatGPT.

### 2.4 External resources and benchmark auditing

CancerMine, DISEASES 2.0 and PubTator 3.0 were downloaded within the same four-day window (2026-09-26 to 2026-09-30, exact identities in `data/external/PROVENANCE.md`), audited and corrected before any number was reported final (Supplementary Table S3, Supplementary Methods S2-S4).

**CancerMine** citation counts were summed across matched disease subtypes and roles (Lever et al., 2019a), checked against a fully deduplicated alternative (negligible difference, Supplementary Methods S2). The audit found a real entity-recognition error. All 1,329 sentences credited to PTGDR were unrelated lncRNA mentions, placing PTGDR in its top 10 for 13/19 cancers, so it was removed from the benchmark and the author notified. **DISEASES 2.0** used the full unfiltered release, each gene's single highest z-score per cancer, with disease-type overrides preventing hybrid-term misrouting (e.g. `Liver lymphoma`). **PubTator 3.0** mapped each cancer to its parent MeSH descriptor plus every matched subtype (Supplementary Table S5), since single-parent mapping silently excluded subtype-tagged relation data, with genes ranked by distinct supporting-PMID count. **ChatGPT practical baseline.** All 19 cancer names were submitted in one prompt to the free, logged-out web interface (27 September 2026), asking for each cancer's top 10 most frequently mutated genes, with prompt, response and access conditions saved. The response stated it used live web searches grounded in TCGA/cBioPortal mutation-frequency data, so as one non-versioned response, it is treated as an informal practical reference, not a competing system.

### 2.5 Literature weight of unrecovered annotations

To determine whether IntOGen driver gene-cancer pairs absent from OncoDigger's top 10/50 (Results 5.3, Supplementary Results S6) reflect weak literature coverage rather than a ranking or vocabulary failure, all 1,740 pairs were scored for rank at any depth and unweighted supporting-paper count, by re-running the search engine directly against the corpus databases at an extended entity limit (`scripts/analysis/intogen_literature_weight_audit.py`).

This supporting-paper count is not fully independent of rank, since it is one of the three multiplied inputs to the relevance score (rank-decayed mention weight, background-enrichment ratio, log(1+supporting papers)). The deep extended ranking run matched the official top 50 exactly for all 430 IntOGen-driver genes present in both, with zero rank discrepancies.

Every IntOGen driver gene never appearing in OncoDigger's ranking at any depth was also checked directly against the gene vocabulary lexicon, since a gene absent from the vocabulary cannot be recognized regardless of literature volume. Genes were not classified by mutation mechanism (no curated field in IntOGen's data), and no external live-PubMed cross-check was performed for the larger never-surfaced set (Results 5.5).

---

## 3 Algorithm

### 3.1 Document ranking

Each abstract was ranked using Okapi BM25 (Robertson and Zaragoza, 2009), with parameters k1=1.5 and b=0.75 held fixed across all 19 cancers and never subsequently modified (Supplementary Results S8 confirms Precision@10 is not sharply sensitive to these values). Inverse-document frequency and average document length were calculated separately within each cancer-specific collection.

The BM25 score was combined with fixed multipliers, namely 1.2 for a title match, 1.15 for a matching MeSH term, a recency multiplier (1.1 for 2020+, 0.95 for pre-2010), 0.65 for reviews, and the "Primary Evidence" content-type multiplier, defined in OncoDigger's own source code (private live-application repository, commit `4e39fc068102`, 2026-05-31) and vendored here for reproducibility (Supplementary Methods S1). The highest-scoring 1,000 papers formed the gene-ranking pool.

### 3.2 Gene enrichment

Gene names were recognized using an HGNC-derived vocabulary of approved gene symbols and curated aliases (Seal et al., 2026). A mention's weight decayed with the inverse square root of its paper's BM25 rank, so a mention in the #1-ranked paper carries more weight than the same mention in the #100-ranked paper.

For each gene, an enrichment score checked whether it appeared in the query-relevant pool much more often than in background literature.

`enrichment = (supporting_papers / pool_size) / (background_papers / corpus_size)`

A gene mentioned everywhere in the background has a large denominator, lowering its score and suppressing generic genes.

`Relevance = weighted_score × log(1 + enrichment) × log(1 + supporting_papers)`

A gene needed >=2 independent supporting papers to be reported, rewarding genes preferential to the most query-relevant papers while suppressing those common throughout the background literature.

### 3.3 Gene-name quality control

The HGNC gene vocabulary was cleaned before benchmarking. 143 aliases colliding with another gene's own approved symbol were removed, aliases doubling as ordinary English words were filtered out, a longest-match rule was applied to 8,203 nested-name cases (e.g. EGF not counted inside `epidermal growth factor receptor`), 7 symbols with no safer alias were excluded entirely (`GENE_EXCLUDE`), and 13 further symbols had only their bare ambiguous symbol suppressed, leaving their full name and other aliases matchable (full symbol list, live-policy synchronization status, and downstream impact, Supplementary Methods S1).

Every gene in the top 30 of any cancer, and all 190 top-10 pairs, was checked by hand against a representative PubMed sample, not every posting (a gene contaminated in only a minority of postings, as MET/REST later turned out, could pass this). Extending through rank 50 (223 genes) surfaced one false alias (MB, cell-line names), prompting a sweep of every short symbol (<=4 characters) in the top 50 (66 candidates), finding eight more (Supplementary Methods S1). All nine were corrected and the pipeline re-run from scratch. Only ELN is an IntOGen driver, and it appeared in no top-50 under a safe alias either, so no detection was lost (full metrics, Supplementary Methods S1). All numbers in this paper reflect this corrected pipeline.

A further audit applied a uniform two-step verification (standalone capitalized match, then a pre-identified alternate meaning, both frozen before any comparison against this paper's statistics) to every confirmed IntOGen driver gene in the top-50 rankings (196 genes). A separate audit flagged 23 as plausible collision candidates, all received the full check. Two required correction, MET (Lung Cancer, rank 10) and REST (Thyroid Cancer, rank 50). Of MET's 3,967 raw postings, 1,675 (42.2%) passed both checks and 2,292 (57.8%) were excluded (the ordinary word "met," or capitalized matches to alternate meanings, namely a process acronym, a bacterial gene, a PET tracer, a drug-dosing term, or a mesothelial cell line). Passing both checks means no recognized counter-evidence was found, not that each posting was individually confirmed genuine (below). REST's 136 postings were all excluded. Both corpora were corrected by redacting the false-positive mentions directly from the title and abstract text of the affected postings, not only removing them from the retrieval index, and the engine was re-run on the corrected text. MET's Lung Cancer supporting-paper count fell from 65 to 56 and it retained its top-10 position (10th), leaving both headline metrics unchanged. REST's Thyroid Cancer supporting-paper count fell to zero, genuinely removing it from the ranking rather than leaving its postings count at zero while its text-based mentions remained, and it was replaced at rank 50 by MAPK3 (not an IntOGen driver), reducing the ranks-31-50 interval from 239/380 to 238/380. The remaining 21 of 23 flagged genes required no correction.

A robustness check addressed whether "passing both checks" might be too permissive, since it only means no counter-evidence was found, not that positive context was present. A random sample of 40 of MET's postings passing the rule but matching no positive gene-context term (Supplementary Methods S1) was read in full. 33/40 (82.5%) were genuine (typically multi-gene panel lists lacking a single-gene qualifying phrase), one ambiguous, and 6 false, revealing four further false-meaning patterns (metformin abbreviation, methionine-codon genotype notation, an abbreviated already-excluded false meaning, and a mesothelial cell line). These patterns were incorporated into the recognition rule, excluding 32 further false matches and reducing the retained MET postings from 1,707 to 1,675. Retention indicates a posting passed the revised recognition criteria, not that it was individually verified as a genuine MET reference. A fifth candidate pattern (metastasis-status abbreviation) was discarded after also matching genuine high-volume usage ("MET-amplified") and the unrelated word "unmet," so this single paper's false contribution is disclosed as unresolved rather than patched with an overly broad rule. This 33/40 estimate is diagnostic, from the same sample used to discover the patterns, not independent validation. As a separate, stricter check free of that circularity, the pipeline was re-run keeping only the 1,300 postings (32.8%) with both a capitalized match and a positive context term. Lung Cancer's top-10 and MET's rank were unchanged under both the 1,675- and 1,300-posting conditions.

### 3.4 Post hoc decomposition of pipeline components

After the vocabulary audits above, a further, explicitly post hoc analysis examined how much of the pipeline's IntOGen agreement depends on document retrieval, vocabulary cleaning and the relevance formula individually, and whether they interact, under a frozen written protocol (`FROZEN_ANALYSIS_SPEC_retrieval_vocab_decomposition.md`, companion repository) fixing conditions, metrics and stopping criteria in advance.

Two experiments were run across all 19 cancers using production code. The first is a 2x2 factorial crossing document retrieval (BM25 pool vs. complete corpus) with vocabulary state (cleaned vs. "true-raw"), scored by flat mention-count ranking. This requires retrieval to be unaffected by vocabulary. `bm25_search`'s signature takes no vocabulary argument, confirmed empirically by comparing exact ordered (PMID, score) pairs under both conditions for all 19 cancers, all byte-identical (`scripts/analysis/retrieval_vocab_independence_check.py`).

The second holds retrieval fixed at the BM25 pool and compares the full relevance formula against flat counting under both vocabulary states. A complete three-way factorial was not run. The enrichment term is forced to exactly 1 for every gene when the pool is the full corpus (in-pool and full-corpus rate become identical), and the remaining rank-weighting term would then depend on arbitrary fetch order, so this condition was not evaluated. A further ablation decomposed the formula's two remaining terms (rank-weighted mention position, enrichment) under true-raw vocabulary. Both metrics (Methods 2.2-2.3) used the same case-bootstrap methodology (2,000 resamples, seed 42).

---

## 4 Implementation

The software environment was locked down for reproducibility (`PYTHONHASHSEED=0`, with per-connection caches cleared between cancers), confirmed empirically by independently re-running the pipeline multiple times and getting byte-identical top-50 lists and rank order every time.

Most data and figures came from a frozen internal codebase snapshot (git tag `manuscript1.1-v1.0-vocab-audit-2026-10-04`, retained only in the author's private development history, not independently resolvable by readers), kept separate from the continuously updated public web application. Two parts postdate this snapshot, the MET/REST/FH vocabulary-collision rules were further refined (Methods 3.3), and the post hoc decomposition (Methods 3.4, Results 5.8) was conducted separately, using the same underlying code. The companion repository's `v1.1` release is the actual publicly accessible frozen code and data, covering both, and was verified by clean-checkout reproduction to reproduce every IntOGen-based number this paper reports (Data and code availability).

The public `oncodigger.com` interface is under active development, with corpora that update monthly, and may return different rankings, since it also defaults to "Overview" mode rather than this paper's "Primary Evidence" mode. Exact reproduction requires the frozen release. Figure 1 summarises the pipeline.

<!-- INSERT FIGURE 1 HERE -->

---

## 5 Results

### 5.1 Cancer-specific collections produce different rankings

The common mutation query surfaced both pan-cancer drivers (TP53, top 10 in 19/19 cancers, and KRAS, PIK3CA and BRAF in 11-13/19) and disease-specific drivers confined to one cancer type (Figure 2, Supplementary Results S1). Mean pairwise top-10 Jaccard similarity across 171 cancer pairs was 0.189 (1.0 = identical lists, 0.0 = no shared genes), confirming cancer-tailored rather than generic rankings (Supplementary Results S1).

A pooled-literature control merging all 19 cancers' abstracts into one index returned the identical top-10 gene list for every cancer (100% any-cancer Precision@10 but zero cancer-type differentiation, Supplementary Results S1). Separating the literature instead recovered a mean of 5.5 additional genes per cancer (4.8 confirmed IntOGen drivers, Supplementary Figure S1), 67.3% of which were drivers for that specific cancer, including lineage-defining genes such as VHL/PBRM1/BAP1 for Kidney Cancer (Cancer Genome Atlas Research Network, 2013, Supplementary Results S1). Scored against IntOGen's cancer-type-specific calls, the pooled control's exact concordance falls to 58.4% (111/190), a 17.4-point drop from the separated pipeline's 75.8% (p=0.0004), since its inflated any-cancer score reflects famous pan-cancer genes (Supplementary Results S1).

<!-- INSERT FIGURE 2 HERE -->

### 5.2 IntOGen calibration and the scoring controls

Across all 190 top-10 slots, 177 were verified IntOGen driver genes, a headline Precision@10 of 93.2% (95% cancer-bootstrap CI 88.9-96.9%, Methods 2.2, with a Wilson CI in Supplementary Results S1). Eleven of 19 cancers scored a perfect 10/10, the rest 70-90%, with no significant association between collection size and Precision@10 despite an 11-fold range (Supplementary Figure S2).

Beyond the top 10, IntOGen drivers ranked above non-drivers through the full top 50 (mean rank AUC 0.681, positive in every cancer, range 0.518-0.798), and the rank-driver Spearman correlation was negative in all 19, significant after Bonferroni correction in 2 (Supplementary Results S1).

Three controls isolate why OncoDigger works. A frequency-only control (curated vocabulary, no BM25/enrichment) reached only 43.2% Precision@10, a 50.0-point drop, lower in all 19 cancers (p=3.2x10^-10, Supplementary Results S1). A genuinely unfiltered vocabulary through the full pipeline produced a closely matching 48.9-point drop (Supplementary Results S7). Per Results 5.8's properly isolated factorial, document selection and vocabulary cleaning jointly account for most of this advantage. A stricter, primary-symbols-only variant reached 53.7% (39.5 points below OncoDigger). Including aliases at all (vs. primary symbols only) accounts for 10.5 of the 50.0-point gap, though this does not separately quantify valid aliases versus ambiguous false matches within that 10.5 points, and the remaining 39.5 points reflect retrieval, enrichment and vocabulary changing together, not separable here (Supplementary Results S1).

Together, these controls show that neither cancer-specific corpus separation nor accurate document selection and vocabulary alone is sufficient (Figure 3). **Cancer-specific corpus separation produces disease differentiation**, evidenced by the pooled control's 100% global precision but zero cancer-type specificity (Results 5.1). **Document selection and vocabulary cleaning, acting jointly rather than additively, produce most of the IntOGen-driver enrichment seen here** (Results 5.8). The scoring formula's own further contribution is modest once the vocabulary is already clean.

<!-- INSERT FIGURE 3 HERE -->

### 5.3 Exact cancer-type driver concordance

Exact cancer-type driver concordance, scored against IntOGen (Methods 2.3), was 75.8% (144/190, 95% CI 71.6-80.0%). Of the 46 non-matching top-10 hits, 33 (71.7%) are genes IntOGen confirms as drivers in a different study cancer, often via secondary-tissue patterns (e.g. BRCA1/BRCA2, primary drivers in Breast/Ovarian, also surfacing in Prostate via germline predisposition, Pritchard et al., 2016). The remaining 13 (28.3%) correspond to nine genes IntOGen does not recognise as a significant driver in any cancer type. Individual examination identified plausible biological and clinical explanations for their literature prominence that do not require a significant somatic point-mutation signal in the evaluated cohorts (full per-gene review, Supplementary Results S2, Table S7).

Cervical and Esophageal Cancer, the two lowest-scoring cancers, reach 80% and 60% exact concordance (Results 5.5). A `driver mutation` query variant gave 95.3% any-cancer precision and 78.4% exact concordance (paired t-test p=0.056), confirming the 75.8% result is not sharply wording-dependent (Supplementary Results S3).

Recall of IntOGen's full 1,740-pair driver catalogue (Methods 2.3) increased from 8.3% at rank 10 to 14.1% at rank 20 and 24.7% at rank 50. This limited recovery reflects several contributing factors, including the fixed top-50 analysis window, uneven literature coverage and vocabulary limitations, examined further in Results 5.5 and Supplementary Results S9-S12. Exact matches concentrate sharply toward the top of the ranking (75.8% of rank-1-10 positions vs. 32.5% at ranks 21-50, Supplementary Results S1).

### 5.4 Comparison with CancerMine, DISEASES 2.0 and PubTator 3.0

After benchmark auditing, CancerMine reached 87.9% Precision@10 against IntOGen (167/190, up from 84.2%, Supplementary Results S4), compared with 93.2% for OncoDigger, a non-significant 5.3-point difference (p=0.086). Exact cancer-type concordance is clearer, 75.8% vs. 60.0%, a significant 15.8-point difference (95% CI 10.0-21.1 pp, p=5.1x10^-5), against tumour-sequencing driver evidence, not another literature-curated resource.

DISEASES 2.0 and PubTator 3.0 reached 60.0% and 67.9% Precision@10 respectively, both significantly lower (p<10^-5 each, Supplementary Results S4), though not like-for-like since both capture broader gene-disease associations. Illustrative non-driver outputs (e.g. PDCD1/CTLA4 in Melanoma) are frequently still biologically meaningful (Supplementary Results S4). PubTator's subtype-entity aggregation has a disclosed limitation, letting broadly co-mentioned genes (TP53 in 18/19 cancers) accumulate large unnormalised totals independent of cancer-type specificity, not corrected since that would score PubTator by a method it does not use. An unrelated ALS gene pattern that contaminated Liver Cancer under an earlier, since-corrected mapping did not recur once queried via its own subtypes (Supplementary Results S4).

Figure 4 summarises this head-to-head benchmark.

<!-- INSERT FIGURE 4 HERE -->

### 5.5 Literature weight of unrecovered annotations

To determine whether OncoDigger's top-50 ranking, recovering less than a quarter of IntOGen's full driver catalogue (Supplementary Results S6), reflects a ranking failure on well-supported genes, all 1,740 pairs were scored for rank and unweighted supporting-paper count (Methods 2.5). Support tracked rank sharply. Median papers fell from 87.5 in the top 10 to 15.0 at ranks 11-50 and 3.0 beyond rank 50 (Spearman rho=-0.942, p<10^-300), holding in all 19 cancers. Only 1.0% of non-top-10 pairs matched or exceeded their cancer's weakest top-10 gene (Supplementary Results S1, Table S6). Cervical and Esophageal Cancer follow this same gradient, not a cancer-specific failure (Results 5.3).

Restricting recall to IntOGen pairs both corpus-present and supported by >=2 of its seven detection methods raised recall at every depth (8.3% to 19.6% at rank 10, 14.1% to 30.8% at rank 20, 24.7% to 48.9% at rank 50, 43.3% to 70.2% at any rank, Supplementary Results S11). Most excluded pairs enter IntOGen via its single-method CGC-rescue route (Supplementary Results S9). Directly sampling corpus text for the 191 pairs still unrecovered under this filter assigns the large majority to an identifiable diagnostic category, leaving four pairs unclassified (Supplementary Results S12, Scope and Limitations). Figure 5 summarises this analysis.

<!-- INSERT FIGURE 5 HERE -->

### 5.6 Practical ChatGPT baseline

On exact cancer-type driver concordance, a single ChatGPT free-tier response reached 76.3% (145/190) against IntOGen, no significant difference from OncoDigger's 75.8% (95% CI -7.4 to 7.4 pp). The saved response stated it grounded its answer in TCGA/cBioPortal mutation-frequency data, the same kind of tumour-sequencing evidence IntOGen is computed from (Supplementary Results S4). The gap reopens on gene-level Precision@10, where ChatGPT reached only 78.9% (95% CI 6.8-21.1 pp, Supplementary Figure S3), skewed toward large passenger genes (TTN, MUC16, Supplementary Results S4).

### 5.7 Ranking depth

IntOGen-driver enrichment persisted below rank 10, declining steadily with depth (mean cumulative Precision@k 96.8% at k=5, 93.2% at 10, 88.8% at 15, 84.7% at 20, 80.0% at 30, 76.7% at 40). Non-overlapping rank intervals show the same gradient (93.2% at ranks 1-10 down to 62.6% at ranks 31-50, Supplementary Figure S4), consistently within each cancer individually (Supplementary Results S5).

### 5.8 Post hoc decomposition of retrieval, vocabulary and scoring-formula contributions

A full 2x2 factorial (Methods 3.4) crossing document retrieval (BM25-selected pool vs. complete corpus) with gene-vocabulary state (cleaned vs. genuinely unfiltered "true-raw") was run across all 19 cancers using flat mention-count ranking, isolating the interaction from the scoring formula itself (Table 2).

<!-- INSERT TABLE 2 HERE -->

Retrieval and vocabulary cleaning interacted positively and substantially, +45.3 pp on any-cancer Precision@10 beyond the two main effects' additive prediction (95% cancer-bootstrap CI 36.8-53.7 pp) and +47.4 pp on exact concordance, positive in all 19 of 19 cancers individually (Supplementary Figure S9). The two effects are asymmetric. Document selection alone improved any-cancer Precision@10 by only 4.2 pp (14.7% vs. 10.5%), while vocabulary cleaning alone improved it by 31.6 pp (42.1% vs. 10.5%). Applied together the two reached 91.6%, well above the 48.4% additive prediction. This should be read as a property of this evaluated pipeline and benchmark scale, not a scale-independent synergy measure (odds-ratio equivalents ~10.2 any-cancer, ~1.9 exact concordance).

The vocabulary effect traces overwhelmingly to a small number of specific alias collisions, not general filtering intensity. 10 ambiguous single-token aliases, including ordinary English function words and biomedical abbreviations, account for most of it, and removing only these recovers 76.3% of the full any-cancer gap (Supplementary Results S7, exploratory, not pre-registered).

A second, separate experiment held retrieval fixed at the BM25 pool and compared the full relevance formula (Methods 3.2) against flat counting under both vocabulary states. The formula's own marginal contribution depended on vocabulary quality, +1.6 pp any-cancer under the cleaned vocabulary but +29.5 pp under the true-raw vocabulary. Decomposing the formula's two remaining terms under true-raw vocabulary found rank-weighting alone contributed essentially nothing over flat counting, while enrichment alone (67.9%) exceeded the combined formula (44.2%), plausibly because rank-weighting rewards mentions in highly-ranked documents regardless of whether they are genuine gene references, reintroducing noise enrichment is designed to suppress. This reduction could reflect positional weighting, the additional frequency weighting embedded in that term, or both, which were not separately isolated (Supplementary Results S7).

A complete three-way factorial was not run (Methods 3.4 details why). These findings are two separate, complete experiments, specific to this mutation-focused query and discovery mode, not a general claim about every OncoDigger mode.

Figure 6 shows the complete precision@k curve.

<!-- INSERT FIGURE 6 HERE -->

---

## 6 Discussion

OncoDigger showed strong agreement with a cohort-computed driver-gene reference while producing rankings that differed substantially among cancers. OncoDigger measures **cancer-specific literature prominence**, not biological importance directly, and IntOGen concordance measures agreement with statistically-confirmed somatic driver signal, not independent proof of biological truth.

The controls clarify where performance comes from, and the properly isolated post hoc decomposition (Results 5.8) sharpens this further. Cancer-specific separation contributes **cancer differentiation** (a pooled index reached 100% Precision@10 but returned an identical list for every cancer), though this comparison characterizes the complete cancer-aware versus cancer-blind strategy rather than isolating index separation alone (Results 5.1). Document selection and vocabulary cleaning, not corpus separation, drive **IntOGen-driver enrichment**, interacting substantially rather than additively (Results 5.8), a property of this evaluated pipeline and benchmark scale, not a claim of general mechanistic synergy. The scoring formula's own contribution is similarly conditional on vocabulary quality, traced specifically to its enrichment term rather than rank-weighting under true-raw vocabulary (Results 5.8).

Exact cancer-type driver concordance (Results 5.1, 5.3) shows this differentiation is real. The cancer-blind pooled control scores lower (58.4%) despite exceeding OncoDigger's any-cancer precision (100% vs. 93.2%). Of the 46 non-matching pairs, most reflect genes IntOGen confirms as drivers in a different cancer or biological relationships outside somatic point-mutation detection, not necessarily ranking errors (Results 5.3).

CancerMine was the most informative external comparator, with significantly lower exact concordance (Results 5.4) despite a non-significant Precision@10 difference, though its methods also differ, since CancerMine captures any driver/oncogene/suppressor evidence rather than specifically mutation-focused literature. DISEASES and PubTator's larger gaps need more caution, since both address broader gene-disease relationships. This exercise found and corrected real problems on both sides (Results 5.4), and the same standard applied to OncoDigger's own output found and corrected multiple vocabulary artefacts, including the short-symbol cases and the MET, REST and FH context errors (Methods 3.3). Further OncoDigger data should be expected to need the same treatment, not assumed clean.

The ChatGPT comparison addresses a different practical question. No significant difference was detected on exact concordance (Results 5.6), consistent with its stated grounding in tumour-sequencing data, but its lower gene-level Precision@10, driven by large passenger genes (TTN, MUC16), shows mutation frequency and driver evidence differ (Lawrence et al., 2013). Its live-web-search grounding means its IntOGen agreement cannot be assumed fully independent, and only one non-reproducible response was tested, so no broader AI-system claim applies.

OncoDigger is therefore best viewed as a **literature-ranking and evidence-prioritization system**, not a clinical variant-interpretation tool or replacement for expert-curated driver-gene references. Reaching this performance from freely-accessible PubMed abstracts alone, with no subscription required, is a practical advantage. The results support separating cancer literature, while query relevance and background enrichment improve recovery over simple counting, providing a reproducible way to prioritize mutation-focused literature without treating literature prominence as biological truth.

---

## 7 Scope and limitations

This section states, in one place, what the results do and do not show.

**Recall-gap caveats (Results 5.5).** Never-surfaced IntOGen pairs concentrate in IntOGen's single-method, lower-significance tier (76.6% vs. 40.2% of surfaced pairs, Supplementary Results S9), not necessarily false driver calls. "Never surfaced" is not "absent from the literature." Only 11.4% of pairs have zero lexicon-detected mentions in the full corpus (reflecting vocabulary detection, not confirmed absence), and three genes (H3-3A, H3C2, H4C9) are missing from the vocabulary under current HGNC symbols (Supplementary Results S10). Four pairs with substantial mutation-framed literature remain unresolved after the strictest filter, reported as an open question (Supplementary Results S12).

**Query, corpus and reference-standard scope.** The study validates one mutation-focused query domain against one primary reference standard, and the `driver mutation` check demonstrates wording robustness within it. Each corpus pools molecular subtypes into one ranking (e.g. triple-negative/HER2-positive breast cancer not distinguished), so top-10 lists are subtype-blended, and the 19 corpora are not fully independent by construction (e.g. 1,388 papers, 5.2% of the Myeloma corpus, are also indexed under Leukemia via shared MeSH hierarchy, 365 of those carrying the Plasma Cell Leukemia term specifically), though this did not inflate the two cancers' ranking similarity (Jaccard 0.111, below the 0.189 mean, Results 5.1). OncoDigger's corpora are PubMed abstracts only, unlike the full-text-indexing comparators (Table 1), and the live OncoDigger website differs from this frozen release.

**CancerMine vintage and dependency.** CancerMine's fixed 2023-03-01 snapshot (Zenodo 7689627) cannot incorporate literature after that date, unlike OncoDigger's through September 2026 (effect not separately quantified). IntOGen's own post-processing also removes candidate non-Tier-1-CGC genes lacking a CancerMine literature reference, so the CancerMine comparison is not fully independent outside CGC Tier 1 (intogen-plus.readthedocs.io, "Post-processing").

**IntOGen's own scope limits.** IntOGen detects somatic point-mutation signal only, so fusion- or germline-driven genes can be structurally undetectable (Methods 2.3, Results 5.3). Its methods are also calibrated using CGC gene membership, with a structurally lower significance bar for CGC genes than non-CGC genes (intogen-plus.readthedocs.io, `bbglab/intogen-plus` tag `v2024`), so IntOGen is independent in its evidence but not fully independent in its method.

**Capitalization filter and article-title conventions (Methods 3.3).** The capitalization-based false-positive filter treats non-sentence-initial capitalization as a sign of genuine gene-symbol usage, but journal titles are often Title Cased regardless of meaning, so a common word capitalized only by title-casing convention can pass as genuine when a human reader would recognize it as the ordinary word. A check across all 19 corpora and 15 symbols found a small effect, 0% for most, up to about 6% for REST (11-12 of 194 retained postings, versus 3-7 of 4,505 for MET). This real, previously undocumented blind spot is not corrected here and does not affect the manually-read samples (Results 5.2, Supplementary Methods S1), since a human reader is not fooled by title casing.

**Post hoc decomposition scope (Results 5.8).** This decomposition is exploratory, not confirmatory. The two full-corpus/real-formula cells were not evaluated (Methods 3.4), so results are two separate experiments, not one three-way design. The 10-alias mechanism was identified by inspecting results rather than pre-registered. Rank-weighting's own positional vs. frequency sub-components were not further decomposed. And all of it used a single query and discovery mode, not a general claim about every OncoDigger mode.

---

## Data and code availability

Most results, figures and rankings in this paper were generated with a frozen internal pipeline snapshot (git tag `manuscript1.1-v1.0-vocab-audit-2026-10-04`, retained only in the author's private development history, not part of the public release below). The refined MET/REST/FH vocabulary-collision rules (Methods 3.3) and the post hoc decomposition (Methods 3.4, Results 5.8) were both produced after this snapshot. The companion repository's `v1.1` release, linked below, covers both and was confirmed by clean-checkout reproduction to reproduce every IntOGen-based number this paper reports (Implementation).

github.com/jfataphd/OncoDigger-mutation-paper1, release `v1.1`, archived at Zenodo: https://doi.org/10.5281/zenodo.23269250

The repository contains the corrected gene vocabulary, quality-control rules, ranking code, settings and derived result tables used here.

The continuously developed web interface is available at the following address.

`oncodigger.com`

The live site may differ from the frozen manuscript pipeline and should not be used to reproduce the publication results unless it has been updated to the same release.

The IntOGen driver gene compendium used throughout this paper, for both Precision@10 and exact cancer-type driver concordance, is release 2024.09.20 (`Compendium_Cancer_Genes.tsv`, Martínez-Jiménez et al., 2020), licensed CC0 (public domain) and included in the companion repository at `data/external/intogen/` with no access restrictions. Every derived result in this paper (rankings, precision scores, cancer-type mappings) is fully available in the repository.

---

## CRediT author contribution statement

**Jimmie E. Fata.** Conceptualization, Methodology, Software, Validation, Formal analysis, Investigation, Resources, Data Curation, Writing – Original Draft, Writing – Review & Editing, Visualization, Supervision, Project administration.

## Declaration of competing interests

The author is the creator and developer of OncoDigger, the software evaluated in this study, and declares no financial competing interests.

## Funding

This work received no specific funding from any funding agency in the public, commercial, or not-for-profit sectors.

## Use of generative AI

OncoDigger and the study reported here were conceived and developed by Jimmie E. Fata. The author defined the biological questions, research objectives and overall direction of the study, selected the reference standard and comparison framework, developed OncoDigger, determined the scientific content and final presentation of the figures, interpreted the results and determined the scientific conclusions. Large language models, principally Claude (Anthropic) and ChatGPT (OpenAI), were used extensively as computational, statistical and writing assistants throughout the project. In addition to Python code generation and debugging, these tools played a substantial role in proposing, developing and refining the statistical and computational analyses used to evaluate OncoDigger. The author directed this iterative process by defining the questions to be addressed, evaluating proposed approaches, requesting additional controls and validation where needed, and deciding which analyses and interpretations were scientifically appropriate. The tools also assisted with implementation of the resulting analysis and validation scripts, figure-generation code, and drafting and editing manuscript passages. The computational and statistical scope of the study in its present form would not have been feasible for the author without substantial AI assistance. All AI-assisted analyses, code, figures, numerical results and manuscript text were reviewed by the author, with key findings checked against the underlying data, source code and primary sources. The analysis, validation, statistical and figure-generation scripts underlying the reported results are provided with the companion repository to support direct inspection and reproducibility. The author made all final scientific and editorial decisions and assumes full responsibility for the accuracy, integrity and reproducibility of the work.

---

## References

Bryant HE, Schultz N, Thomas HD, et al. Specific killing of BRCA2-deficient tumours with inhibitors of poly(ADP-ribose) polymerase. *Nature.* 2005;434(7035):913–917. PMID: 15829966.

Cancer Genome Atlas Research Network. Comprehensive molecular characterization of clear cell renal cell carcinoma. *Nature.* 2013;499(7456):43–49. PMID: 23792563.

Gorlov IP, Pikielny CW, Frost HR, et al. Gene characteristics predicting missense, nonsense and frameshift mutations in tumor samples. *BMC Bioinformatics.* 2018;19(1):430. PMID: 30453881.

Grissa D, Junge A, Oprea TI, Jensen LJ. Diseases 2.0: a weekly updated database of disease–gene associations from text mining and data integration. *Database (Oxford).* 2022;2022:baac019. PMID: 35348648.

Hegi ME, Diserens AC, Gorlia T, et al. MGMT gene silencing and benefit from temozolomide in glioblastoma. *N Engl J Med.* 2005;352(10):997–1003. PMID: 15758010.

Idos G, Hampel H, Valle L. Lynch Syndrome. 2004 Feb 5 [Updated 2026 Sep 17]. In: Adam MP, Feldman J, Mirzaa GM, et al., editors. *GeneReviews.* Seattle (WA): University of Washington; 1993–2026.

Kamihara J, Schultz KA, Rana HQ. FH Tumor Predisposition Syndrome. 2006 Jul 31 [Updated 2025 May 8]. In: Adam MP, Feldman J, Mirzaa GM, et al., editors. *GeneReviews.* Seattle (WA): University of Washington; 1993–2026.

Kuehl WM, Bergsagel PL. Early genetic events provide the basis for a clinical classification of multiple myeloma. *Hematology Am Soc Hematol Educ Program.* 2005:346–352. PMID: 16304402.

Lai PT, Wei CH, Luo L, Chen Q, Lu Z. BioREx: Improving Biomedical Relation Extraction by Leveraging Heterogeneous Datasets. *J Biomed Inform.* 2023;146:104487. PMID: 37673376.

Lawrence MS, Stojanov P, Polak P, et al. Mutational heterogeneity in cancer and the search for new cancer-associated genes. *Nature.* 2013;499(7457):214–218. PMID: 23770567.

Lever J, Zhao EY, Grewal J, Jones MR, Jones SJM. CancerMine: a literature-mined resource for drivers, oncogenes and tumor suppressors in cancer. *Nat Methods.* 2019a;16(6):505–507.

Lever J, Jones MR, Danos AM, et al. Text-mining clinically relevant cancer biomarkers for curation into the CIViC database. *Genome Med.* 2019b;11(1):78.

Martínez-Jiménez F, Muiños F, Sentís I, et al. A compendium of mutational cancer driver genes. *Nat Rev Cancer.* 2020;20(10):555–572. PMID: 32778778.

O'Connell FP, Pinkus JL, Pinkus GS. CD138 (Syndecan-1), a Plasma Cell Marker: Immunohistochemical Profile in Hematopoietic and Nonhematopoietic Neoplasms. *Am J Clin Pathol.* 2004;121(2):254–263.

Piñero J, Bravo À, Queralt-Rosinach N, et al. DisGeNET: a comprehensive platform integrating information on human disease-associated genes and variants. *Nucleic Acids Res.* 2017;45(D1):D833–D839. PMID: 27924018.

Pletscher-Frankild S, Pallejà A, Tsafou K, Binder JX, Jensen LJ. DISEASES: text mining and data integration of disease-gene associations. *Methods.* 2015;74:83–89. PMID: 25484339.

Pritchard CC, Mateo J, Walsh MF, et al. Inherited DNA-repair gene mutations in men with metastatic prostate cancer. *N Engl J Med.* 2016;375(5):443–453. PMID: 27433846.

Robertson S, Zaragoza H. The Probabilistic Relevance Framework: BM25 and Beyond. *Found Trends Inf Retr.* 2009;3(4):333–389. DOI: 10.1561/1500000019.

Seal RL, Braschi B, Gray K, McClay J, Tweedie S, Bruford EA. Genenames.org: the HGNC and PGNC resources in 2026. *Nucleic Acids Res.* 2026;54(D1):D1098–D1107. PMID: 41287213.

Topalian SL, Hodi FS, Brahmer JR, et al. Safety, Activity, and Immune Correlates of Anti–PD-1 Antibody in Cancer. *N Engl J Med.* 2012;366(26):2443–2454. PMID: 22658127.

Wei CH, Allot A, Lai PT, et al. PubTator 3.0: an AI-powered literature resource for unlocking biomedical knowledge. *Nucleic Acids Res.* 2024;52(W1):W540–W546. PMID: 38572754.

---

# SUPPLEMENTARY MATERIAL

<!--
Everything below this line is supplementary material. It preserves the detailed numerical results,
mappings, audit information and figure legends from the source manuscript so they are not lost
when the main article is compressed to Bioinformatics length.
-->

## Supplementary Methods S1. Full collection characteristics

**Supplementary Table S1. Cancer-specific literature collections.** Coverage is the percentage of each collection represented in the 1,000-paper scoring pool. Vocabulary size is the number of distinct word tokens indexed for that collection, and abstract length is the mean number of words per abstract.

| Cancer | Papers | Coverage @ 1,000 | Vocabulary size | Mean abstract length (words) |
|---|---:|---:|---:|---:|
| Breast Cancer | 250,061 | 0.40% | 427,569 | 226 |
| Lung Cancer | 178,908 | 0.56% | 339,571 | 225 |
| Colorectal Cancer | 168,174 | 0.59% | 327,311 | 229 |
| Liver Cancer | 127,252 | 0.79% | 273,114 | 223 |
| Prostate Cancer | 110,937 | 0.90% | 239,258 | 235 |
| Leukemia | 106,922 | 0.94% | 232,482 | 201 |
| Brain Cancer | 103,180 | 0.97% | 220,319 | 216 |
| Lymphoma | 79,022 | 1.27% | 171,981 | 202 |
| Melanoma | 66,447 | 1.50% | 184,880 | 209 |
| Pancreatic Cancer | 65,223 | 1.53% | 168,545 | 218 |
| Stomach Cancer | 64,372 | 1.55% | 163,064 | 225 |
| Ovarian Cancer | 61,538 | 1.63% | 164,033 | 221 |
| Cervical Cancer | 50,875 | 1.97% | 137,776 | 234 |
| Kidney Cancer | 46,739 | 2.14% | 125,170 | 219 |
| Thyroid Cancer | 38,640 | 2.59% | 105,443 | 228 |
| Esophageal Cancer | 37,925 | 2.64% | 108,665 | 228 |
| Bladder Cancer | 35,620 | 2.81% | 110,528 | 226 |
| Myeloma | 26,591 | 3.76% | 88,809 | 209 |
| Endometrial Cancer | 22,079 | 4.53% | 75,563 | 231 |

**Supplementary Table S2. Exact corpus construction queries (Methods 2.1).** Each collection was built via the NCBI PubMed E-utilities API using the query `(<MeSH term>[mh]) AND english[la] AND hasabstract`, restricted to publication dates 2000/01/01-2026/09/27, with the cancer-specific MeSH term below. Unlike the PubTator subtype-mapping limitation corrected in Methods 2.4, this single-term query does not under-count subtype-specific literature. PubMed's `[mh]` tag explodes by default to every narrower term in the MeSH tree (nlm.nih.gov, "Explosion," Foundations of MeSH in MEDLINE, e.g. `liver neoplasms[mh]` already includes `carcinoma, hepatocellular` and other child terms, confirmed via NCBI E-utilities, 217,667 results with explosion versus 199,945 with `[mh:noexp]`), whereas PubTator's bulk relation file carries only the specific disease identifier each paper was tagged with and has no equivalent automatic hierarchy at query time. All 19 collections were retrieved and incrementally refreshed from the same live corpus infrastructure (`data/external/PROVENANCE.md` documents the equivalent identity anchors for external comparator datasets, and the per-corpus retrieval manifests, including PMID lists, deduplication counts and exact timestamps, are retained by the corresponding author and available on request).

| Cancer | MeSH term |
|---|---|
| Bladder Cancer | urinary bladder neoplasms |
| Brain Cancer | brain neoplasms |
| Breast Cancer | breast neoplasms |
| Cervical Cancer | uterine cervical neoplasms |
| Colorectal Cancer | colorectal neoplasms |
| Endometrial Cancer | endometrial neoplasms |
| Esophageal Cancer | esophageal neoplasms |
| Kidney Cancer | kidney neoplasms |
| Leukemia | leukemia |
| Liver Cancer | liver neoplasms |
| Lung Cancer | lung neoplasms |
| Lymphoma | lymphoma |
| Melanoma | melanoma |
| Myeloma | multiple myeloma |
| Ovarian Cancer | ovarian neoplasms |
| Pancreatic Cancer | pancreatic neoplasms |
| Prostate Cancer | prostatic neoplasms |
| Stomach Cancer | stomach neoplasms |
| Thyroid Cancer | thyroid neoplasms |

**Pool-size sensitivity (Methods 2.1).** Relative to the 1,000-paper analysis pool used throughout this paper, mean top-10 Jaccard similarity was 0.722 at a 250-paper pool, 0.943 at 750 papers and 0.852 at 2,000 papers. Adjacent-pool comparisons (250/500, 500/750, 750/1,000, 1,000/1,500, 1,500/2,000) gave 0.835, 0.817, 0.943, 0.914 and 0.923 respectively (mean 0.886). Both views support 1,000 papers as a stable operating point rather than an isolated or cherry-picked choice.

**MET positive gene-context terms (Methods 3.3).** The pre-identified positive gene-context terms used to classify MET's 40-posting robustness sample, and to build the stricter 1,300-posting condition, are frozen in `scripts/analysis/met_disambiguation_check.py`'s `POSITIVE_TERMS` list, namely "c-met," "proto-oncogene," "protooncogene," "receptor tyrosine kinase," "hepatocyte growth factor," "hgf," "amplification," "exon 14," "crizotinib," "capmatinib," "tepotinib," "savolitinib," "kinase inhibitor," "juxtamembrane," and the compounds "met receptor/gene/oncogene/mutation/inhibitor/amplif/kinase/signaling/signalling/pathway/alterations," plus "copy number." A posting was classified "genuine" if it matched one of these terms and none of the negative-phrase patterns (Methods 3.3). This is the same rule used to build the 1,300-posting positive-context-only condition.

**Vocabulary exclusion and term-level suppression list (Methods 3.3).** Seven symbols producing unresolved common-word or abbreviation collisions with no viable safer alias (MICE, DCR, PC, C2, FBN1, SGCG, WDHD1) are excluded from all rankings entirely (`GENE_EXCLUDE`). A separate set of 13 symbols is instead handled by term-level suppression (`suppress_ambiguous_symbols()`), under which only the bare, ambiguous symbol is removed as a matchable term, while the gene's full approved name and every other alias remain matchable under the same identity. This paper's audit was subsequently adopted into the live application, which now shows users side-by-side Raw and Vetted gene results, the live site's own analogue of this paper's true-raw versus cleaned vocabulary comparison. OncoDigger's `gene_quality.py` policy version 2026-10-05-v2 implements this 13-symbol suppression set, the four safer-alias substitutions below, and the MET/FH negative-context corrections (Methods 3.3) without modification, as part of a broader live policy that additionally suppresses further symbols not evaluated in this paper. Its `GENE_EXCLUDE` additionally blanket-excludes REST, whereas this paper corrected REST via a Thyroid-specific postings removal rather than a global exclusion (Results 5.2). Four of the 13 were originally identified and are additionally indexed under a distinct safer alias, namely GC (colliding with "gastric cancer") under VDBP, HR ("hazard ratio") under KDM3D, HCCS ("hepatocellular carcinoma") under CCHL, and BPIFA4P under LATH. The remaining nine were found by a systematic post-hoc sweep of every gene symbol of four characters or fewer appearing anywhere in the top-50 rankings (66 candidates), each verified against live PubMed titles and abstracts. MB (myoglobin) matches MDA-MB-231/MDA-MB-436 cell lines, "mutations per megabase," and "medulloblastoma." ELN (elastin) matches "European LeukemiaNet." CP (ceruloplasmin) matches "chronic phase" and "chronic pancreatitis." UBC (ubiquitin C) matches "urothelial bladder cancer." CAMP (cathelicidin antimicrobial peptide) matches "cyclic AMP." EFS (embryonal Fyn-associated substrate) matches "event-free survival." IVD (isovaleryl-CoA dehydrogenase) matches "in vitro diagnostic." SCT (secretin) matches "solid and cystic tumor" of the pancreas. SON matches the English word "son." Two further candidates surfaced by the same sweep, INS and TPO, were checked against live PubMed but not confirmed as false aliases with comparable clarity, and were left unmodified in the vocabulary rather than acted on with inconclusive evidence. Twenty symbols are therefore handled in total, 7 excluded and 13 term-level-suppressed.

**Downstream impact of the nine newly-found false aliases (Methods 3.3, Results 5.1-5.3).** This vocabulary correction was made before this paper's conversion to IntOGen as the primary reference. All numbers reported anywhere in this paper reflect the final, already-corrected pipeline; no historical pre-correction figures are reported. These nine were first corrected via a blunt, symbol-level exclusion (matching the treatment already used for the original 7), then this was itself refined to the term-level suppression policy described above, since blunt exclusion treated these nine more harshly than the four originally-identified symbols. SCT alone among the nine recovers a genuine match once only its bare symbol is suppressed rather than the whole gene excluded, since "secretin" independently earns Pancreatic Cancer rank 48, displacing SMARCB1 (annotated for malignant rhabdoid tumour, not pancreatic cancer) from that cancer's top 50.

**Vendored engine source (Methods 3.1).** The "Primary Evidence" content-type multiplier is a term-and-entity-weighted function, not a single fixed value, defined in the OncoDigger application source (`src/oncodigger/search/discovery_modes.py`, commit `4e39fc068102`, 2026-05-31, of the author's private live-application repository, not independently resolvable by readers). Because that repository continues to be actively developed past this commit and is not public, an exact copy of this file and the other engine source files this pipeline depends on (BM25 retrieval, entity aggregation, HGNC vocabulary loading) is vendored in this repository at `reference_engine_code/`, the actual citable and reproducible source for this paper's engine code.

## Supplementary Methods S2. Benchmark-audit details

**Supplementary Table S3. How each comparison was audited.**

| Resource | Scoring used here | Audit | Finding | Correction / interpretation | Final result |
|---|---|---|---|---|---|
| CancerMine | Citation counts summed across matched disease terms and three roles | Disease terms checked; supporting sentences inspected for suspicious PTGDR result; cross-resource crosswalk audit (Results 5.4) checked for missed subtype disease terms | 3 disease-term mapping mistakes; all 1,329 PTGDR-supporting sentences were unrelated `AS1` lncRNA mentions; Brain Cancer's mapping was separately found to be missing ependymoma/medulloblastoma terms present in CancerMine's own data | 3 terms reassigned; PTGDR excluded; Brain Cancer mapping extended to include the missing subtype terms | 84.2% before correction; 87.9% after correction (against IntOGen) |
| DISEASES 2.0 | Maximum matched association z-score per cancer | Overlapping disease terms and ontology definitions checked | Gallbladder/bladder substring collision; 32 `[organ] lymphoma`/`[organ] melanoma` routing issues in the benchmark mapping; unresolved MYOM2/LOC102723407 pattern | Mapping rules corrected; unresolved pattern retained as an open flag | 60.0% (114/190) |
| PubTator 3.0 | Summed distinct-PMID Gene-Disease relations across parent + subtype MeSH entities | Single-parent-ID mapping found to miss subtype-tagged relations; Liver Cancer bucket (single ambiguous parent code) dominated by ALS-related genes and papers | Single-parent mapping omitted subtype-tagged evidence and produced an ALS-contaminated Liver Cancer ranking | Parent + subtype mapping restored subtype evidence and eliminated the Liver/ALS artifact; Liver Cancer included without exclusion once queried via its own specific subtypes | 67.9% (129/190) after correction; the exact pre-correction percentage was not separately preserved and is not restated here to avoid citing an unverified figure. A live reconstruction of the original single-parent-MeSH methodology against PubTator's current API (`pubtator_v1_single_parent_reconstruction.py`) returned 67.9% (129/190), indistinguishable from the corrected value rather than confirming the originally reported figure, indicating PubTator's underlying relations data has changed since the original benchmark window; the pre-correction value is omitted from Figure 4C for this reason |
| ChatGPT free tier | Single generated answer | Prompt, response and access conditions saved | Not a bulk/versioned dataset; the exact underlying free-tier model could not be determined | No correction; practical baseline only | 78.9% gene-level Precision@10 (145/190 = 76.3% on exact cancer-type concordance, Results 5.6) |

**CancerMine deduplication check (Methods 2.4).** CancerMine's citation-sum scoring was checked against a fully deduplicated alternative (distinct supporting PMIDs per gene per cancer, from the sentence-level evidence file, Zenodo 7689627). Deduplication changed almost nothing. Top-10 membership changed in only 4/19 cancers, pooled Precision@10 was nearly identical (167/190=87.9% citation-sum vs. 166/190=87.4% unique-PMID, `scripts/analysis/intogen_cancermine_double_counting_check.py`), and exact cancer-type concordance likewise (114/190=60.0% vs. 113/190=59.5%, `scripts/analysis/intogen_cancermine_dedup_exact_concordance_check.py`). Dual-diagnosis terms (e.g. adult T-cell leukemia/lymphoma) were allowed to count toward both relevant cancers (Supplementary Table S4). OncoDigger and CancerMine were compared per cancer using a paired t-test, Wilcoxon signed-rank test, and paired 2,000-resample case bootstrap, with effect sizes and CIs treated as the primary result.

## Supplementary Methods S3. CancerMine disease-name mapping

The complete unabridged 251 cancer-term pairs are stored with the manuscript data, summarised below.

**Supplementary Table S4. CancerMine cancer-to-disease-term mapping summary.**

| Cancer | Name pattern(s) | Matched terms | Audit correction / reassignment |
|---|---|---:|---|
| Breast Cancer | breast | 15 | 0 |
| Lung Cancer | lung | 12 | 0 |
| Colorectal Cancer | colorectal, colon, rectal, rectum | 8 | 0 |
| Prostate Cancer | prostate | 4 | 0 |
| Melanoma | melanoma | 15 | 0 |
| Bladder Cancer | bladder | 9 | 0 |
| Kidney Cancer | kidney, renal | 16 | 0 |
| Pancreatic Cancer | pancrea | 8 | 0 |
| Liver Cancer | liver, hepato | 8 | 1 liver-lymphoma term reassigned to Lymphoma |
| Stomach Cancer | stomach, gastric | 7 | 0 |
| Esophageal Cancer | esophag, gastroesophageal | 7 | 0 |
| Ovarian Cancer | ovar | 20 | 1 ovarian-lymphoma term reassigned to Lymphoma |
| Endometrial Cancer | endometri, uterine corpus | 13 | 1 ovarian-cancer subtype reassigned to Ovarian Cancer |
| Cervical Cancer | cervi | 5 | 0 |
| Thyroid Cancer | thyroid | 9 | 0 |
| Brain Cancer | brain, glioma, glioblastoma, astrocytoma | 16 | 0 |
| Leukemia | leukemia | 44 | 0 |
| Lymphoma | lymphoma | 34 | 0 |
| Myeloma | myeloma | 1 | 0 |

## Supplementary Methods S4. PubTator disease mapping

An initial version of this table mapped each cancer to a single parent MeSH descriptor only (the column repeated below as "Parent"). This was found, during construction of the benchmark, to systematically under-count PubTator's own relation data for any cancer with well-known named subtypes. PubTator 3.0's entity recognition and identifier mapping (Wei et al., 2024), upstream of its BioREx relation-extraction model (Lai et al., 2023), tags subtype mentions (e.g. "glioblastoma") under their own distinct MeSH identifier, separate from the broad parent term, so a parent-only query never retrieves them. Verified directly against PubTator3's public API (2026-09-30), EGFR alone carries 1,604 distinct-PMID relations under Glioblastoma's own identifier (D005909) and FLT3 carries 3,199 under Acute Myeloid Leukemia's (D015470), none reachable through the parent codes Brain Neoplasms (D001932) or Leukemia (D007938) alone. Every cancer below was therefore remapped to its parent descriptor plus every matched subtype entity in PubTator's own disease vocabulary, found via PubTator's entity-autocomplete API using the same subtype terminology already curated for the cancer-type mapping (`cancer_aware_precision.py`'s `CANCER_TERMS`).

**Supplementary Table S5. PubTator 3.0 disease-entity mapping used in the benchmark.**

| Cancer | Parent MeSH | Entities used in the benchmark |
|---|---|---|
| Breast Cancer | D001943 | Breast Neoplasms (D001943) |
| Lung Cancer | D008175 | Lung Neoplasms (D008175), Non-Small Cell Lung Carcinoma (D002289), Lung Adenocarcinoma (D000077192) |
| Colorectal Cancer | D015179 | Colorectal Neoplasms (D015179), Colonic Neoplasms (D003110), Rectal Neoplasms (D012004) |
| Prostate Cancer | D011471 | Prostatic Neoplasms (D011471) |
| Melanoma | D008545 | Melanoma (D008545), Uveal Melanoma (C536494) |
| Bladder Cancer | D001749 | Urinary Bladder Neoplasms (D001749) |
| Kidney Cancer | D007680 | Kidney Neoplasms (D007680), Renal Cell Carcinoma (D002292), Wilms Tumor (D009396) |
| Pancreatic Cancer | D010190 | Pancreatic Neoplasms (D010190) |
| Liver Cancer | D008113 | Hepatocellular Carcinoma (D006528), Cholangiocarcinoma (D018281) (the ambiguous parent D008113 itself was not queried; see Results 5.4) |
| Stomach Cancer | D013274 | Stomach Neoplasms (D013274) |
| Esophageal Cancer | D004938 | Esophageal Neoplasms (D004938) |
| Ovarian Cancer | D010051 | Ovarian Neoplasms (D010051) |
| Endometrial Cancer | D016889 | Endometrial Neoplasms (D016889) |
| Cervical Cancer | D002583 | Uterine Cervical Neoplasms (D002583) |
| Thyroid Cancer | D013964 | Thyroid Neoplasms (D013964) |
| Brain Cancer | D001932 | Brain Neoplasms (D001932), Glioma (D005910), Glioblastoma (D005909), Medulloblastoma (D008527), Astrocytoma (D001254), Ependymoma (D004806), Diffuse Intrinsic Pontine Glioma (D000080443) |
| Leukemia | D007938 | Leukemia (D007938), Acute Myeloid Leukemia (D015470), Precursor Cell Lymphoblastic Leukemia-Lymphoma (D054198), Chronic Myeloid Leukemia (D015464), Chronic Lymphocytic Leukemia (D015451), Myelodysplastic Syndromes (D009190), Acute Promyelocytic Leukemia (D015473) |
| Lymphoma | D008223 | Lymphoma (D008223), Non-Hodgkin Lymphoma (D008228), Diffuse Large B-Cell Lymphoma (D016403), Burkitt Lymphoma (D002051), Mantle Cell Lymphoma (D020522), Marginal Zone B-Cell Lymphoma (D018442) |
| Myeloma | D009101 | Multiple Myeloma (D009101) |

Full per-cancer entity mapping and audit trail are in `data/derived/pubtator_v2_entity_mapping.csv` (`scripts/analysis/pubtator_benchmark_v2_subtype_aggregated.py`).

**Supplementary Table S6. Literature weight of unrecovered annotations (Results 5.5, Methods 2.5).**

*Panel A, all 16 boundary-case pairs (non-top-10 IntOGen driver pairs with supporting-paper count at or above that cancer's weakest top-10 gene).*

| Cancer | Gene | Rank | Supporting papers | Weakest top-10 threshold |
|---|---|---:|---:|---:|
| Leukemia | TET2 | 11 | 86 | 81 |
| Kidney Cancer | MTOR | 16 | 65 | 39 |
| Melanoma | CDK4 | 13 | 61 | 48 |
| Pancreatic Cancer | MEN1 | 11 | 59 | 51 |
| Stomach Cancer | ERBB2 | 11 | 58 | 45 |
| Ovarian Cancer | PTEN | 11 | 57 | 31 |
| Brain Cancer | PTEN | 11 | 55 | 42 |
| Lymphoma | BCL6 | 21 | 54 | 48 |
| Stomach Cancer | ARID1A | 13 | 54 | 45 |
| Melanoma | TP53 | 15 | 48 | 48 |
| Kidney Cancer | MET | 20 | 46 | 39 |
| Kidney Cancer | TSC1 | 15 | 44 | 39 |
| Thyroid Cancer | PTEN | 14 | 40 | 27 |
| Prostate Cancer | FOXA1 | 18 | 38 | 38 |
| Thyroid Cancer | TSHR | 13 | 36 | 27 |
| Ovarian Cancer | PALB2 | 16 | 34 | 31 |

*Panel B, IntOGen driver genes confirmed, by direct inspection of OncoDigger's gene vocabulary lexicon, to have no entry in it at all (structurally unmatchable regardless of literature volume). These are IntOGen's own current HGNC symbols for three genes already identified as vocabulary gaps under their deprecated names in an earlier internal audit of this vocabulary, confirmed here independently via a differently-sourced gene list.*

| Cancer(s) where listed as a driver | Gene (current symbol) | Deprecated symbol (pre-2021 HGNC) | External PubMed records (deprecated symbol) |
|---|---|---|---:|
| Brain Cancer | H3-3A | H3F3A | 196 |
| Brain Cancer | H3C2 | HIST1H3B | 53 |
| Lymphoma | H4C9 | HIST1H4I | 0 |

No Panel C (mutation-mechanism breakdown) is reported here, since IntOGen's own data does not carry a curated mechanism field (Methods 2.5), and no external live-PubMed cross-check was performed for the full never-surfaced set, for the same reason given there.

**Supplementary Table S7. Individual audit of the 13 gene-cancer pairs not detected as a driver anywhere by IntOGen (Results 5.3).** Established-role sources are Lynch syndrome (Idos et al., 2004), FH/HLRCC (Kamihara et al., 2006), IGH translocations in myeloma (Kuehl and Bergsagel, 2005), MGMT/temozolomide (Hegi et al., 2005), PARP1 synthetic lethality (Bryant et al., 2005), PDCD1/PD-1 (Topalian et al., 2012), SDC1/CD138 (O'Connell et al., 2004).

| Cancer | Gene | Mechanism category | Established role |
|---|---|---|---|
| Kidney Cancer | FH | Germline predisposition | Causes hereditary leiomyomatosis and renal cell cancer (HLRCC) |
| Colorectal Cancer | MLH1 | Germline predisposition | Lynch syndrome mismatch-repair gene |
| Endometrial Cancer | MLH1 | Germline predisposition | Lynch syndrome mismatch-repair gene |
| Ovarian Cancer | MLH1 | Germline predisposition | Lynch syndrome mismatch-repair gene |
| Colorectal Cancer | MSH2 | Germline predisposition | Lynch syndrome mismatch-repair gene |
| Endometrial Cancer | MSH2 | Germline predisposition | Lynch syndrome mismatch-repair gene |
| Colorectal Cancer | MSH6 | Germline predisposition | Lynch syndrome mismatch-repair gene |
| Endometrial Cancer | MSH6 | Germline predisposition | Lynch syndrome mismatch-repair gene |
| Myeloma | IGH | Translocation/fusion partner | Recurrent partner locus in myeloma translocations (e.g. t(4;14), t(11;14)), not a point-mutated gene |
| Brain Cancer | MGMT | Epigenetic biomarker | Promoter methylation status predicts temozolomide response |
| Prostate Cancer | PARP1 | Therapeutic biomarker | Synthetic-lethality drug target discussed alongside BRCA mutation status |
| Melanoma | PDCD1 | Therapeutic biomarker | PD-1; target of anti-PD-1 checkpoint immunotherapy |
| Myeloma | SDC1 | Diagnostic biomarker | CD138; immunohistochemical marker of plasmacytic differentiation |

All 9 unique genes (13 gene-cancer pairs) were individually checked against their established clinical or biological role rather than assumed to fit a general pattern, and none has a primary somatic point-mutation driver mechanism, consistent with IntOGen's disclosed scope limit (Methods 2.3).

## Supplementary Methods S5. IntOGen cancer-type mapping

IntOGen's release 2024.09.20 (`Compendium_Cancer_Genes.tsv`, CC0 licence) assigns each driver gene to the specific cohort(s) in which it reached significance, each cohort labelled with one of 86 cancer-type codes, combined across its seven driver-detection methods into one significance call per gene per cohort (Martínez-Jiménez et al., 2020). These codes were mapped to the study's 19 cancer categories conservatively. A code was included only when it corresponded unambiguously to a single study cancer, and codes that could plausibly span more than one study cancer, or that duplicated a more specific code already included, were excluded rather than guessed at (full code list and exclusion reasons in Supplementary Table S8, scoring script `scripts/analysis/intogen_full_rescoring.py`). A gene-cancer pair counted as an exact match only when IntOGen listed the gene as a significant driver specifically within a cohort mapped to that exact cancer. This mapping and matching rule is the single source used for every exact-concordance and recall figure reported for every system in this paper (OncoDigger, the pooled control, CancerMine and the ChatGPT baseline).

**Supplementary Table S8. IntOGen cancer-type codes mapped to each study cancer, and codes deliberately excluded as ambiguous.**

| Study cancer | IntOGen codes included |
|---|---|
| Breast Cancer | BRCA |
| Lung Cancer | LUAD, LUSC, NSCLC, SCLC |
| Colorectal Cancer | COAD, READ, COADREAD |
| Prostate Cancer | PRAD |
| Melanoma | MEL, SKCM, UM |
| Bladder Cancer | BLCA, UTUC |
| Kidney Cancer | CCRCC, CHRCC, PRCC, RCC, WT |
| Pancreatic Cancer | PAAD, PANET |
| Liver Cancer | HCC, CHOL, LIHB |
| Stomach Cancer | STAD |
| Esophageal Cancer | ESCA, ESCC |
| Ovarian Cancer | OVT |
| Endometrial Cancer | UCEC, UCS |
| Cervical Cancer | CESC, CEAD |
| Thyroid Cancer | WDTC |
| Brain Cancer | GB, GBM, HGGNOS, LGGNOS, PAST, MBL, EPM, ATRT |
| Leukemia | ALL, AML, CLLSLL, MDS |
| Lymphoma | BL, DLBCLNOS, NHL, MLYM |
| Myeloma | PCM |

Of the 86 codes present in the data, 50 are included above. The remaining 36 were excluded, in two distinct categories. Eight are ambiguous or duplicative of a more specific code already included, namely EGC (spans Stomach and Esophageal Cancer), LNM (spans Leukemia, Lymphoma and Myeloma), and the generic LUNG, PANCREAS, STOMACH, BLADDER, PROSTATE and SKIN buckets. The remaining 28 denote cancer types with no corresponding study cancer at all, rather than an ambiguity within this study's scope, namely ACC, ACYC, ANGS, ANSC, BCC, CSCC, ES, GBC, GIST, HNSC, LIPO, LMS, MGCT, MT, NBL, NETNOS, NPC, OS, PGNG, PLMESO, RBL, RMS, SACA, SARCNOS, SIC, SOFT_TISSUE, THYM and VULVA (sarcomas, gastrointestinal stromal tumours, head-and-neck and non-melanoma skin cancers, and other rare or paediatric tumour types not among this study's 19 cancers).

## Supplementary Results S1. Additional cancer-specificity statistics

**Headline-precision detail (Results 5.2).** A Wilson CI treating all 190 observations as independent gives 88.6-96.0%, narrower than the primary 95% cancer-bootstrap CI (88.9-96.9%) used throughout this paper, since it ignores that genes cluster within 19 cancers rather than spreading across 190 truly independent trials (Methods 2.2). Mean rank AUC's per-cancer range (0.518-0.798) reflects IntOGen's broader, cohort-computed driver catalogue making the top-50 sorting task somewhat harder than a smaller, more selective reference would. The per-cancer Spearman rank-vs-driver rho ranged from -0.030 to -0.446 in absolute value (`scripts/analysis/secondary_stats.py`, `data/derived/suppl_rank_spearman_intogen.csv`). The frequency-only control also counted SHE as an ordinary pronoun, alongside NODAL in `nodal metastasis`, as gene mentions (Methods 2.2's "curated" refers to the gene symbol/alias list itself, not the separate ambiguous-symbol exclusion step applied downstream, Methods 3.3).

**Exact-concordance and recall detail (Results 5.3).** Among the 33 driver-elsewhere non-matches, CD274 (IntOGen-confirmed in Melanoma) and PALB2/CHEK2 (germline-predisposition relevance in a secondary tissue context with weak somatic signal there) illustrate the secondary-tissue pattern. By rank interval, exact matches fall from 75.8% at ranks 1-10 to 53.2% at ranks 11-20 and 32.5% at ranks 21-50. Cervical Cancer has the highest Recall@50 of any cancer (41.7%), consistent with its top-10 driver genes being recovered directly rather than requiring deep-rank search (Results 5.5).

**Literature-weight detail (Results 5.5).** Supporting-paper IQRs by band were 54-145 for top 10 (n=144), 9-26 for ranks 11-50 (n=286), and 2-5 beyond rank 50 (n=324) (Kruskal-Wallis H=583.7, p=1.8x10^-127, Supplementary Table S6). Only 16 of 1,596 non-top-10 pairs (1.0%) matched or exceeded the weakest top-10 gene's supporting-paper count for their cancer, all narrow boundary cases at ranks 11-21 (e.g. TET2 in Leukemia, 86 papers at rank 11 vs. an 81-paper threshold, full list in Supplementary Table S6). Cervical and Esophageal Cancer recover their confirmed drivers directly within the top 10 (8/48 and 6/79 mapped driver genes respectively) at the rates their 80% and 60% exact-concordance figures already report.

NPM1 and FLT3 reached the top 20 only in Leukemia, and CDK12 only in Prostate Cancer, illustrating the pipeline's capture of lineage-specific as well as pan-cancer drivers. The most similar cancer pair by top-10 Jaccard was Liver and Stomach Cancer (0.67), followed by Cervical/Lung, Esophageal/Lung, Colorectal/Liver and Endometrial/Ovarian (0.54 each).

The pooled top 10 was EGFR, KRAS, TP53, PIK3CA, BRAF, POLE, NRAS, BRCA1, ERBB2 and BRCA2, all genuine IntOGen drivers (100% any-cancer Precision@10) but identical for every cancer (zero differentiation, treating Leukemia the same as Prostate Cancer). Because a pooled index has one ranking, its pairwise cancer Jaccard similarity is 1.0 by construction. OncoDigger's corresponding mean pairwise Jaccard was 0.189.

Separation-recovered genes included VHL, PBRM1 and BAP1 for Kidney Cancer, FLT3, NPM1 and DNMT3A for Leukemia, and MYD88, CD79B and EZH2 for Lymphoma, not obscure findings. Across the 91 IntOGen-driver gene-cancer pairs recovered only through corpus separation, median supporting PubMed count was 368 papers (IQR 179-799). 70 of the 104 separation-recovered genes (67.3%) were also confirmed IntOGen drivers for that specific cancer.

The pooled control's inflated any-cancer score traces to famous pan-cancer genes. IntOGen recognises TP53 and KRAS as significant drivers in 19/19 study cancers each, so the pooled list's cancer-type-specificity score is elevated by the breadth of these particular genes rather than by any cancer-type-aware mechanism. A perfect general any-cancer score is achievable simply by returning famous pan-cancer genes everywhere, while genuine cancer-type-specific accuracy still requires separating the literature by cancer.

CancerMine's corrected mean pairwise top-10 Jaccard was 0.207, down from 0.233 before audit correction. DISEASES was 0.160, and PubTator (subtype-aggregated, all 19 cancers) was 0.184. These values show why cancer specificity cannot be inferred from architecture alone, since DISEASES and PubTator were at least as differentiated as OncoDigger (0.189) by raw cross-cancer overlap even though they address broader tasks.

## Supplementary Results S2. Detailed exact-concordance decomposition

Across 190 OncoDigger gene-cancer observations, scored against IntOGen (Methods 2.3), the breakdown is as follows.

- 144/190 (75.8%) were confirmed IntOGen drivers in the exact cancer OncoDigger ranked them in.
- 33/190 (17.4%) were confirmed IntOGen drivers, but in a different one of the 19 study cancers.
- 13/190 (6.8%) were not detected as significant drivers by IntOGen in any of its 86 cancer-type codes.
- Among the 177 observations confirmed as IntOGen drivers somewhere, 144/177 (81.4%) retained exact cancer-type concordance.

The IntOGen non-matches concentrate in genes with a plausible mechanistic reason to evade point-mutation-based detection, namely TERT and EGFR (5 each), CD274 (4), and the Lynch-syndrome mismatch-repair genes MLH1, MSH2 and MSH6 (3, 2 and 2 respectively, germline-predisposition-related). The 13 observations not detected as a driver anywhere by IntOGen (9 unique genes, namely FH, IGH, MGMT, MLH1, MSH2, MSH6, PARP1, PDCD1, SDC1) were each individually checked against their established clinical/biological role, rather than assumed to fit a general pattern. For all 9, literature prominence in these pairings reflects biological relationships not captured by somatic point-mutation driver detection. FH (Kidney Cancer), MLH1, MSH2 and MSH6 (Colorectal, Endometrial and/or Ovarian Cancer) are germline-predisposition genes. MLH1, MSH2 and MSH6 are Lynch syndrome mismatch-repair genes (Idos et al., 2004), and FH causes hereditary leiomyomatosis and renal cell cancer (Kamihara et al., 2006). Because FH is not an IntOGen-confirmed driver in any cohort, it was not among the 23 genes flagged for vocabulary-collision scrutiny by the IntOGen-driver-guided audit (Methods 3.3), and FH is also a common clinical abbreviation ("family history") unrelated to the fumarate hydratase gene, so its supporting evidence was checked directly rather than assumed correct by analogy with the audited genes. A random sample of 40 of FH's 297 raw postings in the Kidney Cancer corpus was read in full. 36 (90.0%) were genuine fumarate hydratase gene references, and 4 (10.0%) were false, all sharing a single pattern, "favorable histology (FH)," a Wilms tumor clinical-staging abbreviation rather than "family history" as initially anticipated. This pattern matched 31 of FH's 297 raw postings corpus-wide. Redacting them directly from the title and abstract text, not only from the retrieval index, and re-running the full retrieval and scoring engine reduced FH's supporting-paper count from 88 to 86 and left Kidney Cancer's top-10 gene list, including FH's own rank-5 position, unchanged (`scripts/analysis/vocab_collision_pipeline.py`, `scripts/analysis/vocab_corrected_rerun.py`, `vocab_corrected_rescoring.py`). IGH (Myeloma) is a translocation partner locus (e.g. t(4;14), t(11;14)), not a point-mutated gene (Kuehl and Bergsagel, 2005). MGMT (Brain Cancer), PARP1 (Prostate Cancer), PDCD1 (Melanoma) and SDC1 (Myeloma) are biomarker- or epigenetically-driven rather than mutation-driven. MGMT promoter methylation predicts temozolomide response (Hegi et al., 2005), PARP1 is a synthetic-lethality drug target discussed alongside BRCA mutation status (Bryant et al., 2005), PDCD1/PD-1 is the target of anti-PD-1 checkpoint immunotherapy (Topalian et al., 2012), and SDC1/CD138 is an immunohistochemical marker of plasmacytic differentiation (O'Connell et al., 2004, full per-gene classification in Supplementary Table S7). These results do **not** show that IntOGen should have detected all 46 pairings. They show that IntOGen's statistical driver-detection framework is a narrower criterion, restricted to somatic point-mutation signal, than the existence of mutation-focused literature for a gene-cancer pairing, related but not equivalent evidence, with its own disclosed scope limits (Methods 2.3).

## Supplementary Results S3. Sensitivity to `driver mutation`

Replacing the primary query with `driver mutation driver mutations` produced mean top-10 Jaccard similarity of 0.895 to the primary rankings (range 0.818–1.0). IntOGen any-cancer precision changed from 93.2% to 95.3%, and exact cancer-type driver concordance changed from 75.8% to 78.4%, a 2.6-percentage-point difference not clearly distinguishable from no difference (paired t-test p=0.056). This supports robustness to wording within the mutation-focused domain only.

## Supplementary Results S4. Detailed comparator statistics

### CancerMine

- Corrected Precision@10 (IntOGen). 87.9% (167/190)
- Unique-PMID deduplicated alternative Precision@10 (post-correction). 87.4% (166/190)
- OncoDigger-CancerMine difference. 5.3 pp
- 95% bootstrap CI. −0.0–10.5 pp
- Paired t-test. p=0.086
- Wilcoxon signed-rank. p=0.196
- OncoDigger higher / tie / lower. 8 / 9 / 2 cancers
- CancerMine exact cancer-type driver concordance (IntOGen). 60.0% (114/190)
- Exact-concordance difference. 15.8 pp
- 95% bootstrap CI. 10.0–21.1 pp
- Paired t-test. p=5.1×10^-5
- Wilcoxon. p<0.001
- Exact-concordance higher / tie / lower. 17 / 0 / 2 cancers
- Retention among CancerMine's own IntOGen any-cancer hits. 114/167 = 68.3%

CancerMine achieved 10/10 confirmed IntOGen drivers in Cervical Cancer, Lung Cancer, Lymphoma and Thyroid Cancer. Its lowest cancer score was 60% in Esophageal Cancer.

Before correction, the CancerMine benchmark was 84.2% against IntOGen. Most of the change to 87.9% resulted from removing PTGDR after the confirmed 1,329-sentence normalization artifact (Results 5.4), and three disease-term assignments were also corrected. A cross-resource crosswalk audit then found Brain Cancer's disease-term mapping was separately missing ependymoma and medulloblastoma terms that CancerMine's own data already carried. Extending the mapping to include them brought two additional genes (GLI1, SHH) into that cancer's top 10 and moved its score from 90% to 80%. CancerMine had 10/10 confirmed IntOGen drivers in 4 cancers, compared with 11 for OncoDigger.

### DISEASES 2.0

- Precision@10 (IntOGen). 60.0% (114/190)
- OncoDigger higher in 17 cancers, tied in 2
- Paired t-test p=5.6×10^-7
- Wilcoxon p=2.7×10^-4
- 95% bootstrap CI. 24.2–41.1 pp
- Mean pairwise top-10 Jaccard. 0.160

Illustrative Melanoma non-driver genes were CTLA4, PDCD1, CD8A, CD4, IFNG, PMEL and IL2, while the confirmed IntOGen drivers were BRAF, CD274 and GNA11. Illustrative Myeloma non-driver genes were CD38, SDC1, CD19, IL6 and CRBN, with MYOM2 and LOC102723407 treated separately as an unresolved pattern. PDCD1/PD-1 is an established immunotherapy target and SDC1/CD138 an established plasma-cell marker, illustrating why a gene outside the reference catalogue does not mean biologically irrelevant.

### PubTator 3.0

- Precision@10 across all 19 cancers (subtype-aggregated, IntOGen). 67.9% (129/190)
- OncoDigger higher in 17, tied in 2
- Paired t-test p=1.2×10^-7
- Wilcoxon p=2.6×10^-4
- 95% bootstrap CI. 19.5–31.1 pp
- Mean pairwise top-10 Jaccard. 0.184 across all 19 cancers

Thyroid Cancer non-driver examples included TG, SLC5A5 and LGALS3 (TSHR is itself a confirmed IntOGen driver in Thyroid Cancer). Lymphoma example KRT20 was supported by papers genuinely concerning lymphoma. Subtype-entity aggregation's own disclosed limitation (Results 5.4) is that TP53 occupied a top-10 slot in 18 of 19 cancers and AKT1 in 17 of 19, since summed relation counts with no background normalisation reward genes broadly co-mentioned across many disease terms, not only cancer-type-specific ones.

Under the earlier, since-corrected single-parent-MeSH-ID mapping, the ambiguous Liver Cancer bucket (`Liver Neoplasms`, D008113) was dominated by SOD1 (397 supporting PMIDs), TARDBP (250), C9ORF72 (230) and FUS (120), an unrelated ALS/neurodegeneration pattern also present in VAPB, TBK1, UBQLN2 and NEFL. Eight top-gene evidence sets were checked at the title/PMID level and confirmed unrelated to liver cancer. This did not recur once Liver Cancer was queried via its own specific subtype entities (Hepatocellular Carcinoma, Cholangiocarcinoma) instead of the single ambiguous parent code.

### ChatGPT practical baseline

- Gene-level Precision@10 (IntOGen). 78.9%
- OncoDigger higher / tie / lower. 12 / 5 / 2 cancers
- Paired t-test p=0.0017
- Wilcoxon p=0.0067
- 95% bootstrap CI. 6.8–21.1 pp
- Exact cancer-type driver concordance (IntOGen). 76.3% (145/190)
- Difference from OncoDigger. −0.5 pp (ChatGPT marginally higher, not significant)
- 95% bootstrap CI. −7.4 to 7.4 pp
- Paired t-test p=0.90
- Wilcoxon p=0.72
- Exact-concordance higher / tie / lower for OncoDigger. 7 / 2 / 10 cancers
- Retention among each system's own IntOGen any-cancer hits. OncoDigger 81.4%, CancerMine 68.3%, ChatGPT 96.7%

TTN appeared in 15/19 ChatGPT lists, MUC16 in 12/19, and LRP1B, RYR2, SYNE1 and FAT4 in 5–6 each. Leukemia and Lymphoma were the two 10/10 ChatGPT cancer lists.

## Supplementary Results S5. Ranking-depth detail

Cumulative mean IntOGen-driver Precision@k is as follows.
- k=5. 96.8% (95% CI 93.7%–100%)
- k=10. 93.2% (95% CI 88.9%–96.9%)
- k=15. 88.8%
- k=20. 84.7%
- k=30. 80.0%
- k=40. 76.7%
- k=50. 73.1%

Non-overlapping rank intervals are as follows.
- ranks 1–10. 177/190 = 93.2%
- ranks 11–20. 145/190 = 76.3%
- ranks 21–30. 134/190 = 70.5%
- ranks 31–50. 238/380 = 62.6%

Rank-membership correlations ranged from -0.030 to -0.446 across cancers. Cervical, Endometrial, Leukemia, Lymphoma and Melanoma remained 100% IntOGen-driver precision through ranks 11–20. Brain Cancer and Breast Cancer fell to 40% in that interval.

These are descriptive rank intervals, not claims of biological novelty or discovery timing.

## Supplementary Results S6. IntOGen driver recall

Precision (Results 5.3) asks whether OncoDigger's top-10 genes are IntOGen drivers specifically confirmed for that cancer. The reverse question, recall, asks how much of IntOGen's complete cancer-type-resolved driver catalogue OncoDigger's ranking actually contains. This uses the same top-50 rankings as the main analysis (`data/canonical/precision_at_k_detail.csv`) and the full IntOGen compendium, mapped to the 19 study cancers using the same conservative mapping as the primary analysis (Methods 2.3, Supplementary Methods S5).

Across the 19 cancers, IntOGen's conservative mapping yields 1,740 gene-cancer driver pairs. Pooled recall increased from 8.3% at rank 10 (144/1,740) to 14.1% at rank 20 (245/1,740) and 24.7% at rank 50 (430/1,740). The unweighted per-cancer macro-average follows a similar pattern, 9.3% at rank 10, 15.5% at rank 20, 26.1% at rank 50 (median 7.6%/14.7%/26.5%). These absolute recall values are modest mainly because the denominator spans the full top-50 depth across every cancer, not because OncoDigger's ranking recovers few genes in absolute terms.

Recall is far from uniform across cancers (Supplementary Table S9). Cervical Cancer has the *highest* Recall@50 of any cancer (41.7%), and Esophageal Cancer reaches 26.6%, both consistent with Results 5.3 and 5.5's finding that these two cancers' driver genes are recovered directly within the top 10 rather than requiring deep-rank search.

Exact matches are nonetheless strongly concentrated toward the top of the ranking rather than spread evenly. Of the 190 available positions, 75.8% are exact matches at ranks 1–10, 53.2% at ranks 11–20, and 32.5% at ranks 21–50. Match density therefore declines with rank depth, consistent with a prioritization ranking rather than a flat or random one.

**Supplementary Table S9. IntOGen driver recall by cancer, ranks 10/20/50 (sorted by Recall@50).**

| Cancer | IntOGen driver pairs mapped | Recall@10 | Recall@20 | Recall@50 |
|---|---:|---:|---:|---:|
| Lung Cancer | 148 | 6.1% | 7.4% | 16.2% |
| Liver Cancer | 142 | 4.9% | 9.9% | 18.3% |
| Breast Cancer | 118 | 7.6% | 11.0% | 19.5% |
| Colorectal Cancer | 96 | 7.3% | 11.5% | 20.8% |
| Prostate Cancer | 95 | 7.4% | 10.5% | 21.1% |
| Brain Cancer | 115 | 7.0% | 9.6% | 21.7% |
| Pancreatic Cancer | 71 | 9.9% | 15.5% | 22.5% |
| Endometrial Cancer | 79 | 7.6% | 16.5% | 25.3% |
| Stomach Cancer | 74 | 9.5% | 16.2% | 25.7% |
| Bladder Cancer | 102 | 6.9% | 14.7% | 26.5% |
| Esophageal Cancer | 79 | 7.6% | 13.9% | 26.6% |
| Kidney Cancer | 82 | 8.5% | 17.1% | 26.8% |
| Lymphoma | 118 | 7.6% | 14.4% | 27.1% |
| Melanoma | 99 | 8.1% | 18.2% | 28.3% |
| Thyroid Cancer | 42 | 19.0% | 23.8% | 28.6% |
| Leukemia | 125 | 7.2% | 14.4% | 31.2% |
| Ovarian Cancer | 54 | 14.8% | 18.5% | 31.5% |
| Myeloma | 53 | 13.2% | 28.3% | 35.8% |
| Cervical Cancer | 48 | 16.7% | 22.9% | 41.7% |

The complete list of IntOGen driver genes not recovered at any rank up to 50, per cancer, is available in `data/derived/intogen_rescoring/intogen_recall_missed_genes.csv`.

## Supplementary Results S7. Genuinely unfiltered vocabulary sensitivity check on the primary pipeline

An earlier version of this check used a "raw, uncleaned vocabulary" baseline that, on closer inspection (Results 5.8, Methods 3.4), still carried the vocabulary-auditing module's own common-word filter and was therefore not genuinely unfiltered. It is replaced here with the actual unfiltered comparison, built from the same vocabulary-auditing module's raw lookup table with every filtering rule removed (69,228 single-token aliases and 11,150 phrase anchors, versus the much smaller filtered vocabulary this paper otherwise uses), run through the identical BM25 retrieval, matching, and scoring pipeline for all 19 cancers (`scripts/analysis/true_raw_vs_cleaned_all19.py`, `data/derived/intogen_rescoring/true_raw_vs_cleaned_all19.csv`).

Pooled Precision@10 against IntOGen fell from 93.2% (177/190) with this paper's vocabulary QC applied to 44.2% (84/190) with the genuinely unfiltered vocabulary, a 48.9-percentage-point drop, present in all 19 cancers individually (Supplementary Table S10). The drop is driven overwhelmingly by a small number of specific alias collisions, not a general property of filtering intensity. Four genes (BRIP1, CD44, SPI1, WWOX) occupy a top-10 slot in every one of the 19 cancers under the unfiltered vocabulary, and a fifth (WAS) in 18 of 19, traced to 10 ambiguous single-token aliases, ordinary English function words and two biomedical abbreviations ("of," "in," "for," "was," "all," "not," "as," "an," "cml," "phl"). Removing only these 10 entries, scored with the real relevance formula rather than the flat mention-count ranking used in the main factorial, recovers 76.3% of the full any-cancer gap and 76.1% of the exact-concordance gap between true-raw (44.2%/27.4%) and cleaned (93.2%/75.8%) vocabulary (81.6%/64.2% with only this minimal correction, BM25 pool, all 19 cancers, Results 5.8).

**Formula-decomposition detail (Results 5.8).** The formula's full exact-concordance figures are -0.5 pp under cleaned vocabulary (93.2%/75.8% vs. 91.6%/76.3%) and +22.1 pp under true-raw vocabulary (44.2%/27.4% vs. 14.7%/5.3%). Decomposing the formula's two remaining multiplicative terms under true-raw vocabulary (rank-weighted mention position and enrichment, supporting-paper-count term held fixed), rank-weighting alone reached 14.2%/4.7% (vs. flat counting's 14.7%/5.3%, essentially no improvement), while enrichment alone reached 67.9%/49.5%. Whether the combined formula's underperformance relative to enrichment alone is caused by rank-weighting's positional component specifically, or by the additional frequency weighting it implicitly reintroduces, was not separately isolated and is left as disclosed future work.

**Supplementary Table S10. Genuinely unfiltered vocabulary sensitivity against IntOGen, primary pipeline, Precision@10 by cancer.**

| Cancer | Cleaned (this paper's pipeline) | True-raw (unfiltered) |
|---|---:|---:|
| Breast Cancer | 100% | 60% |
| Lung Cancer | 100% | 30% |
| Bladder Cancer | 100% | 50% |
| Pancreatic Cancer | 100% | 50% |
| Liver Cancer | 100% | 40% |
| Stomach Cancer | 100% | 40% |
| Esophageal Cancer | 100% | 30% |
| Cervical Cancer | 100% | 30% |
| Thyroid Cancer | 100% | 60% |
| Leukemia | 100% | 50% |
| Lymphoma | 100% | 40% |
| Prostate Cancer | 90% | 60% |
| Melanoma | 90% | 50% |
| Kidney Cancer | 90% | 30% |
| Ovarian Cancer | 90% | 50% |
| Brain Cancer | 90% | 40% |
| Myeloma | 80% | 40% |
| Colorectal Cancer | 70% | 40% |
| Endometrial Cancer | 70% | 50% |

## Supplementary Results S8. Parameter-tuning sensitivity check

Git history for the scoring parameters (BM25 k1 and b, title/MeSH match multipliers, review down-weight, review multiplier, analysis-pool size) in the author's private live-application repository shows all were set once in a single commit (`eb4da47`, 2026-05-04) and never subsequently modified. This is evidence against iterative tuning against validation performance. To test this more directly, the primary "mutation, mutations" query was re-run for all 19 cancers at a range of alternate values for each parameter, varying one parameter at a time while holding all others at their reported value, and pooled Precision@10 against IntOGen's any-cancer driver set was rescored at each setting (`scripts/analysis/intogen_parameter_sensitivity.py`, `data/derived/intogen_rescoring/intogen_parameter_sensitivity.csv`).

Across all 24 unique configurations tested (23 parameter perturbations plus the shared baseline, which is repeated within each parameter block in the table below, giving 29 displayed rows, k1 = 1.0-2.0, b = 0.5-1.0, title-match multiplier = 1.0-1.5, MeSH-match multiplier = 1.0-1.35, review multiplier = 0.5-1.0, analysis-pool size = 500-2,000 papers), pooled Precision@10 ranged only from 91.6% (174/190) to 93.2% (177/190), a spread of 1.6 percentage points around the reported baseline of 93.2% (177/190), with no configuration producing a sharp peak or a value outside this narrow band (Supplementary Table S11). This is inconsistent with the reported result depending sensitively on the exact parameter values used, since a tuned or overfit configuration would be expected to show larger degradation under at least some of these variations.

**Supplementary Table S11. Parameter-sensitivity sweep, pooled Precision@10 against IntOGen (baseline 177/190 = 93.2%).**

| Parameter | Value | Precision@10 |
|---|---:|---:|
| k1 | 1.0 | 92.6% (176/190) |
| k1 | 1.2 | 92.6% (176/190) |
| k1 | 1.5 (baseline) | 93.2% (177/190) |
| k1 | 1.8 | 93.2% (177/190) |
| k1 | 2.0 | 93.2% (177/190) |
| b | 0.5 | 92.6% (176/190) |
| b | 0.65 | 93.2% (177/190) |
| b | 0.75 (baseline) | 93.2% (177/190) |
| b | 0.85 | 93.2% (177/190) |
| b | 1.0 | 93.2% (177/190) |
| title match multiplier | 1.0 | 93.2% (177/190) |
| title match multiplier | 1.1 | 93.2% (177/190) |
| title match multiplier | 1.2 (baseline) | 93.2% (177/190) |
| title match multiplier | 1.3 | 93.2% (177/190) |
| title match multiplier | 1.5 | 92.6% (176/190) |
| MeSH match multiplier | 1.0 | 93.2% (177/190) |
| MeSH match multiplier | 1.05 | 93.2% (177/190) |
| MeSH match multiplier | 1.15 (baseline) | 93.2% (177/190) |
| MeSH match multiplier | 1.25 | 93.2% (177/190) |
| MeSH match multiplier | 1.35 | 93.2% (177/190) |
| review multiplier | 0.5 | 93.2% (177/190) |
| review multiplier | 0.65 (baseline) | 93.2% (177/190) |
| review multiplier | 0.8 | 93.2% (177/190) |
| review multiplier | 1.0 (no down-weight) | 92.6% (176/190) |
| analysis-pool size | 500 | 91.6% (174/190) |
| analysis-pool size | 750 | 91.6% (174/190) |
| analysis-pool size | 1,000 (baseline) | 93.2% (177/190) |
| analysis-pool size | 1,500 | 93.2% (177/190) |
| analysis-pool size | 2,000 | 93.2% (177/190) |

## Supplementary Results S9. IntOGen's own confidence signals for never-surfaced pairs

Results 5.5 shows that most IntOGen driver pairs never entered OncoDigger's ranking at any evaluated depth. To check whether this reflects a uniform sample of IntOGen's driver catalogue or a concentration in its weaker-evidence calls, three signals already present in IntOGen's own compendium were compared between the pairs OncoDigger's ranking surfaced and the pairs it never surfaced. None of these three signals is used anywhere else in this paper's scoring. They are how many of IntOGen's seven driver-detection methods (dNdScv, OncodriveFML, OncodriveCLUSTL, HotMAPS, CBaSE, MutPanning and smRegions) agreed on that gene's significance in that cohort, the combined q-value, and the percentage of the cohort carrying a mutation in that gene (`scripts/analysis/intogen_confidence_breakdown.py`).

Never-surfaced pairs are concentrated in IntOGen's single-method, lower-confidence tier. 76.6% rest on agreement from only one of the seven methods, versus 40.2% of surfaced pairs. Median q-value was 7.9×10^-4 for never-surfaced pairs versus 5.6×10^-9 for surfaced pairs, roughly five orders of magnitude weaker, and median percentage of the cohort mutated was 3.3% versus 6.6%. All three differences are significant in exploratory pair-level Mann-Whitney U tests (p<10^-35 for each, Supplementary Figure S5, Supplementary Table S12), which treat gene-cancer pairs rather than the 19 cancer types as the observation unit and so do not account for within-cancer clustering, unlike this paper's primary cancer-bootstrap analyses (Methods 2.2). This does not establish that any individual never-surfaced gene is a false positive in IntOGen's own catalogue, since IntOGen's own methods already apply multiple-testing correction within each call. It does show that the recall gap is not spread evenly across IntOGen's driver catalogue. It concentrates in the calls IntOGen itself supports most thinly, which is also where large-gene background-mutation-rate artefacts of the kind already discussed for TTN and MUC16 (Results 5.6, Lawrence et al., 2013, and Gorlov et al., 2018) are most likely to occur.

**Supplementary Table S12. IntOGen confidence signals, surfaced vs never-surfaced driver pairs (Results 5.5, Supplementary Results S9).**

| Metric | Surfaced (n=754) | Never surfaced (n=986) | Mann-Whitney p |
|---|---:|---:|---:|
| Median number of IntOGen methods agreeing | 2 | 1 | 1.6×10^-71 |
| Median q-value | 5.6×10^-9 | 7.9×10^-4 | 2.9×10^-73 |
| Median % of cohort mutated | 6.6% | 3.3% | 3.2×10^-36 |
| Share backed by only 1 method | 40.2% | 76.6% | — |

## Supplementary Results S10. Full-corpus presence of never-surfaced pairs

Results 5.5 defines "never surfaced" as not entering OncoDigger's ranking within the standard 1,000-paper, BM25-selected analysis pool used throughout this paper (Methods 2.5). This is narrower than the full corpus. Each cancer's full collection ranges from 22,079 to 250,061 abstracts (Supplementary Table S1), so the analysis pool this check uses covers well under 5% of most corpora. To separate genuine literature absence from an artefact of this analysis window, each never-surfaced pair's raw mention count was queried directly against the full corpus index, bypassing BM25 retrieval and entity-table extraction entirely (`scripts/analysis/intogen_full_corpus_presence_check.py`).

Only 198 pairs (11.4% of all 1,740) have zero lexicon-detected mentions anywhere in the full corpus, meaning this paper's gene vocabulary found no occurrence of that gene's recognized symbols or aliases in any abstract of that cancer's corpus. This is not the same as confirmed biological or literature absence. Detection depends on vocabulary completeness and specificity, and this paper's own vocabulary audit independently confirms both false negatives (H3-3A, H3C2 and H4C9 missing from the vocabulary entirely under their current HGNC symbols, below) and false positives (WAS and MAX's raw full-corpus counts dominated by common-word collisions, below) in this same detection step. These 198 pairs are more precisely described as having no lexicon-detected mentions, not as confirmed absent from the literature. The remaining 788 (45.3% of all 1,740) have at least one lexicon-detected mention somewhere in the corpus but did not reach the 1,000-paper analysis pool. The same vocabulary caveat applies in the other direction here, since some of these 788 detected mentions may themselves be false-positive alias collisions rather than genuine gene mentions. This presence is typically sparse, with a median of 5 mentions across the whole corpus, so most of these pairs would likely remain unranked even under a substantially larger analysis pool. The practical conclusion is not that recall would approach 100% with more compute, but that the raw "never surfaced" figure should not be read as "no literature support." A smaller group has more substantial full-corpus presence (49 pairs with ≥50 mentions, 20 with ≥100), listed in Supplementary Table S13.

Two false-alias collisions were found by this check itself. WAS (Brain Cancer) returned 68,416 mentions, because "was" is a common English auxiliary verb, not a genuine mention of the WAS (Wiskott-Aldrich syndrome) gene. MAX (Brain Cancer) returned 427 mentions, confirmed by direct inspection of a sample of the underlying titles (Supplementary Results S12) to be dominated by "maximum [tolerated dose/uptake/...]" and the unrelated soybean species *Glycine max*, not the MAX proto-oncogene. Both are the same class of vocabulary bug already documented for MB, ELN and seven other symbols (Methods 3.3), not previously caught because this specific check had not been run before. Both are excluded from the magnitude statistics above (but retained in the pair counts, since they do have a nonzero raw count).

**Supplementary Table S13. Full-corpus presence breakdown of never-surfaced pairs (Results 5.5, Supplementary Results S10).**

| Category | n (of 1,740) | % of 1,740 |
|---|---:|---:|
| Surfaced within the standard analysis window | 754 | 43.3% |
| Never surfaced in window, lexicon-detected mentions elsewhere in full corpus | 788 | 45.3% |
| Never surfaced in window, no lexicon-detected mentions in full corpus | 198 | 11.4% |

## Supplementary Results S11. Recall against IntOGen's higher-confidence driver calls

Supplementary Results S9 and S10 establish two independent reasons a pair may never surface without reflecting a ranking failure. It may be backed by only one of IntOGen's seven methods (overwhelmingly a CGC-rescue case, Supplementary Table S12) rather than agreement from at least two methods, or it may have little to no presence in the analysed corpus at all (Supplementary Table S13). Combining both as explicit filters and recomputing recall shows how much of the raw recall figures reported in Results 5.5 and Supplementary Results S6 is attributable to these two factors (`scripts/analysis/intogen_cleaned_recall.py`).

**Supplementary Table S14. Recall under progressively stricter, independently justified filters.**

| Filter | Pairs | Any rank | Rank 10 | Rank 20 | Rank 50 |
|---|---:|---:|---:|---:|---:|
| None (raw, as reported in Results 5.5) | 1,740 | 43.3% | 8.3% | 14.1% | 24.7% |
| Remove pairs with zero full-corpus mentions | 1,542 | 48.9% | 9.3% | 15.9% | 27.9% |
| Remove zero-mention pairs **and** require ≥2 of IntOGen's seven methods | 642 | 70.2% | 19.6% | 30.8% | 48.9% |

Restricting to the 642 pairs that are both detectable in the corpus and supported by at least two of IntOGen's own driver-detection methods, a higher-evidence subset defined for this sensitivity analysis and not IntOGen's own formal tier classification (most pairs this excludes enter IntOGen's compendium through its single-method CGC-rescue route, Supplementary Results S9), raises recall from 24.7% to 48.9% at rank 50, nearly twice the original value, and from 43.3% to 70.2% at any evaluated rank, approximately 62% higher. This does not mean the raw recall figures are incorrect. It means most of the apparent gap concentrates in pairs that are either undetectable in this analysis pool or only weakly, CGC-dependently supported by IntOGen itself, rather than in pairs representing a genuine failure to rank independently well-evidenced drivers. The remaining gap under the strictest filter, 29.8% of these 642 pairs (191) still never surfaced, is addressed directly in Supplementary Results S12.

## Supplementary Results S12. Root-cause classification of the remaining recall gap

The 191 pairs that are both detectable in the full corpus and backed by ≥2 of IntOGen's seven methods, yet still never entered OncoDigger's ranking at any evaluated depth (Supplementary Results S11), were individually classified by directly sampling the underlying corpus text (bypassing BM25 retrieval and ranking entirely) rather than left as an unexplained residual (`scripts/analysis/intogen_miss_root_cause.py`).

For the 31 of these 191 pairs with substantial full-corpus presence (≥20 mentions), up to 20 paper titles mentioning the gene in that cancer's corpus were sampled directly and checked for mutation-focused framing (a title containing "mutat*"). The remaining 160 pairs have fewer than 20 full-corpus mentions and were classified as low-volume without a title sample, consistent with the sparse-presence pattern already established for the broader never-surfaced set (Supplementary Results S10).

**Supplementary Table S15. Root-cause classification of the 191 pairs unrecovered under the strictest filter (Supplementary Results S11, S12).**

| Category | n | % | Explanation |
|---|---:|---:|---|
| Low corpus volume | 160 | 83.8% | Fewer than 20 full-corpus mentions; consistent with the broader sparse-presence pattern (Supplementary Results S10), unlikely to rank regardless of analysis-pool size. |
| Non-mutation-framed literature | 26 | 13.6% | ≥20 full-corpus mentions, but under 20% of a 20-title sample mention mutation; real, on-topic cancer literature (expression, transcription-factor biology, immune/HLA typing, germline risk polymorphisms) that does not compete well in a mutation-focused query regardless of corpus size. |
| Vocabulary artefact | 1 | 0.5% | MAX (Brain Cancer): confirmed, by direct title inspection, to be a false-alias collision with "maximum" and the soybean species *Glycine max*, not the MAX gene (Supplementary Results S10). |
| Unresolved | 4 | 2.1% | CBFB (Breast), SF3B1 (Breast), ARID2 (Lung) and CREBBP (Lung): ≥20 full-corpus mentions with ≥20% mutation-framed titles in the sample — genuine, mutation-relevant literature that still did not rank. Reported as an open question rather than explained away. |

187 of the 191 pairs remaining under the strictest filter (97.9%) could be assigned to an identifiable diagnostic category associated with limited recoverability. Most have literature too sparse to plausibly rank, some have substantial literature that is not framed around mutation status, and one case reflects a newly-identified vocabulary collision rather than genuine gene literature. These classifications identify plausible limitations on recoverability, but they do not experimentally prove that each factor caused the corresponding ranking failure. Only four pairs (2.1%) remain unclassified after this audit.

---

# Figure legends

## Figure 1. OncoDigger analysis workflow

Nineteen separate cancer-specific PubMed collections are searched with the same mutation query. BM25 ranks papers, enrichment ranks genes, vocabulary quality control reduces ambiguous entity matches, and supporting PMIDs preserve traceability. The resulting rankings are evaluated against IntOGen, internal pooled/frequency controls and external literature-mining systems.

## Figure 2. Shared and cancer-differentiated gene-ranking patterns

Heatmap of genes reaching the top 20 in at least two cancers, split into a ranks-1-10 column block (genes reaching the top 10 in at least one cancer) and a ranks-11-20 block (genes reaching only ranks 11-20). Rows are cancers and columns are genes, and shading is a continuous gradient from darkest (rank 1) to lightest (rank 20), white indicating the gene is outside that cancer's top 20. The key's 1/5/10/15/20 ticks sample this scale, not discrete bins. The figure illustrates a shared pan-cancer backbone together with cancer-differentiated literature signals.

## Figure 3. Cancer separation and the complete pipeline contribute different properties

(A) Mean Precision@k for OncoDigger and the frequency-only baseline. (B) Per-cancer Precision@10 for the same methods. (C) Number of OncoDigger top-10 genes absent from the pooled global top 10. Panels A/B assess the combined effect of the complete pipeline relative to simple publication counting, which changes document retrieval, vocabulary behaviour and scoring together rather than isolating any one of them (Results 5.8 reports the properly isolated decomposition). Panel C assesses cancer differentiation.

## Figure 4. Comparison with published literature-mining resources

(A) Mean IntOGen-driver Precision@10 for OncoDigger, audited CancerMine, DISEASES 2.0 and subtype-aggregated PubTator 3.0 (all 19 cancer categories included). (B) Per-cancer OncoDigger-CancerMine comparison. (C) Effect of benchmark auditing on CancerMine and PubTator. PubTator's pre-correction value is omitted as unverifiable (Supplementary Table S3).

## Figure 5. Literature weight of unrecovered driver annotations

(A) Unweighted supporting-paper count by recovery band (OncoDigger top 10, ranks 11–50, ranked beyond 50) across all 1,740 IntOGen driver gene-cancer pairs. Most never surfaced within the standard analysis window, of which the large majority have a lexicon-detected mention elsewhere in the full corpus and a minority have none (Supplementary Figure S6, Supplementary Table S13). (B) OncoDigger rank versus supporting-paper count (log-log) for the 754 pairs that did surface, coloured by band, with Spearman ρ=−0.942 and the same recovery-band ordering (top 10 > ranks 11–50 > beyond 50) holding in all 19 study cancers individually. Panel B is a mechanistic ranking diagnostic, not an independent validation, since supporting-paper count is one multiplicative component of the OncoDigger score it is plotted against (Methods 2.5).

## Figure 6. IntOGen-driver enrichment across ranking depth

Cumulative IntOGen-driver Precision@k through rank 50 across 19 cancers. The bold line is the mean and the shaded region the bootstrap confidence interval. Background rank ranges are descriptive only. The inset table gives each cancer's own Precision@10 and Precision@50.

## Supplementary Figure S1. Cancer differentiation relative to the pooled control

(A) Composition of each OncoDigger top 10 as genes shared with or absent from the pooled global top 10. (B) Pairwise Jaccard similarity among the 19 cancer top-10 lists. Mean OncoDigger Jaccard=0.189, pooled control=1.0.

## Supplementary Figure S2. Collection size versus Precision@10

Since Precision@10 can only take multiples of 10%, cancers sharing the same value are offset vertically for visibility, so points do not sit exactly on the gridlines. Precision was not associated with collection size (Spearman rho=0.18, permutation p=0.47). Each cancer's own Precision@10 is given in Figure 6's inset table.

## Supplementary Figure S3. Practical ChatGPT baseline

(A) Gene-level IntOGen-driver Precision@10 for OncoDigger and the single ChatGPT response. (B) Per-cancer comparison. ChatGPT achieved 78.9% gene-level Precision@10 and 76.3% exact cancer-type driver concordance against IntOGen, with no significant difference detected from OncoDigger on the latter measure. This is a practical baseline, not a peer-system benchmark.

## Supplementary Figure S4. IntOGen-driver precision across four rank intervals

Heatmap of IntOGen-driver precision for ranks 1–10, 11–20, 21–30 and 31–50 across 19 cancers. The intervals describe ranking depth only.

## Supplementary Figure S5. IntOGen's own confidence signals for never-surfaced driver pairs

Distribution of IntOGen method-count agreement (1–7 of its 7 driver-detection methods), surfaced versus never-surfaced pairs. Never-surfaced pairs concentrate in IntOGen's single-method, lower-significance, lower-mutation-frequency tier (median q-value and median percentage of cohort mutated, Supplementary Results S9, Supplementary Table S12).

## Supplementary Figure S6. Full-corpus presence of never-surfaced driver pairs

(A) Breakdown of all 1,740 driver pairs into surfaced within the standard window, present elsewhere in the full corpus but outside it, and no lexicon-detected mention anywhere in the full corpus. 788 (45.3%) have a lexicon-detected mention elsewhere in the full corpus but did not reach the analysis pool, and 198 (11.4%) have no lexicon-detected mention anywhere in the full corpus. Neither "never surfaced" nor "no lexicon-detected mention" is the same as "absent from the literature," since detection depends on vocabulary completeness and specificity. (B) Full-corpus mention count for the 787 present-but-outside-window pairs (excludes WAS/Brain Cancer, a confirmed false-alias collision), log-binned. Most have only a handful of mentions (median 5), so a larger analysis pool alone would not recover most of them (Supplementary Results S10, Supplementary Table S13).

## Supplementary Figure S7. Recall under progressively stricter, independently justified filters

Recall at ranks 10, 20, 50 and any rank, recomputed after removing pairs with zero full-corpus mentions and, additionally, pairs backed by only one of IntOGen's seven driver-detection methods. Recall rises under the strictest filter at every depth, roughly doubling at rank 50 (24.7% to 48.9%) and rising about 62% at any rank (43.3% to 70.2%), showing that most of the apparent recall gap concentrates in pairs that are either undetectable in the analysed pool or only weakly, CGC-dependently supported by IntOGen itself (Supplementary Results S11, Supplementary Table S14).

## Supplementary Figure S8. Root-cause classification of the remaining recall gap

Classification of the 191 pairs that remain unrecovered even after restricting to detectable, ≥2-method-corroborated IntOGen driver pairs (Supplementary Results S11), based on directly sampled corpus text rather than inference. 83.8% have fewer than 20 full-corpus mentions (low volume), 13.6% have substantial literature that is not framed around mutation status, 0.5% (MAX, Brain Cancer) is a confirmed vocabulary collision with "maximum" and the soybean species *Glycine max*, and 2.1% (4 pairs) remain genuinely unresolved (Supplementary Results S12, Supplementary Table S15).

## Supplementary Figure S9. Per-cancer retrieval x vocabulary interaction

Interaction term (BM25-pool + cleaned vocabulary) minus (BM25-pool + true-raw vocabulary) minus (full-corpus + cleaned vocabulary) plus (full-corpus + true-raw vocabulary), per cancer, flat mention-count ranking against IntOGen any-cancer Precision@10. Positive in all 19 of 19 cancers individually (range +1 to +7 of 10 top-10 slots), pooled interaction +45.3 percentage points (95% cancer-bootstrap CI +36.8 to +53.7 pp, Results 5.8, Table 2).

---

# Table 1. Comparison of literature-mining resources

| Feature | OncoDigger | CancerMine | DISEASES 2.0 | PubTator 3.0 |
|---|---|---|---|---|
| Literature scope | 19 cancer-specific PubMed collections, 2000–2026 | PubMed + PMC across cancers | PubMed + PMC across diseases | PubMed + PMC across biomedical entities |
| Ranking/evidence used here | BM25 retrieval + gene enrichment | Cancer-role extraction + citation counts | Gene-disease z-score | Extracted Gene-Disease relations + distinct PMIDs |
| Separate cancer-specific source corpora | Yes | No | No | No |
| Primary output | Per-cancer ranked gene list | Driver/oncogene/TSG gene-cancer relations | Gene-disease association score | Normalized entities and relations |
| Paper-level traceability | Yes | Yes | Yes; underlying abstracts via viewer | Yes |
| Mutation-focused in this analysis | Query-specific; yes for this study | No | No | No |
| Benchmarked against IntOGen here | Yes | Yes, audited | Yes, audited | Yes, audited |

---

# Table 2. Retrieval x vocabulary factorial against IntOGen, flat mention-count ranking, all 19 cancers (Methods 3.4, Results 5.8)

| Retrieval | Vocabulary | Any-cancer P@10 | Exact concordance |
|---|---|---|---|
| Full corpus | True-raw | 20/190 = 10.5% | 2/190 = 1.1% |
| Full corpus | Cleaned | 80/190 = 42.1% | 47/190 = 24.7% |
| BM25 pool | True-raw | 28/190 = 14.7% | 10/190 = 5.3% |
| BM25 pool | Cleaned | 174/190 = 91.6% | 145/190 = 76.3% |

Interaction (observed combined effect minus the two main effects' additive prediction) is +45.3 percentage points any-cancer Precision@10 (95% cancer-bootstrap CI 36.8-53.7 pp), +47.4 pp exact cancer-type concordance (95% CI 40.0-54.2 pp), positive in all 19 of 19 cancers individually.

