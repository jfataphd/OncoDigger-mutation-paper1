# Frozen analysis specification: retrieval/vocabulary/scoring-formula decomposition

Date frozen: 2026-10-04. This document exists specifically to protect against
benchmark-driven tuning. Every experiment below was designed by reading real
text and reasoning about the pipeline's mechanics, not by watching IntOGen
scores and adjusting until a desired number appeared. Everything in this
document is fixed before any further experimental result is examined. Any
deviation from this spec must be logged in the amendment log at the bottom,
with the reason, before the deviating result is looked at.

## Why this document exists

Over the course of this audit, a large number of iterative adjustments have
been made (vocabulary collision checks, a matched frequency-only baseline, a
2x2 retrieval/vocabulary decomposition that turned out to use an unfaithful
matcher, repeated external review). Each individual step followed the
discipline of freezing a rule before checking its effect on IntOGen agreement,
but the cumulative, iterative nature of the process creates a real risk of
implicitly adapting the pipeline to the benchmark even without any single
dishonest step. This document is the safeguard against that risk going
forward. Everything from here on is explicitly exploratory, post hoc
mechanistic analysis, not confirmatory validation, and will be described that
way in the manuscript regardless of outcome.

## What is being measured

Four components of the OncoDigger pipeline, and whether their contribution to
agreement with IntOGen driver genes can be isolated and attributed:

1. Document retrieval (BM25 + title/MeSH/recency/review multipliers +
   Primary Evidence content-type multiplier, selecting the top 1,000 papers
   per cancer for the `mutation, mutations` query).
2. Gene vocabulary cleaning (the existing 20-symbol GENE_EXCLUDE/
   suppress_ambiguous_symbols policy plus the longest-match nested-name
   correction).
3. The ranking formula (rank-weighted mention position x log(1+enrichment) x
   log(1+supporting papers)).
4. Entity-matching correctness itself (whether a gene's supporting papers
   genuinely discuss that gene, independent of whether the gene is correctly
   ranked) -- the blind spot IntOGen Precision@10 cannot see on its own.

## Conditions to be tested (frozen, will not be added to or removed without
logging an amendment first)

- Full corpus vs. BM25-selected 1,000-paper pool (retrieval axis).
- Less-cleaned vs. fully-cleaned vocabulary (vocabulary axis), using the
  REAL production code path (`load_lexicons -> suppress_ambiguous_symbols ->
  longest_match_lookups`), not a hand-reimplementation.
- Keyword-filtered-only, random-sampled-from-filtered, plain-BM25-only, and
  production-BM25 retrieval variants, to separate "any relevance filter
  helps" from "BM25's specific sophistication helps."
- Flat mention-count ranking vs. the real rank-weighted+enrichment formula,
  holding retrieval and vocabulary fixed at their production state.
- Entity-matching correctness for a stratified sample of top-10 IntOGen-hit
  genes (not a convenience sample), oversampling high collision-risk cases
  (short symbols, verified ambiguous aliases, boundary ranks, single-term-
  dependent rankings), plus a separate, non-overlapping random stratified
  sample for an unbiased population-level error estimate.

## Metrics to be reported (both, always, for every condition)

- Any-cancer Precision@10 against IntOGen's full driver gene set.
- Exact cancer-type driver concordance against IntOGen's cancer-mapped
  driver pairs (Supplementary Table S8 mapping).
- Supplementary, not primary: Precision@5/20/50, rank-membership correlation,
  paired cancer-level differences with case-bootstrap CIs (same methodology
  already used throughout this paper, not a new ad hoc method).

## What counts as a valid retrieval/matching implementation

An implementation is only trusted once it reproduces, against the known-
faithful reference path, ALL of the following for a held-out check, not just
the final percentage:
- Exact retrieved PMIDs, in rank order.
- Exact recognized genes and the specific matching term/alias per paper.
- Exact supporting-paper counts per gene.
- Exact candidate-gene exclusions and minimum-support behavior.
- Exact final ranking and tie-breaking outcome.
Matching the headline Precision@10/exact-concordance percentage alone is
NOT sufficient evidence of fidelity and will not be treated as such.

## What is NOT allowed

- No vocabulary rule, matching correction, or retrieval parameter may be
  adjusted because it changes (improves or worsens) an IntOGen-agreement
  number. Corrections are justified only by independently verified
  entity-recognition errors (reading real text), exactly as MET and REST
  were handled.
- No cherry-picking of which experiments get reported. All results in this
  program, favorable or not, are reported in the manuscript's account of
  this work.
- No treating a single random-sample draw as an adequate baseline where a
  random-sampling condition is used; repeated draws, reported as a
  distribution, are required.
- No assuming additive main effects where an interaction has not been
  explicitly tested and ruled out.

## Stopping criteria

This program stops when either: (a) all conditions above have been run and
reported with their full detail (not just pooled percentages), or (b) a
specific, logged decision is made to scope down remaining phases, with the
reason and what is being left as disclosed future work rather than silently
dropped.

## Possible outcome, stated in advance

It is possible that no single, uniquely attributable decomposition exists --
that retrieval, vocabulary, and the scoring formula interact rather than
contribute additively, since they were designed to work together. That is a
legitimate scientific finding, not a failure of this analysis, and will be
reported as such if it is what the data show.

## Amendment log

**2026-10-04, scope reduction on full-corpus cells (1 and 2 of the 2x2).**
Reason: empirically timed, not assumed. Full-corpus entity matching via the
real production matching logic took 378.9 seconds for Thyroid Cancer's 37,789
documents, the smallest of the 19 corpora. Breast Cancer has 245,265 documents
(6.5x larger); the full 19-cancer x 2-vocabulary-state full-corpus run would
require many hours, infeasible within this session's compute budget. Decision,
made BEFORE looking at whether this affects the retrieval-effect conclusion:
run cells 3/4 (BM25 pool, cheap, ~10-15s/cancer) for all 19 cancers as planned,
and run cells 1/2 (full corpus) for as many of the smallest corpora as the
remaining budget allows, prioritized by document count ascending, reporting
exactly how many were completed and leaving the rest as explicit, disclosed
future work rather than extrapolating from a partial sample to all 19.

**Outcome:** cells 3/4 completed for all 19 cancers. Cells 1/2 (full corpus)
completed for the 5 smallest corpora (Thyroid, Bladder, Esophageal,
Endometrial, Myeloma; 21,657-37,789 documents each), using the real
production engine throughout (verified via exact gene-list/rank match against
the already-faithful matched_frequency_baseline.py reference for the pool
cells). The remaining 14 cancers (up to 245,265 documents for Breast) were
not run for the full-corpus cells; this is disclosed, not extrapolated from.

**2026-10-04, true-raw-vocabulary re-baseline of the vocabulary axis
(corrects an earlier flawed "no QC" baseline, logged as a new amendment
rather than editing the entry above).** Reason: the "no QC" condition used
throughout this document up to this point (`load_lexicons()`'s raw output,
labeled "dirty"/"less-cleaned") still had `useful_gene_term()` filtering
baked into `load_lexicons()` itself, so it was never truly unfiltered. The
correct "no QC" counterfactual is `gene_quality.py`'s genuinely unfiltered
`raw` dict (zero filtering of any kind, 69,228 single tokens / 11,150 phrase
anchors vs. the much smaller "dirty" vocabulary) via a separate,
consistently-versioned copy of the engine at commit `d592370` (the live
site's deployed commit; `reference_engine_code/` stays untouched, pinned to
`4e39fc06`). This amendment re-runs the flat-count 2x2 (`true_raw_flat_count_2x2.py`)
with this corrected vocabulary baseline, reusing the existing cleaned-
vocabulary cells (2 and 4) unchanged since those did not depend on the "no
QC" baseline being wrong.

Before committing compute, Thyroid Cancer's full-corpus matching was timed
with the true-raw vocabulary specifically, since it has far more candidate
terms than the old "dirty" vocabulary and could plausibly be slower:
measured at 25.6s (fetch 0.3s + match 25.3s), NOT slower than the prior
378.9s dirty-vocabulary measurement -- in fact ~15x faster. This is not a
contradiction: the old 378.9s figure timed a block that ran BOTH the dirty
AND clean vocabulary matching passes together, and the clean pass invokes
`longest_match_lookups`'s regex-substitution overlay on every document's
full text (a monkey-patched `full_document_text`), which is the actual
dominant cost, not raw vocabulary size. True-raw vocabulary uses no such
overlay, so a single true-raw pass is fast regardless of candidate-term
count.

Scope decision, made before looking at whether it affects the retrieval-
effect conclusion: given this measured speed, all 19 cancers were feasible
for cell 1 (true-raw, full corpus) within budget. However, cell 2 (cleaned,
full corpus) remains fixed at the 5-cancer sample from the prior amendment
(not re-run here, reused as-is). Extending cell 1 alone to 19 cancers would
produce a cell 1 vs. cell 2 comparison on mismatched cancer sets. Decision:
cell 1 was run for the same 5 smallest cancers as cell 2 (Thyroid, Bladder,
Esophageal, Endometrial, Myeloma), preserving a valid matched comparison,
rather than extending cell 1 alone. Extending both cells 1 and 2 to all 19
cancers is disclosed as feasible future work (cell 1's cost is now known to
be low; cell 2's cost, using the real `longest_match_lookups` overlay, was
not re-measured here and was the likely driver of the original 378.9s
figure).

**Outcome (flat-count 2x2, true-raw-vocabulary baseline, matched 5-cancer
sample for cells 1-2, all 19 cancers for cells 3-4, pooled n=190; a matched
5-cancer n=50 subset of cells 3-4 is also reported for the paired
interaction check):**

| Cell | Condition | Any-cancer P@10 | Exact concordance |
|---|---|---|---|
| 1 | Full corpus, true-raw vocab (n=50, 5 cancers) | 5/50 = 10.0% | 0/50 = 0.0% |
| 2 | Full corpus, cleaned vocab (n=50, same 5 cancers, reused) | 21/50 = 42.0% | 9/50 = 18.0% |
| 3 | BM25 pool, true-raw vocab (n=190, all 19) | 28/190 = 14.7% | 10/190 = 5.3% |
| 3' | BM25 pool, true-raw vocab (n=50, same 5 cancers, matched) | 8/50 = 16.0% | 3/50 = 6.0% |
| 4 | BM25 pool, cleaned vocab (n=190, all 19, reused) | 174/190 = 91.6% | 145/190 = 76.3% |
| 4' | BM25 pool, cleaned vocab (n=50, same 5 cancers, matched, reused) | 44/50 = 88.0% | 35/50 = 70.0% |

Retrieval effect under true-raw vocabulary (matched 5-cancer, cell 3' - cell 1):
+6.0pp any-cancer, +6.0pp exact concordance.
Retrieval effect under cleaned vocabulary (matched 5-cancer, cell 4' - cell 2):
+46.0pp any-cancer, +52.0pp exact concordance.

Vocabulary effect under full corpus (cell 2 - cell 1): +32.0pp any-cancer,
+18.0pp exact concordance.
Vocabulary effect under BM25 pool (matched 5-cancer, cell 4' - cell 3'):
+72.0pp any-cancer, +64.0pp exact concordance (or, pooled 19-cancer,
cell 4 - cell 3: +76.9pp any-cancer, +71.1pp exact concordance).

**Interaction finding: the earlier "no interaction, purely additive" claim
does NOT hold once the vocabulary axis is correctly baselined to true-raw.**
Both retrieval and vocabulary cleaning have small individual effects in
isolation (retrieval alone: +6pp; vocabulary alone: +18 to +32pp) but a
large combined effect when applied together.

**2026-10-05, arithmetic correction (external adversarial review caught a
real error, logged per protocol before further results were examined).**
The interaction-size statement originally given here ("+78 to +91.6pp
combined vs. a naive additive prediction of only +38 to +50pp") incorrectly
mixed the any-cancer and exact-concordance metrics together in a single
number. The correct, metric-separated 2x2 interaction term
(cell4' - cell2 - cell3' + cell1, the standard definition) on the matched
5-cancer sample is:
- Any-cancer P@10: 88.0 - 42.0 - 16.0 + 10.0 = **+40.0pp** interaction
  (naive additive prediction from main effects alone: 10+6+32 = 48%; actual
  88%).
- Exact concordance: 70.0 - 18.0 - 6.0 + 0.0 = **+46.0pp** interaction
  (naive additive prediction: 0+6+18 = 24%; actual 70%).
These replace the earlier mixed-metric figures. The qualitative conclusion
(large positive interaction, not additive) is unchanged; only the earlier
arithmetic was wrong.

**2026-10-05, mechanistic follow-up: the raw-vocabulary deterioration
traces to a small, identifiable set of alias collisions, not a broad
"unfiltered vocabulary" property (external review hole #2).** Gene-level
inspection of all 19 cancers' true-raw top-10 lists (BM25 pool, real
formula; see `true_raw_vs_cleaned_all19.csv`) found that BRIP1, CD44, SPI1,
and WWOX appear in literally every single cancer's top-10 regardless of
cancer type (19/19), and WAS in 18/19 -- together consuming 95/190 (50%) of
all raw-vocabulary top-10 slots pooled across all 19 cancers. Inspection of
`gene_quality.py`'s own audit log (the "Short or common-word alias filter"
category) identified the root cause: these five genes' raw single-token
alias sets include ordinary English function words -- "of" (BRIP1, SPI1),
"in" (CD44), "for" (WWOX), and the gene symbol "was" itself (WAS) -- plus a
further 6 entries for other genes ("all", "cml", "phl" -> BCR; "not" ->
NR4A2; "as" -> UBE3A; "an" -> CCDC102B/DIAPH3/PAX6). A targeted test
(`stopword_alias_sensitivity_check.py`) removing ONLY these 13 term entries from
the true-raw vocabulary (vs. the full cleaning policy's much larger set of
corrections) and re-running the real BM25 pool + real formula for all 19
cancers recovered **76.3% of the any-cancer gap and 76.1% of the exact-
concordance gap** between true-raw (44.2%/27.4%) and fully-cleaned
(93.2%/75.8%): minimal-fix scored 81.6% any-cancer / 64.2% exact. This
means the large vocabulary-QC effect is substantially, though not entirely,
attributable to a handful of catastrophic function-word alias collisions
rather than being evidence of a broad, general property of "any amount of
cleaning helps a lot." The remaining ~24% of the gap is attributable to the
rest of the cleaning policy (143 alias-collision removals, 20-symbol
exclusions, longest-match nested-name correction). Per the external
review's recommendation, this distinction is now explicit and will be
stated as such in any manuscript integration of this work.

**2026-10-05, full 19-cancer matched 2x2 completed (addresses external
review hole #4 directly: is the 5-cancer sample representative?).** Cell 2
(full-corpus, cleaned vocabulary, flat count) was extended to all 19
cancers, timed per-cancer before committing compute (measured clean-vocab
full-corpus matching rate: ~8.2ms/document, consistent across all 19
corpora from Thyroid, 37,789 docs, to Breast, 245,265 docs; total compute
~3h15m). All four cells now cover all 19 cancers with the identical
real-engine, flat-count methodology throughout:

| Cell | Condition | Any-cancer P@10 | Exact concordance |
|---|---|---|---|
| 1 | Full corpus, true-raw vocab (n=190) | 20/190 = 10.5% | 2/190 = 1.1% |
| 2 | Full corpus, cleaned vocab (n=190) | 80/190 = 42.1% | 47/190 = 24.7% |
| 3 | BM25 pool, true-raw vocab (n=190) | 28/190 = 14.7% | 10/190 = 5.3% |
| 4 | BM25 pool, cleaned vocab (n=190) | 174/190 = 91.6% | 145/190 = 76.3% |

**Full 19-cancer interaction (cell4 - cell2 - cell3 + cell1):**
- Any-cancer P@10: **+45.3pp** [95% case-bootstrap CI: +36.8, +53.7pp;
  n=19 cancers, 2000 reps, seed 42]
- Exact concordance: **+47.4pp** [95% CI: +40.0, +54.2pp]

Both CIs exclude zero by a wide margin. These replace the matched 5-cancer
estimates (+40.0pp / +46.0pp) as the primary reported figures; the 5-cancer
estimates are retained above as the compute-constrained interim result that
motivated running the full set, and agree closely with the full-19 result.

**Per-cancer breakdown: the interaction is positive for all 19 of 19
cancers** (range +1 to +7 any-cancer hits out of 10, zero cancers with a
null or negative interaction). This directly answers external review hole
#4's concern that a handful of collision-prone cancers might be driving the
whole effect: they are not. The interaction is a general property of this
pipeline across every cancer type tested, not an artifact of which 5
cancers happened to be computationally cheapest to run first.

**Combined with the mechanistic finding above (13 stopword aliases explain
~76% of the vocabulary-QC main effect), the full picture is:** the
production pipeline's headline performance depends on BM25 retrieval and
vocabulary cleaning jointly, not on either alone; a substantial share of
why vocabulary cleaning matters so much is that a small number of specific
alias collisions (not a general property of "any filtering helps") would
otherwise contaminate every single cancer's gene ranking with the same
handful of non-specific terms, and this contamination is far more damaging
when it also has to survive retrieval-stage competition for rank position
against genuine candidates in a focused 1,000-document pool. All results in
this program have been reported, favorable and unfavorable, per the "no
cherry-picking" rule above.

**2026-10-05, the scoring formula's contribution also depends on vocabulary
state (answers the open "3rd variable" question raised during this
audit).** The formula's effect (real rank-weighted+enrichment formula vs.
flat mention-count, holding retrieval and vocabulary fixed) was previously
measured only under cleaned vocabulary: +1.6pp any-cancer / -0.5pp exact
concordance (`matched_frequency_baseline.py`: 93.2%/75.8% real vs.
91.6%/76.3% flat). `formula_interaction_dirty_vocab.py` repeats this
comparison under true-raw vocabulary instead, computing both rankings from
a single BM25-pool fetch and matching pass per cancer (not two separately-
constructed pools) to avoid cross-script pairing risk; it exactly reproduces
both pre-existing cross-script numbers (84/190 real, 28/190 flat), which
cross-validates the pairing.

Result: under true-raw vocabulary, the formula's effect is
**+29.5pp any-cancer / +22.1pp exact concordance** -- roughly 18x larger on
any-cancer P@10, and a sign flip from slightly negative to strongly
positive on exact concordance, versus the clean-vocabulary condition.

Mechanistic explanation: the enrichment term compares a gene's in-pool
mention rate to its full-corpus background rate. Stopword-collision alias
matches (e.g. "of"/"in"/"for" -> BRIP1/CD44/WWOX) occur at a near-identical
rate in the pool and the background corpus, so their enrichment ratio is
~1 and `log(1+enrichment)` correctly suppresses them; flat counting has no
such mechanism and is flooded by the same handful of stopword-aliased genes
regardless of cancer type (see the BRIP1/CD44/SPI1/WWOX/WAS finding above).
Once vocabulary is already clean, there is little such contamination left
for enrichment to suppress, so its marginal contribution is small.

**This means the scoring formula is not a general refinement with a fixed
small effect size -- its value is itself contingent on vocabulary being
dirty.** All three components (retrieval, vocabulary, formula) interact
with each other rather than contributing fixed, independent amounts; this
is consistent with, and extends, the retrieval x vocabulary interaction
finding above.

**2026-10-05, correction per external review: this is NOT a completed
three-way (2x2x2) factorial.** The two full-corpus + real-formula cells were
never empirically measured. The claim that they are "analytically
degenerate" is only partly correct: the enrichment term IS forced to
exactly 1.0 under full-corpus retrieval (proven algebraically from
`enrichment_score`'s definition -- when the analyzed set equals the
background set, `relevant_fraction` and `background_fraction` are
identical by construction), but the rank-weighting term (`W_g`, a sum of
`1/sqrt(rank+1)` over supporting documents) still depends on whatever
arbitrary order full-corpus documents are fetched in, so those two cells
are NOT provably equal to the flat-count cells and were not measured. This
audit therefore consists of two separate, complete experiments -- a full
2x2 retrieval x vocabulary factorial (flat-count, all 19 cancers), and a
conditional formula comparison (flat vs. real formula, BM25-pool only,
both vocabulary states) -- not one combined three-way design. The
manuscript must describe these as two experiments, not a completed 2x2x2.

**2026-10-05, rank-weighting vs. enrichment ablation (external review's
requested cheap check, confirms which formula component is responsible).**
`rank_vs_enrichment_ablation.py` decomposes the real formula's two
multiplicative terms (beyond the shared `log(1+support_count)` factor)
under true-raw vocabulary, same frozen BM25 pool, same matching evidence,
re-ranking four ways. Validation: the "flat" and "full" modes exactly
reproduce the already-established reference numbers (28/190, 10/190 and
84/190, 52/190 respectively), confirming this reimplementation is faithful.

| Scoring mode | Any-cancer P@10 | Exact concordance |
|---|---|---|
| Flat count | 28/190 = 14.7% | 10/190 = 5.3% |
| Rank-weighting only (enrichment forced =1) | 27/190 = 14.2% | 9/190 = 4.7% |
| Enrichment only (rank-weighting forced =1) | 129/190 = 67.9% | 94/190 = 49.5% |
| Full formula (both terms) | 84/190 = 44.2% | 52/190 = 27.4% |

**Finding: rank-weighting alone contributes essentially nothing (14.2% vs.
flat's 14.7%, within noise). Enrichment alone is not just responsible for
the formula's protective effect -- it outperforms the full production
formula** (67.9%/49.5% vs. 44.2%/27.4%) under dirty vocabulary. Combining
rank-weighting with enrichment makes performance WORSE than enrichment
alone in this regime: stopword-collision terms ("of"/"in"/"for" etc.)
appear throughout essentially every document regardless of true relevance,
including the most highly BM25-ranked ones, so rewarding mentions in
top-ranked documents does not discriminate against them and reintroduces
exactly the noise enrichment suppresses. This resolves the open question
from the external review (hole: "you haven't isolated which component is
responsible") with a direct, measured answer rather than an inference from
the formula's definition alone.

Not yet tested, left as disclosed future work: this ablation was only run
under true-raw vocabulary and BM25-pool retrieval; whether enrichment's
dominance over rank-weighting holds under cleaned vocabulary or full-corpus
retrieval was not checked.

**2026-10-05, final caveat on the rank-weighting vs. enrichment ablation
(external review, not yet independently isolated).** `W` (the rank-weighted
mention score) confounds two effects: how many papers mention the gene,
and how highly those papers rank in the BM25 pool (`W = n * mean positional
weight per supporting paper`). The ablation above establishes that adding
the full `W` term on top of enrichment reduces performance substantially
under true-raw vocabulary (67.9%/49.5% enrichment-only -> 44.2%/27.4% full
formula), but it has NOT been independently established whether this drop
is caused specifically by rewarding high-BM25-rank documents, or simply by
`W`'s implicit extra frequency weighting (beyond the already-present
`log(1+n)` term) overwhelming enrichment's correction. Decomposing `W` into
`n` and mean positional weight to isolate these two explanations was not
done and is left as explicit future work. What IS established without this
further decomposition: enrichment-only outperforms the full formula under
true-raw vocabulary; the full formula adds little after vocabulary
cleaning; and these findings are specific to this mutation-focused
benchmark and discovery mode, not a general claim about all of
OncoDigger's scoring behavior.
