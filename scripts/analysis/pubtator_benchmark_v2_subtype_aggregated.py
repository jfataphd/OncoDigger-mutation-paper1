"""
pubtator_benchmark_v2_subtype_aggregated.py -- redo of the PubTator 3.0 head-to-head
benchmark (pubtator_benchmark.py), addressing red-team audit item 6.

CONFIRMED PROBLEM with the original benchmark: it maps each of the 19 cancers to a
single, broad parent MeSH descriptor (e.g. Brain Cancer -> D001932 "Brain Neoplasms")
and searches PubTator's Gene<->Disease relation data for that ID only. Verified directly
against PubTator3's own live API (2026-09-30): PubTator's NER tags well-known cancer
subtypes under their OWN distinct MeSH ID, not the broad parent -- e.g. glioblastoma is
tagged D005909, separate from D001932, and EGFR alone has 1,604 PMIDs under D005909 that
the original benchmark's single-ID search never sees. Confirmed on a second cancer too:
AML is tagged D015470, separate from Leukemia's D007938, with FLT3 (AML's single most
famous driver) carrying 3,199 PMIDs under D015470 alone.

This script aggregates across every clinically recognized subtype per cancer, using the
SAME subtype vocabulary already curated and audited for the cancer-to-subtype mapping
(cancer_aware_precision.py's CANCER_TERMS), converted to full disease-name phrases for
PubTator's entity autocomplete. For each matched disease entity (parent + subtypes), gene
publication counts are pulled from PubTator3's public relations API and SUMMED across all
of a cancer's matched entities.

Retired 2026-10-08: this script originally also scored its own top-10 against COSMIC
Cancer Gene Census Tier 1 (`tier1` column, `cosmic_census_2026-05-23.csv`, not included in
this repository). That COSMIC-based score is not used by manuscript1.1's IntOGen-calibrated
analysis and has been removed; only the ranking itself (gene, rank, cancer,
summed_publications) is still used, scored against IntOGen by
intogen_diseases_pubtator_rescoring.py.

LIMITATION, stated explicitly rather than hidden: the relations API returns per-entity
publication COUNTS, not PMID lists, so a paper mentioning two subtypes of the same cancer
(e.g. "AML" and "leukemia" in the same abstract) could be counted twice toward that gene's
total. This is a documented approximation, not silently corrected; the original single-ID
benchmark has the opposite bias (missing subtype-tagged evidence entirely), so this is not
assumed to be a strict superset, but the direction and scale of the fix is independently
motivated by the two confirmed cases above.

Run:  python pubtator_benchmark_v2_subtype_aggregated.py
Output: data/derived/pubtator_v2_top10.csv, data/derived/pubtator_v2_precision_at_10.csv,
        data/derived/pubtator_v2_entity_mapping.csv (audit trail: every matched entity per cancer)
"""
from __future__ import annotations

import csv
import re
import time
from collections import defaultdict
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[2]
OUT_TOP10 = ROOT / "data" / "derived" / "pubtator_v2_top10.csv"
OUT_MAPPING = ROOT / "data" / "derived" / "pubtator_v2_entity_mapping.csv"

API_BASE = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api"
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "oncodigger-manuscript-audit/1.0 (academic use)"})
SLEEP = 0.4  # polite rate limiting

# Parent MeSH descriptor per cancer, same as the original benchmark (pubtator_benchmark.py).
PARENT_MESH = {
    "Breast Cancer": "D001943", "Lung Cancer": "D008175", "Colorectal Cancer": "D015179",
    "Prostate Cancer": "D011471", "Melanoma": "D008545", "Bladder Cancer": "D001749",
    "Kidney Cancer": "D007680", "Pancreatic Cancer": "D010190", "Liver Cancer": "D008113",
    "Stomach Cancer": "D013274", "Esophageal Cancer": "D004938", "Ovarian Cancer": "D010051",
    "Endometrial Cancer": "D016889", "Cervical Cancer": "D002583", "Thyroid Cancer": "D013964",
    "Brain Cancer": "D001932", "Leukemia": "D007938", "Lymphoma": "D008223", "Myeloma": "D009101",
}

# Full-name search phrases for PubTator's disease-entity autocomplete, derived from the
# same subtype vocabulary already curated in cancer_aware_precision.py's CANCER_TERMS
# (there as COSMIC-tumour-type regex fragments; here as full clinical names for entity
# lookup). The cancer's own base name is always included as a query too.
SUBTYPE_QUERIES = {
    "Breast Cancer": ["breast cancer"],
    "Lung Cancer": ["lung cancer", "non-small cell lung cancer", "small cell lung cancer",
                     "lung adenocarcinoma", "lung squamous cell carcinoma"],
    "Colorectal Cancer": ["colorectal cancer", "colon cancer", "rectal cancer"],
    "Prostate Cancer": ["prostate cancer"],
    "Melanoma": ["melanoma", "uveal melanoma", "acral melanoma"],
    "Bladder Cancer": ["bladder cancer", "urothelial carcinoma"],
    "Kidney Cancer": ["kidney cancer", "renal cell carcinoma", "wilms tumor"],
    "Pancreatic Cancer": ["pancreatic cancer"],
    "Liver Cancer": ["liver cancer", "hepatocellular carcinoma", "cholangiocarcinoma"],
    "Stomach Cancer": ["gastric cancer", "stomach cancer"],
    "Esophageal Cancer": ["esophageal cancer"],
    "Ovarian Cancer": ["ovarian cancer"],
    "Endometrial Cancer": ["endometrial cancer", "endometrial carcinoma", "uterine corpus cancer"],
    "Cervical Cancer": ["cervical cancer"],
    "Thyroid Cancer": ["thyroid cancer"],
    "Brain Cancer": ["glioma", "glioblastoma", "medulloblastoma", "astrocytoma",
                      "ependymoma", "diffuse intrinsic pontine glioma"],
    "Leukemia": ["acute myeloid leukemia", "acute lymphoblastic leukemia",
                  "chronic myeloid leukemia", "chronic lymphocytic leukemia",
                  "myelodysplastic syndrome", "acute promyelocytic leukemia"],
    "Lymphoma": ["non-hodgkin lymphoma", "hodgkin lymphoma", "diffuse large b-cell lymphoma",
                  "burkitt lymphoma", "mantle cell lymphoma", "marginal zone lymphoma"],
    "Myeloma": ["multiple myeloma"],
}

GENE_ENTITY_RE = re.compile(r"^@GENE_(.+)$")


def autocomplete_disease(query: str) -> list[dict]:
    """Return candidate disease entities for a query phrase, each with entity id + mesh id."""
    r = SESSION.get(f"{API_BASE}/entity/autocomplete/", params={"query": query, "concept": "disease"}, timeout=20)
    r.raise_for_status()
    return r.json()


def get_relations(entity_id: str) -> list[dict]:
    r = SESSION.get(f"{API_BASE}/relations", params={"e1": entity_id, "type": "associate"}, timeout=30)
    r.raise_for_status()
    return r.json()


def main() -> None:
    mapping_rows = []
    top10_rows = []
    summary_rows = []
    n_total = 0

    for cancer, parent_mesh in PARENT_MESH.items():
        queries = [cancer.lower()] + SUBTYPE_QUERIES.get(cancer, [])
        entity_ids = {}  # entity_id -> mesh_id (dedup)
        for q in queries:
            time.sleep(SLEEP)
            try:
                candidates = autocomplete_disease(q)
            except Exception as e:
                print(f"  [{cancer}] autocomplete failed for {q!r}: {e}", flush=True)
                continue
            # Only the top-ranked candidate, and only when it matched on the name itself
            # (not a synonym) or the name is a clear substring/superset of the query --
            # avoids pulling in unrelated synonym collisions (e.g. "glioblastoma" also
            # autocompleting to "Retinoblastoma" via a "Glioblastoma, Retinal" synonym).
            if not candidates:
                continue
            c = candidates[0]
            eid = c.get("_id")
            mesh = c.get("db_id")
            name = (c.get("name") or "").strip()
            # Reject generic, non-neoplastic organ-disease buckets (e.g. "Brain Diseases",
            # which includes stroke/epilepsy/etc., not cancer specifically) while accepting
            # genuine cancer/neoplasm entities even when autocomplete reports "Multiple
            # matches" rather than "Matched on name".
            generic_disease_bucket = "disease" in name.lower() and not any(
                tok in name.lower() for tok in ("neoplasm", "cancer", "carcinoma", "oma", "leukemia", "lymphoma"))
            if eid and eid.startswith("@DISEASE_") and not generic_disease_bucket:
                entity_ids[eid] = (mesh, name)
            elif eid:
                print(f"    [{cancer}] skipping unrelated autocomplete match for {q!r}: {name!r} ({c.get('match','')})", flush=True)

        mapping_rows.append({"cancer": cancer, "parent_mesh": parent_mesh,
                              "matched_entities": "; ".join(f"{e}({m})" for e, (m, n) in entity_ids.items())})
        print(f"[{cancer}] matched {len(entity_ids)} disease entities: {list(entity_ids.keys())}", flush=True)

        gene_counts = defaultdict(int)
        for eid in entity_ids:
            time.sleep(SLEEP)
            try:
                rels = get_relations(eid)
            except Exception as e:
                print(f"  relations fetch failed for {eid}: {e}", flush=True)
                continue
            for rel in rels:
                target = rel.get("target", "")
                source = rel.get("source", "")
                gene_entity = target if target.startswith("@GENE_") else (source if source.startswith("@GENE_") else None)
                if not gene_entity:
                    continue
                symbol = gene_entity[len("@GENE_"):].upper()
                gene_counts[symbol] += rel.get("publications", 0)

        ranked = sorted(gene_counts.items(), key=lambda kv: kv[1], reverse=True)[:10]
        n_total += len(ranked)
        for i, (gene, count) in enumerate(ranked, 1):
            top10_rows.append({"cancer": cancer, "rank": i, "gene": gene, "summed_publications": count})
        print(f"  -> top10: {[g for g,_ in ranked]}", flush=True)

    print(f"\nWrote {n_total} ranked (cancer, gene) pairs; scored against IntOGen by "
          "intogen_diseases_pubtator_rescoring.py", flush=True)

    for path, rows, fields in [
        (OUT_MAPPING, mapping_rows, ["cancer", "parent_mesh", "matched_entities"]),
        (OUT_TOP10, top10_rows, ["cancer", "rank", "gene", "summed_publications"]),
    ]:
        with path.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)
        print(f"Saved {path}", flush=True)


if __name__ == "__main__":
    main()
