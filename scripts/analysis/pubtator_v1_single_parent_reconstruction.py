"""
pubtator_v1_single_parent_reconstruction.py -- reconstructs the ORIGINAL (pre-correction,
single-parent-MeSH-only) PubTator 3.0 benchmark live from PubTator3's public API, to verify
whether the 65.8% (125/190) IntOGen-any-cancer Precision@10 figure shown in Figure 4C as
PubTator's "before" value is reproducible. The source file this number was originally
computed from (data/derived/pubtator_top10.csv) no longer exists in this repository and has
no known generator script, so Supplementary Table S3 currently disclaims the number as
unverified. This script closes that gap if its output matches.

Uses the SAME parent MeSH descriptors already documented and audited for the corrected
(v2, subtype-aggregated) benchmark (pubtator_benchmark_v2_subtype_aggregated.py's
PARENT_MESH dict, identical to Supplementary Table S5's "Parent" column) -- single entity
per cancer, no subtype aggregation, exactly reproducing the original benchmark's
(confirmed-flawed) single-parent-MeSH mapping. Genes ranked by distinct supporting-PMID
count (the "publications" field on each gene-disease relation), per Methods 2.4's
description of the original methodology.

Run:  python -X utf8 scripts/analysis/pubtator_v1_single_parent_reconstruction.py
Output: data/derived/intogen_rescoring/pubtator_v1_single_parent_top10.csv
        data/derived/intogen_rescoring/pubtator_v1_single_parent_intogen_detail.csv
"""
from __future__ import annotations

import csv
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "data" / "derived" / "intogen_rescoring"
INTOGEN_DRIVERS = ROOT / "data/external/intogen/drivers/2024-06-18_IntOGen-Drivers/Compendium_Cancer_Genes.tsv"

API_BASE = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api"
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "oncodigger-manuscript-audit/1.0 (academic use)"})
SLEEP = 0.4

# Identical to pubtator_benchmark_v2_subtype_aggregated.py's PARENT_MESH (= Supplementary
# Table S5's "Parent" column) -- the single broad parent descriptor per cancer used by the
# ORIGINAL, pre-correction benchmark.
PARENT_MESH = {
    "Breast Cancer": "D001943", "Lung Cancer": "D008175", "Colorectal Cancer": "D015179",
    "Prostate Cancer": "D011471", "Melanoma": "D008545", "Bladder Cancer": "D001749",
    "Kidney Cancer": "D007680", "Pancreatic Cancer": "D010190", "Liver Cancer": "D008113",
    "Stomach Cancer": "D013274", "Esophageal Cancer": "D004938", "Ovarian Cancer": "D010051",
    "Endometrial Cancer": "D016889", "Cervical Cancer": "D002583", "Thyroid Cancer": "D013964",
    "Brain Cancer": "D001932", "Leukemia": "D007938", "Lymphoma": "D008223", "Myeloma": "D009101",
}


def load_genes_any() -> set:
    genes = set()
    with open(INTOGEN_DRIVERS, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            genes.add(row["SYMBOL"].upper())
    return genes


def autocomplete_disease(query: str) -> list[dict]:
    r = SESSION.get(f"{API_BASE}/entity/autocomplete/", params={"query": query, "concept": "disease"}, timeout=20)
    r.raise_for_status()
    return r.json()


def get_relations(entity_id: str) -> list[dict]:
    r = SESSION.get(f"{API_BASE}/relations", params={"e1": entity_id, "type": "associate"}, timeout=30)
    r.raise_for_status()
    return r.json()


def main() -> None:
    genes_any = load_genes_any()

    top10_rows = []
    detail_rows = []
    hits_total, n_total = 0, 0

    for cancer, parent_mesh in PARENT_MESH.items():
        time.sleep(SLEEP)
        candidates = autocomplete_disease(cancer.lower())
        match = next((c for c in candidates if c.get("db_id") == parent_mesh), None)
        if match is None:
            print(f"  [{cancer}] WARNING: autocomplete did not return parent MeSH {parent_mesh}; "
                  f"candidates were {[(c.get('db_id'), c.get('name')) for c in candidates[:5]]}")
            continue
        eid = match["_id"]

        time.sleep(SLEEP)
        rels = get_relations(eid)
        gene_counts = {}
        for rel in rels:
            target, source = rel.get("target", ""), rel.get("source", "")
            gene_entity = target if target.startswith("@GENE_") else (source if source.startswith("@GENE_") else None)
            if not gene_entity:
                continue
            symbol = gene_entity[len("@GENE_"):].upper()
            gene_counts[symbol] = gene_counts.get(symbol, 0) + rel.get("publications", 0)

        ranked = sorted(gene_counts.items(), key=lambda kv: kv[1], reverse=True)[:10]
        n_hit = sum(1 for gene, _ in ranked if gene in genes_any)
        hits_total += n_hit
        n_total += len(ranked)
        for i, (gene, count) in enumerate(ranked, 1):
            hit = gene in genes_any
            top10_rows.append({"cancer": cancer, "rank": i, "gene": gene, "distinct_pmids": count,
                                "intogen_any": hit})
            detail_rows.append({"cancer": cancer, "rank": i, "gene": gene, "distinct_pmids": count,
                                 "intogen_any": hit})
        print(f"[{cancer}] parent={parent_mesh} entity={eid} top10={[g for g,_ in ranked]} "
              f"IntOGen hits: {n_hit}/{len(ranked)}", flush=True)

    print(f"\nReconstructed PubTator v1 (single-parent-MeSH, pre-correction) vs IntOGen any-cancer: "
          f"{hits_total}/{n_total} = {hits_total/n_total:.1%}")

    for path, rows in [
        (OUT_DIR / "pubtator_v1_single_parent_top10.csv", top10_rows),
        (OUT_DIR / "pubtator_v1_single_parent_intogen_detail.csv", detail_rows),
    ]:
        with path.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["cancer", "rank", "gene", "distinct_pmids", "intogen_any"])
            w.writeheader()
            w.writerows(rows)
        print(f"Saved {path}")


if __name__ == "__main__":
    main()
