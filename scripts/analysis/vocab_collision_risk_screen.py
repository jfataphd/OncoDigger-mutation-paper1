"""
vocab_collision_risk_screen.py -- Stage 1 of the vocabulary-collision audit. For every
gene appearing anywhere in the top-50 rankings of any of the 19 study cancers, computes
two generic, gene-agnostic diagnostic metrics using the exact same formula for every gene
(no per-gene keyword tuning at this stage -- that would be exactly the kind of cherry-
picking this audit is designed to avoid):

  C_g = fraction of the gene's full-corpus, lowercase-token-matched papers (queried
        directly from the real postings table) whose RAW, case-preserved title+abstract
        text (from pubmed_raw/<corpus>/abstracts.jsonl, before this project's own
        preprocessing lowercased everything) contains the gene's symbol as a standalone,
        whole-word, ALL-CAPS match (\\bSYMBOL\\b against the unmodified-case text).

  G_g = fraction of the gene's matched papers where the standalone capitalized symbol
        match appears within a ~3-word window of one of a small, GENERIC, gene-domain-
        agnostic qualifier set (gene, mutation, mutations, amplification, expression,
        protein, receptor, kinase, inhibitor, oncogene, pathway, signaling, variant,
        allele, locus, promoter, deletion, overexpression, copy number). This list is the
        same for every one of the 327 genes screened -- no gene-specific terms -- which is
        what keeps it a fair, uniform screen rather than per-gene tuning. (An earlier
        version of this script used "contains the gene's full HGNC Approved Name" as G_g;
        that was discarded before Stage 2 began, since established genes are almost always
        referred to by symbol only, never their spelled-out name, which made that
        definition read as near-zero for essentially every gene -- including clean ones --
        and useless for separating risk. Full-name presence is kept only as a weak,
        non-primary secondary signal below.)

Both metrics are capped at a 500-paper random sample per gene (seeded, reproducible) when
a gene's full-corpus postings exceed 500, since C_g/G_g are estimated proportions and this
keeps runtime bounded without biasing the estimate.

Thresholds (predeclared before running, not fit to the data):
  RED:   C_g < 0.50  (the majority of matches don't even contain the capitalized symbol --
         the exact pattern confirmed for MET this session: 2,203/3,967 = 55.5% had no
         standalone capitalized "MET" at all)
  AMBER: 0.50 <= C_g < 0.85  OR  G_g < 0.30 even when C_g is high (capitalization alone
         does not rule out same-case collisions, e.g. MET's own "mesenchymal-epithelial
         transition" collision was also properly capitalized)
  GREEN: C_g >= 0.85 AND G_g >= 0.30

Output: data/derived/vocab_rescoring/vocab_collision_risk_screen.csv (one row per
gene-in-top50, with every cancer it ranks in, raw postings count, sample size, C_g, G_g,
bucket).

Run: python -X utf8 scripts/analysis/vocab_collision_risk_screen.py
"""
from __future__ import annotations

import os

import csv
import json
import random
import re
import sqlite3
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ONCODIGGER_ROOT = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET"))
PUBMED_RAW = ONCODIGGER_ROOT / "data" / "pubmed_raw"
PUBMED_PROCESSED = ONCODIGGER_ROOT / "data" / "pubmed_processed"
HUMAN_GENES_CSV = ONCODIGGER_ROOT / "Human Genes.csv"
CANONICAL = ROOT / "data" / "canonical" / "precision_at_k_detail.csv"
OUT_DIR = ROOT / "data" / "derived" / "vocab_rescoring"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SAMPLE_CAP = 500
SEED = 42

# Generic, gene-domain-agnostic qualifier words for the corrected G_g definition. Same set
# applied to every gene -- no per-gene tuning.
GENERIC_QUALIFIERS = [
    "gene", "mutation", "mutations", "amplification", "expression", "protein", "receptor",
    "kinase", "inhibitor", "oncogene", "pathway", "signaling", "signalling", "variant",
    "allele", "locus", "promoter", "deletion", "overexpression", "copy number",
]
WINDOW_CHARS = 25  # roughly 3-4 words either side of the symbol match

CANCER_CORPUS = {
    "Breast Cancer": "breast_cancer_inclusive", "Lung Cancer": "lung_cancer_inclusive",
    "Colorectal Cancer": "colorectal_cancer_inclusive", "Prostate Cancer": "prostate_cancer_inclusive",
    "Melanoma": "melanoma_cancer_inclusive", "Bladder Cancer": "bladder_cancer_inclusive",
    "Kidney Cancer": "kidney_cancer_inclusive", "Pancreatic Cancer": "pancreas_cancer_inclusive",
    "Liver Cancer": "liver_cancer_inclusive", "Stomach Cancer": "stomach_cancer_inclusive",
    "Esophageal Cancer": "esophageal_cancer_inclusive", "Ovarian Cancer": "ovarian_cancer_inclusive",
    "Endometrial Cancer": "endometrial_cancer_inclusive", "Cervical Cancer": "cervical_cancer_inclusive",
    "Thyroid Cancer": "thyroid_cancer_inclusive", "Brain Cancer": "brain_cancer_inclusive",
    "Leukemia": "leukemia_cancer_inclusive", "Lymphoma": "lymphoma_cancer_inclusive",
    "Myeloma": "myeloma_cancer_inclusive",
}


def load_approved_names() -> dict[str, str]:
    names = {}
    with open(HUMAN_GENES_CSV, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            sym = str(row.get("symbol2") or row.get("symbol") or "").strip().upper()
            name = str(row.get("Approved name") or "").strip()
            if sym and name:
                names[sym] = name
    return names


def load_top50_genes() -> dict[str, set[str]]:
    """gene -> set of cancers it ranks in (rank <= 50)."""
    gene_cancers = defaultdict(set)
    with open(CANONICAL, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if int(row["rank"]) <= 50:
                gene_cancers[row["gene"].strip().upper()].add(row["cancer"])
    return dict(gene_cancers)


def raw_text_index(corpus_key: str) -> dict[str, str]:
    """pmid -> title+abstract (raw, case-preserved) for one corpus. Built once per corpus,
    reused across every gene that ranks in that cancer (avoids re-reading the JSONL file
    once per gene)."""
    path = PUBMED_RAW / corpus_key / "abstracts.jsonl"
    idx = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            rec = json.loads(line)
            pmid = str(rec.get("pmid", ""))
            idx[pmid] = (rec.get("title") or "") + " " + (rec.get("abstract") or "")
    return idx


def postings_for_gene(corpus_key: str, gene: str) -> list[str]:
    db = PUBMED_PROCESSED / f"{corpus_key}_all" / "lexical_index.db"
    con = sqlite3.connect(str(db))
    pmids = [str(p[0]) for p in con.execute(
        "SELECT DISTINCT pmid FROM postings WHERE token = ?", (gene.lower(),)).fetchall()]
    con.close()
    return pmids


def classify_bucket(c_g: float, g_g: float) -> str:
    if c_g < 0.50:
        return "RED"
    if c_g < 0.85 or g_g < 0.30:
        return "AMBER"
    return "GREEN"


def main() -> None:
    approved_names = load_approved_names()
    gene_cancers = load_top50_genes()
    print(f"genes to screen: {len(gene_cancers)}")

    rng = random.Random(SEED)
    raw_cache: dict[str, dict[str, str]] = {}
    rows = []

    for i, (gene, cancers) in enumerate(sorted(gene_cancers.items())):
        approved_name = approved_names.get(gene, "")
        caps_pattern = re.compile(r"\b" + re.escape(gene) + r"\b")
        name_lower = approved_name.lower()
        has_window_qualifier = re.compile(
            r"(" + "|".join(re.escape(q) for q in GENERIC_QUALIFIERS) + r")", re.IGNORECASE)

        # Use the first cancer this gene ranks in as the representative corpus for the
        # full-corpus collision check (a gene's symbol collision risk is a property of the
        # symbol and English/biomedical usage generally, not specific to which cancer corpus
        # it happens to rank in -- checking every corpus it ranks in would be redundant for
        # this generic screen and is unnecessary given Stage 2/3 do per-cancer tracing for
        # anything flagged here).
        corpus_key = CANCER_CORPUS[sorted(cancers)[0]]
        if corpus_key not in raw_cache:
            raw_cache[corpus_key] = raw_text_index(corpus_key)
        raw_idx = raw_cache[corpus_key]

        pmids = postings_for_gene(corpus_key, gene)
        total_postings = len(pmids)
        if total_postings == 0:
            rows.append({"gene": gene, "cancers": ";".join(sorted(cancers)),
                         "screened_corpus": corpus_key, "approved_name": approved_name,
                         "total_postings": 0, "sample_size": 0, "C_g": "", "G_g": "", "bucket": "NO_POSTINGS"})
            continue

        sample = pmids if total_postings <= SAMPLE_CAP else rng.sample(pmids, SAMPLE_CAP)
        n_caps = n_window_qual = n_name = n_checked = 0
        for pmid in sample:
            text = raw_idx.get(pmid)
            if text is None:
                continue
            n_checked += 1
            has_caps = False
            for m in caps_pattern.finditer(text):
                has_caps = True
                window = text[max(0, m.start() - WINDOW_CHARS):m.end() + WINDOW_CHARS]
                if has_window_qualifier.search(window):
                    n_window_qual += 1
                    break
            if has_caps:
                n_caps += 1
            if name_lower and name_lower in text.lower():
                n_name += 1

        c_g = n_caps / n_checked if n_checked else 0.0
        g_g = n_window_qual / n_checked if n_checked else 0.0
        name_rate = n_name / n_checked if (n_checked and name_lower) else float("nan")
        bucket = classify_bucket(c_g, g_g)

        rows.append({
            "gene": gene, "cancers": ";".join(sorted(cancers)), "screened_corpus": corpus_key,
            "approved_name": approved_name, "total_postings": total_postings, "sample_size": n_checked,
            "C_g": round(c_g, 4), "G_g": round(g_g, 4),
            "full_name_rate_secondary": (round(name_rate, 4) if name_lower else ""), "bucket": bucket,
        })
        if (i + 1) % 25 == 0:
            print(f"  screened {i + 1}/{len(gene_cancers)}...")

    out_path = OUT_DIR / "vocab_collision_risk_screen.csv"
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["gene", "cancers", "screened_corpus", "approved_name",
                                                "total_postings", "sample_size", "C_g", "G_g", "bucket"])
        writer.writeheader()
        writer.writerows(rows)

    from collections import Counter
    counts = Counter(r["bucket"] for r in rows)
    print(f"\nWrote {out_path}")
    print(f"Bucket counts: {dict(counts)}")
    print("\nRED genes (require Stage 2 rule derivation):")
    for r in rows:
        if r["bucket"] == "RED":
            print(f"  {r['gene']}: C_g={r['C_g']} G_g={r['G_g']} postings={r['total_postings']} cancers={r['cancers']}")
    print("\nAMBER genes (candidates for Stage 2, lower priority than RED):")
    for r in rows:
        if r["bucket"] == "AMBER":
            print(f"  {r['gene']}: C_g={r['C_g']} G_g={r['G_g']} postings={r['total_postings']} cancers={r['cancers']}")


if __name__ == "__main__":
    main()
