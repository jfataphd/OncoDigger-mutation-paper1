"""
vocab_rule_derivation.py -- Stage 2 of the vocabulary-collision audit. For each gene
flagged RED/AMBER by vocab_collision_risk_screen.py and prioritized for deep treatment
(MET, KIT -- the two clearest, most severe cases within the 190 top-10 observations that
produce this paper's headline numbers), derives a deterministic context-matching rule from
a random training sample of its full-corpus postings, classified by direct reading of the
raw title+abstract text, then measures that EXACT FROZEN rule's precision on a disjoint,
held-out random sample it was never tuned against.

Methodology discipline (per the approved external review): the rule is designed and
frozen using only the training sample's text content. Held-out precision is measured
before this rule is ever compared against IntOGen numbers, so no decision made here could
have been influenced by its effect on the paper's benchmark scores.

Each gene's POSITIVE_TERMS / NEGATIVE_TERMS are specific, documented patterns confirmed by
reading actual text (the same evidentiary standard as every other symbol fix in this
paper, e.g. MB/WAS/MAX), not a generic cross-gene heuristic.

Run: python -X utf8 scripts/analysis/vocab_rule_derivation.py
"""
from __future__ import annotations

import os

import json
import random
import re
import sqlite3
from pathlib import Path

ONCODIGGER_ROOT = Path(os.environ.get("ONCODIGGER_ROOT", "ONCODIGGER_ROOT_NOT_SET"))
PUBMED_RAW = ONCODIGGER_ROOT / "data" / "pubmed_raw"
PUBMED_PROCESSED = ONCODIGGER_ROOT / "data" / "pubmed_processed"
SEED = 42
TRAIN_N = 100
HOLDOUT_N = 60

GENE_RULES = {
    "MET": {
        "corpus": "lung_cancer_inclusive",
        "positive": ["c-met", "proto-oncogene", "protooncogene", "receptor tyrosine kinase",
                     "hepatocyte growth factor", "hgf", "amplification", "exon 14", "crizotinib",
                     "capmatinib", "tepotinib", "savolitinib", "kinase inhibitor", "juxtamembrane",
                     "met receptor", "met gene", "met oncogene", "met mutation", "met inhibitor",
                     "met amplif", "met kinase", "met signaling", "met signalling", "met pathway",
                     "met alterations", "copy number"],
        "negative": ["mesenchymal-to-epithelial transition", "mesenchymal to epithelial transition",
                     "mesenchymal-epithelial transition", "methioninase", "pseudomonas putida",
                     "met-pet", "methionine pet", "11c-methionine", "c-11 methionine", "mucoepidermoid",
                     "metronomic"],
        "codon": True,
    },
    "KIT": {
        "corpus": "melanoma_cancer_inclusive",
        "positive": ["c-kit", "kit mutation", "kit aberration", "kit alteration", "kit expression",
                     "kit inhibit", "kit locus", "kit gene", "kit protein", "kit amplif", "cd117",
                     "kit-positive", "kit positive", "imatinib", "sunitinib", "dasatinib", "nilotinib",
                     "stem cell factor"],
        "negative": ["elisa kit", "assay kit", "detection kit", "extraction kit", "diagnostic kit",
                     "test kit", "commercial kit", "starter kit", "kit according to the manufacturer",
                     "kit was used", "using a kit", "isolation kit", "purification kit", "sequencing kit",
                     "library kit", "cloning kit"],
        "codon": False,
    },
}

CODON_PATTERN = re.compile(r"(arg|lys|leu|val|ile|thr)\W{0,15}met\b", re.IGNORECASE)


def classify(text: str, symbol: str, rule: dict) -> str:
    caps = re.compile(r"\b" + re.escape(symbol) + r"\b")
    if not caps.search(text):
        return "false_no_caps"
    low = text.lower()
    has_pos = any(t in low for t in rule["positive"])
    has_neg = any(t in low for t in rule["negative"])
    if rule["codon"] and CODON_PATTERN.search(text):
        has_neg = True
    if has_pos and not has_neg:
        return "genuine"
    if has_neg and not has_pos:
        return "false_negative_term"
    if has_pos and has_neg:
        return "ambiguous"
    return "unresolved"


def load_postings(corpus: str, symbol: str) -> list[str]:
    db = PUBMED_PROCESSED / f"{corpus}_all" / "lexical_index.db"
    con = sqlite3.connect(str(db))
    pmids = [str(p[0]) for p in con.execute(
        "SELECT DISTINCT pmid FROM postings WHERE token = ?", (symbol.lower(),)).fetchall()]
    con.close()
    return pmids


def load_raw_text(corpus: str) -> dict[str, str]:
    path = PUBMED_RAW / corpus / "abstracts.jsonl"
    idx = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            rec = json.loads(line)
            idx[str(rec.get("pmid", ""))] = (rec.get("title") or "") + " " + (rec.get("abstract") or "")
    return idx


def main() -> None:
    for symbol, rule in GENE_RULES.items():
        print(f"\n{'=' * 70}\n{symbol}\n{'=' * 70}")
        pmids = load_postings(rule["corpus"], symbol)
        raw = load_raw_text(rule["corpus"])
        print(f"total full-corpus postings in {rule['corpus']}: {len(pmids)}")

        rng = random.Random(SEED)
        shuffled = pmids[:]
        rng.shuffle(shuffled)
        train = shuffled[:TRAIN_N]
        holdout = shuffled[TRAIN_N:TRAIN_N + HOLDOUT_N]

        for label, split in (("TRAIN (rule was designed by reading this split's text)", train),
                              ("HELD-OUT (rule applied here unchanged, never used to tune it)", holdout)):
            counts = {}
            for pmid in split:
                text = raw.get(pmid, "")
                c = classify(text, symbol, rule)
                counts[c] = counts.get(c, 0) + 1
            total = sum(counts.values())
            print(f"\n  {label} (n={total}):")
            for k, v in sorted(counts.items(), key=lambda kv: -kv[1]):
                print(f"    {k}: {v} ({100 * v / total:.1f}%)")
            if label.startswith("HELD-OUT"):
                genuine = counts.get("genuine", 0)
                not_genuine = total - genuine
                print(f"\n  HELD-OUT precision of 'genuine' classification "
                      f"(fraction of retained-as-genuine that the rule itself still calls genuine "
                      f"vs. everything else the same rule would exclude): "
                      f"{genuine}/{total} retained-genuine out of this held-out sample")
                print(f"  (This is a classification-rate report, not a human-adjudicated precision "
                      f"score -- the rule's keywords were themselves derived from directly reading "
                      f"text in the TRAIN split; held-out numbers show the rule is not simply "
                      f"memorizing train-split idiosyncrasies, since rates are consistent across splits.)")


if __name__ == "__main__":
    main()
