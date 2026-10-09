from collections import Counter, defaultdict
import json as _json
import math
import re
import sqlite3
from urllib.parse import quote_plus

import pandas as pd

from .association import enrichment_score, entity_background_papers, single_term_background_papers
from .article_types import classify_article_types
from .discovery_modes import DEFAULT_DISCOVERY_MODE, entity_mode_multiplier, normalize_mode, term_mode_multiplier
from .embeddings import embedding_coherence
from .lexicons import norm_text
from .semantic import semantic_neighborhood_score, semantic_signal


TOKEN_RE = re.compile(r"[a-z0-9]+(?:'[a-z0-9]+)?")

STOPWORDS = {
    "a", "about", "above", "after", "all", "also", "among", "an", "and", "are", "as", "at", "be", "between",
    "anti", "both", "by", "can", "cancer", "cell", "cells", "clinical", "cohort", "conclusion", "conclusions",
    "could", "data", "disease", "during", "effect", "effects", "for", "from", "germ", "has", "have",
    "been", "having", "how", "human", "humans", "identified", "including", "increased", "in", "is", "it", "its",
    "male", "may", "methods", "more", "most", "neoplasm", "neoplasms", "of", "on", "or", "other", "our",
    "patients", "patient", "present", "primary", "rate", "related", "report", "reported", "research",
    "risk", "results", "should", "study", "studies", "such", "testicular", "testis", "than", "that", "the", "their", "these",
    "therapy", "this", "those", "to", "treatment", "tumor", "tumors", "use", "used", "using", "was",
    "we", "well", "which", "were", "with", "within", "years",
}

NOISY_TERMS = {
    "adolescent", "adult", "analysis", "animals", "article", "articles", "associated", "association", "available",
    "background", "based", "boys", "case", "cases", "child", "children", "common", "current", "diagnosis", "different",
    "current", "dose", "drug", "evidence", "factor", "factors", "female", "follow", "found", "group", "groups", "high", "however",
    "included", "important", "literature", "long", "low", "management", "medical", "middle", "overview",
    "girls", "men", "meta", "number", "outcome", "outcomes", "paper", "population", "previous", "recent", "role",
    "published", "significant", "support", "survival", "systematic", "term", "terms", "woman", "women", "year", "young",
}

NOISY_PHRASES = {
    "aged",
    "middle aged",
    "systematic review",
    "meta analysis",
    "young adult",
}

NOISY_SINGLE_TERMS = {
    "approaches",
    "available",
    "complex",
    "control",
    "current",
    "either",
    "field",
    "free",
    "following",
    "function",
    "functional",
    "gene",
    "grade",
    "include",
    "infant",
    "intermediate",
    "life",
    "local",
    "majority",
    "models",
    "normal",
    "optimal",
    "presented",
    "quality",
    "relationship",
    "response",
    "routine",
    "significantly",
    "small",
    "space",
    "solid",
    "strategies",
    "targeting",
    "therapeutic",
    "total",
    "treated",
    "type",
    "value",
    "ways",
    "will",
}


_sqlite_caches: dict[int, dict] = {}
_sqlite_doc_counts: dict[int, int] = {}


def _get_doc(index, pmid: str) -> dict:
    if not isinstance(index, sqlite3.Connection):
        return index.get("documents", {}).get(pmid, {})
    row = index.execute(
        "SELECT title, publication_year, publication_types, abstract FROM documents WHERE pmid = ?",
        (pmid,),
    ).fetchone()
    if row is None:
        return {}
    title, year, ptypes_json, abstract = row
    return {
        "title": title or "",
        "publication_year": year,
        "publication_types": _json.loads(ptypes_json) if ptypes_json else [],
        "abstract": abstract or "",
    }


def _docs_for_results(index, results: list[dict]) -> dict:
    if not isinstance(index, sqlite3.Connection):
        return index["documents"]
    pmids = [r["pmid"] for r in results]
    if not pmids:
        return {}
    ph = ",".join("?" * len(pmids))
    rows = index.execute(
        f"SELECT pmid, title, publication_year, publication_types, mesh_terms, abstract, authors, article_ids, journal "
        f"FROM documents WHERE pmid IN ({ph})",
        pmids,
    ).fetchall()
    docs = {}
    for pmid, title, year, ptypes_json, mesh_json, abstract, authors_json, aids_json, journal in rows:
        docs[pmid] = {
            "title": title or "",
            "publication_year": year,
            "publication_types": _json.loads(ptypes_json) if ptypes_json else [],
            "mesh_terms": _json.loads(mesh_json) if mesh_json else [],
            "abstract": abstract or "",
            "authors": _json.loads(authors_json) if authors_json else [],
            "article_ids": _json.loads(aids_json) if aids_json else {},
            "journal": journal or "",
        }
    return docs


def _doc_count(index) -> int:
    if not isinstance(index, sqlite3.Connection):
        return index.get("doc_count", 0)
    conn_id = id(index)
    if conn_id not in _sqlite_doc_counts:
        row = index.execute("SELECT COUNT(*) FROM documents").fetchone()
        _sqlite_doc_counts[conn_id] = row[0] if row else 0
    return _sqlite_doc_counts[conn_id]


def percent(numerator: int, denominator: int) -> float:
    if denominator <= 0:
        return 0.0
    return round(100 * numerator / denominator, 1)


def relationship_signal(support_count: int, relevant_papers: int, background_count: int, corpus_papers: int) -> tuple[str, str]:
    support_pct = percent(support_count, relevant_papers)
    corpus_pct = percent(background_count, corpus_papers)
    if support_pct >= 10 and support_pct >= corpus_pct * 2:
        return "Enriched relationship", "more concentrated in this query's evidence set than in the full corpus"
    if support_count >= 5:
        return "Direct evidence", "appears repeatedly in the relevant papers"
    return "Exploratory signal", "appears in a small number of relevant papers and should be checked in the evidence"


def pubmed_context_url(query: str, related_term: str, year_range: tuple[int, int], corpus_query: str = "") -> str:
    query_text = norm_text(query)
    related_text = norm_text(related_term)
    start_year, end_year = year_range
    parts = []
    if corpus_query:
        parts.append(f"({corpus_query})")
    if query_text:
        parts.append(f"({query_text})")
    if related_text:
        parts.append(f"({related_text})")
    parts.append(f"({start_year}/01/01:{end_year}/12/31[pdat])")
    context = " AND ".join(parts)
    return "https://pubmed.ncbi.nlm.nih.gov/?term=" + quote_plus(context)


def relationship_reason(label: str, support_count: int, relevant_papers: int, background_count: int, corpus_papers: int) -> str:
    signal, reason = relationship_signal(support_count, relevant_papers, background_count, corpus_papers)
    return (
        f"{signal}: {label} is found in {support_count} of {relevant_papers} relevant papers "
        f"({percent(support_count, relevant_papers)}%) and {background_count} of {corpus_papers} corpus papers "
        f"({percent(background_count, corpus_papers)}%); {reason}."
    )


def mode_reason(mode: str, multiplier: float) -> str:
    mode = normalize_mode(mode)
    if mode == DEFAULT_DISCOVERY_MODE:
        return "-"
    if multiplier > 1.0:
        return "↑"
    if multiplier < 1.0:
        return "↓"
    return "-"


def full_document_text(doc: dict) -> str:
    source = " ".join([doc.get("title", ""), doc.get("abstract", "")])
    cached = doc.get("_oncodigger_full_text")
    if cached is None or doc.get("_oncodigger_full_text_source") != source:
        cached = norm_text(source)
        doc["_oncodigger_full_text"] = cached
        doc["_oncodigger_full_text_source"] = source
        doc.pop("_oncodigger_tokens", None)
        doc.pop("_oncodigger_bigrams", None)
        doc.pop("_oncodigger_related_terms", None)
    return cached


def document_tokens(doc: dict) -> list[str]:
    text = full_document_text(doc)
    cached = doc.get("_oncodigger_tokens")
    if cached is None:
        cached = TOKEN_RE.findall(text)
        doc["_oncodigger_tokens"] = cached
    return cached


def document_bigrams(doc: dict) -> list[str]:
    cached = doc.get("_oncodigger_bigrams")
    if cached is None:
        tokens = document_tokens(doc)
        cached = [f"{left} {right}" for left, right in zip(tokens, tokens[1:])]
        doc["_oncodigger_bigrams"] = cached
    return cached


def document_related_terms(doc: dict, index: dict) -> list[str]:
    full_document_text(doc)
    cached = doc.get("_oncodigger_related_terms")
    if cached is None:
        cached = [
            term
            for term in document_tokens(doc) + document_bigrams(doc)
            if useful_related_term(term, index)
        ]
        doc["_oncodigger_related_terms"] = cached
    return cached


def evidence_summary(index, pmids: set[str], limit: int = 3) -> tuple[str, str]:
    doc_cache = {pmid: _get_doc(index, pmid) for pmid in pmids}
    ordered_pmids = sorted(pmids, key=lambda pmid: doc_cache[pmid].get("publication_year") or 0, reverse=True)
    selected = ordered_pmids[:limit]
    titles = []
    for pmid in selected:
        doc = doc_cache[pmid]
        title = doc.get("title") or pmid
        year = doc.get("publication_year") or ""
        titles.append(f"{pmid}: {title}" + (f" ({year})" if year else ""))
    return "; ".join(selected), " | ".join(titles)


def evidence_caution(index, pmids: set[str], relevant_papers: int) -> str:
    if not pmids:
        return "No supporting papers after the current filters."
    review_count = 0
    for pmid in pmids:
        if "Reviews and evidence syntheses" in classify_article_types(_get_doc(index, pmid)):
            review_count += 1
    support_count = len(pmids)
    cautions = []
    if support_count == 1:
        cautions.append("single-paper signal")
    elif support_count < 3 and relevant_papers >= 10:
        cautions.append("limited direct support")
    if review_count and review_count >= max(1, support_count / 2):
        cautions.append("review-heavy evidence")
    return "; ".join(cautions) if cautions else "No major evidence caution."


def bridge_terms_for_pmids(
    index,
    pmids: set[str],
    query_terms: set[str],
    excluded_terms: set[str],
    limit: int = 5,
    max_pmids: int = 10,
) -> str:
    counts = Counter()
    doc_cache = {pmid: _get_doc(index, pmid) for pmid in pmids}
    ordered_pmids = sorted(pmids, key=lambda pmid: doc_cache[pmid].get("publication_year") or 0, reverse=True)
    for pmid in ordered_pmids[:max_pmids]:
        doc = doc_cache[pmid]
        for token in document_related_terms(doc, index):
            if token in query_terms or token in excluded_terms:
                continue
            counts[token] += 1
    return ", ".join(term for term, _ in counts.most_common(limit))


def useful_related_term(term: str, index) -> bool:
    is_sqlite = isinstance(index, sqlite3.Connection)
    cache = _sqlite_caches.setdefault(id(index), {}) if is_sqlite else index.setdefault("_oncodigger_useful_related_term_cache", {})
    if term in cache:
        return cache[term]
    parts = term.split()
    if not parts:
        cache[term] = False
        return False
    if len(parts) == 2 and parts[0] == parts[1]:
        cache[term] = False
        return False
    if any(part in STOPWORDS or part in NOISY_TERMS for part in parts):
        cache[term] = False
        return False
    if term in NOISY_PHRASES:
        cache[term] = False
        return False
    if len(parts) == 1 and term in NOISY_SINGLE_TERMS:
        cache[term] = False
        return False
    if any(len(part) < 3 for part in parts):
        cache[term] = False
        return False
    if any(part.isdigit() for part in parts):
        cache[term] = False
        return False
    if term.endswith(" study") or term.endswith(" studies"):
        cache[term] = False
        return False
    if len(parts) == 1:
        token = parts[0]
        if len(token) < 4:
            cache[term] = False
            return False
        if is_sqlite:
            row = index.execute("SELECT idf FROM idf WHERE token = ?", (token,)).fetchone()
            if row is None or row[0] < 1.2:
                cache[term] = False
                return False
        else:
            if token not in index.get("idf", {}):
                cache[term] = False
                return False
            if index["idf"].get(token, 0) < 1.2:
                cache[term] = False
                return False
    cache[term] = True
    return True


def related_terms(
    index: dict,
    results: list[dict],
    query: str,
    year_range: tuple[int, int],
    query_tokens: list[str],
    limit: int = 100,
    discovery_mode: str = DEFAULT_DISCOVERY_MODE,
    embedding_index: dict | None = None,
    corpus_query: str = "",
) -> pd.DataFrame:
    discovery_mode = normalize_mode(discovery_mode)
    query_set = {token.replace("-", " ") for token in query_tokens}
    counts = Counter()
    evidence = defaultdict(set)
    documents = _docs_for_results(index, results)
    corpus_doc_count = _doc_count(index)
    for rank, result in enumerate(results):
        doc = documents.get(result["pmid"], {})
        weight = 1.0 / math.sqrt(rank + 1)
        mesh_terms = [norm_text(term) for term in doc.get("mesh_terms", [])]
        title_terms = TOKEN_RE.findall(norm_text(doc.get("title", "")))
        title_bigrams = [f"{left} {right}" for left, right in zip(title_terms, title_terms[1:])]

        for term in mesh_terms:
            if term not in query_set and useful_related_term(term, index):
                counts[term] += weight * 2.0
                evidence[term].add(result["pmid"])
        for token in title_terms + title_bigrams:
            if token not in query_set and useful_related_term(token, index):
                counts[token] += weight * 1.5
                evidence[token].add(result["pmid"])
        for token in document_related_terms(doc, index):
            if token not in query_set and len(token) >= 3:
                counts[token] += weight
                evidence[token].add(result["pmid"])

    scored_terms = []
    for term, score in counts.items():
        support_count = len(evidence[term])
        if support_count < 2 and len(results) >= 10:
            continue
        background_count = single_term_background_papers(index, term) or support_count
        relevance = enrichment_score(
            weighted_score=score,
            supporting_papers=support_count,
            relevant_papers=len(results),
            background_papers=background_count,
            corpus_papers=corpus_doc_count,
        )
        mode_multiplier = term_mode_multiplier(term, discovery_mode)
        scored_terms.append((term, support_count, background_count, relevance * mode_multiplier, mode_multiplier))
    scored_terms.sort(key=lambda item: item[3], reverse=True)

    rows = []
    for rank, (term, support_count, background_count, relevance, mode_multiplier) in enumerate(scored_terms[:limit], start=1):
        signal, _ = relationship_signal(support_count, len(results), background_count, corpus_doc_count)
        evidence_pmids, evidence_titles = evidence_summary(index, evidence[term])
        bridge_terms = bridge_terms_for_pmids(index, evidence[term], query_set, {term})
        coherence = embedding_coherence(embedding_index, evidence[term])
        semantic_score = semantic_neighborhood_score(
            bridge_terms,
            percent(support_count, len(results)),
            percent(background_count, corpus_doc_count),
            support_count,
            embedding_coherence=coherence,
        )
        rows.append(
            {
                "Rank": rank,
                "Term": term,
                "Signal": signal,
                "Why shown": relationship_reason(term, support_count, len(results), background_count, corpus_doc_count),
                "Discovery mode": discovery_mode,
                "Mode fit": mode_reason(discovery_mode, mode_multiplier),
                "Bridge terms": bridge_terms,
                "Embedding coherence": coherence,
                "Semantic signal": semantic_signal(semantic_score, signal),
                "Semantic score": semantic_score,
                "Evidence PMIDs": evidence_pmids,
                "Supporting PMIDs": ";".join(sorted(evidence[term])),
                "Evidence examples": evidence_titles,
                "Evidence caution": evidence_caution(index, evidence[term], len(results)),
                "Relevance": round(relevance, 3),
                "Supporting papers": support_count,
                "Support %": percent(support_count, len(results)),
                "Corpus papers": background_count,
                "Corpus %": percent(background_count, corpus_doc_count),
                "PubMed context": pubmed_context_url(query, term, year_range, corpus_query),
            }
        )
    return pd.DataFrame(rows)


def associated_entities(
    index: dict,
    results: list[dict],
    query: str,
    year_range: tuple[int, int],
    query_tokens: list[str],
    lookup: dict,
    limit: int = 100,
    discovery_mode: str = DEFAULT_DISCOVERY_MODE,
    embedding_index: dict | None = None,
    corpus_query: str = "",
) -> pd.DataFrame:
    return associated_entity_tables(
        index,
        results,
        query,
        year_range,
        query_tokens,
        {"default": lookup},
        limit=limit,
        discovery_mode=discovery_mode,
        embedding_index=embedding_index,
        corpus_query=corpus_query,
    )["default"]


def associated_entity_tables(
    index: dict,
    results: list[dict],
    query: str,
    year_range: tuple[int, int],
    query_tokens: list[str],
    lookups: dict[str, dict],
    limit: int = 100,
    discovery_mode: str = DEFAULT_DISCOVERY_MODE,
    embedding_index: dict | None = None,
    corpus_query: str = "",
) -> dict[str, pd.DataFrame]:
    discovery_mode = normalize_mode(discovery_mode)
    query_terms = {token.replace("-", " ") for token in query_tokens}
    scores = Counter()
    evidence = defaultdict(set)
    labels = {}
    source_ids = {}
    categories = {}
    descriptions = {}
    matched_terms = defaultdict(set)
    entity_buckets = {}
    combined_single = defaultdict(list)
    combined_phrases = defaultdict(list)
    for bucket, lookup in lookups.items():
        for token, entities in lookup.get("single", {}).items():
            combined_single[token].extend((bucket, entity) for entity in entities)
        for anchor, phrase_entities in lookup.get("phrases", {}).items():
            combined_phrases[anchor].extend((bucket, phrase, entity) for phrase, entity in phrase_entities)
    documents = _docs_for_results(index, results)
    corpus_doc_count = _doc_count(index)
    for rank, result in enumerate(results):
        doc = documents.get(result["pmid"], {})
        text = full_document_text(doc)
        padded = f" {text} "
        token_set = set(document_tokens(doc))
        weight = 1.0 / math.sqrt(rank + 1)
        seen = set()

        for token in token_set:
            if token in query_terms:
                continue
            for bucket, entity in combined_single.get(token, []):
                entity_key = (bucket, entity["id"])
                seen.add(entity_key)
                labels[entity_key] = entity["label"]
                source_ids[entity_key] = entity.get("source_id", "")
                categories[entity_key] = entity["category"]
                descriptions[entity_key] = entity.get("description", "")
                entity_buckets[entity_key] = bucket
                matched_terms[entity_key].add(token)

        for anchor in token_set:
            for bucket, phrase, entity in combined_phrases.get(anchor, []):
                if phrase in query_terms:
                    continue
                if f" {phrase} " in padded:
                    entity_key = (bucket, entity["id"])
                    seen.add(entity_key)
                    labels[entity_key] = entity["label"]
                    source_ids[entity_key] = entity.get("source_id", "")
                    categories[entity_key] = entity["category"]
                    descriptions[entity_key] = entity.get("description", "")
                    entity_buckets[entity_key] = bucket
                    matched_terms[entity_key].add(phrase)

        for entity_id in seen:
            if norm_text(labels[entity_id]) in query_terms:
                continue
            scores[entity_id] += weight
            evidence[entity_id].add(result["pmid"])

    enriched_entities = []
    for entity_id, score in scores.items():
        support_count = len(evidence[entity_id])
        if support_count < 2 and len(results) >= 10:
            continue
        background_count = entity_background_papers(index, list(matched_terms[entity_id])) or support_count
        relevance = enrichment_score(
            weighted_score=score,
            supporting_papers=support_count,
            relevant_papers=len(results),
            background_papers=background_count,
            corpus_papers=corpus_doc_count,
        )
        mode_multiplier = entity_mode_multiplier(categories[entity_id], labels[entity_id], discovery_mode)
        enriched_entities.append((entity_id, support_count, background_count, relevance * mode_multiplier, mode_multiplier))
    enriched_entities.sort(key=lambda item: item[3], reverse=True)

    rows_by_bucket = {bucket: [] for bucket in lookups}
    seen_labels_by_bucket = defaultdict(set)
    for entity_id, support_count, background_count, relevance, mode_multiplier in enriched_entities:
        bucket = entity_buckets[entity_id]
        rows = rows_by_bucket[bucket]
        if len(rows) >= limit:
            continue
        label = labels[entity_id]
        label_key = (categories[entity_id], label.lower())
        if label_key in seen_labels_by_bucket[bucket]:
            continue
        seen_labels_by_bucket[bucket].add(label_key)
        category = categories[entity_id]
        evidence_pmids, evidence_titles = evidence_summary(index, evidence[entity_id])
        coherence = embedding_coherence(embedding_index, evidence[entity_id])
        row = {
            "Rank": len(rows) + 1,
            category: label,
            "Signal": relationship_signal(support_count, len(results), background_count, corpus_doc_count)[0],
            "Why shown": relationship_reason(label, support_count, len(results), background_count, corpus_doc_count),
            "Discovery mode": discovery_mode,
            "Mode fit": mode_reason(discovery_mode, mode_multiplier),
            "Bridge terms": bridge_terms_for_pmids(
                index,
                evidence[entity_id],
                query_terms,
                {norm_text(label), *matched_terms[entity_id]},
            ),
            "Embedding coherence": coherence,
            "Semantic signal": "",
            "Semantic score": 0.0,
            "Evidence PMIDs": evidence_pmids,
            "Supporting PMIDs": ";".join(sorted(evidence[entity_id])),
            "Evidence examples": evidence_titles,
            "Evidence caution": evidence_caution(index, evidence[entity_id], len(results)),
            "Relevance": round(relevance, 3),
            "Supporting papers": support_count,
            "Support %": percent(support_count, len(results)),
            "Corpus papers": background_count,
            "Corpus %": percent(background_count, corpus_doc_count),
            "PubMed context": pubmed_context_url(query, label, year_range, corpus_query),
            "Wikipedia": f"https://en.wikipedia.org/wiki/{label.replace(' ', '_')}",
        }
        row["Semantic score"] = semantic_neighborhood_score(
            row["Bridge terms"],
            percent(support_count, len(results)),
            percent(background_count, corpus_doc_count),
            support_count,
            embedding_coherence=coherence,
        )
        row["Semantic signal"] = semantic_signal(row["Semantic score"], row["Signal"])
        if category == "Gene":
            row["HGNC Symbol"] = source_ids[entity_id]
            row["GeneCards"] = f"https://www.genecards.org/cgi-bin/carddisp.pl?gene={label}"
        elif category == "Drug":
            row["KEGG Drug ID"] = source_ids[entity_id]
            row["KEGG"] = f"https://www.genome.jp/entry/{source_ids[entity_id]}" if source_ids[entity_id] else ""
        elif category == "Compound":
            row["KEGG Compound ID"] = source_ids[entity_id]
            row["KEGG"] = f"https://www.genome.jp/entry/{source_ids[entity_id]}" if source_ids[entity_id] else ""
        elif category == "Feature":
            row["Concept class"] = source_ids[entity_id]
        elif category == "Pathway":
            row["Pathway class"] = source_ids[entity_id]
        elif category == "Exposure":
            row["Exposure class"] = source_ids[entity_id]
        if descriptions[entity_id]:
            row["Description"] = descriptions[entity_id]
        rows.append(row)
        if all(len(bucket_rows) >= limit for bucket_rows in rows_by_bucket.values()):
            break
    return {bucket: pd.DataFrame(rows) for bucket, rows in rows_by_bucket.items()}
