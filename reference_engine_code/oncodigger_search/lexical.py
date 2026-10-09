from collections import Counter, defaultdict
import heapq
import json
import re
import sqlite3

from .discovery_modes import document_mode_multiplier


TOKEN_RE = re.compile(r"[a-z0-9]+(?:[-'][a-z0-9]+)*")


def normalize_query(query: str) -> list[str]:
    tokens = []
    for segment in re.split(r"[,;|/]+", query.lower().strip()):
        segment = re.sub("[.?!&*:]", "", segment)
        segment = re.sub(r"\s+", "-", segment.strip())
        if not segment:
            continue
        segment_tokens = TOKEN_RE.findall(segment)
        tokens.extend(segment_tokens)
        if len(segment_tokens) == 1 and "-" in segment:
            tokens.extend(part for part in segment.split("-") if part)
    return list(dict.fromkeys(tokens))


def review_like(publication_types: list[str]) -> bool:
    return any("review" in publication_type.lower() for publication_type in publication_types)


def metadata_multiplier(doc: dict, query_tokens: list[str], options: dict) -> tuple[float, list[str]]:
    multiplier = 1.0
    reasons = []
    publication_types = doc.get("publication_types", [])
    if options.get("exclude_reviews") and review_like(publication_types):
        return 0.0, ["excluded_review_like"]
    if options.get("downrank_reviews", True) and review_like(publication_types):
        multiplier *= options.get("review_multiplier", 0.65)
        reasons.append("review_downrank")

    title = doc.get("title", "")
    mesh_terms = " ".join(doc.get("mesh_terms", []))
    for token in query_tokens:
        loose_token = token.replace("-", " ")
        if token in title or loose_token in title:
            multiplier *= options.get("title_match_multiplier", 1.2)
            reasons.append(f"title:{token}")
            break
    for token in query_tokens:
        loose_token = token.replace("-", " ")
        if token in mesh_terms or loose_token in mesh_terms:
            multiplier *= options.get("mesh_match_multiplier", 1.15)
            reasons.append(f"mesh:{token}")
            break

    year = doc.get("publication_year")
    if options.get("recency_boost") and year:
        if year >= 2020:
            multiplier *= 1.1
            reasons.append("recent_2020_plus")
        elif year < 2010:
            multiplier *= 0.95
            reasons.append("older_pre_2010")
    mode_multiplier, mode_reasons = document_mode_multiplier(doc, options.get("discovery_mode"))
    multiplier *= mode_multiplier
    reasons.extend(mode_reasons)
    return multiplier, reasons


def snippet(abstract: str, query_tokens: list[str], radius: int = 160) -> str:
    if not abstract:
        return ""
    lowered = abstract.lower()
    positions = [lowered.find(token.replace("-", " ")) for token in query_tokens]
    positions.extend(lowered.find(token) for token in query_tokens)
    positions = [position for position in positions if position >= 0]
    start = max(min(positions) - radius, 0) if positions else 0
    end = min(start + radius * 2, len(abstract))
    return abstract[start:end].strip()


def bm25_search(
    index,
    query_tokens: list[str],
    topn: int,
    options: dict,
    include_total: bool = False,
) -> list[dict] | tuple[list[dict], int]:
    if isinstance(index, sqlite3.Connection):
        return _bm25_sqlite(index, query_tokens, topn, options, include_total)
    return _bm25_dict(index, query_tokens, topn, options, include_total)


def _bm25_dict(
    index: dict,
    query_tokens: list[str],
    topn: int,
    options: dict,
    include_total: bool,
) -> list[dict] | tuple[list[dict], int]:
    k1 = 1.5
    b = 0.75
    avgdl = index["avg_doc_length"] or 1
    documents = index["documents"]
    scores = defaultdict(float)
    score_reasons = defaultdict(list)
    query_counts = Counter(query_tokens)
    year_min = options.get("year_min")
    year_max = options.get("year_max")

    for token, query_count in query_counts.items():
        idf = index["idf"].get(token)
        if idf is None:
            continue
        for pmid, term_frequency in index["postings"].get(token, []):
            doc = documents[pmid]
            year = doc.get("publication_year")
            if year_min and year and year < year_min:
                continue
            if year_max and year and year > year_max:
                continue
            doc_length = doc.get("length") or avgdl
            denominator = term_frequency + k1 * (1 - b + b * doc_length / avgdl)
            scores[pmid] += idf * (term_frequency * (k1 + 1)) / denominator * query_count

    for pmid in list(scores):
        multiplier, reasons = metadata_multiplier(documents[pmid], query_tokens, options)
        if multiplier == 0:
            del scores[pmid]
            score_reasons[pmid].extend(reasons)
            continue
        scores[pmid] *= multiplier
        score_reasons[pmid].extend(reasons)

    total_matches = len(scores)
    ranked = heapq.nlargest(topn, scores.items(), key=lambda item: item[1])
    results = []
    for pmid, score in ranked:
        doc = documents[pmid]
        results.append(
            {
                "pmid": pmid,
                "score": round(score, 4),
                "title": doc.get("title", ""),
                "publication_year": doc.get("publication_year"),
                "journal": doc.get("journal", ""),
                "publication_types": doc.get("publication_types", []),
                "ranking_reasons": score_reasons.get(pmid, []),
                "mesh_terms": doc.get("mesh_terms", [])[:10],
                "authors": doc.get("authors", [])[:5],
                "doi": doc.get("article_ids", {}).get("doi", ""),
                "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                "snippet": snippet(doc.get("abstract", ""), query_tokens),
            }
        )
    if include_total:
        return results, total_matches
    return results


def _bm25_sqlite(
    conn: sqlite3.Connection,
    query_tokens: list[str],
    topn: int,
    options: dict,
    include_total: bool,
) -> list[dict] | tuple[list[dict], int]:
    k1 = 1.5
    b = 0.75

    row = conn.execute("SELECT avg_doc_length FROM corpus_stats").fetchone()
    avgdl = row[0] or 1

    query_counts = Counter(query_tokens)
    unique_tokens = list(query_counts.keys())

    placeholders = ",".join("?" * len(unique_tokens))
    idf_map = dict(
        conn.execute(
            f"SELECT token, idf FROM idf WHERE token IN ({placeholders})", unique_tokens
        ).fetchall()
    )

    tokens_with_idf = [t for t in unique_tokens if t in idf_map]
    if not tokens_with_idf:
        return ([], 0) if include_total else []

    year_min = options.get("year_min")
    year_max = options.get("year_max")
    year_clauses = []
    year_params = []
    if year_min:
        year_clauses.append("d.publication_year >= ?")
        year_params.append(year_min)
    if year_max:
        year_clauses.append("d.publication_year <= ?")
        year_params.append(year_max)
    year_filter = (" AND " + " AND ".join(year_clauses)) if year_clauses else ""

    ph = ",".join("?" * len(tokens_with_idf))
    rows = conn.execute(
        f"""
        SELECT p.token, p.pmid, p.tf,
               d.doc_length, d.publication_year, d.title,
               d.mesh_terms, d.publication_types, d.abstract,
               d.authors, d.article_ids, d.journal
        FROM postings p
        JOIN documents d ON p.pmid = d.pmid
        WHERE p.token IN ({ph}){year_filter}
        """,
        tokens_with_idf + year_params,
    ).fetchall()

    scores: dict[str, float] = defaultdict(float)
    score_reasons: dict[str, list] = defaultdict(list)
    doc_cache: dict[str, dict] = {}

    for token, pmid, tf, doc_length, year, title, mesh_json, ptypes_json, abstract, authors_json, article_ids_json, journal in rows:
        idf = idf_map[token]
        dl = doc_length or avgdl
        query_count = query_counts[token]
        denominator = tf + k1 * (1 - b + b * dl / avgdl)
        scores[pmid] += idf * (tf * (k1 + 1)) / denominator * query_count

        if pmid not in doc_cache:
            doc_cache[pmid] = {
                "length": dl,
                "publication_year": year,
                "title": title or "",
                "mesh_terms": json.loads(mesh_json) if mesh_json else [],
                "publication_types": json.loads(ptypes_json) if ptypes_json else [],
                "abstract": abstract or "",
                "authors": json.loads(authors_json) if authors_json else [],
                "article_ids": json.loads(article_ids_json) if article_ids_json else {},
                "journal": journal or "",
            }

    for pmid in list(scores):
        multiplier, reasons = metadata_multiplier(doc_cache[pmid], query_tokens, options)
        if multiplier == 0:
            del scores[pmid]
            score_reasons[pmid].extend(reasons)
            continue
        scores[pmid] *= multiplier
        score_reasons[pmid].extend(reasons)

    total_matches = len(scores)
    ranked = heapq.nlargest(topn, scores.items(), key=lambda item: item[1])
    results = []
    for pmid, score in ranked:
        doc = doc_cache[pmid]
        results.append(
            {
                "pmid": pmid,
                "score": round(score, 4),
                "title": doc.get("title", ""),
                "publication_year": doc.get("publication_year"),
                "journal": doc.get("journal", ""),
                "publication_types": doc.get("publication_types", []),
                "ranking_reasons": score_reasons.get(pmid, []),
                "mesh_terms": doc.get("mesh_terms", [])[:10],
                "authors": doc.get("authors", [])[:5],
                "doi": doc.get("article_ids", {}).get("doi", ""),
                "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                "snippet": snippet(doc.get("abstract", ""), query_tokens),
            }
        )
    if include_total:
        return results, total_matches
    return results
