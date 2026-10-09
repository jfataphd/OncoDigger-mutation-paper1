DISCOVERY_MODES = [
    "Overview",
    "Mechanisms",
    "Treatments",
    "Targets",
    "Risk & Exposure",
    "Primary Evidence",
]

DEFAULT_DISCOVERY_MODE = "Overview"

MODE_HELP = {
    "Overview": "Balanced relationship discovery across concepts, entities, and evidence.",
    "Mechanisms": "Favors molecular biology, signaling, resistance, pathways, and phenotype terms.",
    "Treatments": "Favors drugs, regimens, surgery, radiation, response, toxicity, and survival terms.",
    "Targets": "Favors genes, biomarkers, receptors, mutations, amplifications, and expression terms.",
    "Risk & Exposure": "Favors exposure, infection, epidemiology, inherited risk, and environment terms.",
    "Primary Evidence": "Favors original research and clinical evidence over broad review-style summaries.",
}

MODE_KEYWORDS = {
    "Mechanisms": {
        "apoptosis",
        "biomarker",
        "damage",
        "differentiation",
        "expression",
        "mechanism",
        "mechanisms",
        "metabolism",
        "methylation",
        "pathway",
        "pathways",
        "phenotype",
        "proliferation",
        "resistance",
        "signaling",
        "transcription",
    },
    "Treatments": {
        "antibody",
        "chemoradiation",
        "chemoradiotherapy",
        "chemotherapy",
        "immunotherapy",
        "induction",
        "inhibitor",
        "laryngectomy",
        "radiation",
        "radiotherapy",
        "regimen",
        "surgery",
        "surgical",
        "targeted",
        "toxicity",
        "transplantation",
    },
    "Targets": {
        "amplification",
        "biomarker",
        "expression",
        "gene",
        "kinase",
        "marker",
        "mutation",
        "mutations",
        "oncogene",
        "receptor",
        "target",
        "targets",
        "translocation",
    },
    "Risk & Exposure": {
        "alcohol",
        "asbestos",
        "exposure",
        "familial",
        "genotype",
        "infection",
        "incidence",
        "inherited",
        "occupational",
        "papillomavirus",
        "polymorphism",
        "prevalence",
        "risk",
        "smoking",
        "tobacco",
        "virus",
    },
    "Primary Evidence": {
        "assay",
        "clinical",
        "cohort",
        "expression",
        "mutation",
        "outcome",
        "phase",
        "trial",
        "xenograft",
    },
}

MODE_ENTITY_CATEGORY_WEIGHTS = {
    "Mechanisms": {"Gene": 1.25, "Compound": 1.1},
    "Treatments": {"Drug": 1.6, "Compound": 1.25, "Gene": 0.85, "Phytochemical": 0.85},
    "Targets": {"Gene": 1.8, "Drug": 0.9, "Compound": 0.9, "Phytochemical": 0.8},
    "Risk & Exposure": {"Exposure": 1.8, "Compound": 1.05, "Phytochemical": 1.1, "Drug": 0.75},
    "Primary Evidence": {"Gene": 1.1, "Drug": 1.1, "Compound": 1.05},
}


def normalize_mode(mode: str | None) -> str:
    return mode if mode in DISCOVERY_MODES else DEFAULT_DISCOVERY_MODE


def term_mode_multiplier(term: str, mode: str | None) -> float:
    mode = normalize_mode(mode)
    if mode == "Overview":
        return 1.0
    parts = set(str(term or "").lower().replace("-", " ").split())
    if not parts:
        return 1.0
    hits = parts & MODE_KEYWORDS.get(mode, set())
    if not hits:
        return 0.72 if mode != "Primary Evidence" else 0.9
    return min(2.4, 1.55 + (0.25 * len(hits)))


def entity_mode_multiplier(category: str, label: str, mode: str | None) -> float:
    mode = normalize_mode(mode)
    category_weight = MODE_ENTITY_CATEGORY_WEIGHTS.get(mode, {}).get(category, 1.0)
    return round(category_weight * term_mode_multiplier(label, mode), 3)


def document_mode_multiplier(doc: dict, mode: str | None) -> tuple[float, list[str]]:
    mode = normalize_mode(mode)
    if mode == "Overview":
        return 1.0, []
    keywords = MODE_KEYWORDS.get(mode, set())
    if not keywords:
        return 1.0, []

    title = str(doc.get("title", "")).lower().replace("-", " ")
    abstract = str(doc.get("abstract", "")).lower().replace("-", " ")
    mesh = " ".join(str(term).lower().replace("-", " ") for term in doc.get("mesh_terms", []))
    title_hits = {keyword for keyword in keywords if keyword in title}
    mesh_hits = {keyword for keyword in keywords if keyword in mesh}
    abstract_hits = {keyword for keyword in keywords if keyword in abstract}
    hits = title_hits | mesh_hits | abstract_hits
    if not hits:
        return (0.7 if mode != "Primary Evidence" else 0.9), [f"mode_miss:{mode}"]

    multiplier = 1.0 + (0.22 * len(title_hits)) + (0.16 * len(mesh_hits)) + (0.08 * len(abstract_hits))
    return min(multiplier, 2.6), [f"mode:{mode}", "mode_terms:" + ",".join(sorted(hits)[:5])]


def mode_explanation(mode: str | None) -> str:
    mode = normalize_mode(mode)
    return MODE_HELP[mode]
