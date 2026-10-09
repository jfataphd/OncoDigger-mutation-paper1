from collections import defaultdict
import re
from pathlib import Path

import pandas as pd


def norm_text(value: str) -> str:
    value = str(value or "").lower().replace("-", " ")
    value = re.sub(r"[^a-z0-9']+", " ", value)
    return re.sub(r"\s+", " ", value).strip()


def split_terms(value: str) -> list[str]:
    terms = []
    for term in re.split(r"[;,|]", str(value or "")):
        cleaned = norm_text(term)
        if cleaned:
            terms.append(cleaned)
    return terms


def add_entity(
    entities: dict,
    entity_id: str,
    label: str,
    category: str,
    terms: list[str],
    source_id: str = "",
    description: str = "",
) -> None:
    cleaned_terms = list(dict.fromkeys(term for term in terms if term and len(term) > 1))
    if not cleaned_terms:
        return
    entities[entity_id] = {
        "id": entity_id,
        "label": label,
        "category": category,
        "terms": cleaned_terms,
        "source_id": source_id,
        "description": description,
    }


EXPOSURE_ENTITIES = {
    "alcohol": {
        "label": "Alcohol exposure",
        "class": "Lifestyle and diet",
        "terms": ["alcohol", "ethanol", "alcohol consumption", "alcohol drinking", "drinking"],
        "description": "Risk and exposure signal, not a treatment compound.",
    },
    "asbestos": {
        "label": "Asbestos exposure",
        "class": "Occupational exposure",
        "terms": ["asbestos", "asbestos exposure"],
        "description": "Occupational/environmental exposure signal.",
    },
    "silica": {
        "label": "Silica exposure",
        "class": "Occupational exposure",
        "terms": ["silica", "crystalline silica", "silica exposure"],
        "description": "Occupational/environmental exposure signal.",
    },
    "smoking": {
        "label": "Smoking",
        "class": "Lifestyle and diet",
        "terms": ["smoking", "tobacco smoking", "cigarette smoking"],
        "description": "Tobacco exposure signal.",
    },
    "tobacco": {
        "label": "Tobacco exposure",
        "class": "Lifestyle and diet",
        "terms": ["tobacco", "tobacco exposure", "tobacco consumption"],
        "description": "Tobacco exposure signal.",
    },
    "air_pollution": {
        "label": "Air pollution",
        "class": "Environmental pollutant",
        "terms": ["air pollution", "ambient air pollution", "pollutants", "pollution"],
        "description": "Environmental exposure signal.",
    },
    "hpv": {
        "label": "HPV infection",
        "class": "Infectious agent",
        "terms": ["hpv", "hpv infection", "human papillomavirus", "papillomavirus infection"],
        "description": "Viral infection and cancer-risk signal.",
    },
    "hbv": {
        "label": "HBV infection",
        "class": "Infectious agent",
        "terms": ["hbv", "hepatitis b", "hepatitis b virus", "hepatitis b infection"],
        "description": "Hepatitis B viral infection and liver-cancer risk signal.",
    },
    "hcv": {
        "label": "HCV infection",
        "class": "Infectious agent",
        "terms": ["hcv", "hepatitis c", "hepatitis c virus", "hepatitis c infection"],
        "description": "Hepatitis C viral infection and liver-cancer risk signal.",
    },
    "ebv": {
        "label": "EBV infection",
        "class": "Infectious agent",
        "terms": ["ebv", "epstein barr virus", "epstein-barr virus"],
        "description": "Epstein-Barr virus infection and cancer-risk signal.",
    },
    "h_pylori": {
        "label": "Helicobacter pylori infection",
        "class": "Infectious agent",
        "terms": ["helicobacter pylori", "h pylori", "h. pylori", "pylori infection"],
        "description": "Bacterial infection associated with gastric cancer biology and risk.",
    },
    "obesity": {
        "label": "Obesity",
        "class": "Lifestyle and metabolic risk",
        "terms": ["obesity", "obese", "body mass index", "bmi", "adiposity", "overweight"],
        "description": "Body weight, adiposity, and metabolic risk signal.",
    },
    "diabetes": {
        "label": "Diabetes",
        "class": "Lifestyle and metabolic risk",
        "terms": ["diabetes", "diabetes mellitus", "type 2 diabetes", "insulin resistance", "hyperglycemia"],
        "description": "Metabolic disease and insulin-resistance risk signal.",
    },
    "diet": {
        "label": "Dietary exposure",
        "class": "Lifestyle and diet",
        "terms": ["diet", "dietary", "red meat", "processed meat", "fiber intake", "dietary fiber", "nutrition"],
        "description": "Dietary pattern, food exposure, and nutrition signal.",
    },
    "physical_activity": {
        "label": "Physical activity",
        "class": "Lifestyle and diet",
        "terms": ["physical activity", "exercise", "sedentary", "sedentary behavior", "sedentary behaviour"],
        "description": "Activity, exercise, and sedentary-behavior signal.",
    },
    "radiation_exposure": {
        "label": "Radiation exposure",
        "class": "Radiation",
        "terms": ["radiation exposure", "ionizing radiation", "radiation-induced", "radiation induced", "ultraviolet radiation", "uv radiation"],
        "description": "Ionizing or ultraviolet radiation exposure signal.",
    },
    "estrogen_hormones": {
        "label": "Hormonal exposure",
        "class": "Hormonal/endocrine factor",
        "terms": ["estrogen", "oestrogen", "hormone replacement therapy", "oral contraceptives", "endogenous hormones", "sex hormones"],
        "description": "Hormonal and endocrine exposure signal.",
    },
    "pesticides": {
        "label": "Pesticide exposure",
        "class": "Environmental pollutant",
        "terms": ["pesticide", "pesticides", "herbicide", "herbicides", "insecticide", "insecticides"],
        "description": "Agricultural or environmental pesticide exposure signal.",
    },
    "heavy_metals": {
        "label": "Heavy metal exposure",
        "class": "Environmental pollutant",
        "terms": ["heavy metals", "arsenic", "cadmium", "chromium", "nickel", "lead exposure"],
        "description": "Heavy-metal and metal-pollutant exposure signal.",
    },
    "aflatoxin": {
        "label": "Aflatoxin exposure",
        "class": "Environmental/toxin exposure",
        "terms": ["aflatoxin", "aflatoxin b1", "mycotoxin", "mycotoxins"],
        "description": "Foodborne mycotoxin exposure signal, especially relevant to liver cancer.",
    },
}


BIOLOGY_FEATURE_ENTITIES = {
    "metastasis_invasion": {
        "label": "Metastasis and invasion",
        "class": "Tumor behavior",
        "terms": [
            "metastasis",
            "metastatic",
            "invasion",
            "invasive",
            "tumor invasion",
            "tumour invasion",
            "lymph node metastasis",
            "distant metastasis",
            "peritoneal metastasis",
            "brain metastasis",
            "bone metastasis",
            "liver metastasis",
            "lung metastasis",
        ],
        "description": "Signals spread, invasion, and metastatic behavior in the query-relevant evidence.",
    },
    "therapy_resistance": {
        "label": "Therapy resistance",
        "class": "Treatment response",
        "terms": [
            "resistance",
            "drug resistance",
            "chemoresistance",
            "chemo resistance",
            "chemotherapy resistance",
            "radioresistance",
            "radiation resistance",
            "endocrine resistance",
            "immunotherapy resistance",
            "acquired resistance",
            "treatment resistance",
            "therapeutic resistance",
        ],
        "description": "Signals reduced sensitivity, acquired resistance, or failure of treatment response.",
    },
    "treatment_response": {
        "label": "Treatment response",
        "class": "Treatment response",
        "terms": [
            "response",
            "treatment response",
            "therapeutic response",
            "complete response",
            "partial response",
            "objective response",
            "pathologic complete response",
            "pathological complete response",
            "response rate",
            "disease control",
        ],
        "description": "Signals response or sensitivity to treatment in the evidence set.",
    },
    "prognosis_survival": {
        "label": "Prognosis and survival",
        "class": "Clinical outcome",
        "terms": [
            "prognosis",
            "prognostic",
            "survival",
            "overall survival",
            "disease free survival",
            "progression free survival",
            "relapse free survival",
            "recurrence free survival",
            "mortality",
            "poor prognosis",
            "survival rate",
        ],
        "description": "Signals outcome, survival, mortality, or prognostic association.",
    },
    "recurrence_relapse": {
        "label": "Recurrence and relapse",
        "class": "Clinical outcome",
        "terms": [
            "recurrence",
            "recurrent",
            "relapse",
            "relapsed",
            "disease recurrence",
            "tumor recurrence",
            "tumour recurrence",
            "local recurrence",
            "biochemical recurrence",
        ],
        "description": "Signals cancer recurrence or relapse after treatment.",
    },
    "immune_evasion": {
        "label": "Immune evasion",
        "class": "Cancer hallmark",
        "terms": [
            "immune evasion",
            "immune escape",
            "immune suppression",
            "immunosuppression",
            "immune checkpoint",
            "checkpoint blockade",
            "pd 1",
            "pd l1",
            "ctla 4",
            "tumor microenvironment",
            "tumour microenvironment",
            "immune microenvironment",
            "t cell exhaustion",
        ],
        "description": "Signals immune escape, checkpoint biology, or suppressive tumor microenvironment.",
    },
    "immune_infiltration": {
        "label": "Immune infiltration",
        "class": "Tumor microenvironment",
        "terms": [
            "immune infiltration",
            "tumor infiltrating lymphocytes",
            "tumour infiltrating lymphocytes",
            "tils",
            "cd8 t cells",
            "t cells",
            "macrophages",
            "tumor associated macrophages",
            "tumour associated macrophages",
            "myeloid derived suppressor cells",
            "mdsc",
        ],
        "description": "Signals immune-cell infiltration or immune composition of the cancer environment.",
    },
    "angiogenesis_hypoxia": {
        "label": "Angiogenesis and hypoxia",
        "class": "Cancer hallmark",
        "terms": [
            "angiogenesis",
            "neovascularization",
            "vascularization",
            "vegf",
            "hypoxia",
            "hypoxic",
            "hif 1",
            "hif1",
            "microvessel density",
        ],
        "description": "Signals blood-vessel formation, oxygen stress, or hypoxia-driven biology.",
    },
    "dna_repair_genomic_instability": {
        "label": "DNA repair and genomic instability",
        "class": "Cancer hallmark",
        "terms": [
            "dna repair",
            "dna damage",
            "dna damage response",
            "genomic instability",
            "chromosomal instability",
            "homologous recombination",
            "homologous recombination deficiency",
            "hrd",
            "mismatch repair",
            "microsatellite instability",
            "msi",
            "tumor mutation burden",
            "tumour mutation burden",
            "tmb",
        ],
        "description": "Signals DNA damage response, repair deficiency, mutation burden, or genome instability.",
    },
    "emt_plasticity": {
        "label": "EMT and cell plasticity",
        "class": "Cancer hallmark",
        "terms": [
            "epithelial mesenchymal transition",
            "emt",
            "mesenchymal transition",
            "cell plasticity",
            "phenotypic plasticity",
            "cell migration",
            "migration",
            "motility",
        ],
        "description": "Signals epithelial-mesenchymal transition, migration, or phenotypic plasticity.",
    },
    "stemness": {
        "label": "Cancer stemness",
        "class": "Cancer hallmark",
        "terms": [
            "stemness",
            "cancer stem cell",
            "cancer stem cells",
            "tumor initiating cells",
            "tumour initiating cells",
            "self renewal",
            "sphere formation",
            "tumorsphere",
            "tumoursphere",
        ],
        "description": "Signals stem-like cancer-cell behavior or tumor-initiating cell biology.",
    },
    "metabolic_reprogramming": {
        "label": "Metabolic reprogramming",
        "class": "Cancer hallmark",
        "terms": [
            "metabolic reprogramming",
            "metabolism",
            "glycolysis",
            "aerobic glycolysis",
            "warburg effect",
            "oxidative phosphorylation",
            "oxphos",
            "lipid metabolism",
            "glutamine metabolism",
            "metabolic pathway",
        ],
        "description": "Signals altered cancer metabolism, glycolysis, oxidative phosphorylation, or nutrient use.",
    },
    "inflammation": {
        "label": "Inflammation",
        "class": "Tumor microenvironment",
        "terms": [
            "inflammation",
            "inflammatory",
            "chronic inflammation",
            "cytokine",
            "cytokines",
            "chemokine",
            "chemokines",
            "nf kappa b",
            "nfkb",
            "il 6",
            "tnf alpha",
        ],
        "description": "Signals inflammatory pathways, cytokines, or inflammation-associated cancer biology.",
    },
    "apoptosis_cell_death": {
        "label": "Apoptosis and cell death",
        "class": "Cancer hallmark",
        "terms": [
            "apoptosis",
            "apoptotic",
            "cell death",
            "programmed cell death",
            "necrosis",
            "ferroptosis",
            "pyroptosis",
            "autophagy",
            "senescence",
        ],
        "description": "Signals cell death, survival programs, autophagy, ferroptosis, or senescence.",
    },
    "proliferation_cell_cycle": {
        "label": "Proliferation and cell cycle",
        "class": "Cancer hallmark",
        "terms": [
            "proliferation",
            "cell proliferation",
            "cell cycle",
            "cell cycle progression",
            "mitosis",
            "ki 67",
            "ki67",
            "cyclin",
            "cdk",
        ],
        "description": "Signals growth, cell-cycle progression, or proliferative cancer biology.",
    },
    "epigenetic_regulation": {
        "label": "Epigenetic regulation",
        "class": "Mechanism",
        "terms": [
            "epigenetic",
            "epigenetics",
            "dna methylation",
            "methylation",
            "histone modification",
            "histone acetylation",
            "chromatin remodeling",
            "chromatin remodelling",
            "non coding rna",
            "microrna",
            "mirna",
            "lncrna",
        ],
        "description": "Signals epigenetic, chromatin, methylation, or non-coding RNA regulation.",
    },
    "liquid_biopsy_biomarkers": {
        "label": "Liquid biopsy and biomarkers",
        "class": "Translational marker",
        "terms": [
            "biomarker",
            "biomarkers",
            "liquid biopsy",
            "circulating tumor dna",
            "circulating tumour dna",
            "ctdna",
            "circulating tumor cells",
            "circulating tumour cells",
            "ctcs",
            "exosome",
            "exosomes",
            "early detection",
        ],
        "description": "Signals biomarkers, liquid biopsy, early detection, or circulating tumor material.",
    },
}


PATHWAY_ENTITIES = {
    "egfr_erbb": {
        "label": "EGFR/ERBB signaling",
        "class": "Growth factor signaling",
        "terms": ["egfr", "erbb", "erbb2", "her2", "her 2", "erbb3", "erbb4", "egfr pathway", "egfr signaling"],
        "description": "Growth-factor receptor signaling involving EGFR, HER2/ERBB2, and related ERBB-family biology.",
    },
    "pi3k_akt_mtor": {
        "label": "PI3K/AKT/mTOR signaling",
        "class": "Growth and survival signaling",
        "terms": ["pi3k", "pik3ca", "akt", "akt1", "mtor", "pten", "pi3k akt", "pi3k akt mtor", "mtor pathway"],
        "description": "Growth, survival, metabolism, and resistance pathway involving PI3K, AKT, mTOR, and PTEN.",
    },
    "mapk_erk": {
        "label": "RAS/MAPK/ERK signaling",
        "class": "Growth factor signaling",
        "terms": ["ras", "kras", "nras", "hras", "braf", "raf", "mek", "mapk", "erk", "map kinase", "mapk pathway"],
        "description": "RAS-RAF-MEK-ERK growth signaling pathway.",
    },
    "tp53_cell_stress": {
        "label": "p53 cell-stress signaling",
        "class": "Tumor suppressor pathway",
        "terms": ["tp53", "p53", "mdm2", "p21", "cdkn1a", "p53 pathway", "tumor suppressor p53", "tumour suppressor p53"],
        "description": "Tumor-suppressor pathway connecting DNA damage, cell-cycle arrest, apoptosis, and stress response.",
    },
    "wnt_beta_catenin": {
        "label": "WNT/beta-catenin signaling",
        "class": "Developmental signaling",
        "terms": ["wnt", "beta catenin", "ctnnb1", "apc", "wnt pathway", "wnt signaling", "catenin signaling"],
        "description": "Developmental and stemness pathway centered on WNT and beta-catenin signaling.",
    },
    "tgf_beta": {
        "label": "TGF-beta signaling",
        "class": "Developmental signaling",
        "terms": ["tgf beta", "tgfb", "tgfb1", "smad", "smad2", "smad3", "smad4", "tgf beta pathway", "tgf beta signaling"],
        "description": "TGF-beta and SMAD signaling linked to invasion, immune suppression, and plasticity.",
    },
    "jak_stat": {
        "label": "JAK/STAT signaling",
        "class": "Cytokine signaling",
        "terms": ["jak", "jak1", "jak2", "stat", "stat3", "stat5", "jak stat", "jak stat pathway", "stat3 signaling"],
        "description": "Cytokine and inflammatory signaling through JAK and STAT proteins.",
    },
    "nfkb_inflammation": {
        "label": "NF-kB inflammatory signaling",
        "class": "Inflammatory signaling",
        "terms": ["nf kb", "nfkb", "nf kappa b", "rela", "ikk", "nuclear factor kappa b", "nf kb pathway"],
        "description": "Inflammatory and survival signaling through NF-kB pathway activity.",
    },
    "notch": {
        "label": "NOTCH signaling",
        "class": "Developmental signaling",
        "terms": ["notch", "notch1", "notch2", "notch3", "notch4", "jagged", "jag1", "dll4", "notch pathway"],
        "description": "Developmental signaling pathway associated with cell fate, stemness, and tumor microenvironment biology.",
    },
    "hedgehog": {
        "label": "Hedgehog signaling",
        "class": "Developmental signaling",
        "terms": ["hedgehog", "sonic hedgehog", "shh", "smoothened", "smo", "gli1", "gli2", "hedgehog pathway"],
        "description": "Developmental pathway involving SHH, SMO, and GLI transcription factors.",
    },
    "dna_repair_hrd": {
        "label": "DNA repair / HRD pathway",
        "class": "Genome maintenance",
        "terms": ["brca1", "brca2", "parp", "parp1", "atm", "atr", "rad51", "homologous recombination", "hrd", "dna repair pathway"],
        "description": "DNA repair and homologous recombination deficiency pathway.",
    },
    "cell_cycle_cdk": {
        "label": "Cell-cycle / CDK signaling",
        "class": "Cell-cycle control",
        "terms": ["cell cycle", "cdk", "cdk4", "cdk6", "cyclin", "ccnd1", "rb1", "e2f", "cell cycle pathway"],
        "description": "Cell-cycle control involving cyclins, CDKs, RB, and E2F.",
    },
    "apoptosis_bcl2": {
        "label": "Apoptosis / BCL-2 signaling",
        "class": "Cell death pathway",
        "terms": ["bcl2", "bcl 2", "bax", "bak", "caspase", "caspases", "apoptosis pathway", "intrinsic apoptosis"],
        "description": "Cell-death pathway involving BCL-2 family members and caspases.",
    },
    "vegf_angiogenesis": {
        "label": "VEGF/angiogenesis signaling",
        "class": "Angiogenesis",
        "terms": ["vegf", "vegfa", "vegfr", "kdr", "flt1", "angiogenesis", "vegf pathway", "angiogenic signaling"],
        "description": "Angiogenic signaling involving VEGF and vascular growth programs.",
    },
    "immune_checkpoint": {
        "label": "Immune checkpoint signaling",
        "class": "Immune oncology",
        "terms": ["pd 1", "pd1", "pd l1", "pdl1", "cd274", "ctla4", "ctla 4", "lag3", "tigit", "immune checkpoint"],
        "description": "Immune checkpoint biology and immunotherapy-related signaling.",
    },
    "androgen_receptor": {
        "label": "Androgen receptor signaling",
        "class": "Hormone receptor signaling",
        "terms": ["androgen receptor", "ar", "androgen signaling", "androgen deprivation", "castration resistant", "castration resistance"],
        "description": "Androgen receptor and hormone signaling, especially relevant to prostate cancer.",
    },
    "estrogen_receptor": {
        "label": "Estrogen receptor signaling",
        "class": "Hormone receptor signaling",
        "terms": ["estrogen receptor", "oestrogen receptor", "er alpha", "esr1", "progesterone receptor", "pgr", "endocrine therapy"],
        "description": "Estrogen/progesterone receptor signaling and endocrine therapy biology.",
    },
}


def useful_compound_term(term: str) -> bool:
    if len(term) < 3:
        return False
    if len(term) <= 3 and any(char.isdigit() for char in term):
        return False
    if len(term) <= 4 and len(re.sub(r"[^0-9]", "", term)) >= 2:
        return False
    if re.fullmatch(r"[0-9a-z]{1,3}", term) and any(char.isdigit() for char in term):
        return False
    structural_fragments = {"10e", "11z", "12r", "14e", "14z", "17z", "18e", "20r", "20z", "23z", "28r"}
    if term in structural_fragments:
        return False
    ambiguous_terms = {
        "acid",
        "amino acid",
        "alpha",
        "base",
        "beta",
        "carbon",
        "carbon dioxide",
        "cis",
        "alcohol",
        "asbestos",
        "crystalline silica",
        "dna",
        "d glucose",
        "acetate",
        "fatty acid",
        "glycolipid",
        "glucose",
        "lead",
        "light",
        "l tyrosine",
        "methionine",
        "methyl ethyl ketone",
        "nitrogen",
        "nucleotide",
        "oxygen",
        "phosphate",
        "peptide",
        "protein",
        "reduced acceptor",
        "rna",
        "silica",
        "starch",
        "sulfate",
        "thyroid",
        "trans",
        "tyrosine",
        "water",
    }
    noisy_chemical_terms = {
        "adp ribose",
        "amp",
        "atp",
        "cholecystokinin 8",
        "formaldehyde",
        "l tyrosine",
        "methyl ethyl ketone",
    }
    if term in noisy_chemical_terms:
        return False
    if term in ambiguous_terms:
        return False
    return True


def useful_gene_term(term: str, symbol: str) -> bool:
    term = norm_text(term)
    if not term:
        return False
    symbol_norm = norm_text(symbol)
    ambiguous_symbols = {
        "cs",
        "dcr",
        "ell",
        "hccs",  # "HCCS" used as hepatocellular carcinoma abbreviation in liver cancer literature
        "impact",
        "mice",
        "set",
        "was",
    }
    if term == symbol_norm and term in ambiguous_symbols:
        return False
    if term == symbol_norm:
        return True
    if len(term) < 4:
        return False
    ambiguous_aliases = {
        "air",
        "alpha",
        "base",
        "beta",
        "ins",
    }
    if term in ambiguous_aliases:
        return False
    return True


def useful_drug_term(term: str) -> bool:
    term = norm_text(term)
    if len(term) < 4:
        return False
    if ";" in term or "," in term:
        return False
    ambiguous_drugs = {
        "ammonia solution",
        "air medical",
        "albumin human recombinant",
        "alcohol",
        "alcohol rubbing",
        "basic",
        "blood",
        "blood whole",
        "carbon dioxide",
        "estrogens conjugated",
        "fat hard",
        "globulin immune",
        "glucose",
        "glucose liquid",
        "insulin dalanated",
        "insulin human zinc extended",
        "insulin neutral",
        "insulin zinc extended",
        "medical",
        "nitrogen",
        "oxygen",
        "rose water",
        "rose water strong",
        "strong",
        "thyroid",
        "tyrosine",
        "water",
        "whole",
    }
    if term in ambiguous_drugs:
        return False
    return True


def useful_drug_label(label: str) -> bool:
    label_norm = norm_text(label)
    if ";" in label:
        return False
    if label.count(",") > 1:
        return False
    noisy_prefixes = (
        "attapulgite",
        "charcoal",
        "insulin",
        "paraffin",
        "silica",
    )
    if label_norm.startswith(noisy_prefixes):
        return False
    noisy_labels = {
        "aluminum carbonate basic",
        "fuchsin basic",
        "globulin immune",
        "lanolin modified",
        "ointment white",
    }
    return label_norm not in noisy_labels


def useful_phytochemical_term(term: str) -> bool:
    term = norm_text(term)
    if len(term) < 4:
        return False
    ambiguous_phyto = {
        "acid",
        "alpha",
        "beta",
        "cis",
        "coffee",
        "methyl",
        "phosphate",
        "pigment",
        "rose",
        "sulfate",
        "trans",
        "yellow",
    }
    return term not in ambiguous_phyto


def canonical_compound_label(terms: list[str], fallback: str) -> str:
    term_set = set(terms)
    if term_set & {"5 fu", "5 fluorouracil", "fluorouracil"}:
        return "Fluorouracil"
    if term_set & {"13 cis retinoic acid", "all trans retinoic acid", "isotretinoin", "tretinoin", "retinoic acid"}:
        return "Retinoic Acid"
    return fallback


def build_lookup(entities: dict) -> dict:
    single = defaultdict(list)
    phrases = defaultdict(list)
    for entity in entities.values():
        for term in entity["terms"]:
            parts = term.split()
            if len(parts) == 1:
                single[parts[0]].append(entity)
            else:
                phrases[parts[0]].append((term, entity))
    return {"single": dict(single), "phrases": dict(phrases)}


def load_lexicons(root: Path) -> dict:
    lexicons = {}

    genes = {}
    gene_df = pd.read_csv(root / "Human Genes.csv").fillna("")
    for row in gene_df.itertuples(index=False):
        symbol = str(getattr(row, "symbol", "")).strip()
        display = str(getattr(row, "symbol2", "") or symbol).strip().upper()
        approved_name = str(getattr(row, "_2", "") or "").strip()
        aliases = str(getattr(row, "_3", "") or "")
        preferred_alias = str(getattr(row, "preferred_alias", "") or "").strip().upper()

        if preferred_alias:
            # Primary symbol is ambiguous (abbreviation collision); index via alias instead.
            # The alias becomes the entity key and display label; source_id retains the
            # original HGNC symbol so GeneCards and HGNC Symbol columns still resolve correctly.
            alias_terms = [preferred_alias, approved_name]
            alias_terms.extend(
                a.strip() for a in aliases.split(",")
                if a.strip().upper() not in {symbol.upper(), display, preferred_alias}
            )
            alias_terms = [t for t in alias_terms if useful_gene_term(t, preferred_alias)]
            add_entity(
                genes,
                preferred_alias,
                preferred_alias,
                "Gene",
                [norm_text(t) for t in alias_terms],
                source_id=display,
                description=f"{approved_name} [primary symbol: {display}]",
            )
        else:
            terms = [symbol, display, approved_name]
            terms.extend(alias.strip() for alias in aliases.split(","))
            terms = [term for term in terms if useful_gene_term(term, display)]
            add_entity(
                genes,
                display,
                display,
                "Gene",
                [norm_text(term) for term in terms],
                source_id=display,
                description=approved_name,
            )
    lexicons["Genes"] = build_lookup(genes)

    drugs = {}
    drug_df = pd.read_csv(root / "kegg_drug_list_lowercase.csv").fillna("")
    for row in drug_df.itertuples(index=False):
        source_id = str(getattr(row, "_0", "") or getattr(row, "Drug_ID", "") or "").strip()
        drug = str(getattr(row, "drugs", "")).strip()
        label = drug.replace("-", " ").title()
        if not useful_drug_label(label):
            continue
        terms = [term for term in split_terms(drug) if useful_drug_term(term)]
        if not terms:
            continue
        add_entity(drugs, drug, label, "Drug", terms, source_id=source_id.upper())
    lexicons["Drugs"] = build_lookup(drugs)

    phytochemicals = {}
    phyto_df = pd.read_csv(root / "phytochemicals.csv").fillna("")
    for row in phyto_df.itertuples(index=False):
        phyto = str(getattr(row, "phyto", "")).strip()
        terms = [term for term in split_terms(phyto) if useful_phytochemical_term(term)]
        if not terms:
            continue
        add_entity(
            phytochemicals,
            phyto,
            phyto.replace("-", " ").title(),
            "Phytochemical",
            terms,
            source_id="PhytoHub list",
        )
    lexicons["Phytochemicals"] = build_lookup(phytochemicals)

    exposures = {}
    for exposure_id, exposure in EXPOSURE_ENTITIES.items():
        add_entity(
            exposures,
            exposure_id,
            exposure["label"],
            "Exposure",
            [norm_text(term) for term in exposure["terms"]],
            source_id=exposure.get("class", "Exposure vocabulary"),
            description=exposure["description"],
        )
    lexicons["Exposures"] = build_lookup(exposures)

    biology_features = {}
    for feature_id, feature in BIOLOGY_FEATURE_ENTITIES.items():
        add_entity(
            biology_features,
            feature_id,
            feature["label"],
            "Feature",
            [norm_text(term) for term in feature["terms"]],
            source_id=feature["class"],
            description=feature["description"],
        )
    lexicons["Biology & Clinical Features"] = build_lookup(biology_features)

    pathways = {}
    for pathway_id, pathway in PATHWAY_ENTITIES.items():
        add_entity(
            pathways,
            pathway_id,
            pathway["label"],
            "Pathway",
            [norm_text(term) for term in pathway["terms"]],
            source_id=pathway["class"],
            description=pathway["description"],
        )
    lexicons["Pathways"] = build_lookup(pathways)

    compounds = {}
    compound_df = pd.read_csv(root / "kegg_compounds_lowercase.csv").fillna("")
    for row in compound_df.itertuples(index=False):
        source_id = str(getattr(row, "_0", "") or getattr(row, "Compound_ID", "") or "").strip()
        raw_names = str(getattr(row, "compound", "")).strip()
        terms = [term for term in split_terms(raw_names) if useful_compound_term(term)]
        if not terms:
            continue
        label = canonical_compound_label(terms, terms[0].replace("-", " ").title())
        aliases = ", ".join(term for term in terms[1:6] if term != terms[0])
        add_entity(
            compounds,
            source_id or terms[0],
            label,
            "Compound",
            terms,
            source_id=source_id.upper(),
            description=f"Also known as: {aliases}" if aliases else "",
        )
    lexicons["Compounds"] = build_lookup(compounds)
    return lexicons
