"""mito_taxa.py - MitoM's taxonomic groups: the --clade choices and their genetic codes, phylum <-> clade, the metazoan
phyla and the class names that stand for them, the usual mitochondrial genetic code of a phylum, and the per-clade
options of the boundary curation (mito_curate.py)."""
from __future__ import annotations

from typing import List, Optional


# --clade choices: taxonomic groups that set the genetic code and the curation rules of the group (names kept
# from earlier MitoM versions so that existing commands keep working)
BASE_CLADES = ("Chordata", "Arthropoda", "Echinodermata", "Annelida-segmented-worms", "Bryozoa", "Mollusca",
                "Nematoda", "Nemertea-ribbon-worms", "Porifera-sponges")
# genetic code per clade when --genetic_code is not given
CLADE_CODES = {"Chordata": 2, "Echinodermata": 9, "Porifera-sponges": 4}
# further groups that need --genetic_code; "Animal" = any other phylum
EXTRA_CLADES = ("Cnidaria", "Platyhelminthes-flatworms", "Chaetognatha-arrow-worms", "Animal")
ALL_CLADES = BASE_CLADES + EXTRA_CLADES
# --clade auto: phylum (from the GenBank taxonomy) -> clade; any other metazoan phylum -> "Animal"
PHYLUM_CLADE = {"Chordata": "Chordata", "Arthropoda": "Arthropoda", "Echinodermata": "Echinodermata",
                "Annelida": "Annelida-segmented-worms", "Bryozoa": "Bryozoa", "Mollusca": "Mollusca",
                "Nematoda": "Nematoda", "Nemertea": "Nemertea-ribbon-worms", "Porifera": "Porifera-sponges",
                "Cnidaria": "Cnidaria", "Platyhelminthes": "Platyhelminthes-flatworms",
                "Chaetognatha": "Chaetognatha-arrow-worms"}
CLADE_PHYLUM = {clade: phylum for phylum, clade in PHYLUM_CLADE.items()}      # inverse: explicit --clade -> phylum
METAZOAN_PHYLA = set(PHYLUM_CLADE) | {
    "Hemichordata", "Xenacoelomorpha", "Ctenophora", "Placozoa", "Tardigrada", "Onychophora", "Nematomorpha",
    "Priapulida", "Kinorhyncha", "Loricifera", "Brachiopoda", "Phoronida", "Entoprocta", "Cycliophora", "Rotifera",
    "Acanthocephala", "Gastrotricha", "Gnathostomulida", "Orthonectida", "Dicyemida", "Micrognathozoa"}
# usual NCBI mitochondrial code when the input has no CDS transl_table (benchmark/mito_sample.py DEFAULT_CODE)
PHYLUM_CODE = {"Chordata": 2, "Echinodermata": 9, "Hemichordata": 9, "Platyhelminthes": 9, "Cnidaria": 4,
               "Porifera": 4, "Placozoa": 4, "Ctenophora": 4}
# class / subphylum names that stand for a phylum when the lineage has no phylum token (Geneious, MitoFish and
# hand-edited headers: "Vertebrata", "Hexapoda", "Gastropoda"...); used only when METAZOAN_PHYLA matches nothing
PHYLUM_SYNONYMS = {
    **dict.fromkeys(("Vertebrata", "Craniata", "Tunicata", "Cephalochordata"), "Chordata"),
    **dict.fromkeys(("Hexapoda", "Insecta", "Crustacea", "Chelicerata", "Myriapoda"), "Arthropoda"),
    **dict.fromkeys(("Gastropoda", "Bivalvia", "Cephalopoda", "Polyplacophora", "Scaphopoda"), "Mollusca"),
    **dict.fromkeys(("Polychaeta", "Clitellata", "Sipuncula"), "Annelida"),
    **dict.fromkeys(("Anthozoa", "Hydrozoa", "Scyphozoa", "Cubozoa"), "Cnidaria"),
    **dict.fromkeys(("Demospongiae", "Hexactinellida", "Calcarea", "Homoscleromorpha"), "Porifera"),
    **dict.fromkeys(("Turbellaria", "Trematoda", "Cestoda", "Monogenea"), "Platyhelminthes"),
    **dict.fromkeys(("Chromadorea", "Enoplea"), "Nematoda"),
}


def aggressive_start_genes(clade: Optional[str]) -> str:
    """mito_curate.py --aggressive_start_genes for a clade: nad5 / nad4L start at the most upstream in-frame
    start codon in NCBI annelid records (7 / 7, annelid benchmark 2026-10) but not in other phyla (173 metazoa
    records: 17 right, 52 wrong)."""
    return "ND5,ND4L" if (clade or "").lower().startswith("annelida") else ""


def abbreviated_overlap(clade: Optional[str]) -> List[str]:
    """mito_curate.py --abbreviated_overlap for the Chordata clade: vertebrate RefSeq records end a gene whose complete
    stop would run 1-2 nt into the next tRNA / gene with T / TA (leave-out benchmark: +14 exact ends in Chordata;
    other phyla keep the complete stop)."""
    return ["--abbreviated_overlap"] if (clade or "").lower().startswith("chordata") else []


def no_upstream_atg(clade: Optional[str]) -> List[str]:
    """mito_curate.py --no_upstream_atg for nematodes and arthropods, whose genes often start with TTG / ATT / ATA
    (moving such starts to an upstream ATG: 6 right, 14 wrong in the leave-out benchmark; other phyla +37)."""
    return ["--no_upstream_atg"] if (clade or "").lower().startswith(("nematoda", "arthropoda")) else []
