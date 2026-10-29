"""mito_names.py - gene names in one place, for MitoM and its benchmark scripts (benchmark/mito_names.py loads this
module): the 13 protein genes and the 2 rRNAs with their aliases and products, the names the annotation core writes,
the NCBI names of the curation (mito_curate N1), amino-acid codes, tRNA gene symbols (trnL1 / trnL2, trnS1 / trnS2 by
anticodon), the orientation genes and their aliases, and the normalised name of an annotated feature (gene_name,
trna_name)."""
import re

PCG = {"cox1", "cox2", "cox3", "cob", "atp6", "atp8", "nad1", "nad2", "nad3", "nad4", "nad4l", "nad5", "nad6"}
AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H",
       "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
       "Tyr": "Y", "Val": "V"}                                            # three-letter -> one-letter amino acid
AA1_TO_3 = {v: k for k, v in AA3.items()}                                  # one-letter -> three-letter


def name_key(name: str) -> str:
    """Lower case, letters and digits only: 'COX-1' and 'cox1' are one key."""
    return re.sub(r"[^a-z0-9]", "", (name or "").lower())


# ── the names the annotation core writes (mito_core.py; MITOS2 gene name -> the core's /gene and /product) ──────────
CORE_NAMES = {
    "atp6": ("ATP6", "ATP synthase F0 subunit 6"), "atp8": ("ATP8", "ATP synthase F0 subunit 8"),
    "cox1": ("COX1", "cytochrome c oxidase subunit I"), "cox2": ("COX2", "cytochrome c oxidase subunit II"),
    "cox3": ("COX3", "cytochrome c oxidase subunit III"), "cob": ("CYTB", "cytochrome b"),
    "nad1": ("ND1", "NADH dehydrogenase subunit 1"), "nad2": ("ND2", "NADH dehydrogenase subunit 2"),
    "nad3": ("ND3", "NADH dehydrogenase subunit 3"), "nad4": ("ND4", "NADH dehydrogenase subunit 4"),
    "nad4l": ("ND4L", "NADH dehydrogenase subunit 4L"), "nad5": ("ND5", "NADH dehydrogenase subunit 5"),
    "nad6": ("ND6", "NADH dehydrogenase subunit 6"),
    "rrnS": ("s-rRNA", "12S ribosomal RNA"), "rrnL": ("l-rRNA", "16S ribosomal RNA"),
}
# product of a gene by its MITOS2 name (atp9: sponges, cnidarians; mito_intron_ref writes it)
PRODUCTS = {**{k: product for k, (_, product) in CORE_NAMES.items()}, "atp9": "ATP synthase F0 subunit 9"}

# ── NCBI gene symbols and products (mito_curate N1, as in RefSeq, e.g. NC_012920.1) and their synonyms (name_key) ──
_PCG_SYNONYMS = {
    "COX1": ("cox1", "coi", "co1", "coxi", "cytochromecoxidasesubuniti"),
    "COX2": ("cox2", "coii", "co2", "coxii", "cytochromecoxidasesubunitii"),
    "COX3": ("cox3", "coiii", "co3", "coxiii", "cytochromecoxidasesubunitiii"),
    "ATP6": ("atp6", "atpase6", "atp", "atpsynthasef0subunit6"),
    "ATP8": ("atp8", "atpase8", "atpsynthasef0subunit8"),
    "CYTB": ("cytb", "cob", "ctb", "cyb", "cytochromeb"),
    "ND1": ("nd1", "nad1", "nadh1", "nadhdehydrogenasesubunit1"),
    "ND2": ("nd2", "nad2", "nadh2", "nadhdehydrogenasesubunit2"),
    "ND3": ("nd3", "nad3", "nadh3", "nadhdehydrogenasesubunit3"),
    "ND4": ("nd4", "nad4", "nadh4", "nadhdehydrogenasesubunit4"),
    "ND4L": ("nd4l", "nad4l", "nadh4l", "nadhdehydrogenasesubunit4l"),
    "ND5": ("nd5", "nad5", "nadh5", "nadhdehydrogenasesubunit5"),
    "ND6": ("nd6", "nad6", "nadh6", "nadhdehydrogenasesubunit6"),
}
_CORE_BY_SYMBOL = {symbol: product for k, (symbol, product) in CORE_NAMES.items() if k in PCG}
PCG_NAMES = {symbol: (_CORE_BY_SYMBOL[symbol], syns) for symbol, syns in _PCG_SYNONYMS.items()}
PCG_SYN = {syn: name for name, (_, syns) in PCG_NAMES.items() for syn in syns}      # synonym key -> NCBI symbol
RRNA_PRODUCT = {g: CORE_NAMES[g][1] for g in ("rrnS", "rrnL")}                    # rrnS / rrnL -> NCBI product

# ── normalised names of annotated features (mito_check, mito_compare, mito_intron_ref, benchmark) ────────────────
RRNA_ALIAS = {"rrns": "rrnS", "srrna": "rrnS", "12s": "rrnS", "12srrna": "rrnS", "rns": "rrnS", "rrn12": "rrnS",
              "rrnl": "rrnL", "lrrna": "rrnL", "16s": "rrnL", "16srrna": "rrnL", "rnl": "rrnL", "rrn16": "rrnL"}
ALIAS = {"coi": "cox1", "coii": "cox2", "coiii": "cox3", "co1": "cox1", "co2": "cox2", "co3": "cox3",
         "cytb": "cob", "cyb": "cob", "nd1": "nad1", "nd2": "nad2", "nd3": "nad3", "nd4": "nad4", "nd4l": "nad4l",
         "nd5": "nad5", "nd6": "nad6", "atpase6": "atp6", "atpase8": "atp8", **RRNA_ALIAS}
PRODUCT = {"12s ribosomal rna": "rrnS", "small subunit ribosomal rna": "rrnS", "s-rrna": "rrnS",
           "16s ribosomal rna": "rrnL", "large subunit ribosomal rna": "rrnL", "l-rrna": "rrnL"}
ROMAN = {"i": "1", "ii": "2", "iii": "3"}
PCG_PRODUCT = [(r"cytochrome c oxidase subunit (i{1,3}|[123])\b", "cox"), (r"cytochrome b\b", "cob"),
               (r"nadh dehydrogenase subunit (4l|[1-6])\b", "nad"), (r"atp synthase (?:f0 )?subunit (6|8)\b", "atp")]
_PCG_PRODUCT_ATP9 = PCG_PRODUCT[:3] + [(r"atp synthase (?:f0 )?subunit (6|8|9)\b", "atp")]


def gene_name(qualifiers: dict, ftype: str, atp9: bool = False) -> str:
    """Normalised PCG / rRNA name (cox1, nad4l, rrnS ...); '' when not a standard one.
    Split genes such as COX1a / COX1b map to their gene (cox1). atp9=True also accepts atp9 (intron sites)."""
    pcg = PCG | {"atp9"} if atp9 else PCG
    g = re.sub(r"[^a-z0-9]", "", (qualifiers.get("gene", [""])[0] or "").lower())
    g = ALIAS.get(g, g)
    m = re.fullmatch(r"(cox[123]|nad[1-6]|nad4l|atp[689]|cob)[abc]" if atp9 else
                     r"(cox[123]|nad[1-6]|nad4l|atp[68]|cob)[abc]", g)
    if m:
        g = m.group(1)
    prod = (qualifiers.get("product", [""])[0] or "").lower()
    if ftype == "rRNA" and g not in ("rrnS", "rrnL"):
        g = PRODUCT.get(prod) or ALIAS.get(re.sub(r"[^a-z0-9]", "", prod), g)      # MitoFinder: /product="rrnS"
    if ftype == "CDS" and g not in pcg:
        for pat, base in (_PCG_PRODUCT_ATP9 if atp9 else PCG_PRODUCT):
            m = re.search(pat, prod)
            if m:
                g = base + ROMAN.get(m.group(1), m.group(1)) if base != "cob" else "cob"
                break
    return g if g in pcg or g in ("rrnS", "rrnL") else ""


def trna_name(qualifiers: dict) -> tuple:
    """(amino acid one-letter, anticodon or '') of a tRNA feature; ('', '') if unknown."""
    prod = qualifiers.get("product", [""])[0] or ""
    gene = qualifiers.get("gene", [""])[0] or ""
    m = re.search(r"tRNA-([A-Z][a-z]{2})", prod) or re.search(r"tRNA-([A-Z][a-z]{2})", gene)
    aa = AA3.get(m.group(1), "") if m else ""
    if not aa:
        m = re.match(r"trn([A-Z])", gene)
        aa = m.group(1) if m else ""
    text = " ".join(qualifiers.get("anticodon", []))
    ac = re.search(r"seq:([acgtu]{3})", text, re.I) or re.fullmatch(r"([acgtu]{3})", text.strip(), re.I)
    if not ac:
        ac = re.search(r"\(([acgtu]{3})\)", gene + " " + prod, re.I) or \
            re.search(r"anticodon[:= ]+([acgtu]{3})", " ".join(qualifiers.get("note", [])), re.I)
    return aa, (ac.group(1).lower().replace("u", "t") if ac else "")


def trna_number(qualifiers: dict) -> str:
    """'1' / '2' when the Leu / Ser gene is numbered (trnL1, trnS2, tRNA-Leu2, tRNA-Ser1 in gene / product), else ''
    (not /note: free text such as 'similar to tRNA-Ser 2 of ...' is not the gene's number)."""
    text = " ".join(qualifiers.get(k, [""])[0] for k in ("gene", "product"))
    m = re.search(r"(?:trn[LS]|tRNA-(?:Leu|Ser))\s*([12])\b", text)
    return m.group(1) if m else ""


# ── tRNA gene symbols ────────────────────────────────────────────────────────────────────────────────────────────────
def trna_symbol(aa: str, ac: str) -> str:
    """trnL1/trnL2 and trnS1/trnS2 by anticodon (invertebrate convention: L1 = CUN (tag), L2 = UUR (taa);
    S1 = AGN (tct/gct), S2 = UCN (tga))."""
    if aa == "L":
        return "trnL1" if ac.endswith("ag") else "trnL2" if ac.endswith("aa") else "trnL"
    if aa == "S":
        return "trnS1" if ac.endswith("ct") else "trnS2" if ac.endswith("ga") else "trnS"
    return f"trn{aa}"


# ── orientation genes: a record starts with the first of these present (mitom_orient) ────────────────────────────────
ORIENTATION_ALIAS_GROUPS = [
    ["cox1", "coi", "coxi", "co1", "cytochromecoxidasesubuniti", "cytochromecoxidasesubunit1"],
    ["cox2", "coii", "coxii", "co2", "cytochromecoxidasesubunitii", "cytochromecoxidasesubunit2"],
    ["cox3", "coiii", "coxiii", "co3", "cytochromecoxidasesubunitiii", "cytochromecoxidasesubunit3"],
    ["cob", "cytb", "cyb", "cytochromeb"],
    ["nad1", "nd1", "nadh1", "nadhdehydrogenasesubunit1"],
    ["atp6", "atpase6", "atpsynthasef0subunit6", "atpsynthasesubunit6"],
]


def gene_aliases(gene: str) -> list:
    """Every name key of the orientation gene `gene` (cox1 -> coi, co1, cytochromecoxidasesubuniti ...)."""
    g = name_key(gene)
    for group in ORIENTATION_ALIAS_GROUPS:
        if g in group:
            return list(group)
    return [g]
