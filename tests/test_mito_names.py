"""mito_names.py is the one place for gene names, aliases, products, amino-acid codes and tRNA symbols. Every call site
gives the name it gave before: tests/data/name_expected.tsv.gz holds the result of each call site, computed by the code
of 4855378 (before mito_names became the single source), for
  - the 2,368 gene / product / tRNA qualifier combinations of 10,947 local GenBank files (NCBI records of the benchmark
    sets, MitoZ / MITOS2 / MitoM outputs, the golden inputs and outputs, tests/data; numbers in anticodon / note text
    replaced by 1),
  - every key of every old name table (aliases, synonyms, products, orientation aliases) in six spellings, as /gene and
    as /product of a CDS, rRNA, tRNA and gene feature.
A deliberate change of a name rule regenerates the file and shows the changed names in its diff (2026-10-31: RNR1 / RNR2,
l-rRNA / small subunit ribosomal RNA products, curate's PCG synonyms NADH4L / COXI in gene_name)."""
import argparse
import csv
import gzip
import json
import sys
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_core  # noqa: E402
import mito_curate  # noqa: E402
import mito_intron  # noqa: E402
import mito_intron_ref  # noqa: E402
import mito_names  # noqa: E402
import mito_plot  # noqa: E402
import mito_rescue  # noqa: E402
import mito_rrna_hmm  # noqa: E402
import mito_trna  # noqa: E402
import mitom_orient  # noqa: E402

COLS = ["type", "gene", "product", "label", "anticodon", "note"]


def qualifiers(row):
    q = {k: [row[k]] for k in ("gene", "product", "label", "anticodon") if row[k]}
    if row["note"]:
        q["note"] = row["note"].split(" | ")
    return q


def results(row):
    """{call site: result as text} for one row (the columns of name_expected.tsv.gz)."""
    q = qualifiers(row)
    t = row["type"]
    r = {"names_gene_name": mito_names.gene_name(q, t), "names_trna_name": "|".join(mito_names.trna_name(q)),
         "names_trna_number": mito_names.trna_number(q), "intron_ref_gene_name": mito_intron_ref.gene_name(q, t)}
    f = SeqFeature(FeatureLocation(10, 100, 1), type=t, qualifiers={k: list(v) for k, v in q.items()})
    r["rrna_hmm_gene_of"] = mito_rrna_hmm.gene_of(f) or ""
    rec = SeqRecord(Seq("A" * 300), id="X")
    rec.features = [f]
    r["trna_from_input"] = "|".join(f"{p.aa}:{p.anticodon}" for p in mito_trna.from_input(rec))
    r["curate_n1"] = ""
    if t in ("CDS", "rRNA", "tRNA"):                     # rule N1 on a record with this feature and its gene feature
        g = SeqFeature(FeatureLocation(10, 100, 1), type="gene", qualifiers={"gene": list(q.get("gene", [""]))})
        rec = SeqRecord(Seq("A" * 300), id="X")
        ff = SeqFeature(FeatureLocation(10, 100, 1), type=t, qualifiers={k: list(v) for k, v in q.items()})
        rec.features = [g, ff]
        c = mito_curate.Curator(rec, 5, argparse.Namespace(aggressive_start_genes=""))
        c.names(ff)
        r["curate_n1"] = json.dumps([ff.qualifiers.get(k) for k in ("gene", "product", "note")] +
                                    [g.qualifiers.get("gene"), [x["rule"] for x in c.rows],
                                     sorted(x for x, s in (("weak", c.weak), ("copy", c.copies)) if id(ff) in s)])
    for k in ("gene", "product", "label"):
        r[f"orient_aliases_{k}"] = ",".join(mitom_orient.gene_aliases(row[k])) if row[k] else ""
    r["orient_feature_names"] = ",".join(mitom_orient.feature_gene_names(f))
    cat = mito_plot.category(f)
    r["plot"] = f"{cat}:{mito_plot.label(f, cat)}" if cat else ""
    r["intron_norm_gene"] = mito_intron.norm_gene(row["gene"]) if row["gene"] else ""
    r["intron_product_for"] = mito_intron.product_for(row["gene"]) if row["gene"] else ""
    r["intron_ref_product"] = mito_intron_ref.PRODUCTS.get(row["gene"], row["gene"]) if row["gene"] else ""
    return r


ROWS = list(csv.DictReader(gzip.open(HERE / "data" / "name_expected.tsv.gz", "rt"), delimiter="\t"))


@pytest.fixture(scope="module")
def got():
    return [results(r) for r in ROWS]


@pytest.mark.parametrize("site", [c for c in ROWS[0] if c not in COLS])
def test_every_call_site_gives_the_old_name(site, got):
    bad = [({c: r[c] for c in COLS if r[c]}, r[site], g[site]) for r, g in zip(ROWS, got) if g[site] != r[site]]
    assert not bad, f"{len(bad)} of {len(ROWS)} rows differ at {site}, e.g. {bad[:3]}"


def test_tables_shared_with_the_core_are_one():
    """The core's /gene and /product names and amino acid table are mito_names' own (no copy)."""
    assert mito_core.NAMES is mito_names.CORE_NAMES and mito_core.AA3 is mito_names.AA1_TO_3
    assert mito_rescue.CORE_NAMES is mito_names.CORE_NAMES and mito_rrna_hmm.GENES == {
        g: mito_core.NAMES[g] for g in ("rrnS", "rrnL")}
    assert mito_trna.AA3 == mito_names.AA3 == mito_plot.AA and mito_trna.ONE2THREE == mito_names.AA1_TO_3


def test_trna_symbol_is_the_old_mito_trna_rule():
    def old(aa, ac):
        if aa == "L":
            return "trnL1" if ac.endswith("ag") else "trnL2" if ac.endswith("aa") else "trnL"
        if aa == "S":
            return "trnS1" if ac.endswith("ct") else "trnS2" if ac.endswith("ga") else "trnS"
        return f"trn{aa}"
    codons = [a + b + c for a in "acgt" for b in "acgt" for c in "acgt"] + [""]
    for aa in list(mito_names.AA1_TO_3) + [""]:
        for ac in codons:
            assert mito_names.trna_symbol(aa, ac) == mito_trna.trna_name(aa, ac) == old(aa, ac)


def test_leu_ser_numbers_follow_mitos2_for_every_taxon():
    """User's decision 2026-10-31: follow MITOS2 (2.1.10: L1 = CUN, L2 = UUR, S1 = AGN, S2 = UCN codons) for every taxon,
    vertebrates included, although NCBI's human record numbers them the other way."""
    t = mito_names.trna_symbol
    assert (t("L", "tag"), t("L", "taa"), t("S", "gct"), t("S", "tct"), t("S", "tga")) == \
        ("trnL1", "trnL2", "trnS1", "trnS1", "trnS2")
