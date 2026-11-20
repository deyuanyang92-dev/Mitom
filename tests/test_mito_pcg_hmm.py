"""mito_pcg_hmm.py: nearest references, the relatives' prior, and finding a gene with a profile HMM."""
import random

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import mito_orf
import mito_pcg_hmm as P

LIN = lambda *x: ["Eukaryota", "Metazoa", *x]          # noqa: E731


def test_nearest_takes_the_longest_shared_lineage_and_leaves_out_genera():
    refs = [("A.1", LIN("Annelida", "Sabellida", "Serpulidae", "Serpula_x"), "MA"),
            ("B.1", LIN("Annelida", "Sabellida", "Sabellidae", "Sabella_y"), "MB"),
            ("C.1", LIN("Mollusca", "Bivalvia", "Mytilus_z"), "MC")]
    q = LIN("Annelida", "Sabellida", "Serpulidae", "Hydroides_q")
    assert [r[0] for r in P.nearest(refs, q, 2)] == ["A.1", "B.1"]
    assert [r[0] for r in P.nearest(refs, q, 2, exclude_genus={"Serpula"})] == ["B.1", "C.1"]
    assert [r[0] for r in P.nearest(refs, q, 3, exclude_acc={"B"})] == ["A.1", "C.1"]


def test_prior_is_low_for_a_gene_the_relatives_lack():
    recs = [(f"N{i}.1", LIN("Nematoda", "Chromadorea", f"G{i}_s"), set(P.PCG) - {"atp8"}) for i in range(20)]
    recs += [(f"S{i}.1", LIN("Annelida", "Sabellida", f"H{i}_s"), set(P.PCG)) for i in range(20)]
    refs = {"_records": recs}
    p_nem, have, of = P.prior(refs, "atp8", LIN("Nematoda", "Chromadorea", "Q_q"), set(), set(), k=10)
    assert have == 0 and p_nem < 0.05 and p_nem / (1 - p_nem) < 0.06      # E threshold x <= 0.06
    p_sab, have, _ = P.prior(refs, "atp8", LIN("Annelida", "Sabellida", "Q_q"), set(), set(), k=10)
    assert have == 10 and p_sab > 0.9


def _protein(n, seed):
    rnd = random.Random(seed)
    return "M" + "".join(rnd.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(n - 1))


def test_a_divergent_gene_is_found_with_the_profile_of_its_relatives(tmp_path):
    """atp8-like protein: 12 relatives mutated at 25 % of the positions, the record's copy at another 25 %; the
    record holds it as the only ORF in unannotated sequence; found at its place."""
    pytest.importorskip("Bio")
    if not (mito_orf.find_tool("mafft") and mito_orf.find_tool("hmmbuild") and mito_orf.find_tool("hmmsearch")):
        pytest.skip("mafft / HMMER not installed")
    rnd = random.Random(3)
    core = _protein(60, 1)

    def mutate(p, f):
        return "".join(rnd.choice("ACDEFGHIKLMNPQRSTVWY") if i and rnd.random() < f else c for i, c in enumerate(p))
    refs = {"atp8": [(f"R{i}.1", LIN("Annelida", "Sabellida", f"G{i}_x"), mutate(core, 0.25)) for i in range(12)],
            "_records": [(f"R{i}.1", LIN("Annelida", "Sabellida", f"G{i}_x"), set(P.PCG)) for i in range(12)]}
    back = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT",
            "K": "AAA", "L": "CTT", "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT", "S": "TCT",
            "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT"}
    gene = "".join(back[a] for a in mutate(core, 0.25)) + "TAA"
    fill = "".join(rnd.choice(["TAA", "TAG", "ATAG", "CTAA"]) for _ in range(150))
    seq = fill + gene + fill
    rec = SeqRecord(Seq(seq), id="Q.1", annotations={"topology": "linear", "molecule_type": "DNA",
                                                       "taxonomy": LIN("Annelida", "Sabellida"), "organism": "Hq x"})
    rec.features = [SeqFeature(FeatureLocation(0, len(seq), 1), type="source"),
                    SeqFeature(FeatureLocation(0, 30, 1), type="CDS", qualifiers={"gene": ["COX1"],
                                                                                   "transl_table": ["5"]})]
    rows = P.search_record(rec, refs, 5, ["atp8"], k=12)
    assert rows[0]["status"] == "found" and rows[0]["start"] == len(fill) + 1 and rows[0]["end"] == len(fill) + len(gene)
    cds = [f for f in rec.features if f.type == "CDS" and f.qualifiers["gene"] == ["ATP8"]]
    assert cds and P.NOTE in cds[0].qualifiers["note"][0]
    tr = cds[0].qualifiers.get("translation", [""])[0]                  # 2026-11-20 review: it was never written
    assert tr and "*" not in tr and len(tr) == (len(gene) - 3) // 3


def test_overlap_across_the_origin():
    """A gene found across the origin (end > n) overlaps one at the start of the sequence."""
    n = 1000
    assert P.overlaps(950, 1030, 10, 40, n)             # 951..1030 runs over the origin onto 11..30
    assert P.overlaps(10, 40, 950, 1030, n)
    assert not P.overlaps(950, 1030, 100, 200, n)
    assert not P.overlaps(100, 200, 300, 400, n) and P.overlaps(100, 300, 250, 400, n)
