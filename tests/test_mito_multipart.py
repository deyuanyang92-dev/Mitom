"""mito_multipart.py: records of one assembly whose genes are complete only together (multipartite mtDNA)."""
import random

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import mito_multipart as M


def rec(rid, genes, seq_seed, comment="The assembly nxTest1.1 is based on PacBio HiFi"):
    rnd = random.Random(seq_seed)
    n = 16000
    r = SeqRecord(Seq("".join(rnd.choice("ACGT") for _ in range(n))), id=rid,
                  annotations={"comment": comment, "organism": "Testus testus", "molecule_type": "DNA"})
    pos = 0
    for g, ln in genes:
        r.features.append(SeqFeature(FeatureLocation(pos, pos + ln, 1), type="CDS", qualifiers={"gene": [g]}))
        pos += ln + 10
    return r


ALL = [(g, M.USUAL[g]) for g in sorted(M.PCG)]


def test_two_circles_with_complementary_genes_are_parts_of_one_genome():
    a = rec("A.1", [x for x in ALL if x[0] in ("cox1", "cox3", "nad4", "nad4l", "nad6")] + [("nad5", 400)], 1)
    b = rec("B.1", [x for x in ALL if x[0] not in ("cox1", "cox3", "atp8")], 2)
    rows = M.analyse([a, b])
    assert len(rows) == 1 and rows[0]["verdict"].startswith("parts of one genome")
    assert rows[0]["together"] == 12 and rows[0]["missing_together"] == "atp8"
    assert "nad5 (A.1)" in rows[0]["fragments"]


def test_a_refseq_copy_is_not_a_part():
    a = rec("OZ1.1", ALL, 3)
    b = SeqRecord(a.seq.reverse_complement(), id="NW_1.1", annotations=dict(a.annotations))
    rows = M.analyse([a, b])
    assert rows[0]["distinct"] == 1 and rows[0]["verdict"] == "copies of one sequence"


def test_records_of_different_assemblies_are_not_grouped():
    a = rec("A.1", ALL[:6], 4, "The assembly x1 is based")
    b = rec("B.1", ALL[6:], 5, "The assembly x2 is based")
    assert M.analyse([a, b]) == []
