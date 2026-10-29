"""Orientation: fallback genes when cox1 is absent, and the 5'-most piece of a split gene.

Run: python3 -m pytest -q tests/test_orientation.py
"""
import argparse
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import mitom_pipeline as bm  # noqa: E402


def rec_with(features):
    r = SeqRecord(Seq("ACGT" * 500), id="R1", name="R1", description="test")
    r.annotations.update(molecule_type="DNA", topology="circular")
    r.features = [SeqFeature(FeatureLocation(a, b, strand=st), type=t, qualifiers={"gene": [g]})
                  for t, g, a, b, st in features]
    return r


def test_aliases():
    assert "coii" in bm.gene_aliases("cox2") and "cytb" in bm.gene_aliases("cob")
    assert bm.gene_aliases("COX1") == bm.gene_aliases("COI")


def test_orientation_genes_from_args():
    a = argparse.Namespace(which_gene_first="cox1", fallback_genes="cox2,cob")
    assert bm.orientation_genes(a) == "cox1,cox2,cob"
    assert bm.orientation_genes(argparse.Namespace(which_gene_first="cox1", fallback_genes="none")) == "cox1"


def test_fallback_when_cox1_absent():
    r = rec_with([("CDS", "COII", 300, 900, 1), ("CDS", "CYTB", 1000, 1500, 1)])
    feat, gene = bm.find_orientation_feature(r, "cox1,cox2,cob")
    assert gene == "cox2" and int(feat.location.start) == 300


def test_split_gene_takes_five_prime_piece():
    plus = rec_with([("CDS", "cox1", 1200, 1700, 1), ("CDS", "cox1", 200, 700, 1)])   # exon 2 listed first
    assert int(bm.find_orientation_feature(plus, "cox1")[0].location.start) == 200
    minus = rec_with([("CDS", "cox1", 200, 700, -1), ("CDS", "cox1", 1200, 1700, -1)])
    assert int(bm.find_orientation_feature(minus, "cox1")[0].location.end) == 1700


def test_nothing_found():
    assert bm.find_orientation_feature(rec_with([("CDS", "nad5", 10, 400, 1)]), "cox1,cox2") == (None, "")


def test_reorient_uses_fallback(tmp_path):
    bm.ORIENTATION_FALLBACK.clear()
    r = rec_with([("CDS", "cox2", 300, 900, 1), ("CDS", "nad5", 1000, 1500, 1)])
    gb_in = tmp_path / "in.gbf"
    SeqIO.write(r, str(gb_in), "genbank")
    n_ok, n_no, _ = bm.reorient_one_gb(gb_in, tmp_path / "out.gbf", tmp_path / "no.gbf", "cox1,cox2", True)
    assert (n_ok, n_no) == (1, 0)
    out = SeqIO.read(str(tmp_path / "out.gbf"), "genbank")
    assert str(out.seq)[:20] == str(r.seq)[300:320]               # sequence now starts at cox2
    assert bm.ORIENTATION_FALLBACK == {"R1": "cox2"}


def test_reorient_after_introns_starts_with_the_first_exon(tmp_path):
    """2026-10-29: the intron step joined cox1 across the origin (the core had annotated exon 2, so the record started
    there and exon 1 + the intron lay at its end); the record is rotated again to start with exon 1."""
    import argparse
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    L = 3000
    rec = SeqRecord(Seq("".join("ACGT"[(i * 7) % 4] for i in range(L))), id="S1", name="S1", description="")
    rec.annotations.update(molecule_type="DNA", topology="circular")
    loc = CompoundLocation([FeatureLocation(2700, 3000, 1), FeatureLocation(0, 900, 1)])     # exon 1 at the end
    rec.features = [SeqFeature(FeatureLocation(0, L, 1), type="source"),
                    SeqFeature(loc, type="gene", qualifiers={"gene": ["COX1"]}),
                    SeqFeature(loc, type="CDS", qualifiers={"gene": ["COX1"]}),
                    SeqFeature(FeatureLocation(1000, 1600, 1), type="CDS", qualifiers={"gene": ["COX2"]})]
    (tmp_path / "f").mkdir()
    SeqIO.write(rec, str(tmp_path / "f" / "S1.gbf"), "genbank")
    a = argparse.Namespace(which_gene_first="cox1", fallback_genes="cox2,cox3,cob,nad1")
    assert bm.reorient_after_introns(a, tmp_path / "f") == 1
    out = SeqIO.read(str(tmp_path / "f" / "S1.gbf"), "genbank")
    cox1 = next(f for f in out.features if f.type == "CDS" and f.qualifiers["gene"] == ["COX1"])
    cox2 = next(f for f in out.features if f.type == "CDS" and f.qualifiers["gene"] == ["COX2"])
    assert int(cox1.location.parts[0].start) == 0 and str(cox1.extract(out.seq)) == str(loc.extract(rec.seq))
    assert int(cox2.location.start) == 1300
    assert bm.reorient_after_introns(a, tmp_path / "f") == 0                  # already starts with exon 1


def test_rotation_keeps_partial_ends():
    """2026-10-29: rotation dropped < / > of the trans-spliced placozoan cox1 pieces (NC_015309.1: table2asn
    StartCodon / NoStop / BadProteinStart ERRORs); a partial end stays partial, also when the piece is cut by the
    new origin."""
    from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation
    r = rec_with([])
    five = CompoundLocation([FeatureLocation(BeforePosition(100), 200, 1), FeatureLocation(500, 700, 1)])
    three = FeatureLocation(1500, AfterPosition(1800), strand=1)
    r.features = [SeqFeature(five, type="CDS", qualifiers={"gene": ["COX1"]}),
                  SeqFeature(three, type="CDS", qualifiers={"gene": ["COX1"]})]
    out = bm.rotate_record(r, 600)
    a, b = (next(f for f in out.features if int(f.location.parts[0].start) == x) for x in (1500, 900))
    assert isinstance(a.location.parts[0].start, BeforePosition) and int(a.location.parts[-1].end) == 100
    assert isinstance(b.location.end, AfterPosition) and int(b.location.end) == 1200
    cut = bm.rotate_record(r, 1600)                                            # the 3'-partial piece spans the origin
    c = next(f for f in cut.features if int(f.location.parts[0].start) == 1900)
    assert isinstance(c.location.parts[-1].end, AfterPosition) and int(c.location.parts[-1].end) == 200


def test_reverse_complement_keeps_and_swaps_partial_ends():
    """2026-10-30 (integration audit): turning a record to the plus strand (cox1 on the minus strand) dropped the
    partial ends < / > of genes cut by the ends of a linear record; the 5' end of a gene stays its 5' end."""
    from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation
    r = rec_with([])
    L = len(r.seq)
    five = FeatureLocation(BeforePosition(0), 200, strand=1)                      # 5'-partial, plus strand
    three = FeatureLocation(1500, AfterPosition(L), strand=-1)                   # 5'-partial on the minus strand
    joined = CompoundLocation([FeatureLocation(BeforePosition(300), 400, 1), FeatureLocation(500, AfterPosition(700), 1)])
    r.features = [SeqFeature(loc, type="CDS", qualifiers={"gene": ["COX1"]}) for loc in (five, three, joined)]
    rc = bm.reverse_complement_record(r)
    a, b, c = (f.location for f in rc.features)
    assert (int(a.start), int(a.end), a.strand) == (L - 200, L, -1) and isinstance(a.end, AfterPosition)
    assert (int(b.start), int(b.end), b.strand) == (0, L - 1500, 1) and isinstance(b.start, BeforePosition)
    # minus-strand join in transcript order: the 3'-partial piece first (its 3' end is its start), the 5'-partial last
    assert isinstance(c.parts[0].start, BeforePosition) and isinstance(c.parts[-1].end, AfterPosition)
    back = bm.reverse_complement_record(rc)
    assert [str(f.location) for f in back.features] == [str(f.location) for f in r.features]


def _triplet(rec, text):
    """The 3 nt named by an /anticodon (pos:...) value, read on its strand."""
    import re as _re
    m = _re.search(r"pos:(complement\()?(?:join\()?(\d+(?:\.\.\d+)?(?:,\d+(?:\.\.\d+)?)*)", text)
    nums = [int(x) for x in _re.findall(r"\d+", m.group(2))]
    ranges = [tuple(int(x) for x in r.split("..")) if ".." in r else (int(r), int(r)) for r in m.group(2).split(",")]
    s = "".join(str(rec.seq[a - 1:b]) for a, b in ranges)
    from Bio.Seq import Seq as _Seq
    return str(_Seq(s).reverse_complement()) if m.group(1) else s, nums


def test_rotation_and_reverse_complement_move_transl_except_and_anticodon():
    """2026-10-30 (integration audit): /transl_except and /anticodon (pos:...) kept their old positions when a record
    was rotated (re-orientation, the rotation after the intron step, `mitom import` of NCBI records) or turned to the
    plus strand: the NCBI submission files then named the wrong bases."""
    from Bio.SeqFeature import CompoundLocation
    seq = "".join("ACGT"[(i * 7 + i // 3) % 4] for i in range(3000))
    seq = seq[:2994] + "ACTAAT"                                          # the CDS ends with a T at the record end
    r = SeqRecord(Seq(seq), id="R1", name="R1", description="")
    r.annotations.update(molecule_type="DNA", topology="circular")
    L = len(seq)
    trna = FeatureLocation(1000, 1070, -1)
    anti = str(Seq(seq[1030:1033]).reverse_complement()).lower()
    r.features = [SeqFeature(FeatureLocation(0, L, 1), type="source"),
                  SeqFeature(FeatureLocation(2500, L, 1), type="CDS", qualifiers={
                      "gene": ["ND5"], "transl_except": [f"(pos:{L},aa:TERM)"]}),
                  SeqFeature(trna, type="tRNA", qualifiers={
                      "gene": ["trnL1"], "anticodon": [f"(pos:complement(1031..1033),aa:Leu,seq:{anti})"]}),
                  SeqFeature(FeatureLocation(200, 1500, 1), type="CDS", qualifiers={"gene": ["COX1"]})]

    def check(rec):
        cds = next(f for f in rec.features if f.qualifiers.get("gene") == ["ND5"])
        t = next(f for f in rec.features if f.type == "tRNA")
        last = cds.location.parts[-1]
        pos = int(last.end) if last.strand == 1 else int(last.start) + 1
        te = cds.qualifiers["transl_except"][0]
        assert te in (f"(pos:{pos},aa:TERM)", f"(pos:complement({pos}),aa:TERM)"), te
        trip, nums = _triplet(rec, t.qualifiers["anticodon"][0])
        assert trip.lower() == anti and all(int(t.location.start) < n <= int(t.location.end) for n in nums)

    check(r)
    rot = bm.rotate_record(r, 200)                                     # cox1 first
    check(rot)
    check(bm.rotate_record(r, 1032))                                   # the anticodon is cut by the new origin
    check(bm.reverse_complement_record(rot))
    check(bm.rotate_record(bm.reverse_complement_record(r), 1500))
