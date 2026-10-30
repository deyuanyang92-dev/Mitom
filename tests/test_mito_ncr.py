"""mito_ncr.py: non-coding regions (stretches no feature covers) written as misc_feature."""
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import mito_ncr as mn  # noqa: E402


def record(topology="circular"):
    """6,000 bp: genes 300..2000, 2050..4000 (with an intron 2500..3000), 4400..5800; free: 0..300, 2000..2050,
    4000..4400, 5800..6000 (with 0..300 one stretch of 500 nt across the origin when circular)."""
    r = SeqRecord(Seq("ACGT" * 1500), id="T1", description="")
    r.annotations.update(molecule_type="DNA", topology=topology)
    r.features = [SeqFeature(FeatureLocation(0, 6000, 1), type="source"),
                  SeqFeature(FeatureLocation(300, 2000, 1), type="CDS", qualifiers={"gene": ["COX1"]}),
                  SeqFeature(CompoundLocation([FeatureLocation(2050, 2500, 1), FeatureLocation(3000, 4000, 1)]),
                             type="CDS", qualifiers={"gene": ["ND5"]}),
                  SeqFeature(FeatureLocation(4400, 5800, -1), type="rRNA", qualifiers={"gene": ["rrnL"]})]
    return r


def test_noncoding_stretches_introns_excluded():
    assert mn.noncoding(record()) == [(4000, 4400), (5800, 6300)]          # one stretch across the origin
    assert mn.noncoding(record("linear")) == [(0, 300), (4000, 4400), (5800, 6000)]
    assert mn.noncoding(record(), 401) == [(5800, 6300)]
    assert mn.noncoding(record(), 0) == []


def test_add_ncr_features_numbered_and_replaced(tmp_path):
    r = record()
    assert mn.add_ncr(r) == 2
    ncr = [f for f in r.features if mn.is_ncr(f)]
    assert [f.qualifiers["note"][0] for f in ncr] == [
        "non-coding region NCR1, 400 bp: no gene annotated (MitoM)",
        "non-coding region NCR2, 500 bp: no gene annotated (MitoM)"]
    assert str(ncr[1].location) == "join{[5800:6000], [0:300]}"            # across the origin
    assert mn.add_ncr(r, 450) == 1 and len([f for f in r.features if mn.is_ncr(f)]) == 1   # old NCRs replaced
    SeqIO.write(r, str(tmp_path / "x.gbf"), "genbank")                      # a valid GenBank record
    back = SeqIO.read(str(tmp_path / "x.gbf"), "genbank")
    assert [mn.is_ncr(f) for f in back.features].count(True) == 1


def test_cli_folder(tmp_path):
    (tmp_path / "in").mkdir()
    SeqIO.write(record(), str(tmp_path / "in" / "a.gbf"), "genbank")
    assert mn.main(["-i", str(tmp_path / "in"), "-o", str(tmp_path / "out"), "--min_len", "100"]) == 0
    out = SeqIO.read(str(tmp_path / "out" / "a.gbf"), "genbank")
    assert sum(mn.is_ncr(f) for f in out.features) == 2
