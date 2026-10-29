"""Tests for ncbi_validate.py: the .tbl written for table2asn.

Run: python3 -m pytest -q tests/test_ncbi_validate.py
"""
import random
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import AfterPosition, BeforePosition, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_env  # noqa: E402
import ncbi_validate as nv  # noqa: E402


def _rec(seq_len=600, seed=5):
    rnd = random.Random(seed)
    r = SeqRecord(Seq("".join(rnd.choice("ACGT") for _ in range(seq_len))), id="T1")
    r.annotations.update(molecule_type="DNA", topology="linear", organism="Test sp.")
    return r


def _cds(loc, gene="x"):
    return SeqFeature(loc, type="CDS", qualifiers={"gene": [gene], "product": ["x"], "transl_table": ["5"],
                                                   "codon_start": ["1"]})


def test_tbl_partial_markers_first_column_lt_second_gt_on_minus_strand():
    """Review C finding 1: in a 5-column .tbl '<' always goes on the FIRST column and '>' on the SECOND, whatever the
    strand (table2asn here: '<150<TAB>40' -> complement(40..>150) OK; '>150<TAB>40' -> 'Bad feature interval', feature
    dropped). mito_curate P5 writes a minus-strand 5'-partial CDS as AfterPosition(end), P6 a 3'-partial as
    BeforePosition(start)."""
    r = _rec()
    r.features = [_cds(FeatureLocation(39, AfterPosition(150), -1)), _cds(FeatureLocation(BeforePosition(200), 350, -1))]
    lines = nv.feature_lines(r)
    assert "<150\t40\tCDS" in lines
    assert "350\t>201\tCDS" in lines


def test_tbl_partial_markers_plus_strand_unchanged():
    r = _rec()
    r.features = [_cds(FeatureLocation(BeforePosition(400), 520, 1)), _cds(FeatureLocation(530, AfterPosition(590), 1))]
    lines = nv.feature_lines(r)
    assert "<401\t520\tCDS" in lines
    assert "531\t>590\tCDS" in lines


def test_minus_strand_partial_cds_survives_table2asn(tmp_path):
    """Review C finding 1: a minus-strand partial CDS was dropped by table2asn ('Bad feature interval') and the
    validation reported nothing about the lost feature."""
    t2a = mito_env.find("table2asn")
    if not t2a:
        pytest.skip("table2asn not available")
    r = _rec()
    r.features = [_cds(FeatureLocation(39, AfterPosition(150), -1), "m5"),
                  _cds(FeatureLocation(BeforePosition(200), 350, -1), "m3")]
    nv.validate_record(r, tmp_path, t2a)
    out = SeqIO.read(str(tmp_path / "T1" / "rec.gbf"), "genbank")
    assert sorted(f.qualifiers["gene"][0] for f in out.features if f.type == "CDS") == ["m3", "m5"]


def test_tbl_skips_internal_review_notes():
    """Review C finding 15: tool-internal hints ('... check boundaries manually', 'mito_rescue.py', 'mito_intron.py')
    must not reach the submitted .tbl; ordinary notes stay."""
    r = _rec()
    f = _cds(FeatureLocation(10, 100, 1))
    f.qualifiers["note"] = ["missed by MitoZ; added from MITOS2 by mito_rescue.py - check boundaries manually",
                            "exon-intron structure: Exonerate ... (mito_intron.py)",
                            "start codon not determined"]
    r.features = [f]
    notes = [ln for ln in nv.feature_lines(r) if "\tnote\t" in ln]
    assert notes == ["\t\t\tnote\tstart codon not determined"]


def test_messages_on_features_of_the_intron_steps_are_marked_ours(tmp_path):
    """2026-10-30 (integration audit): the 'ours' column (--only_ours) looked for 'mito_intron' in the notes, which the
    de novo step no longer writes ('... MFannot method') and the site transfer writes only on its CDS: messages on
    their introns were not counted as MitoM's."""
    import sys as _sys
    r = _rec()
    r.features = [SeqFeature(FeatureLocation(100, 200, 1), type="intron", qualifiers={
                      "note": ["group II intron (Rfam RF00029, E 1e-10); MFannot method"]}),
                  SeqFeature(FeatureLocation(300, 400, 1), type="intron", qualifiers={
                      "note": ["group I intron; site transferred from NC_000933.1"]}),
                  SeqFeature(FeatureLocation(450, 500, 1), type="misc_feature", qualifiers={"note": ["submitter's"]})]
    fake = tmp_path / "table2asn"
    fake.write_text(f"#!{_sys.executable}\nimport sys\nout = sys.argv[sys.argv.index('-o') + 1]\n"
                    "open(out[:-4] + '.val', 'w').write(''.join(f'ERROR: valid [SEQ_FEAT.Test] intron {iv} odd\\n' "
                    "for iv in ('101-200', '301-400', '451-500')))\n")
    fake.chmod(0o755)
    rows = nv.validate_record(r, tmp_path / "out", str(fake))
    assert [(row["message"].split()[1], row["ours"]) for row in rows] == [("101-200", "yes"), ("301-400", "yes"),
                                                                          ("451-500", "no")]
