"""Tests for benchmark/mito_bench_eval.py (agreement with the NCBI annotation, obvious errors). No external tools.

Run: python3 -m pytest -q tests/test_mito_bench_eval.py
"""
import sys
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "benchmark"))
import mito_bench_eval as me  # noqa: E402


def f(kind, s, e, strand=1, **q):
    return SeqFeature(FeatureLocation(s, e, strand), type=kind, qualifiers={k: [v] for k, v in q.items()})


def rec(feats, seq="A" * 5000):
    r = SeqRecord(Seq(seq), id="X")
    r.features = feats
    return r


def classes(truth_feats, pred_feats):
    rows = me.match(me.norm_features(rec(truth_feats)), me.norm_features(rec(pred_feats)))
    return sorted((r["name"], r["class"]) for r in rows)


def test_norm_names():
    fs = me.norm_features(rec([f("CDS", 0, 99, gene="COX1a"), f("CDS", 200, 299, product="cytochrome b"),
                               f("rRNA", 400, 999, product="16S ribosomal RNA"),
                               f("tRNA", 1000, 1066, product="tRNA-Leu", anticodon="(pos:1030..1032,aa:Leu,seq:tag)")]))
    assert [(x.kind, x.name, x.anticodon) for x in fs] == [("PCG", "cox1", ""), ("PCG", "cob", ""),
                                                           ("rRNA", "rrnL", ""), ("tRNA", "L", "1")]


def test_exact_and_le3():
    t = [f("CDS", 100, 1600, gene="cox1")]
    assert classes(t, [f("CDS", 100, 1600, gene="COX1")]) == [("cox1", "exact")]
    assert classes(t, [f("CDS", 102, 1600, gene="cox1")]) == [("cox1", "le3")]
    assert classes(t, [f("CDS", 110, 1600, gene="cox1")]) == [("cox1", "gt3")]


def test_strand_error():
    assert classes([f("CDS", 100, 1600, gene="cox1")], [f("CDS", 100, 1600, -1, gene="cox1")]) == [("cox1", "strand_error")]


def test_missed_and_fp():
    got = classes([f("CDS", 100, 1600, gene="cox1")], [f("CDS", 3000, 3500, gene="nad6")])
    assert got == [("cox1", "missed"), ("nad6", "fp")]


def test_trnL_name_error():
    t = [f("tRNA", 1000, 1066, product="tRNA-Leu", anticodon="(pos:1030..1032,aa:Leu,seq:tag)")]
    p = [f("tRNA", 1000, 1066, gene="trnL2", product="tRNA-Leu", note="anticodon: taa")]
    assert classes(t, p) == [("L", "name_error")]


def test_no_anticodon_unverifiable():
    t = [f("tRNA", 1000, 1066, product="tRNA-Ser")]
    p = [f("tRNA", 1000, 1066, gene="trnS2(tga)")]
    assert classes(t, p) == [("S", "name_unverifiable")]


CODE5 = 5


def cds_rec(seq, **q):
    return rec([f("CDS", 0, len(seq), gene="cox3", **q)], seq=seq + "A" * 100)


def test_internal_stop():
    checks = me.pcg_checks(cds_rec("ATG" + "TAA" + "GCC" * 10 + "TAA"), CODE5)
    assert [c["check"] for c in checks] == ["internal_stop"]


def test_incomplete_stop_ok():
    assert me.pcg_checks(cds_rec("ATG" + "GCC" * 10 + "T"), CODE5) == []
    assert me.pcg_checks(cds_rec("ATG" + "GCC" * 10 + "TA"), CODE5) == []


def test_bad_start():
    assert [c["check"] for c in me.pcg_checks(cds_rec("CCC" + "GCC" * 10 + "TAA"), CODE5)] == ["bad_start"]


def test_unsupported_excluded():
    obv = me.obvious_errors([], [], {"status": "unsupported"})
    assert obv == []


def test_failed_is_obvious_error():
    obv = me.obvious_errors([], [], {"status": "failed", "error": "boom"})
    assert [o["error"] for o in obv] == ["failed"]


def test_obvious_from_rows():
    rows = [{"kind": "PCG", "name": "nad5", "class": "missed", "delta": ""},
            {"kind": "tRNA", "name": "L", "class": "name_error", "delta": ""},
            {"kind": "PCG", "name": "cox1", "class": "gt3", "delta": 45},
            {"kind": "PCG", "name": "cox2", "class": "gt3", "delta": 12},
            {"kind": "tRNA", "name": "W", "class": "missed", "delta": ""}]
    obv = me.obvious_errors(rows, [{"gene": "cox3", "check": "internal_stop"}], {"status": "ok"})
    assert sorted(o["error"] for o in obv) == ["boundary_gt30", "internal_stop", "missing_gene", "name_error"]


def test_duplicate_gene():
    p = me.norm_features(rec([f("CDS", 100, 600, gene="nad1"), f("CDS", 2000, 2500, gene="nad1")]))
    assert me.duplicates(p) == ["nad1"]


def test_ser_same_codon_box_matches():
    """NC_009732.1: NCBI anticodon 'act' vs tools 'gct' - both read Ser AGN codons -> same gene, not a name error."""
    t = [f("tRNA", 1000, 1062, product="tRNA-Ser", anticodon="(pos:1020..1022,aa:Ser,seq:act)")]
    p = [f("tRNA", 1000, 1062, gene="trnS1", anticodon="gct")]
    assert classes(t, p) == [("S", "exact")]
    p = [f("tRNA", 1000, 1062, gene="trnS2", anticodon="tga")]
    assert classes(t, p) == [("S", "name_error")]


def test_wrapped_truth_boundary():
    """Review C2: genes across the origin keep circular coordinates (no min..max span over the genome)."""
    from Bio.SeqFeature import CompoundLocation
    def w(a, b):
        return SeqFeature(CompoundLocation([FeatureLocation(a, 5000, 1), FeatureLocation(0, b, 1)]), type="CDS",
                          qualifiers={"gene": ["cox1"]})
    t = [w(4000, 100)]
    assert classes(t, [w(4040, 100)]) == [("cox1", "gt3")]
    assert classes(t, [w(4000, 100)]) == [("cox1", "exact")]
    rows = me.match(me.norm_features(rec(t)), me.norm_features(rec([f("CDS", 4000, 5000, gene="cox1")])))
    assert [(r["class"], r["delta"]) for r in rows] == [("gt3", 100)]


def test_failed_record_truth_counts_missed():
    """Review I3: a failed run counts its NCBI features as missed (recall), unsupported counts nothing."""
    truth = me.norm_features(rec([f("CDS", 100, 600, gene="nad1")]))
    assert [r["class"] for r in me.record_rows(truth, None, {"status": "failed"})] == ["missed"]
    assert me.record_rows(truth, None, {"status": "unsupported"}) == []


def test_default_cutoff_is_refseq89():
    """Review I4: post_cutoff stratum on by default (refseq89m = RefSeq release 89, 2018-07-09)."""
    assert me.DEFAULT_CUTOFF == "2018-07-09"


def test_table2asn_no_output_raises(tmp_path, monkeypatch):
    """Review I7: a validator crash must not look like zero errors."""
    import pytest
    monkeypatch.setattr(me.subprocess, "run", lambda *a, **k: None)
    p = tmp_path / "x.gbf"
    p.write_text("x")
    with pytest.raises(RuntimeError):
        me.table2asn_errors({"NC_1.1": p}, tmp_path / "w")


def test_ls_class_codon_or_anticodon_or_number():
    assert me.ls_class("L", "tag", "", 5) == "1" and me.ls_class("L", "taa", "", 5) == "2"       # anticodons
    assert me.ls_class("L", "cta", "", 5) == "1" and me.ls_class("L", "tta", "", 5) == "2"       # NCBI trnL1(cta) = codon
    assert me.ls_class("S", "tct", "", 5) == "1" and me.ls_class("S", "gct", "", 5) == "1"       # AGN anticodons
    assert me.ls_class("S", "tga", "", 5) == "2" and me.ls_class("S", "tca", "", 5) == "2"       # anticodon / codon UCN
    assert me.ls_class("L", "xyz", "", 5) == "" and me.ls_class("L", "", "2", 5) == "2" and me.ls_class("F", "gaa", "", 5) == ""
    t = [f("tRNA", 1000, 1066, gene="trnL1(cta)", product="tRNA-Leu")]
    p = [f("tRNA", 1000, 1066, gene="trnL1", product="tRNA-Leu", note="anticodon: UAG")]
    assert classes(t, p) == [("L", "exact")]
    t = [f("tRNA", 1000, 1066, gene="trnS(tca)", product="tRNA-Ser")]
    p = [f("tRNA", 1000, 1066, gene="trnS2", product="tRNA-Ser", note="anticodon: UGA")]
    assert classes(t, p) == [("S", "exact")]


def test_ls_class_triplet_before_number():
    """The anticodon / codon decides L1 / L2, S1 / S2; the submitter's number only when there is no triplet
    (vertebrate records number L1 / L2 the other way round)."""
    assert me.ls_class("L", "taa", "1", 5) == "2"          # anticodon UAA = codon UUA (UUR) -> L2 despite 'trnL1'
    assert me.ls_class("L", "", "1", 5) == "1"
    assert me.ls_class("S", "xyz", "2", 5) == "2"
    assert me.ls_class("L", "", "", 5) == ""
