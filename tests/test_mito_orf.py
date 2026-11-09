"""mito_orf.py: ORFs of unknown function - chance length, kept open by selection, tandem repeats, shadows, ends."""
import random

import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import mito_orf as m

NONSTOP = [c for c in (a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT")
           if c not in ("TAA", "TAG", "TTA", "CTA", "TCA")]   # no stop on this strand nor (table 5) the other


def filler(n, seed=1):
    """AT-rich sequence with stop codons in every frame of both strands every few codons."""
    rnd = random.Random(seed)
    out = []
    while len(out) < n:
        out += list(rnd.choice(["TAA", "TTA", "ATAG", "CTAA", "TTAG", "A", "T"]))
    return "".join(out[:n])


def orf_nt(codons, seed=2, stop="TAA"):
    rnd = random.Random(seed)
    return "ATG" + "".join(rnd.choice(NONSTOP) for _ in range(codons - 1)) + stop


def record(seq, feats=(), rid="R1"):
    r = SeqRecord(Seq(seq), id=rid, annotations={"topology": "linear", "molecule_type": "DNA"})
    r.features = [SeqFeature(FeatureLocation(0, len(seq), 1), type="source")] + [
        SeqFeature(FeatureLocation(s, e, st), type=t, qualifiers={"gene": [g]}) for t, s, e, st, g in feats]
    return r


def test_stop_share_counts_single_changes_to_a_stop():
    assert m.stop_share("TAC", 5) == pytest.approx(2 / 9)         # TAC -> TAA, TAG
    assert m.stop_share("GGG", 5) == 0.0


def test_third_position_test_tells_coding_from_even_differences():
    assert m.third_position_p([10, 5, 45]) < 1e-6                # differences mostly at position 3
    assert m.third_position_p([20, 20, 20]) > 0.3                 # spread evenly: non-coding


def test_tandem_mask_marks_repeat_arrays_only():
    rep = "ATTGTTTGTTT" * 12
    assert m.tandem_mask(rep).mean() > 0.95
    assert m.tandem_mask(orf_nt(80)).mean() < 0.2


def test_candidate_ends_with_an_incomplete_stop_before_a_trna():
    """The frame runs on into a tRNA without a stop: the ORF ends before it with T / TA (mitochondrial punctuation)."""
    body = orf_nt(80, stop="") + "TA"                            # 80 codons, then TA
    seq = filler(300) + body + "CCCGGGCCCGGGAAAGGGCCC" * 3 + filler(300, 3)   # (TA + C: no stop)
    s = 300
    e = s + len(body)
    rec = record(seq, [("tRNA", e, e + 63, 1, "trnA")])
    cs = [o for o in m.candidates(rec, 5) if o.strand == 1 and o.start == s]
    assert cs and cs[0].writable and cs[0].end == e and cs[0].nt.endswith("TA")
    o = cs[0]
    o.cls, o.best = "conserved", {"name": "x", "src": "run", "evalue": 1e-30, "identity": 50.0, "subs": 200,
                                  "p_open": 1e-4, "nt_diffs": 150, "pos3": 0.6, "p_pos3": 1e-8}
    m.write_orf(rec, o, 5)
    cds = [f for f in rec.features if f.type == "CDS"][0]
    assert cds.qualifiers["product"] == ["hypothetical protein"] and "transl_except" in cds.qualifiers
    assert m.is_orf(cds)


def test_repeat_array_is_not_a_candidate():
    seq = filler(300) + "ATG" + "GTTTGTTTATT" * 30 + filler(300, 4)
    rec = record(seq)
    assert not [o for o in m.candidates(rec, 2) if o.codons > 60]


def _hit(name="h", src="run", subs=200, pos3=0.6, p_pos3=1e-8, evalue=1e-30):
    return {"name": name, "src": src, "evalue": evalue, "identity": 50.0, "length": 100, "mismatches": 50,
            "subs": subs, "nt_diffs": 150, "pos3": pos3, "p_pos3": p_pos3}


def test_classify_needs_both_open_and_third_position():
    o = m.Orf("R", 0, 300, 1, "Y" * 99, "TAC" * 99, "Y" * 99, "TAC" * 99 + "TAA", 1.0, True)
    o.hits = [_hit(subs=300)]                                       # many differences, mostly at position 3
    m.classify(o, 5)
    assert o.cls == "conserved"
    o2 = m.Orf("R", 0, 300, 1, "Y" * 99, "TAC" * 99, "Y" * 99, "TAC" * 99 + "TAA", 1.0, True)
    o2.hits = [_hit(subs=300, pos3=0.34, p_pos3=0.4)]               # open, but differences spread evenly
    m.classify(o2, 5)
    assert o2.cls == ""
    o3 = m.Orf("R", 0, 300, 1, "Y" * 99, "TAC" * 99, "Y" * 99, "TAC" * 99 + "TAA", 1e-5, True)
    o3.hits = [_hit(subs=3)]                                        # close relative: no proof, but long
    m.classify(o3, 5)
    assert o3.cls == "long_unconfirmed"
    o4 = m.Orf("R", 0, 300, 1, "Y" * 99, "TAC" * 99, "Y" * 99, "TAC" * 99 + "TAA", 1.0, True)
    o4.hits = [_hit(name="cox1", src="family")]
    m.classify(o4, 5)
    assert o4.cls == "standard_fragment"


def test_overlapping_orfs_keep_the_strongest_and_report_mode_writes_none(monkeypatch):
    body = orf_nt(120)
    seq = filler(400) + body + filler(400, 5)
    def fake_search(batch, *a, **k):
        for o in batch:
            o.hits = [_hit(subs=400 if o.strand == 1 else 250)]
    monkeypatch.setattr(m, "search", fake_search)
    rec = record(seq)
    orfs = m.annotate([rec], 5, "yes", library=None)
    written = [o for o in orfs if o.written]
    assert len(written) == 1 and written[0].strand == 1 and written[0].start == 400
    assert not [o for o in orfs if not o.written and o.start < written[0].end and written[0].start < o.end]
    rec2 = record(seq)
    assert not any(o.written for o in m.annotate([rec2], 5, "report", library=None))
    assert m.warnings(m.annotate([record(seq)], 5, "report", library=None))      # not silent


def test_known_family_with_a_product_gets_its_gene_name(monkeypatch):
    seq = filler(400) + orf_nt(75) + filler(400, 6)
    monkeypatch.setattr(m, "search", lambda batch, *a, **k: [setattr(o, "hits", [_hit(name="atp9", src="family")])
                                                             for o in batch])
    rec = record(seq)
    m.annotate([rec], 4, "yes", library=None)
    cds = [f for f in rec.features if f.type == "CDS"]
    assert cds and cds[0].qualifiers["gene"] == ["atp9"]
    assert cds[0].qualifiers["product"] == ["ATP synthase F0 subunit 9"]


def test_frame_running_into_a_gene_in_frame_is_named():
    """No stop before a protein gene in the same frame: not written, the WARNING names the gene."""
    body = orf_nt(80, stop="")
    gene = orf_nt(100)
    seq = filler(300) + body + "GCT" * 30 + gene + filler(300, 7)
    gs = 300 + len(body)
    rec = record(seq, [("CDS", gs, gs + 90 + len(gene), 1, "nad5")])
    o = [o for o in m.candidates(rec, 5) if o.strand == 1 and o.start == 300][0]
    assert not o.writable and "nad5 (in frame)" in o.runs_into


def test_frame_running_into_a_trna_ends_there_with_an_open_end():
    """tRNA punctuation (the Polydora ORFans end at tRNAs without a stop): the ORF ends at the last whole codon
    before the tRNA, its 3' end open."""
    body = orf_nt(80, stop="") + "G"                              # no T / TA tail
    seq = filler(300) + body + "GCT" * 30 + filler(300, 12)
    ts = 300 + len(body)
    rec = record(seq, [("tRNA", ts, ts + 66, 1, "trnR")])
    o = [o for o in m.candidates(rec, 5) if o.strand == 1 and o.start == 300][0]
    assert o.writable and o.trna_end == "trnR" and o.end == ts - 1
    loc = m.location(o, len(seq))
    assert str(loc).startswith("[300:>")


def test_ncr_note_names_an_orf_not_written():
    from mito_ncr import add_ncr
    rec = record(filler(1000))
    add_ncr(rec, 100, [(200, 500, "contains ORF 201..500 (+, 99 codons): long_unconfirmed, check")])
    notes = [f.qualifiers["note"][0] for f in rec.features if f.type == "misc_feature"]
    assert any("contains ORF 201..500" in n for n in notes)


def test_circular_orf_across_the_origin_is_found():
    body = orf_nt(90)
    seq = body[150:] + filler(600, 8) + body[:150]
    rec = record(seq)
    rec.annotations["topology"] = "circular"
    o = [o for o in m.candidates(rec, 5) if o.strand == 1 and o.codons >= 89]
    assert o and o[0].end > len(seq) and m.location(o[0], len(seq)).parts[1].start == 0
    assert np.isfinite(o[0].e_len)


def test_codon_pattern_rejects_even_differences_and_the_other_strand_frame():
    assert m.third_position_p([71, 35, 117]) < 1e-8               # coding: 2 least, 3 most
    assert m.third_position_p([27, 22, 38]) > 0.05                # nearly even (a fish control region)
    assert m.third_position_p([35, 71, 117]) == 1.0               # position 1 least: the frame on the other strand


def test_fisher_combines_two_tests():
    assert m.fisher(0.01, 0.01) == pytest.approx(1e-4 * (1 - np.log(1e-4)))
    assert m.fisher(1.0, 1.0) == 1.0


def test_mst_bound_adds_up_many_close_relatives():
    """Ten copies, each 5 differences from the query at its own positions: the MST is 50 long, the bound half of it
    (combined_evidence takes the larger of this and one substitution per extra base kind and position)."""
    q = {i: "A" for i in range(100)}
    copies = [q] + [{**q, **{10 * k + j: "G" for j in range(5)}} for k in range(10)]
    assert m.mst_bound(copies) == 25
    assert m.mst_bound([q, dict(q)]) == 0


def test_introns_are_open_to_the_orf_search():
    """A homing endonuclease ORF inside an intron: the joined CDS covers its exons only."""
    from Bio.SeqFeature import CompoundLocation
    body = orf_nt(120)
    seq = filler(300) + "ATG" + "GCT" * 40 + filler(200, 9) + body + filler(200, 10) + "GCT" * 40 + "TAA" + filler(300, 11)
    e1 = (300, 300 + 123)
    i0 = e1[1] + 200
    e2s = i0 + len(body) + 200
    rec = record(seq)
    rec.features.append(SeqFeature(CompoundLocation([FeatureLocation(*e1, 1), FeatureLocation(e2s, e2s + 123, 1)]),
                                   type="CDS", qualifiers={"gene": ["cox1"]}))
    rec.features.append(SeqFeature(FeatureLocation(300, e2s + 123, 1), type="gene", qualifiers={"gene": ["cox1"]}))
    assert [o for o in m.candidates(rec, 5) if o.strand == 1 and o.start == i0]


def test_frame_rank_prefers_the_frame_with_position_2_least_changed():
    a = m.Orf("R", 0, 300, 1, "Y" * 99, "TAC" * 99, "Y" * 99, "TAC" * 100, 1.0, True)
    b = m.Orf("R", 0, 300, -1, "Y" * 99, "TAC" * 99, "Y" * 99, "TAC" * 100, 1.0, True)
    a.cls = b.cls = "conserved"
    a.best = {"d": [71, 35, 117], "pos3": 0.52, "p_joint": 1e-6}
    b.best = {"d": [35, 71, 117], "pos3": 0.52, "p_joint": 1e-9}          # stronger p, wrong pattern
    assert m.frame_rank(a) < m.frame_rank(b)


def test_missing_hmmsearch_is_not_silent(monkeypatch, capsys):
    o = m.Orf("R", 0, 300, 1, "Y" * 99, "TAC" * 99, "Y" * 99, "TAC" * 100, 1.0, True)
    monkeypatch.setattr(m, "find_tool", lambda name, bin_dir=None: None)
    m.pfam_search([o])
    assert "hmmsearch (HMMER) not found" in capsys.readouterr().err and o.hits == []


def test_pfam_family_names_the_product():
    o = m.Orf("R", 0, 300, 1, "Y" * 99, "TAC" * 99, "Y" * 99, "TAC" * 100, 1.0, True)
    o.hits = [{"name": "LAGLIDADG_1", "src": "pfam", "evalue": 1e-20, "identity": 0.0, "length": 0, "mismatches": 0}]
    m.classify(o, 4)
    rec = record("A" * 600)
    m.write_orf(rec, o, 4)
    cds = [f for f in rec.features if f.type == "CDS"][0]
    assert o.cls == "known_family" and cds.qualifiers["product"] == ["LAGLIDADG homing endonuclease"]
