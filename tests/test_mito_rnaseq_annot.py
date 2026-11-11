"""mito_rnaseq_annot: split-read junctions, the overlap cut rule, RNA indels, frameshift joins."""
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import mito_rnaseq_annot as m


def _cds(rec, name, loc, table="5"):
    rec.features.append(SeqFeature(loc, type="CDS", qualifiers={"gene": [name], "product": ["x"],
                                                                 "transl_table": [table], "codon_start": ["1"]}))


def test_tap_counts_split_read_gap_and_indels():
    n = 400
    lines = ["r1\t0\tc\t11\t60\t40M60S\t*\t0\t0\t" + "A" * 100 + "\t*\tSA:Z:c,201,+,40S60M,60,0;",
             "r2\t0\tc\t21\t60\t10M1I20M\t*\t0\t0\t" + "A" * 31 + "\t*"]
    t = m.Tap(lines, n)
    list(t)
    assert t.junctions == {(49, 200): 1}                    # last base before the gap 0-based 49, first after 200
    assert t.ins == {(29, 1): 1}
    assert m.merge_junctions({(49, 200): 8, (50, 201): 3}) == [(49, 200, 11)]


def test_overlap_cut_moves_stop_to_trna_start_when_polya_says_so():
    cds = "ATG" + "GCT" * 19 + "TAC" + "TAA"                # 66 nt; the tRNA starts at 0-based 62 (overlap 4 nt)
    rec = SeqRecord(Seq(cds + "G" * 60), id="t")
    _cds(rec, "COX3", FeatureLocation(0, 66, 1))
    rec.features.append(SeqFeature(FeatureLocation(62, 125, 1), type="tRNA", qualifiers={"gene": ["trnG"]}))
    sites = [{"position": "62", "strand": "+", "reads": "30", "internal_priming": ""}]    # mRNA ends ...TA
    rows = []
    m.step_overlap_cut(rec, sites, stranded=False, run="SRR1", rows=rows)
    f = rec.features[0]
    assert rows[0]["verdict"] == "corrected" and int(f.location.end) == 62
    assert f.qualifiers["transl_except"] == ["(pos:61..62,aa:TERM)"]


def test_overlap_of_one_or_two_nt_left_to_step3():
    cds = "ATG" + "GCT" * 20 + "TAA"
    rec = SeqRecord(Seq(cds + "G" * 60), id="t")
    _cds(rec, "COX3", FeatureLocation(0, 66, 1))
    rec.features.append(SeqFeature(FeatureLocation(65, 125, 1), type="tRNA", qualifiers={"gene": ["trnG"]}))
    rows = []
    m.step_overlap_cut(rec, [{"position": "65", "strand": "+", "reads": "30", "internal_priming": ""}], False, "", rows)
    assert rows == [] and int(rec.features[0].location.end) == 66


def test_single_frame_orf_found_for_frameshift_join():
    # true gene: ATG + 30 codons + TAA in frame 0 at 10; a join annotated from 13 skipping one base
    g = "CCCCCCCCCC" + "ATG" + "GCAGCC" * 15 + "TAA" + "CCCCCCCCCC"
    rec = SeqRecord(Seq(g), id="t")
    _cds(rec, "ATP8", CompoundLocation([FeatureLocation(13, 40, 1), FeatureLocation(41, 106, 1)]))
    st, en, aa = m.single_frame_orf(rec, rec.features[0])
    assert (st, en, aa) == (10, 106, 31)
