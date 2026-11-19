"""EXPERIMENTAL circular-aware MitoM-ML (trna_ml.CIRC_FIX, 2026-11-19): overlaps across the origin and end pairs that
start before position 1. Off by default; the old behaviour (L = 0, circ False) stays as it was."""
import trna_ml


def hit(s, e, model="trnP", bits=50.0, ev=1e-10):
    return {"model": model, "clade": "pan", "depth": 0, "family": "P", "start": s, "end": e, "strand": 1,
            "mfrom": 1, "mto": 60, "clen": 60, "trunc": "no", "bits": bits, "evalue": ev, "gc": 0.4,
            "anticodon": "tgg"}


def test_overlap_across_the_origin():
    a, b = hit(1, 68), hit(1000, 1067)          # b = the same tRNA at L..L+67 of a 999-nt record shifted by L
    assert trna_ml._ov(a, b) < 0                 # old: no overlap
    assert trna_ml._ov(a, b, 999) == 68


def test_loci_merge_pan_and_clade_hit_across_the_origin():
    hits = [hit(1, 68), hit(1000, 1067, model="clade__trnP", bits=60)]
    assert len(trna_ml._loci(hits)) == 2
    assert len(trna_ml._loci(hits, 999)) == 1


def test_end_pairs_start_before_position_one():
    gen = "ACGT" * 250                           # 1,000 nt
    c = {"strand": 1, "start": 1, "_opts": [hit(1, 68)]}
    old = {(p["s"], p["e"]) for p in trna_ml.end_pairs(c, gen)}
    new = {(p["s"], p["e"]) for p in trna_ml.end_pairs(c, gen, circ=True)}
    assert not any(s == 1000 for s, _ in old)    # start 0 (= position L) was lost
    assert (1000, 1068) in new and all(1 <= s <= 1000 for s, _ in new)
    assert old <= new
