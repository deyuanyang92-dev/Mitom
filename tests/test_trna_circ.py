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


def test_mitfi_skipped_only_when_it_cannot_change_the_result(monkeypatch):
    """core_trnas runs MiTFi only when a kept MitoM-ML call has p < 0.9, a family is missing or the record is a
    vertebrate; otherwise the result with MiTFi is the same and its scan is skipped (2026-11-19)."""
    import logging
    from pathlib import Path
    import mito_trna
    from mito_common import CoreFeature
    full = [CoreFeature("tRNA", n, 100 * i, 100 * i + 65, 1, "MitoM-ML", score=0.95)
            for i, n in enumerate(mito_trna.TRNA_NAMES, 1)]
    calls = []
    monkeypatch.setattr(mito_trna, "annotate_trnas", lambda *a, **k: calls.append(1) or [])
    log = logging.getLogger("t")

    def run(ml, lineage="Metazoa;Annelida"):
        calls.clear()
        monkeypatch.setattr(mito_trna, "annotate_trnas_ml", lambda *a: list(ml))
        out = mito_trna.core_trnas("A" * 16000, "A" * 16000, Path("g.fa"), Path("."), 5, 1, {"mitfi": ("j", "c")},
                                   Path("."), log, 16000, None, lineage, True, True, [])
        return len(calls), [f.name for f in out]
    assert run(full) == (0, list(mito_trna.TRNA_NAMES))
    assert run(full, "Metazoa;Chordata;Craniata;Vertebrata")[0] == 1
    assert run(full[:-1])[0] == 1
    weak = full[:-1] + [CoreFeature("tRNA", "trnY", 5000, 5065, 1, "MitoM-ML", score=0.6)]
    assert run(weak) == (1, list(mito_trna.TRNA_NAMES[:-1]))
    weak_copy = full + [CoreFeature("tRNA", "trnY", 6000, 6065, 1, "MitoM-ML", score=0.6)]
    assert run(weak_copy) == (0, list(mito_trna.TRNA_NAMES))      # a weaker copy below 0.9 is dropped anyway
    monkeypatch.setattr(mito_trna, "ALWAYS_MITFI", True)
    assert run(full)[0] == 1


def test_call_merges_one_trna_seen_at_both_ends_of_a_circular_sequence():
    """A pan-model hit at 1..68 and a clade-model hit at L..L+67 are one tRNA: the final filter is circular-aware."""
    import numpy as np

    class Clf:
        def predict_proba(self, X):
            return np.array([[0.1, 0.9]] * len(X))

    L = 1000
    base = {k: 0 for k in trna_ml.FEATS}
    a = dict(base, locus=0, family="M", strand=1, start=1, end=68, p=0)
    b = dict(base, locus=1, family="M", strand=1, start=L, end=L + 67, p=0)
    assert len(trna_ml.call([dict(a), dict(b)], Clf(), 0.5, L)) == 1
    assert len(trna_ml.call([dict(a), dict(b)], Clf(), 0.5, 0)) == 2   # linear: two places, as before
