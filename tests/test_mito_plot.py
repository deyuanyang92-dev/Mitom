"""Tests for mito_plot.py (circular genome map: GC / skew windows, labels, read depth with minimap2 or bwa mem, CLI).

Run: python3 -m pytest -q tests/test_mito_plot.py   (the depth tests need minimap2 / bwa)
"""
import gzip
import random
import sys
from pathlib import Path

import numpy as np
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_env  # noqa: E402
import mito_plot as mp  # noqa: E402

rnd = random.Random(7)


def rand(n):
    return "".join(rnd.choice("ACGT") for _ in range(n))


def feat(kind, loc, gene=None, product=None, note=None):
    q = {k: [v] for k, v in (("gene", gene), ("product", product), ("note", note)) if v}
    return SeqFeature(loc, type=kind, qualifiers=q)


def genome(n=6000, topology="circular"):
    """Both strands, a crowded tRNA cluster, CYTB across the origin (circular), ND1 with an intron, a D-loop."""
    r = SeqRecord(Seq(rand(n)), id="TEST1.1", description="")
    r.annotations.update(molecule_type="DNA", topology=topology, organism="Testus syntheticus")
    fs = [feat("source", FeatureLocation(0, n, 1)),
          feat("gene", FeatureLocation(100, 1600, 1), "COX1"),
          feat("CDS", FeatureLocation(100, 1600, 1), "COX1", "cytochrome c oxidase subunit I"),
          feat("tRNA", FeatureLocation(1600, 1668, -1), "trnS2", "tRNA-Ser"),
          feat("tRNA", FeatureLocation(1668, 1736, -1), "trnD", "tRNA-Asp"),
          feat("tRNA", FeatureLocation(1736, 1804, -1), "trnK", "tRNA-Lys"),
          feat("tRNA", FeatureLocation(1804, 1872, -1), product="tRNA-Phe"),
          feat("CDS", FeatureLocation(1880, 2040, 1), "ATP8"),
          feat("CDS", FeatureLocation(2030, 2700, 1), "ATP6"),
          feat("rRNA", FeatureLocation(2750, 3600, -1), "rrnL", "16S ribosomal RNA"),
          feat("CDS", CompoundLocation([FeatureLocation(3650, 3900, 1), FeatureLocation(4300, 4700, 1)]), "ND1"),
          feat("intron", FeatureLocation(3900, 4300, 1), "ND1"),
          feat("D-loop", FeatureLocation(4800, 5500, None)),
          feat("tRNA", FeatureLocation(5500, 5568, 1), "trnL1", "tRNA-Leu")]
    if topology == "circular":
        fs.append(feat("CDS", CompoundLocation([FeatureLocation(5600, n, 1), FeatureLocation(0, 60, 1)]), "CYTB"))
    r.features = fs
    return r


def test_category_and_label():
    def f(kind, gene=None, product=None, note=None):
        return feat(kind, FeatureLocation(0, 10, 1), gene, product, note)
    assert [mp.category(f("CDS", g)) for g in ("COX1", "cox3", "COI", "ND4L", "nad5", "ATP8", "CYTB", "cob")] == \
        ["COX", "COX", "COX", "NADH", "NADH", "ATP", "CYTB", "CYTB"]
    assert mp.category(f("CDS", "orf314")) == "other" and mp.category(f("D-loop")) == "other"
    assert mp.category(f("gene", "COX1")) is None and mp.category(f("intron")) is None
    assert [mp.label(f("tRNA", g, p), "tRNA") for g, p in
            (("trnL1", "tRNA-Leu"), ("trnS2", None), ("TRNA", None), ("tRNA-Leu", None), (None, "tRNA-Phe"))] == \
        ["L1", "S2", "A", "L", "F"]
    assert mp.label(f("misc_feature", note="putative control region"), "other") == "CR"
    assert mp.label(f("rRNA", "rrnS"), "rRNA") == "rrnS"


def test_gc_windows_values():
    seq = "G" * 100 + "A" * 100                                # 200 nt, windows of 100 every 20
    x, gc, skew = mp.gc_windows(seq, circular=True)
    assert len(x) == len(gc) == len(skew) == 10                 # ceil(200 / 20), wrapping across the origin
    assert x[0] == 50 and gc[0] == 1.0 and skew[0] == 1.0
    assert gc[5] == 0.0 and skew[5] == 0.0                      # all A: no G / C
    assert gc[8] == pytest.approx(0.6)                          # 160..259 = 40 A + 60 G across the origin
    x, gc, skew = mp.gc_windows(seq, circular=False)
    assert len(gc) == 6 and gc[-1] == 0.0                      # full windows only: starts 0..100
    _, gc, skew = mp.gc_windows("G" * 30 + "C" * 70 + "N" * 50 + "A" * 50, circular=False, window=100, step=50)
    assert gc[0] == 1.0 and skew[0] == pytest.approx(-0.4)
    assert gc[1] == 1.0 and gc[2] == 0.0                        # 50..149: 50 C + 50 N, N not counted
    assert len(mp.gc_windows("ACGT" * 10, circular=False)[1]) == 1   # shorter than one window


def test_span_and_fan():
    assert mp.span([(5600, 6000), (0, 60)], 6000, True) == (5600, 460)      # across the origin
    assert mp.span([(3650, 3900), (4300, 4700)], 6000, True) == (3650, 1050)  # exons around an intron
    assert mp.span([(10, 20)], 6000, False) == (10, 10)
    placed = mp.fan([0.10, 0.11, 0.12, 2.0], [0.05] * 4)
    assert placed[3] == pytest.approx(2.0)                      # alone: not moved
    assert np.diff(placed[:3]) == pytest.approx([0.05, 0.05])  # crowded: no overlap, order kept ...
    assert np.mean(placed[:3]) == pytest.approx(0.11)          # ... centred on the targets
    placed = mp.fan([6.27, 0.01], [0.05, 0.05])                # neighbours across angle 0
    assert (placed[1] - placed[0]) % (2 * np.pi) == pytest.approx(0.05)


@pytest.mark.parametrize("topology", ["circular", "linear"])
def test_plot_record_writes_png_and_pdf(tmp_path, topology):
    rec = genome(topology=topology)
    depth = np.linspace(10, 50, len(rec.seq)) if topology == "linear" else None
    paths = mp.plot_record(rec, tmp_path, depth=depth)
    assert [p.name for p in paths] == ["TEST1.1.png", "TEST1.1.pdf"]
    assert paths[0].read_bytes()[:4] == b"\x89PNG" and paths[1].read_bytes()[:4] == b"%PDF"
    assert all(p.stat().st_size > 10000 for p in paths)
    with pytest.raises(ValueError):
        mp.plot_record(rec, tmp_path, depth=np.ones(10))


def write_fastq(path, reads):
    with gzip.open(path, "wt") as fh:
        for i, s in enumerate(reads):
            fh.write(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n")


@pytest.mark.parametrize("aligner", ["minimap2", "bwa"])
def test_depth_from_reads_folds_the_origin(tmp_path, monkeypatch, aligner):
    """Both aligners give the exact depth (minimap2 when installed, else bwa mem)."""
    if mito_env.find(aligner) is None:
        pytest.skip(f"{aligner} not found")
    find = mito_env.find
    monkeypatch.setattr(mp.mito_env, "find", lambda name, *a: None if name in ("minimap2", "bwa") and name != aligner
                        else find(name, *a))
    n = 3000
    rec = SeqRecord(Seq(rand(n)), id="T2", annotations={"topology": "circular"})
    s = str(rec.seq)
    frags = [s[i:i + 300] for i in range(1000, 1500, 50)] + [s[-150:] + s[:150]] * 3   # 10 + 3 across the origin
    write_fastq(tmp_path / "r1.fq.gz", [f[:150] for f in frags])
    write_fastq(tmp_path / "r2.fq.gz", [str(Seq(f[-150:]).reverse_complement()) for f in frags])
    depth = mp.depth_from_reads(rec, [tmp_path / "r1.fq.gz", tmp_path / "r2.fq.gz"], threads=2)
    expected = np.zeros(n, dtype=int)
    for i in range(1000, 1500, 50):
        expected[i:i + 150] += 1                                # r1
        expected[i + 150:i + 300] += 1                          # r2
    expected[n - 150:] += 3                                     # origin pairs: r1 before the end ...
    expected[:150] += 3                                         # ... r2 after it, folded onto the start
    assert depth.shape == (n,) and (depth == expected).all()
    rec.annotations["topology"] = "linear"                      # no folding: the origin pairs map to the ends only
    lin = mp.depth_from_reads(rec, [tmp_path / "r1.fq.gz"])
    assert lin[n - 150:].min() == 3 and lin[:1000].sum() == 0


def test_cli_writes_files(tmp_path, capsys):
    rec = genome()
    SeqIO.write([rec], tmp_path / "x.gb", "genbank")
    other = genome(topology="linear")
    other.id = other.name = "TEST2.1"
    SeqIO.write([other], tmp_path / "y.gbf", "genbank")
    assert mp.main(["-i", str(tmp_path / "x.gb"), "-o", str(tmp_path / "one"), "--formats", "pdf"]) == 0
    assert (tmp_path / "one" / "TEST1.1.pdf").stat().st_size > 0
    assert mp.main(["-i", str(tmp_path), "-o", str(tmp_path / "all"), "--formats", "png"]) == 0
    assert sorted(p.name for p in (tmp_path / "all").iterdir()) == ["TEST1.1.png", "TEST2.1.png"]
    with pytest.raises(SystemExit):                            # reads need a single record
        mp.main(["-i", str(tmp_path), "-o", str(tmp_path / "all"), "--fq1", str(tmp_path / "r.fq")])
    assert "one record" in capsys.readouterr().err


# ── batch_mitoz.run_plots: the genome-map step of `run` / `annotate` ─────────
def final_dir(tmp_path):
    d = tmp_path / "final"
    d.mkdir()
    a, b = genome(), genome(topology="linear")
    b.id = b.name = "TEST2.1"
    SeqIO.write([a], d / "S1.gbf", "genbank")
    SeqIO.write([b], d / "S2.gbf", "genbank")
    return d, a


def test_run_plots_one_map_per_record_and_switch(tmp_path):
    import argparse
    import batch_mitoz as bm
    d, _ = final_dir(tmp_path)
    args = argparse.Namespace(plot="no", infq=None, threads=1, max_tasks=2)
    assert bm.run_plots(args, d, tmp_path) == "skipped (--plot no)" and not (tmp_path / "plots").exists()
    args.plot = "yes"
    assert bm.run_plots(args, d, tmp_path) == "2 map(s) in plots/"           # worker processes
    assert sorted(p.name for p in (tmp_path / "plots").iterdir()) == \
        ["TEST1.1.pdf", "TEST1.1.png", "TEST2.1.pdf", "TEST2.1.png"]


@pytest.mark.skipif(not (mito_env.find("minimap2") or mito_env.find("bwa")), reason="minimap2 / bwa not found")
def test_run_plots_depth_ring_from_the_samples_reads(tmp_path):
    """--infq: reads found by sample name (here S1 only, single file) add the depth ring to that sample's map."""
    import argparse
    import batch_mitoz as bm
    d, a = final_dir(tmp_path)
    (tmp_path / "fq").mkdir()
    s = str(a.seq)
    write_fastq(tmp_path / "fq" / "S1_1.fq.gz", [s[i:i + 150] for i in range(0, len(s) - 150, 50)])
    args = argparse.Namespace(plot="yes", infq=str(tmp_path / "fq"), fq_position="flat",
                              suffix_fq="_1.fq.gz,_2.fq.gz", threads=1, max_tasks=2)
    assert bm.run_plots(args, d, tmp_path / "o") == "2 map(s) in plots/, 1 with read depth"
