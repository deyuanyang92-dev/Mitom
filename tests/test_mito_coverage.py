"""mito_coverage: depth from SAM lines, folding of the circular extension, the samtools-depth reader, summary,
tables and the flat map."""
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import mito_coverage as mc


def sam(rname, pos, cigar, flag=0, mapq=60):
    return f"r\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\tACGT\tIIII\n"


def test_count_sam_spans_and_flags():
    lines = ["@HD\tVN:1.6\n", sam("a", 3, "4M"), sam("a", 9, "2M2D2M"), sam("a", 1, "4M", flag=4),
             sam("a", 1, "4M", flag=256), sam("a", 1, "4M", flag=2048), sam("a", 1, "2S4M"), sam("b", 1, "4M"),
             sam("a", 1, "4M", mapq=0)]
    d = mc.count_sam(lines, {"a": 12}, min_mapq=1)["a"]
    assert d.tolist() == [1, 1, 2, 2, 1, 1, 0, 0, 1, 1, 1, 1]      # 4M@3, 6M-span@9, 2S4M@1 (soft clip not counted)


def test_depth_file_and_summary(tmp_path):
    rec = SeqRecord(Seq("A" * 6), id="x")
    f = tmp_path / "d.tsv"
    f.write_text("x\t1\t5\nx\t2\t5\nx\t4\t10\ny\t1\t9\n")
    d = mc.read_depth_file(str(f), [rec])["x"]
    assert d.tolist() == [5, 5, 0, 10, 0, 0]
    s = mc.summarise(d)
    assert s["max"] == 10 and abs(s["pct_ge1"] - 100 * 3 / 6) < 1e-9


def test_tables_and_flat_map(tmp_path):
    rec = SeqRecord(Seq("A" * 50), id="x")
    d = np.r_[np.zeros(5), np.full(45, 20)].astype(int)
    paths = mc.write_tables({"x": d}, tmp_path / "p")
    assert paths[0].read_text().splitlines()[0] == "x\t1\t0"
    assert len(paths[0].read_text().splitlines()) == 50
    fig = mc.flat_figure(rec, d)
    assert mc.save_figure(fig, tmp_path / "flat", ["png"])[0].stat().st_size > 1000
    assert mc._runs(d == 0) == [(0, 5)]
