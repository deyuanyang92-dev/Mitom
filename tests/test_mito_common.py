"""mito_common.py: the helpers that replaced the modules' own copies give the old results."""
import os
import re
import stat
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_common as mc  # noqa: E402
import mito_engines  # noqa: E402


def test_iter_gb_files(tmp_path):
    (tmp_path / "a" / "b").mkdir(parents=True)
    for name in ("a/x.gb", "a/b/y.GBK", "a/b/z.gbff", "a/note.txt", "w.genbank"):
        (tmp_path / name).write_text("")
    (tmp_path / "dir.gb").mkdir()                               # a folder with a GenBank suffix is not a file
    assert [p.relative_to(tmp_path).as_posix() for p in mc.iter_gb_files(tmp_path)] == \
        ["a/b/y.GBK", "a/b/z.gbff", "a/x.gb", "w.genbank"]
    assert mc.iter_gb_files(tmp_path / "a" / "note.txt") == [tmp_path / "a" / "note.txt"]      # a file as given
    assert mc.iter_gb_files(tmp_path / "missing") == []


def rec(seq="ATGAAATAG" * 20, rid="R1"):
    r = SeqRecord(Seq(seq), id=rid, name=rid, description="")
    r.annotations.update(molecule_type="DNA", topology="circular")
    r.features = [SeqFeature(FeatureLocation(0, 9, 1), type="CDS", qualifiers={"gene": ["ND1"]})]
    return r


def test_write_gb_is_atomic_and_keeps_the_mode(tmp_path, monkeypatch):
    f = tmp_path / "x.gbf"
    mc.write_gb([rec()], f)
    assert SeqIO.read(str(f), "genbank").id == "R1"
    umask = os.umask(0)
    os.umask(umask)
    assert stat.S_IMODE(f.stat().st_mode) == 0o666 & ~umask                 # a new file: as SeqIO.write makes it
    f.chmod(0o640)
    before = f.read_text()

    def broken(*a, **k):
        raise RuntimeError("disk full")
    monkeypatch.setattr(SeqIO, "write", broken)
    with pytest.raises(RuntimeError):
        mc.write_gb([rec(rid="R2")], f)
    assert f.read_text() == before and [p.name for p in tmp_path.iterdir()] == ["x.gbf"]   # old file, no temp left
    monkeypatch.undo()
    mc.write_gb([rec(rid="R2")], f)
    assert SeqIO.read(str(f), "genbank").id == "R2" and stat.S_IMODE(f.stat().st_mode) == 0o640


def test_read_fasta_matches_seqio_and_the_old_engines_reader(tmp_path):
    fa = tmp_path / "in.fa"
    fa.write_text(">c1 topology=linear sample one\nACGT acgt\nNNNN\n>c2\nGGGG\n\n>c3  topology = Linear\nTTT\n")
    new = [(r.id, r.description, str(r.seq)) for r in mc.read_fasta(fa)]
    assert new == [(r.id, r.description, str(r.seq)) for r in SeqIO.parse(str(fa), "fasta")]
    old = [(r.id, str(r.seq).upper(), "linear" if re.search(r"topology\s*=\s*linear", r.description, re.I)
            else "circular") for r in SeqIO.parse(str(fa), "fasta")]
    assert mito_engines.read_fasta(fa) == old == [("c1", "ACGTACGTNNNN", "linear"), ("c2", "GGGG", "circular"),
                                                  ("c3", "TTT", "linear")]


def test_contig_record_id():
    assert mc.contig_record_id("S00001", "contig 1", 0, 1) == "S00001"
    assert mc.contig_record_id("S00001", "NODE_1_length_15000+x", 0, 2) == "NODE_1_length_15"
    assert mc.contig_record_id("S00001", "", 2, 3) == "S00001_3"


def old_stop_qualifiers(rec, loc):
    """mito_intron_ref.stop_qualifiers before it moved here."""
    nt = str(loc.extract(rec.seq)).upper()
    tail = len(nt) % 3
    if tail == 0 or nt[-tail:] not in ("T", "TA"):
        return {}
    last = loc.parts[-1]
    if last.strand == -1:
        a = int(last.start) + 1
        pos = f"complement({a}..{a + tail - 1})" if tail > 1 else f"complement({a})"
    else:
        b = int(last.end)
        pos = f"{b - tail + 1}..{b}" if tail > 1 else f"{b}"
    return {"transl_except": [f"(pos:{pos},aa:TERM)"],
            "note": ["TAA stop codon is completed by the addition of 3' A residues to the mRNA"]}


@pytest.mark.parametrize("loc", [FeatureLocation(10, 41, 1), FeatureLocation(10, 42, 1), FeatureLocation(10, 43, 1),
                                 FeatureLocation(10, 41, -1), FeatureLocation(10, 42, -1),
                                 CompoundLocation([FeatureLocation(10, 30, 1), FeatureLocation(50, 62, 1)]),
                                 CompoundLocation([FeatureLocation(50, 62, -1), FeatureLocation(10, 30, -1)])])
def test_stop_qualifiers_as_before(loc):
    r = SeqRecord(Seq("".join("ACGT"[(i * 7 + i // 5) % 4] for i in range(120))), id="X")
    assert mc.stop_qualifiers(r, loc, 5) == old_stop_qualifiers(r, loc)


def test_transl_except_term_is_curates_formula():
    for strand in (1, -1):
        for tail in (1, 2):
            gs, ge = 101, 400
            old = (f"{ge - tail + 1}..{ge}" if tail > 1 else f"{ge}") if strand == 1 else \
                (f"complement({gs}..{gs + tail - 1})" if tail > 1 else f"complement({gs})")
            assert mc.transl_except_term(gs, ge, strand, tail) == f"(pos:{old},aa:TERM)"


def test_codon_translate():
    for nt in ("ATGAAATAGT", "ATGAAATAGTA", "ATG", "AT", ""):
        assert mc.codon_translate(nt, 5) == str(Seq(nt[: len(nt) - len(nt) % 3]).translate(table=5))
        assert mc.codon_translate(Seq(nt), 2) == str(Seq(nt)[: len(nt) - len(nt) % 3].translate(table=2))


def test_run_tool_logs_command_exit_code_time_and_stderr(tmp_path, caplog):
    import logging
    import subprocess
    with caplog.at_level(logging.DEBUG, logger="mitom.tools"):
        p = mc.run_tool([sys.executable, "-c", "import sys; print('out'); sys.exit(0)"], cwd=tmp_path)
        assert p.returncode == 0 and p.stdout == "out\n"
        q = mc.run_tool([sys.executable, "-c", "import sys; print('bad input', file=sys.stderr); sys.exit(3)"])
        assert q.returncode == 3 and "bad input" in q.stderr
        with pytest.raises(subprocess.CalledProcessError):
            mc.run_tool([sys.executable, "-c", "import sys; sys.exit(2)"], check=True)
        with pytest.raises(OSError):
            mc.run_tool([str(tmp_path / "no_such_program")])
    msgs = [r.getMessage() for r in caplog.records if r.name == "mitom.tools"]
    assert any(m.startswith("[CMD rc=0 ") and f"(in {tmp_path})" in m for m in msgs)
    assert any(m.startswith("[CMD rc=3 ") for m in msgs) and any(m == "[CMD stderr] bad input" for m in msgs)
    assert any(m.startswith("[CMD failed ") and "no_such_program" in m for m in msgs)


def test_abbreviated_stop_and_ends_with_stop():
    stops = {"TAA", "TAG"}
    cases = {"ATGAAAT": (True, True), "ATGAAATA": (True, True), "ATGAAATG": (False, False), "ATGAAATAA": (False, True),
             "ATGAAAAAA": (False, False), "": (False, False), "T": (True, True)}
    for nt, (abbr, stop) in cases.items():
        old_tail = len(nt) % 3
        assert mc.abbreviated_stop(nt) == abbr == bool(old_tail and nt[-old_tail:] in ("T", "TA"))
        assert mc.ends_with_stop(nt, stops) == stop == (nt[-3:] in stops if not old_tail else nt[-old_tail:] in ("T", "TA"))


@pytest.mark.parametrize("text,L,shift,rc,want", [
    ("(pos:100..102,aa:Leu,seq:tag)", 1000, 50, False, "(pos:50..52,aa:Leu,seq:tag)"),
    ("(pos:complement(10..12),aa:Leu,seq:tag)", 1000, 11, False, "(pos:complement(join(999..1000,1)),aa:Leu,seq:tag)"),
    ("(pos:complement(join(999..1000,1)),aa:TERM)", 1000, -5, False, "(pos:complement(4..6),aa:TERM)"),
    ("(pos:1000,aa:TERM)", 1000, 0, True, "(pos:complement(1),aa:TERM)"),
    ("(pos:join(999..1000,1),aa:TERM)", 1000, 0, True, "(pos:complement(join(1000,1..2)),aa:TERM)"),
    ("(pos:2..3,aa:TERM)", 1000, -998, False, "(pos:join(1000,1),aa:TERM)"),
    ("no position here", 1000, 7, True, "no position here"),
])
def test_move_positions(text, L, shift, rc, want):
    assert mc.move_positions(text, L, shift, rc) == want


def test_move_position_qualifiers_only_touches_positions():
    q = {"gene": ["ND5"], "transl_except": ["(pos:10,aa:TERM)"], "note": ["pos:1 is not a position qualifier"]}
    out = mc.move_position_qualifiers(q, 100, 5)
    assert out == {"gene": ["ND5"], "transl_except": ["(pos:5,aa:TERM)"], "note": q["note"]} and q["transl_except"] == [
        "(pos:10,aa:TERM)"]
