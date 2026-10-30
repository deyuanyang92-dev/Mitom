"""Tests for mitom_pipeline.py deduplication (RefSeq NC_ kept over the INSDC submission it was made from).

Run: python3 -m pytest -q tests/test_dedup.py
"""
import copy
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mitom_pipeline as bm  # noqa: E402

DATA = HERE / "data"


def insdc_copy(refseq, acc: str, comment: str = ""):
    """The INSDC submission a RefSeq record was derived from (same sequence, own accession)."""
    r = copy.deepcopy(refseq)
    r.id, r.name = f"{acc}.1", acc
    r.annotations["accessions"] = [acc]
    r.annotations["comment"] = comment
    return r


def test_refseq_source_parsed():
    rec = SeqIO.read(str(DATA / "NC_000933.gb"), "genbank")
    assert bm.is_refseq_acc("NC_000933") and not bm.is_refseq_acc("AF000023")
    assert bm.refseq_source_accessions(rec) == ["AF000023"]      # "derived from AF000023"


def test_refseq_kept_source_dropped_both_orders(tmp_path):
    nc = SeqIO.read(str(DATA / "NC_000933.gb"), "genbank")
    src = insdc_copy(nc, "AF000023")
    other = insdc_copy(nc, "ZZ999999")
    other.seq = Seq(str(nc.seq)[::-1])                            # unrelated, different sequence
    for order in ([src, nc, other], [nc, other, src]):
        f = tmp_path / "in.gb"
        SeqIO.write(order, str(f), "genbank")
        kept, _, stats = bm.dedup_gb_records([f], dedup_report_path=tmp_path / "rep.tsv")
        ids = sorted(m.acc_versioned for m in kept)
        assert ids == ["NC_000933.1", "ZZ999999.1"], ids
        assert stats["dup_refseq"] == 1
        rows = {r.split("\t")[0]: r.split("\t") for r in (tmp_path / "rep.tsv").read_text().splitlines()
                if r and not r.startswith("#")}
        assert rows["AF000023.1"][1:4] == ["dropped", "INSDC submission that RefSeq record was derived from (its COMMENT)", "NC_000933.1"]
        assert rows["NC_000933.1"][1:4] == ["kept", "", "AF000023.1"]
        assert rows["ZZ999999.1"][1:4] == ["kept", "", ""]


def test_identical_sequence_prefers_refseq_without_comment(tmp_path):
    nc = SeqIO.read(str(DATA / "NC_000933.gb"), "genbank")
    nc.annotations["comment"] = ""                                 # no "derived from" link
    src = insdc_copy(nc, "AB123456")
    f = tmp_path / "in.gb"
    SeqIO.write([src, nc], str(f), "genbank")
    kept, _, stats = bm.dedup_gb_records([f])
    assert [m.acc_versioned for m in kept] == ["NC_000933.1"]
    assert stats["dup_seq"] == 1


def test_dedup_subcommand_writes_gb_and_fasta(tmp_path):
    nc = SeqIO.read(str(DATA / "NC_000933.gb"), "genbank")
    f = tmp_path / "in.gb"
    SeqIO.write([insdc_copy(nc, "AF000023"), nc], str(f), "genbank")
    subprocess.run([sys.executable, "-W", "ignore", str(HERE.parent / "mitom_pipeline.py"), "dedup",
                    "-i", str(f), "-o", str(tmp_path / "out")], check=True, capture_output=True)
    assert [r.id for r in SeqIO.parse(str(tmp_path / "out" / "dedup.gb"), "genbank")] == ["NC_000933.1"]
    assert [r.id for r in SeqIO.parse(str(tmp_path / "out" / "dedup.fasta"), "fasta")] == ["NC_000933.1"]


def _mitom(*args) -> int:
    """mitom_pipeline.py with these arguments, as the `mitom` command runs it; its exit code."""
    return subprocess.run([sys.executable, str(HERE.parent / "mitom_pipeline.py"), *map(str, args)],
                          capture_output=True, text=True).returncode


def _fa(path, records):
    """A FASTA file from (header, sequence) pairs."""
    path.write_text("".join(f">{h}\n{s}\n" for h, s in records))
    return path


def test_fasta_dedup_versions_identical_sequences_and_local_names(tmp_path):
    """2026-10-30: `mitom dedup` takes FASTA too - versions of one accession and identical sequences go, local
    names (NODE_1 of two assemblies) are never merged on their name."""
    a, b, c = "ACGT" * 50, "TTGCA" * 40, "GGCCA" * 40
    f1 = _fa(tmp_path / "x.fasta", [("MK123456.1 old", a), ("NODE_1 assembly x", b), ("sample_7", c)])
    f2 = _fa(tmp_path / "y.fa", [("MK123456.2 new", a + "A"), ("NODE_1 assembly y", c), ("OR555555.1 copy", b)])
    kept, _, stats = bm.dedup_gb_records([f1, f2], dedup_report_path=tmp_path / "rep.tsv")
    ids = sorted(m.acc_versioned for m in kept)
    # MK123456.1 -> .2; NODE_1 of x = OR555555.1 (first kept); sample_7 = NODE_1 of y (first kept); the two NODE_1
    # differ and are not merged on their name
    assert ids == ["MK123456.2", "NODE_1", "sample_7"], ids
    assert stats["dup_version"] + stats["dup_exact"] == 1 and stats["dup_seq"] == 2
    rows = [r.split("\t") for r in (tmp_path / "rep.tsv").read_text().splitlines() if r and not r.startswith("#")][1:]
    assert sorted(r[2] for r in rows if r[1] == "dropped") == ["identical sequence", "identical sequence", "older version"]


def test_dedup_command_fasta_and_genbank_outputs(tmp_path):
    nc = SeqIO.read(str(DATA / "NC_000933.gb"), "genbank")
    gb_dir = tmp_path / "in"
    gb_dir.mkdir()
    SeqIO.write([nc], str(gb_dir / "a.gb"), "genbank")
    _fa(gb_dir / "b.fasta", [("my_sample topology=circular", str(nc.seq)), ("other", "ACGTT" * 30)])
    assert _mitom("dedup", "-i", gb_dir, "-o", tmp_path / "out") == 0
    assert [r.id for r in SeqIO.parse(str(tmp_path / "out" / "dedup.gb"), "genbank")] == ["NC_000933.1"]
    heads = [line[1:].strip() for line in (tmp_path / "out" / "dedup.fasta").read_text().splitlines() if line.startswith(">")]
    assert heads[0].startswith("NC_000933.1 ") and heads[0].endswith("topology=circular")   # GenBank: tag added
    assert heads[1:] == ["other"]                         # FASTA: unchanged header; my_sample = NC_000933.1, dropped
    fa_only = tmp_path / "fa_only"
    fa_only.mkdir()
    _fa(fa_only / "c.fasta", [("s1", "ACGT" * 30), ("s2", "ACGT" * 30)])
    assert _mitom("dedup", "-i", fa_only, "-o", tmp_path / "out2") == 0
    assert not (tmp_path / "out2" / "dedup.gb").exists()
    assert (tmp_path / "out2" / "dedup.fasta").read_text().count(">") == 1


def test_gb2fasta_headers_split_and_mitom_run_topology(tmp_path):
    """2026-10-30: `mitom gb2fasta` - NCBI GenBank and MitoZ .gbf to FASTA with topology= in the header, so that
    `mitom run` keeps the genomes of a multi-record FASTA circular."""
    nc = SeqIO.read(str(DATA / "NC_000933.gb"), "genbank")
    other = copy.deepcopy(nc)
    other.id, other.name, other.annotations["topology"] = "XX000001.1", "XX000001", "linear"
    src = tmp_path / "gb"
    src.mkdir()
    SeqIO.write([nc, other], str(src / "two.gb"), "genbank")
    (src / "mitoz.gbf").write_text((DATA / "KY753833_mitoz.gbf").read_text())
    assert _mitom("gb2fasta", "-i", src, "-o", tmp_path / "out") == 0
    recs = list(SeqIO.parse(str(tmp_path / "out" / "two.fasta"), "fasta"))
    assert [r.id for r in recs] == ["NC_000933.1", "XX000001.1"]
    assert recs[0].description.endswith("topology=circular") and recs[1].description.endswith("topology=linear")
    assert str(recs[0].seq) == str(nc.seq)
    assert bm.infer_topology_from_description(recs[0].description) == "circular"
    assert (tmp_path / "out" / "mitoz.fasta").exists()
    assert _mitom("gb2fasta", "-i", src / "two.gb", "-o", tmp_path / "split", "--split") == 0
    assert sorted(p.name for p in (tmp_path / "split").iterdir()) == ["NC_000933.1.fasta", "XX000001.1.fasta"]
    assert _mitom("gb2fasta", "-i", tmp_path / "nothing", "-o", tmp_path / "x") == 2
