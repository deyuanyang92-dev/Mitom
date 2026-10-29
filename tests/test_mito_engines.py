"""mito_engines.py and its use in mitom_pipeline.py: MITOS2 results as the core's record, existing results found under a
folder (mitom import), MitoZ / MITOS2 run as engines (fake programs), the import workflow end to end, the old name
batch_mitoz.py, mitom setup (user 2026-10-28: MitoZ / MITOS2 annotation optional, results importable, modern install).

Run: python3 -m pytest -q tests/test_mito_engines.py
"""
import hashlib
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_engines as me  # noqa: E402
import mito_env  # noqa: E402
import mitom_pipeline as mp  # noqa: E402

SEQ = "ACGT" * 1000                                            # 4,000 nt
BED = ("X\t0\t64\ttrnM(cat)\t1e-12\t+\n"
       "X\t100\t1600\tcox1\t1000\t+\n"
       "X\t1600\t1700\ttrnL1(tag)\t1e-9\t-\n"
       "X\t1700\t2500\trrnS\t1e-30\t-\n"
       "X\t2500\t2600\tOH\t1e-3\t+\n"
       "X\t2700\t2900\tIntron_gpII\t1e-5\t+\n"
       "X\t3900\t60\tnad6\t500\t+\n")                            # across the origin


def mitos_dir(tmp_path, circular=True):
    d = tmp_path / "mitos"
    d.mkdir()
    (d / "result.bed").write_text(BED)
    (d / "result.gff").write_text(f"X\tmitos\tregion\t1\t4000\t.\t+\t.\tID=X:1..4000;Is_circular={circular}\n")
    return d


def test_mitos2_record_has_the_core_names_and_qualifiers(tmp_path):
    rec = me.mitos2_record(mitos_dir(tmp_path), SEQ, "S1", 5, "circular")
    f = {x.qualifiers["gene"][0]: x for x in rec.features if x.type in ("CDS", "rRNA", "tRNA")}
    assert set(f) == {"trnM", "COX1", "trnL1", "s-rRNA", "ND6"}                 # OH and the intron left out
    assert f["COX1"].qualifiers["product"] == ["cytochrome c oxidase subunit I"]
    assert (f["COX1"].qualifiers["codon_start"], f["COX1"].qualifiers["transl_table"]) == (["1"], ["5"])
    assert f["trnL1"].qualifiers["product"] == ["tRNA-Leu"] and "anticodon: UAG" in f["trnL1"].qualifiers["note"]
    assert f["trnL1"].location.strand == -1 and f["s-rRNA"].qualifiers["product"] == ["12S ribosomal RNA"]
    assert len(f["ND6"].location.parts) == 2                                       # across the origin: a join
    assert sum(x.type == "gene" for x in rec.features) == 5 and rec.annotations["topology"] == "circular"


def mitoz_record(rid="S1"):
    r = SeqRecord(Seq(SEQ), id=rid, name=rid, description="topology=linear.")    # as MitoZ writes every record
    r.annotations.update(molecule_type="DNA", topology="linear")
    r.features = [SeqFeature(FeatureLocation(100, 1600, 1), type="CDS", qualifiers={"gene": ["COX1"]})]
    return r


def test_import_takes_mitoz_result_folders_not_its_drafts(tmp_path):
    w = tmp_path / "run" / "w"
    (w / "S1.in.fa.result").mkdir(parents=True)
    (w / "tmp_S1").mkdir()
    SeqIO.write(mitoz_record(), str(w / "S1.in.fa.result" / "S1.gbf"), "genbank")
    SeqIO.write(mitoz_record(), str(w / "tmp_S1" / "draft.gbf"), "genbank")
    gbs, mitos = me.find_results(tmp_path / "run")
    assert [p.name for p in gbs] == ["S1.gbf"] and mitos == []
    recs, notes = me.import_results(tmp_path / "run")
    assert recs[0].annotations["topology"] == "circular" and recs[0].description == "" and "MitoZ" in notes[0]
    recs, _ = me.import_results(tmp_path / "run", topology="linear")
    assert recs[0].annotations["topology"] == "linear"


def test_import_mitos2_needs_the_genetic_code_and_the_sequence(tmp_path):
    d = mitos_dir(tmp_path, circular=False)
    with pytest.raises(ValueError, match="genetic_code"):
        me.import_results(d)
    with pytest.raises(ValueError, match="--fasta"):
        me.import_results(d, code=5)
    (tmp_path / "genome.fa").write_text(f">X\n{SEQ}\n")                           # in the folder above the result
    recs, notes = me.import_results(d, code=5)
    assert [r.id for r in recs] == ["X"] and recs[0].annotations["topology"] == "linear" and "MITOS2" in notes[0]


def test_engine_cmd_per_engine():
    def cmd(*extra):
        a = mp.build_parser().parse_args(["run", "-i", "x", "-o", "y", "--clade", "Mollusca", "--genetic_code", "5",
                                          *extra])
        return mp.engine_cmd(Path("s.fa"), Path("w"), "S00001", a)
    native = cmd()
    assert Path(native[3]).name == "mito_core.py" and native[4] == "annotate"
    mz = cmd("--engine", "mitoz", "--mitoz_bin", "/opt/mz/bin")
    assert Path(mz[3]).name == "mito_engines.py" and mz[4] == "mitoz"
    assert mz[mz.index("--clade") + 1] == "Mollusca" and mz[mz.index("--mitoz_bin") + 1] == "/opt/mz/bin"
    m2 = cmd("--engine", "mitos2", "--mitos_bin", "/opt/m2/bin")
    assert m2[4] == "mitos2" and m2[m2.index("--mitos_bin") + 1] == "/opt/m2/bin" and "--clade" not in m2


FAKE_MITOZ = """#!/bin/sh
# fake MitoZ: the prepared record where MitoZ writes its result
while [ $# -gt 0 ]; do case "$1" in --workdir) wd="$2"; shift;; --outprefix) p="$2"; shift;; esac; shift; done
mkdir -p "$wd/$p.in.fa.result" "$wd/tmp_$p" && cp "$(dirname "$0")/prepared.gbf" "$wd/$p.in.fa.result/$p.gbf"
"""


def test_run_mitoz_with_a_fake_program(tmp_path):
    b = tmp_path / "bin"
    b.mkdir()
    (b / "mitoz").write_text(FAKE_MITOZ)
    (b / "mitoz").chmod(0o755)
    SeqIO.write(mitoz_record(), str(b / "prepared.gbf"), "genbank")
    fa = tmp_path / "S1.fasta"
    fa.write_text(f">S1 topology=circular\n{SEQ}\n")
    recs = me.run_mitoz(fa, tmp_path / "w", "S1", "Mollusca", 5, bin_dir=str(b))
    assert [r.id for r in recs] == ["S1"] and recs[0].annotations["topology"] == "circular"
    assert recs[0].description == ""
    with pytest.raises(RuntimeError, match="MitoZ not found"):
        me.run_mitoz(fa, tmp_path / "w2", "S1", "Mollusca", 5, bin_dir=str(tmp_path / "none"))


def test_resolve_tools_checks_the_engine(tmp_path, monkeypatch):
    def parse(*extra):
        return mp.build_parser().parse_args(["run", "-i", str(tmp_path), "-o", str(tmp_path / "o"), *extra])
    a = parse("--clade", "Cnidaria", "--genetic_code", "4", "--engine", "mitoz")
    assert "MitoZ has no protein database" in mp.resolve_tools(a)
    monkeypatch.setattr(mito_env, "mitos_bin", lambda explicit=None: None)
    a = parse("--clade", "Mollusca", "--engine", "mitos2")
    assert "--engine mitos2" in mp.resolve_tools(a)


def test_import_workflow_end_to_end(tmp_path):
    """MitoZ's annotation of Aphrodita KY753833 -> reoriented (cox1 first), curated, NCBI names; round 2 off."""
    out = tmp_path / "out"
    r = subprocess.run([sys.executable, "-W", "ignore", str(HERE.parent / "mitom_pipeline.py"), "import",
                        "-i", str(HERE / "data" / "KY753833_mitoz.gbf"), "-o", str(out), "--clade",
                        "Annelida-segmented-worms", "--genetic_code", "5", "--plot", "no", "--rrna_hmm", "no",
                        "--intron_fix", "no", "--threads", "2"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-2000:]
    rec = SeqIO.read(str(out / "final_annotation" / "KY753833.gb"), "genbank")
    first = next(f for f in rec.features if f.type == "CDS")
    assert first.qualifiers["gene"] == ["COX1"] and int(first.location.start) == 0
    assert rec.annotations["topology"] == "circular"                     # MitoZ wrote linear
    assert "round 2" not in (out / "summary.tsv").read_text()
    assert (out / "original_vs_mitom.tsv").is_file() and not (out / "00.imported").exists()


def test_old_name_batch_mitoz_still_works():
    import batch_mitoz
    assert batch_mitoz is mp and batch_mitoz.workflow_run is mp.workflow_run
    r = subprocess.run([sys.executable, str(HERE.parent / "batch_mitoz.py"), "import", "-h"], capture_output=True, text=True)
    assert r.returncode == 0 and "--topology" in r.stdout


def test_setup_downloads_checks_and_writes_the_config(tmp_path, monkeypatch):
    (tmp_path / "src" / "refseq89m").mkdir(parents=True)
    (tmp_path / "src" / "refseq89m" / "x.txt").write_text("x")
    tgz = tmp_path / "refseq89m.tar.bz2"
    with tarfile.open(tgz, "w:bz2") as tf:
        tf.add(tmp_path / "src" / "refseq89m", arcname="refseq89m")
    monkeypatch.setattr(mito_env.urllib.request, "urlretrieve", lambda url, dest: shutil.copy(tgz, dest))
    monkeypatch.setattr(mito_env, "REF_MD5", hashlib.md5(tgz.read_bytes()).hexdigest())
    monkeypatch.setenv("YZZ_MITO_CONFIG", str(tmp_path / "cfg.json"))
    monkeypatch.delenv("MITOS_REFDIR", raising=False)
    mito_env.config.cache_clear()
    try:
        mito_env.setup(str(tmp_path / "refs"))
        assert (tmp_path / "refs" / "refseq89m" / "x.txt").is_file() and not (tmp_path / "refs" / tgz.name).exists()
        assert str(tmp_path / "refs") in (tmp_path / "cfg.json").read_text()
        monkeypatch.setattr(mito_env, "REF_MD5", "0" * 32)                # present: nothing downloaded again
        mito_env.setup(str(tmp_path / "refs"))
        monkeypatch.setattr(mito_env.urllib.request, "urlretrieve", lambda url, dest: Path(dest).write_text("bad"))
        assert mito_env.setup(str(tmp_path / "other")) == 1               # checksum mismatch
    finally:
        mito_env.config.cache_clear()


def test_one_mitos2_command_for_engine_rescue_and_trna():
    """mito_engines.mitos2_command gives the command lines the engine, mito_rescue and mito_trna used to build
    themselves (they ran it as `env PATH=<bin>:/usr/bin:/bin ...`; the environment is now engine_env(bin, isolated))."""
    b, d, fa, out = "/m/bin", "/refs", Path("/w/in.fa"), Path("/w/res")
    assert me.mitos2_command(fa, out, 5, b, d, linear=True) == [
        "/m/bin/runmitos", "-i", "/w/in.fa", "-c", "5", "-o", "/w/res", "-r", "refseq89m", "--refdir", "/refs",
        "--noplots", "--linear"]
    fast = str(me.RUNMITOS_FAST)
    assert fast.endswith("scripts/runmitos_fast.py")
    assert me.mitos2_command(fa, out, 5, b, d, mode="sensitive", skip=("trna", "intron", "oril", "orih")) == [
        "/m/bin/python", fast, "-i", "/w/in.fa", "-c", "5", "-o", "/w/res", "-r", "refseq89m", "--refdir", "/refs",
        "--trna", "0", "--intron", "0", "--oril", "0", "--orih", "0", "--noplots"]
    assert me.mitos2_command(fa, out, 9, b, d, ref="r1", mode="fast", skip=("prot", "rrna", "intron", "oril"),
                             noplots=False) == [
        "/m/bin/python", fast, "--no_sensitive", "-i", "/w/in.fa", "-c", "9", "-o", "/w/res", "-r", "r1", "--refdir",
        "/refs", "--prot", "0", "--rrna", "0", "--intron", "0", "--oril", "0"]
    assert me.engine_env(b, isolated=True)["PATH"] == "/m/bin:/usr/bin:/bin"


def test_result_bed_reader(tmp_path):
    (tmp_path / "result.bed").write_text("s1\t10\t80\ttrnL1(tag)\t1e-9\t+\nbad line\ns1\t15000\t20\tcox1\t0\t-\n")
    rows = me.read_result_bed(tmp_path)
    assert rows == [me.BedRow("s1", 10, 80, "trnL1(tag)", "1e-9", 1), me.BedRow("s1", 15000, 20, "cox1", "0", -1)]
    assert me.mitos2_seq_ids(tmp_path) == ["s1"]
