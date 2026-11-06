"""Fixes of the 2026-11-04 / 11-05 model reviews (Mitom_review/VERIFICATION_REPORT.md, NEW_BUGS_2026-11-04.md):
not validated is never 0 ERROR (R2 / R4 / N6), REJECT / Critical / FATAL count as ERROR (R3), value-less qualifiers
in the .tbl (N1 / N2 / D8), --topology linear never trimmed (R5), import after the lock (R6), unpacked-name
collisions (R8) and gzip of gzip (E19), the previous run's results (R9 / N14), tunicate code of `import --clade
auto` (D1), exonerate code of the intron step (D2), install.sh on macOS (D3 / D13 / D23), failed group runs (N3),
_concat_tsv headers (D17), lock of another user's run (N8), rescue / reorient gene names (E15 / E16), rotate
(E22), import names (E5 / E12), metadata (E13 / E14), NCR note (E3).

Run: python3 -m pytest -q tests/test_review_2026_11_05.py
"""
import gzip
import io
import os
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))
import mitom_pipeline as bm  # noqa: E402
import mito_common  # noqa: E402
import mito_engines  # noqa: E402
import mito_intron  # noqa: E402
import mito_rescue  # noqa: E402
import mito_tbl  # noqa: E402
import mitom_groups  # noqa: E402
import mitom_inputs  # noqa: E402
import mitom_metadata  # noqa: E402
import mitom_output  # noqa: E402
import mitom_steps  # noqa: E402
import ncbi_validate as nv  # noqa: E402


def _rec(rid="S1", n=300):
    r = SeqRecord(Seq("ACGT" * (n // 4)), id=rid, name=rid, description="test")
    r.annotations.update(molecule_type="DNA", topology="linear", organism="Foo bar")
    r.features = [SeqFeature(FeatureLocation(0, n, 1), type="source",
                             qualifiers={"organism": ["Foo bar"], "organelle": ["mitochondrion"]})]
    return r


def _fake_table2asn(tmp_path, body: str) -> str:
    """A table2asn stand-in: `body` is Python run with out = the -o path (…/rec.sqn)."""
    fake = tmp_path / "fake_table2asn"
    fake.write_text(f"#!{sys.executable}\nimport sys\nout = sys.argv[sys.argv.index('-o') + 1]\n{body}\n")
    fake.chmod(0o755)
    return str(fake)


# ---- .tbl: value-less qualifiers, NCR note ---------------------------------------------------------------------
@pytest.mark.parametrize("value", [None, ""])
def test_value_less_qualifier_is_key_only(value):
    r = _rec()
    r.features.append(SeqFeature(FeatureLocation(10, 100, 1), type="gene",
                                 qualifiers={"gene": ["nad1"], "pseudo": [value], "ribosomal_slippage": [value]}))
    lines = nv.feature_lines(r)
    assert "\t\t\tpseudo" in lines and "\t\t\tribosomal_slippage" in lines
    assert not any("None" in ln or ln.endswith("pseudo\t") for ln in lines)


def test_ncr_note_without_tool_name():
    r = _rec()
    r.features.append(SeqFeature(FeatureLocation(10, 100, 1), type="misc_feature", qualifiers={
        "note": ["non-coding region NCR1, 90 bp: no gene annotated (MitoM)"]}))
    assert "\t\t\tnote\tnon-coding region NCR1, 90 bp: no gene annotated" in nv.feature_lines(r)


# ---- table2asn failures and severities ---------------------------------------------------------------------------
def test_crashed_table2asn_is_not_zero_error(tmp_path):
    rows = nv.validate_record(_rec(), tmp_path / "out", _fake_table2asn(tmp_path, "sys.exit(3)"))
    assert [(r["severity"], r["code"]) for r in rows] == [("ERROR", "MITOM.Table2asnFailed")]
    rows = nv.validate_record(_rec(), tmp_path / "out", "/bin/false")
    assert rows and rows[0]["code"] == "MITOM.Table2asnFailed"
    # exit 2 with a .val: table2asn 1.28 on records with validation errors - the messages count, not the exit code
    t2a = _fake_table2asn(tmp_path, "open(out[:-4] + '.val', 'w').write('ERROR: valid [SEQ_FEAT.NoStop] x\\n'); sys.exit(2)")
    assert [r["code"] for r in nv.validate_record(_rec(), tmp_path / "out", t2a)] == ["SEQ_FEAT.NoStop"]


def test_missing_table2asn_program_is_not_zero_error(tmp_path):
    rows = nv.validate_record(_rec(), tmp_path / "out", str(tmp_path / "no_such_program"))
    assert rows[0]["code"] == "MITOM.Table2asnFailed"


def test_critical_message_counts_as_reject(tmp_path):
    t2a = _fake_table2asn(tmp_path, "open(out[:-4] + '.val', 'w').write("
                                    "'Critical: valid [SEQ_FEAT.Range] Location: SeqLoc [lcl|S1:50-60] out of range\\n')")
    rows = nv.validate_record(_rec(), tmp_path / "out", t2a)
    assert [(r["severity"], r["code"]) for r in rows] == [("REJECT", "SEQ_FEAT.Range")]


def _gbf_dir(tmp_path, recs):
    d = tmp_path / "gb"
    d.mkdir()
    SeqIO.write(recs, str(d / "a.gbf"), "genbank")
    return d


def test_submission_counts_reject_and_fails_when_not_validated(tmp_path, monkeypatch):
    import mito_env
    import argparse
    gb = _gbf_dir(tmp_path, [_rec("S1"), _rec("S2")])
    t2a = _fake_table2asn(tmp_path, "open(out[:-4] + '.val', 'w').write("
                                    "'Critical: valid [SEQ_FEAT.Range] x [lcl|S1:1-2]\\n'"
                                    " + 'FATAL: valid [SEQ_INST.Bad] y\\n')")
    monkeypatch.setattr(mito_env, "find", lambda name, *a, **k: t2a)
    res = mitom_steps.run_submission(argparse.Namespace(), gb, tmp_path / "o1")
    rows = (tmp_path / "o1" / "06.ncbi_submission" / "submission_summary.tsv").read_text().splitlines()[1:]
    assert [r.split("\t")[1] for r in rows] == ["2", "2"] and res.ok
    monkeypatch.setattr(mito_env, "find", lambda name, *a, **k: "/bin/false")
    res = mitom_steps.run_submission(argparse.Namespace(), gb, tmp_path / "o2")
    assert not res.ok
    rows = (tmp_path / "o2" / "06.ncbi_submission" / "submission_summary.tsv").read_text().splitlines()[1:]
    assert all(r.split("\t")[1] == "1" and "MITOM.Table2asnFailed" in r for r in rows)


def test_gb2tbl_failed_table2asn_marks_every_seqid(tmp_path):
    out = tmp_path / "o"
    out.mkdir()
    per = mito_tbl.validate(out, "sub", "/bin/false", None, ["A", "B"])
    assert per["A"]["ERROR"] == 1 and per["B"]["ERROR"] == 1
    t2a = _fake_table2asn(tmp_path, "open(out[:-4] + '.val', 'w').write("
                                    "'Critical: valid [SEQ_FEAT.Range] x [lcl|A:1-2]\\n')")
    per = mito_tbl.validate(out, "sub", t2a, None, ["A", "B"])
    assert per["A"]["ERROR"] == 1 and per["B"]["ERROR"] == 0


# ---- inputs ------------------------------------------------------------------------------------------------------
def _overlap_fasta(path, tag=""):
    import random
    rnd = random.Random(5)
    s = "".join(rnd.choice("ACGT") for _ in range(12000))
    path.write_text(f">c1{tag}\n{s + s[:100]}\n")


def test_topology_linear_keeps_the_terminal_repeat(tmp_path):
    fa = tmp_path / "a.fa"
    _overlap_fasta(fa)
    for topo, n in ((None, 12000), ("linear", 12100)):
        samples, _ = mitom_inputs.split_or_sanitize_fastas([fa], tmp_path / f"s{topo}", True, 20, topology=topo)
        rec = next(SeqIO.parse(str(samples[0].internal_fasta), "fasta"))
        assert len(rec.seq) == n


def test_unpacked_name_collision_is_an_error(tmp_path):
    d = tmp_path / "in"
    d.mkdir()
    (d / "same.fa").write_text(">a\nACGT\n")
    (d / "same.fa.gz").write_bytes(gzip.compress(b">b\nACGT\n"))
    with pytest.raises(mito_common.InputError, match="both become"):
        mitom_inputs.normalize_input(d, tmp_path / "out")


def test_gzip_of_gzip_is_unpacked(tmp_path):
    f = tmp_path / "a.fa.gz"
    f.write_bytes(gzip.compress(gzip.compress(b">a\nACGT\n")))
    new, _ = mitom_inputs.normalize_input(f, tmp_path / "out")
    assert Path(new).read_text().startswith(">a")


def test_genbank_without_topology_short_is_linear(tmp_path):
    r = _rec("X1", 2000)
    del r.annotations["topology"]
    p = tmp_path / "x.gb"
    SeqIO.write(r, str(p), "genbank")
    txt = p.read_text().replace(" linear ", "        ").replace(" circular ", "          ")
    p.write_text(txt)
    kept, _, _ = mitom_inputs.dedup_gb_records([p])
    assert kept[0].topology == "linear"


def test_auto_clade_tunicate_without_cds_uses_lineage_code(tmp_path):
    r = _rec("T1", 400)
    r.annotations["taxonomy"] = ["Eukaryota", "Metazoa", "Chordata", "Tunicata", "Ascidiacea"]
    SeqIO.write(r, str(tmp_path / "t.gb"), "genbank")
    clade, code, phylum, msg = mitom_inputs.auto_clade(tmp_path / "t.gb")
    assert (phylum, code, msg) == ("Chordata", 13, "")


# ---- output folder -----------------------------------------------------------------------------------------------
def test_previous_results_removed_but_never_the_input(tmp_path):
    out = tmp_path / "out"
    (out / "final_annotation").mkdir(parents=True)
    (out / "final_annotation" / "a.gb").write_text("x")
    (out / "plots").mkdir()
    (out / "original_vs_mitom.tsv").write_text("old")
    (out / "my_notes.txt").write_text("user's")
    mitom_output.clear_previous_work(out, protect=[out / "final_annotation"])
    assert (out / "final_annotation" / "a.gb").exists()          # the input of this run
    assert not (out / "plots").exists() and not (out / "original_vs_mitom.tsv").exists()
    assert (out / "my_notes.txt").exists()
    mitom_output.clear_previous_work(out)
    assert not (out / "final_annotation").exists()


def test_lock_of_another_users_run_is_alive(tmp_path, monkeypatch):
    out = tmp_path / "o"
    out.mkdir()
    (out / bm.LOCK_NAME).write_text("12345\n")

    def kill(pid, sig):
        raise PermissionError(1, "Operation not permitted")
    monkeypatch.setattr(os, "kill", kill)
    assert "another MitoM run" in bm.lock_output(out)


def test_import_is_prepared_after_the_lock(tmp_path, monkeypatch):
    out = tmp_path / "o"
    out.mkdir()
    (out / bm.LOCK_NAME).write_text(f"{os.getppid()}\n")
    called = []
    monkeypatch.setattr(bm, "prepare_import", lambda a: called.append(1) or "")
    monkeypatch.setattr(sys, "argv", ["mitom", "import", "-i", str(tmp_path), "-o", str(out)])
    assert bm.main() == 2 and not called


# ---- groups ------------------------------------------------------------------------------------------------------
def test_concat_tsv_with_different_headers(tmp_path):
    a, b, out = tmp_path / "a.tsv", tmp_path / "b.tsv", tmp_path / "o.tsv"
    a.write_text("x\ty\n1\t2\n")
    b.write_text("x\tz\n3\t4\n")
    mitom_groups._concat_tsv([a, b], out)
    assert out.read_text().splitlines() == ["x\ty\tz", "1\t2\t", "3\t\t4"]


def test_failed_group_is_not_merged(tmp_path):
    outs = {}
    for g in ("ok", "bad"):
        o = tmp_path / g / "out"
        (o / "final_annotation").mkdir(parents=True)
        (o / "final_annotation" / f"{g}.gb").write_text(g)
        (o / "reports").mkdir()
        (o / "reports" / "x.txt").write_text(g + "\n")
        outs[g] = o
    root = tmp_path / "root"
    root.mkdir()
    mitom_groups.merge_outputs(root, outs, failed=["bad"])
    assert sorted(p.name for p in (root / "final_annotation").iterdir()) == ["ok.gb"]
    assert (root / "reports" / "x.txt").read_text() == "ok\nbad\n"       # the failed group's reports are kept


# ---- annotation steps --------------------------------------------------------------------------------------------
def test_intron_exonerate_gets_a_table_for_codes_it_lacks(tmp_path, monkeypatch):
    seen = []
    monkeypatch.setattr(mito_intron, "run_tool", lambda cmd, *a, **k: seen.append(cmd) or "")
    try:
        mito_intron.run_exonerate("exonerate", tmp_path / "p.faa", tmp_path / "g.fna", 33, 5000)
    except Exception:
        pass
    cmd = seen[0]
    code = cmd[cmd.index("--geneticcode") + 1]
    assert len(code) == 64                                     # code 33: exonerate 2.4 has no such table
    seen.clear()
    try:
        mito_intron.run_exonerate("exonerate", tmp_path / "p.faa", tmp_path / "g.fna", 5, 5000)
    except Exception:
        pass
    assert seen[0][seen[0].index("--geneticcode") + 1] == "5"


def test_rescue_knows_genes_named_in_another_style():
    r = _rec()
    r.features += [SeqFeature(FeatureLocation(0, 90, 1), type="CDS", qualifiers={"gene": ["cox1"]}),
                   SeqFeature(FeatureLocation(100, 190, 1), type="rRNA", qualifiers={"gene": ["rrnS"]})]
    have = mito_rescue.present_genes(r)
    assert {"COX1", "s-rRNA"} <= have


def test_rotate_leaves_out_a_feature_it_cannot_move(monkeypatch):
    import mitom_orient
    r = _rec(n=400)
    r.annotations["topology"] = "circular"
    r.features.append(SeqFeature(FeatureLocation(10, 100, 1), type="CDS", qualifiers={"gene": ["COX1"]}))
    real = mitom_orient.shift_location

    def shift(loc, start, n):
        if int(loc.start) == 10:
            raise ValueError("cannot")
        return real(loc, start, n)
    monkeypatch.setattr(mitom_orient, "shift_location", shift)
    out = mitom_orient.rotate_record(r, 50)
    assert [f.type for f in out.features] == ["source"]


def test_rescue_shares_cpus_between_mitos2_runs(monkeypatch, tmp_path):
    got = {}
    monkeypatch.setattr(mito_rescue, "run_mitos", lambda *a, **k: got.update(k) or (_ for _ in ()).throw(RuntimeError("x")))
    r = _rec(n=400)
    SeqIO.write(r, str(tmp_path / "a.gbf"), "genbank")
    monkeypatch.setattr(os, "cpu_count", lambda: 16)
    mito_rescue.rescue_dir(tmp_path, 5, jobs=4)
    assert got["jobs"] == 4


def test_mitos2_command_passes_jobs():
    cmd = mito_engines.mitos2_command(Path("a.fa"), Path("o"), 5, "/b", "/r", mode="sensitive", jobs=3)
    assert cmd[cmd.index("--jobs") + 1] == "3"
    assert "--jobs" not in mito_engines.mitos2_command(Path("a.fa"), Path("o"), 5, "/b", "/r", mode="full")


# ---- import / metadata -------------------------------------------------------------------------------------------
def test_find_results_keeps_genbank_outside_mitoz_folders(tmp_path):
    (tmp_path / "mz" / "S.result").mkdir(parents=True)
    (tmp_path / "mz" / "S.result" / "S.gbf").write_text("x")
    (tmp_path / "mz" / "tmp").mkdir()
    (tmp_path / "mz" / "tmp" / "draft.gbf").write_text("x")
    (tmp_path / "other.gb").write_text("x")
    gbs, _ = mito_engines.find_results(tmp_path)
    assert sorted(p.name for p in gbs) == ["S.gbf", "other.gb"]


def test_accession_from_the_full_id():
    r = _rec("SAMPLE_WITH_A_LONG_NAME_1")
    r.name = "SAMPLE_WITH_A_LO"
    r.annotations.pop("accessions", None)
    assert mitom_metadata.gb_record_accession(r) == "SAMPLE_WITH_A_LONG_NAME_1"


def test_metadata_length_follows_the_new_record():
    meta = "LOCUS       NC_1                   16000 bp    DNA     circular INV 01-JAN-2020\n"
    rec = "LOCUS       NC_1                   15990 bp    DNA     circular INV 01-JAN-2020\n"
    assert "15990 bp" in mitom_metadata._with_length(meta, rec)


# ---- install.sh --------------------------------------------------------------------------------------------------
def test_install_sh_macos_asset_and_mitotools():
    txt = (ROOT / "install.sh").read_text()
    assert 'MF_OS="MacOSX"' in txt
    assert 'ENVS="mitotools $ENVS"' in txt


def test_lock_of_a_reused_pid_is_stale(tmp_path):
    out = tmp_path / "o"
    out.mkdir()
    (out / bm.LOCK_NAME).write_text(f"{os.getppid()} 1\n")        # a running pid, but started at another time
    assert bm.lock_output(out) == ""
    pid, started = (out / bm.LOCK_NAME).read_text().split()
    assert pid == str(os.getpid()) and started == bm._process_start(os.getpid())


# ---- GLM 5.3 review (verified 2026-11-05) ------------------------------------------------------------------------
def test_input_trna_across_the_origin_is_not_the_whole_record():
    import mito_trna
    from Bio.SeqFeature import CompoundLocation
    r = _rec(n=16000)
    r.features.append(SeqFeature(CompoundLocation([FeatureLocation(15970, 16000, 1), FeatureLocation(0, 35, 1)]),
                                 type="tRNA", qualifiers={"gene": ["trnY"], "product": ["tRNA-Tyr"]}))
    p = mito_trna.from_input(r)[0]
    assert (p.start, p.end) == (15971, 16035)


def test_trna_call_between_two_others_joins_their_group():
    import mito_trna
    ps = [mito_trna.Pred(t, s, e, 1, "F", "") for t, s, e in (("arwen", 1, 70), ("trnascan", 30, 100),
                                                               ("mitos2", 60, 130))]
    assert len(mito_trna.group(ps)) == 1


def test_ncr_command_keeps_pseudo_value_less(tmp_path):
    import subprocess
    r = _rec(n=2000)
    r.features.append(SeqFeature(FeatureLocation(10, 400, 1), type="gene", qualifiers={"gene": ["nad1"], "pseudo": [None]}))
    SeqIO.write(r, str(tmp_path / "a.gb"), "genbank")
    subprocess.run([sys.executable, "-I", str(ROOT / "mito_ncr.py"), "-i", str(tmp_path / "a.gb"), "-o",
                    str(tmp_path / "b.gb")], check=True, capture_output=True)
    txt = (tmp_path / "b.gb").read_text()
    assert "/pseudo\n" in txt and '/pseudo=""' not in txt


def test_explicit_chordata_uses_the_tunicate_code(tmp_path):
    r = _rec("T1", 400)
    r.annotations["taxonomy"] = ["Eukaryota", "Metazoa", "Chordata", "Tunicata", "Ascidiacea"]
    SeqIO.write(r, str(tmp_path / "t.gb"), "genbank")
    args = bm.build_parser().parse_args(["run", "-i", str(tmp_path / "t.gb"), "-o", str(tmp_path / "o"),
                                         "--clade", "Chordata"])
    bm.resolve_tools(args)
    assert args.genetic_code == 13


# ---- Grok 4.7 review (verified 2026-11-05) -----------------------------------------------------------------------
def test_arabic_numeral_cox_product_and_bare_atp():
    import mito_names
    assert mito_names.gene_name({"product": ["cytochrome c oxidase subunit 1"]}, "CDS") == "cox1"
    assert mito_names.gene_name({"gene": ["atp"]}, "CDS") == ""             # ambiguous: atp6 or atp8


def test_pterobranch_and_urochordata_codes():
    for tax, code in ((["Metazoa", "Hemichordata", "Pterobranchia", "Rhabdopleuridae"], 24),
                      (["Metazoa", "Hemichordata", "Pterobranchia", "Cephalodiscidae"], 33),
                      (["Metazoa", "Chordata", "Urochordata", "Ascidiacea"], 13),
                      (["Metazoa", "Hemichordata", "Enteropneusta"], 9)):
        r = _rec()
        r.annotations["taxonomy"] = tax
        r.features.append(SeqFeature(FeatureLocation(0, 90, 1), type="CDS", qualifiers={"gene": ["COX1"]}))
        assert mito_tbl.record_code(r, None)[0] == code


def test_fetch_query_takes_complete_mitochondrial_genome_titles():
    import mito_fetch
    assert '"complete mitochondrial genome"[Title]' in mito_fetch.MT_CORE


def test_curate_without_code_uses_the_taxonomy_and_writes_it(tmp_path):
    import mito_curate
    r = next(SeqIO.parse(str(HERE / "data" / "NC_000933.gb"), "genbank"))
    for f in r.features:
        f.qualifiers.pop("transl_table", None)
    SeqIO.write(r, str(tmp_path / "a.gb"), "genbank")
    mito_curate.main(["-i", str(tmp_path / "a.gb"), "-o", str(tmp_path / "o"), "--no_validate"])
    out = next((tmp_path / "o").rglob("*.gb*"))
    codes = {f.qualifiers.get("transl_table", [""])[0] for f in SeqIO.read(str(out), "genbank").features if f.type == "CDS"}
    assert codes == {"4"}                                                   # Cnidaria (NC_000933: Metridium)


# ---- Cursor review of the fixes (2026-11-06) ---------------------------------------------------------------------
def test_explicit_chordata_with_two_codes_is_an_error(tmp_path):
    recs = []
    for rid, sub_ in (("V1", ["Vertebrata", "Mammalia"]), ("T1", ["Tunicata", "Ascidiacea"])):
        r = _rec(rid, 400)
        r.annotations["taxonomy"] = ["Eukaryota", "Metazoa", "Chordata"] + sub_
        recs.append(r)
    SeqIO.write(recs, str(tmp_path / "m.gb"), "genbank")
    args = bm.build_parser().parse_args(["run", "-i", str(tmp_path / "m.gb"), "-o", str(tmp_path / "o"),
                                         "--clade", "Chordata"])
    err = bm.resolve_tools(args)
    assert err.startswith("ERROR") and "--clade auto" in err


def test_auto_clade_reads_gzip_genbank(tmp_path):
    import gzip
    r = _rec("T1", 400)
    r.annotations["taxonomy"] = ["Eukaryota", "Metazoa", "Chordata", "Tunicata", "Ascidiacea"]
    buf = io.StringIO()
    SeqIO.write(r, buf, "genbank")
    (tmp_path / "t.gb.gz").write_bytes(gzip.compress(buf.getvalue().encode()))
    from mitom_inputs import auto_clade
    assert auto_clade(tmp_path / "t.gb.gz")[:2] == ("Chordata", 13)


def test_double_gzip_loses_every_gz_suffix(tmp_path):
    import gzip
    from mitom_inputs import normalize_input
    d = tmp_path / "in"
    d.mkdir()
    (d / "a.fa.gz.gz").write_bytes(gzip.compress(gzip.compress(b">a\nACGT\n")))
    out, _ = normalize_input(d, tmp_path / "o")
    assert [p.name for p in out.iterdir()] == ["a.fa"]
    out, _ = normalize_input(d / "a.fa.gz.gz", tmp_path / "o2")
    assert out.name == "a.fa" and out.read_text() == ">a\nACGT\n"


def test_failed_group_tables_are_not_merged(tmp_path):
    outs = {}
    for g in ("ok", "bad"):
        o = tmp_path / g / "out"
        (o / "reports").mkdir(parents=True)
        (o / "reports" / "t.tsv").write_text(f"sample\n{g}\n")
        (o / "reports" / "only_bad.tsv").write_text("sample\nx\n") if g == "bad" else None
        outs[g] = o
    root = tmp_path / "root"
    root.mkdir()
    mitom_groups.merge_outputs(root, outs, failed=["bad"])
    assert (root / "reports" / "t.tsv").read_text().splitlines() == ["sample", "ok"]
    assert not (root / "reports" / "only_bad.tsv").exists()


def test_call_between_two_separate_trnas_does_not_merge_them():
    import mito_trna
    ps = [mito_trna.Pred(t, s, e, 1, "F", "") for t, s, e in (("arwen", 1, 70), ("trnascan", 36, 105),
                                                               ("mitos2", 71, 140))]
    assert len(mito_trna.group(ps)) == 2


def test_site_at_the_origin_of_a_circular_record():
    r = _rec("C1", 400)
    r.annotations["topology"] = "circular"
    r.features.append(SeqFeature(FeatureLocation(400, 400, 1), type="misc_feature", qualifiers={"note": ["x"]}))
    loc = mito_tbl.clean_record(r)[0].features[-1].location
    assert [(int(p.start), int(p.end)) for p in loc.parts] == [(399, 400), (0, 1)]
    r.annotations["topology"] = "linear"
    loc = mito_tbl.clean_record(r)[0].features[-1].location
    assert (int(loc.start), int(loc.end)) == (399, 400)


def test_trna_conflict_found_through_a_later_member():
    import mito_trna
    f = [mito_trna.Pred("arwen", 1, 70, 1, "F", ""), mito_trna.Pred("mitos2", 30, 100, 1, "F", "")]
    v = [mito_trna.Pred("trnascan", 60, 125, 1, "V", "")]
    groups = mito_trna.group(f + v)
    assert len(groups) == 2 and mito_trna.conflicts(groups)       # F's 2nd call overlaps V; first calls do not


def test_call_without_anticodon_does_not_chain_two_leu_trnas():
    import mito_trna
    ps = [mito_trna.Pred("arwen", 1, 70, 1, "L", "TAG"), mito_trna.Pred("trnascan", 30, 100, 1, "L", ""),
          mito_trna.Pred("mitos2", 60, 130, 1, "L", "TAA")]
    assert len(mito_trna.group(ps)) == 2


def test_curate_keeps_the_code_when_a_record_fails(tmp_path, monkeypatch):
    import mito_curate
    r = next(SeqIO.parse(str(HERE / "data" / "NC_000933.gb"), "genbank"))
    for f in r.features:
        f.qualifiers.pop("transl_table", None)
    SeqIO.write(r, str(tmp_path / "a.gb"), "genbank")

    def boom(self):
        raise RuntimeError("x")
    monkeypatch.setattr(mito_curate.Curator, "run", boom)
    mito_curate.main(["-i", str(tmp_path / "a.gb"), "-o", str(tmp_path / "o"), "--no_validate"])
    out = next((tmp_path / "o").rglob("*.gb*"))
    codes = {f.qualifiers.get("transl_table", [""])[0] for f in SeqIO.read(str(out), "genbank").features if f.type == "CDS"}
    assert codes == {"4"}
