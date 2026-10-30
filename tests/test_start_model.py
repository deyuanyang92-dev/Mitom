"""P11, the learned start codon (mito_start.py, mito_curate.py): candidates and features, the JSON trees, the rule in
a curation run and the training dump (2026-10-30).

Run: python3 -m pytest -q tests/test_start_model.py
"""
import argparse
import csv
import json
import sys
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_curate  # noqa: E402
import mito_start  # noqa: E402

STARTS, STOPS = {"ATG", "ATA", "ATT", "ATC", "GTG", "TTG"}, {"TAA", "TAG"}


def test_candidates_in_frame_open_and_featured():
    #            1   4   7   10  13  16  19  22
    s = "GCC" * 3 + "TAA" + "ATA" + "GCC" + "ATG" + "GCC" * 40 + "TAA"   # stop at 10, ATA at 13, ATG at 19 (P1)
    rows = mito_start.candidates(s, 19, len(s), STARTS, STOPS, core_a=13, aln=(5, 200, 220), up=(11, "tRNA", True),
                                 gene="ND2", phylum="Annelida", code=5)
    pos = [r["pos"] for r in rows]
    assert 13 in pos and 19 in pos and 1 not in pos                      # the stop at 10 closes the frame upstream
    ata = next(r for r in rows if r["pos"] == 13)
    assert (ata["codon"], ata["cand"], ata["rank"], ata["is_p1"], ata["core_off"]) == ("ATA", -6, 2, 0, 0)
    assert (ata["aln_missing5"], ata["aln_cover"], ata["up_gap"], ata["first_after_up"]) == (4, 0.891, 1, 1)
    p1 = next(r for r in rows if r["is_p1"])
    assert p1["pos"] == 19 and p1["first_after_up"] == 0
    none = mito_start.candidates(s, 19, len(s), STARTS, STOPS, None, None, None, "ND2", "", 5)
    assert none[0]["core_off"] == mito_start.MISSING and none[0]["up_type"] == "none" and none[0]["phylum"] == "?"


def tiny_model(tmp_path, prefer_codon="ATA"):
    """One tree: codon=<prefer_codon> -> +1, else 0 (columns as the training writes them)."""
    cats = {"codon": ["ATA", "ATG"], "gene": ["ND2"], "phylum": ["Annelida"], "up_type": ["none", "tRNA"]}
    cols = mito_start.columns(cats)
    j = cols.index(f"codon={prefer_codon}")
    tree = {"left": [1, -1, -1], "right": [2, -1, -1], "feature": [j, 0, 0], "threshold": [0.5, 0, 0],
            "value": [0.0, 0.0, 1.0]}
    path = tmp_path / "m.json"
    path.write_text(json.dumps({"columns": cols, "categories": cats, "trees": [tree]}))
    return path


def test_json_trees_score_and_choose(tmp_path):
    m = mito_start.load(str(tiny_model(tmp_path)))
    rows = [{"codon": "ATG", "gene": "ND2", "phylum": "Annelida", "up_type": "tRNA", "is_p1": 1,
             **{k: 0 for k in mito_start.NUMERIC if k != "is_p1"}},
            {"codon": "ATA", "gene": "ND2", "phylum": "Annelida", "up_type": "tRNA", "is_p1": 0,
             **{k: 0 for k in mito_start.NUMERIC if k != "is_p1"}},
            {"codon": "GTG", "gene": "COX9", "phylum": "Unknown", "up_type": "x", "is_p1": 0,      # unseen values
             **{k: 0 for k in mito_start.NUMERIC if k != "is_p1"}}]
    i, scores = m.best(rows)
    assert scores == [0.0, 1.0, 0.0] and i == 1
    assert mito_start.load(str(tmp_path / "missing.json")) is None


def record():
    """ND2 whose core start (ATG, P1 keeps it) lies 6 nt after an in-frame ATA right behind a tRNA."""
    body = "ATA" + "GCC" + "ATG" + "GCC" * 60 + "TAA"
    seq = "C" * 60 + "G" * 66 + "CC" + body + "C" * 200
    rec = SeqRecord(Seq(seq), id="X1", description="")
    rec.annotations.update(molecule_type="DNA", topology="linear", taxonomy=["Eukaryota", "Metazoa", "Annelida"])
    a = 60 + 66 + 2
    rec.features = [SeqFeature(FeatureLocation(0, len(seq), 1), type="source"),
                    SeqFeature(FeatureLocation(60, 126, 1), type="gene", qualifiers={"gene": ["trnW"]}),
                    SeqFeature(FeatureLocation(60, 126, 1), type="tRNA", qualifiers={"gene": ["trnW"], "product": ["tRNA-Trp"]}),
                    SeqFeature(FeatureLocation(a + 6, a + len(body), 1), type="gene", qualifiers={"gene": ["nad2"]}),
                    SeqFeature(FeatureLocation(a + 6, a + len(body), 1), type="CDS",
                               qualifiers={"gene": ["nad2"], "codon_start": ["1"], "transl_table": ["5"],
                                           "note": ["ref NC_1:nad2 score 300 residues 3-60/62"]})]
    return rec, a


def curate(rec, **kw):
    a = argparse.Namespace(genetic_code=5, rrna_bounds="fill", max_fill=300, max_start_shift=30, aggressive_start=False,
                           aggressive_start_genes="", dump_starts=None, **kw)
    cur = mito_curate.Curator(rec, 5, a)
    return cur.run(), cur


def cds(rec):
    return next(f for f in rec.features if f.type == "CDS")


def test_p11_moves_the_start_and_dumps_candidates(tmp_path):
    rec, a = record()
    rows, cur = curate(rec, start_model_obj=mito_start.load(str(tiny_model(tmp_path))))
    assert int(cds(rec).location.start) == a                              # the ATA right after the tRNA
    assert any(r["rule"] == "P11" for r in rows)
    rec, a = record()
    rows, cur = curate(rec, start_model_obj=None)                         # no model: P1's start stays
    assert int(cds(rec).location.start) == a + 6 and not any(r["rule"] == "P11" for r in rows)
    rec, a = record()
    args = dict(start_model_obj=None)
    rows, cur = curate(rec, **args)
    cur2 = mito_curate.Curator(record()[0], 5, argparse.Namespace(**{**vars(cur.a), "dump_starts": "x"}))
    cur2.run()
    dumped = {r["codon"]: r for r in cur2.start_rows}
    assert set(dumped) >= {"ATA", "ATG"} and dumped["ATG"]["is_p1"] == 1
    assert dumped["ATA"]["aln_missing5"] == 2 and dumped["ATA"]["phylum"] == "Annelida"
    assert dumped["ATA"]["genome_pos"] == a + 1 and dumped["ATA"]["up_gap"] == 2
    rec, a = record()                                                     # another annotator's CDS (no core note): kept
    cds(rec).qualifiers["note"] = ["exon 1"]
    rows, cur = curate(rec, start_model_obj=mito_start.load(str(tiny_model(tmp_path))))
    assert int(cds(rec).location.start) == a + 6 and not any(r["rule"] == "P11" for r in rows)


def test_shipped_model_loads_and_ranks():
    m = mito_start.load()
    assert m is not None and m.meta["cross_validation_by_family"]["all"]["model_exact"] > \
        m.meta["cross_validation_by_family"]["all"]["p1_exact"]
    s = "ATA" + "GCC" + "ATG" + "GCC" * 60 + "TAA"
    rows = mito_start.candidates(s, 7, len(s), STARTS, STOPS, 7, (1, 60, 60), (0, "tRNA", True), "ND2", "Annelida", 5)
    i, scores = m.best(rows)
    assert len(scores) == len(rows) >= 2 and 0 <= i < len(rows)


def test_curate_cli_dump_and_model_options(tmp_path):
    rec, _ = record()
    from Bio import SeqIO
    src = tmp_path / "in"
    src.mkdir()
    SeqIO.write(rec, str(src / "X1.gbf"), "genbank")
    assert mito_curate.main(["-i", str(src), "-o", str(tmp_path / "out"), "--no_validate", "--start_model", "no",
                             "--dump_starts", str(tmp_path / "starts.tsv")]) == 0
    rows = list(csv.DictReader(open(tmp_path / "starts.tsv"), delimiter="\t"))
    assert rows and {"record", "strand", "genome_pos", "codon", "is_p1"} <= set(rows[0])
