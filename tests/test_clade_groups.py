"""--clade auto per record (mitom_groups.py) and the phylum / genetic code of sequences without a taxonomy
(mito_classify.py), 2026-10-31: a GenBank file of several phyla was refused ("run each phylum separately") and FASTA
input fell back to Annelida.

Run: python3 -m pytest -q tests/test_clade_groups.py
"""
import argparse
import sys
from collections import Counter
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_classify as mc  # noqa: E402
import mito_env  # noqa: E402
import mitom_groups as mg  # noqa: E402

DATA = HERE / "data"


def rec(acc, lineage, code=None, seq="ATG" * 200):
    r = SeqRecord(Seq(seq), id=acc, name=acc.split(".")[0], description=f"{acc} test")
    r.annotations.update(molecule_type="DNA", topology="circular", organism="Genus species",
                         taxonomy=["Eukaryota", "Metazoa", *lineage] if lineage is not None else [])
    if code:
        r.features.append(SeqFeature(FeatureLocation(0, 99, 1), type="CDS", qualifiers={"transl_table": [str(code)]}))
    return r


def args_for(path, **kw):
    return argparse.Namespace(input=str(path), genetic_code=kw.get("genetic_code"), threads=1, max_tasks=1,
                              mitos_refdir=None)


def test_mixed_phyla_genbank_gives_one_group_per_phylum_and_code(tmp_path):
    p = tmp_path / "mixed.gb"
    SeqIO.write([rec("NC_1.1", ["Chordata", "Craniata", "Vertebrata"], 2),
                 rec("NC_2.1", ["Cnidaria", "Scyphozoa"], 4),
                 rec("NC_3.1", ["Chordata", "Tunicata", "Ascidiacea"], 13),
                 rec("NC_4.1", ["Cnidaria", "Anthozoa"], 4)], str(p), "genbank")
    plan = mg.decide(args_for(p))
    assert not plan.error
    got = {k: sorted(r.record for r in v) for k, v in plan.groups().items()}
    assert got == {("Cnidaria", 4, "Cnidaria"): ["NC_2.1", "NC_4.1"], ("Chordata", 2, "Chordata"): ["NC_1.1"],
                   ("Chordata", 13, "Chordata"): ["NC_3.1"]}
    assert {r.record: r.group for r in plan.rows}["NC_3.1"] == "Chordata_gc13"


def test_code_from_lineage_when_there_is_no_cds():
    """A record with a taxonomy but no annotation: the code of the reference genomes of its group (tunicates 13)."""
    assert mc.code_of_lineage("Eukaryota;Metazoa;Chordata;Tunicata;Ascidiacea") == 13
    assert mc.code_of_lineage("Eukaryota;Metazoa;Chordata;Craniata;Vertebrata;Mammalia") == 2
    assert mc.code_of_lineage("Eukaryota;Metazoa;Echinodermata") == 9
    assert mc.code_of_lineage("Eukaryota;Metazoa;Bilateria") == 0          # above a phylum: no code


def test_explicit_genetic_code_is_used_for_every_record(tmp_path):
    p = tmp_path / "two.gb"
    SeqIO.write([rec("NC_1.1", ["Mollusca"], 5), rec("NC_2.1", ["Mollusca"], 2)], str(p), "genbank")
    plan = mg.decide(args_for(p, genetic_code=5))
    assert list(plan.groups()) == [("Mollusca", 5, "Mollusca")]


def test_not_an_animal_is_an_error(tmp_path):
    p = tmp_path / "plant.gb"
    SeqIO.write([rec("NC_1.1", None)], str(p), "genbank")
    r = SeqIO.read(str(p), "genbank")
    r.annotations["taxonomy"] = ["Eukaryota", "Viridiplantae", "Streptophyta"]
    SeqIO.write([r], str(p), "genbank")
    assert "no animal phylum" in mg.decide(args_for(p)).error


def test_split_input_keeps_the_genbank_text_of_each_record(tmp_path):
    p = tmp_path / "mixed.gb"
    SeqIO.write([rec("NC_1.1", ["Chordata"], 2), rec("NC_2.1", ["Cnidaria"], 4), rec("NC_3.1", ["Chordata"], 2)],
                str(p), "genbank")
    plan = mg.decide(args_for(p))
    inputs = mg._split_input(plan, args_for(p), tmp_path / "_groups")
    ids = {g: [r.id for r in SeqIO.parse(str(f), "genbank")] for g, f in inputs.items()}
    assert ids == {"Chordata_gc2": ["NC_1.1", "NC_3.1"], "Cnidaria_gc4": ["NC_2.1"]}
    text = p.read_text()
    block = text[text.index("LOCUS       NC_2"):text.index("//", text.index("LOCUS       NC_2")) + 2]
    assert inputs["Cnidaria_gc4"].read_text().strip() == block.strip()


def test_merge_outputs_concatenates_tables_once_per_header(tmp_path):
    outs = {}
    for g, recs in (("A_gc5", ["r1", "r2"]), ("B_gc2", ["r3"])):
        o = tmp_path / g / "out"
        (o / "final_annotation").mkdir(parents=True)
        (o / "ncbi_submission").mkdir()
        (o / "reports").mkdir()
        for r in recs:
            (o / "final_annotation" / f"{r}.gb").write_text(f"LOCUS {r}\n//\n")
            (o / "ncbi_submission" / f"{r}.tbl").write_text(">Feature\n")
        (o / "final_annotation_all.gb").write_text("".join(f"LOCUS {r}\n//\n" for r in recs))
        (o / "summary.tsv").write_text("record\tprotein_genes\n" + "".join(f"{r}\t13\n" for r in recs))
        (o / "ncbi_submission" / "submission_summary.tsv").write_text("record\terrors\n" + "".join(f"{r}\t0\n" for r in recs))
        (o / "reports" / "dedup_report.tsv").write_text("record\tstatus\n" + "".join(f"{r}\tkept\n" for r in recs))
        (o / "reports" / "run_params.json").write_text("{}")
        (o / "README.txt").write_text("readme\n")
        outs[g] = o
    dest = tmp_path / "merged"
    dest.mkdir()
    (dest / "summary.tsv").write_text("stale\n")                 # a previous run's result is replaced
    mg.merge_outputs(dest, outs)
    assert (dest / "summary.tsv").read_text().splitlines() == ["record\tprotein_genes", "r1\t13", "r2\t13", "r3\t13"]
    assert sorted(p.name for p in (dest / "final_annotation").iterdir()) == ["r1.gb", "r2.gb", "r3.gb"]
    assert (dest / "final_annotation_all.gb").read_text().count("LOCUS") == 3
    assert len((dest / "ncbi_submission" / "submission_summary.tsv").read_text().splitlines()) == 4
    assert len((dest / "reports" / "dedup_report.tsv").read_text().splitlines()) == 4
    assert sorted(p.name for p in (dest / "reports" / "groups").glob("*.json")) == ["A_gc5.run_params.json",
                                                                                    "B_gc2.run_params.json"]
    assert "several phyla" in (dest / "README.txt").read_text()


def test_codon_scores_pick_the_code_of_the_reference_residues():
    """AGR read as Ser, ATA as Met, AAA as Lys by the references: code 5 fits; AGR = Arg, ATA = Ile: code 4."""
    obs = {"AGA": Counter(S=20, R=2), "AGG": Counter(S=10), "ATA": Counter(M=12, I=3), "AAA": Counter(K=30),
           "TGA": Counter(W=40), "TAA": Counter()}
    sc = mc.codon_scores(obs)
    assert max(sc, key=sc.get) == 5
    obs4 = {"AGA": Counter(R=20), "AGG": Counter(R=10), "ATA": Counter(I=15, M=1), "AAA": Counter(K=30),
            "TGA": Counter(W=40), "TAA": Counter()}
    sc4 = mc.codon_scores(obs4)
    assert max(sc4, key=sc4.get) == 4


def test_codons_walk_both_strands():
    genome = "AAATGAAGA" + "CCC"
    assert list(mc._codons(genome, 1, 9, "KWS", "KWS")) == [("AAA", "K"), ("TGA", "W"), ("AGA", "S")]
    rc = str(Seq(genome).reverse_complement())                    # the same codons on the minus strand
    assert list(mc._codons(rc, 12, 4, "KWS", "KWS")) == [("AAA", "K"), ("TGA", "W"), ("AGA", "S")]


REFDIR = mito_env.mitos_refdir(None)
BLAST_OK = bool(REFDIR and mito_env.find("blastx"))


@pytest.mark.skipif(not BLAST_OK, reason="refseq89m / blastx not available")
def test_fasta_without_taxonomy_is_classified(tmp_path):
    """A FASTA genome (Aphrodita, Annelida, code 5) without any taxonomy: Annelida, 5, from blastx."""
    src = next(SeqIO.parse(str(DATA / "KY753833_mitoz.gbf"), "genbank"))
    fa = tmp_path / "x.fasta"
    fa.write_text(f">sample1 topology=circular\n{src.seq}\n")
    plan = mg.decide(argparse.Namespace(input=str(fa), genetic_code=None, threads=8, max_tasks=1, mitos_refdir=None))
    (key, rows), = plan.groups().items()
    assert key == ("Annelida-segmented-worms", 5, "Annelida") and rows[0].source.startswith("blastx")


def test_codons_seen_no_more_often_than_tag_are_noise():
    """2026-10-31 stress test: a vertebrate genome with 10 % random substitutions has AGR and TAA inside its genes about
    as often as TAG (a stop in every code); they must not exclude code 2 or pick 13 / 33."""
    obs = {"AGA": Counter(G=4, S=3, R=3), "AGG": Counter(G=3, S=2), "ATA": Counter(M=30, I=10), "AAA": Counter(K=40),
           "TGA": Counter(W=50), "TAA": Counter(Y=4, L=3), "TAG": Counter(L=3, F=2)}
    sc = mc.codon_scores(obs, n_codons=8000)
    assert 2 in sc and max(sc, key=sc.get) in (2, 5)       # AGR / TAA not informative: 2 and 5 tie on ATA AAA TGA


def test_merge_never_overwrites_a_name_of_another_group(tmp_path):
    outs = {}
    for g in ("A_gc5", "B_gc4"):
        o = tmp_path / g / "out"
        (o / "final_annotation").mkdir(parents=True)
        (o / "final_annotation" / "NODE_1.gb").write_text(f"LOCUS {g}\n//\n")
        outs[g] = o
    dest = tmp_path / "m"
    dest.mkdir()
    mg.merge_outputs(dest, outs)
    got = sorted(p.name for p in (dest / "final_annotation").iterdir())
    assert got == ["NODE_1.B_gc4.gb", "NODE_1.gb"]


def test_fasta_names_unique_over_groups(tmp_path):
    """2026-10-31 review: SPAdes 'NODE_1' in two files of two phyla: one result replaced the other."""
    p = tmp_path / "in"
    p.mkdir()
    (p / "a.fasta").write_text(">NODE_1 x\n" + "ACGT" * 300 + "\n")
    (p / "b.fasta").write_text(">NODE_1 y\n" + "ACGT" * 300 + "\n")
    plan = mg.Plan(rows=[mg.Row("NODE_1", "a.fasta", group="A_gc5"), mg.Row("NODE_1", "b.fasta", group="B_gc4")])
    ins = mg._split_input(plan, argparse.Namespace(input=str(p)), tmp_path / "_groups")
    heads = sorted(line for f in ins.values() for x in (f.iterdir() if f.is_dir() else [f])
                   for line in x.read_text().splitlines() if line.startswith(">"))
    assert heads == [">NODE_1 x", ">NODE_1_2 y"]


def test_one_bad_record_does_not_stop_the_classification_of_the_others(monkeypatch):
    def fake(rid, seq, *a, **k):
        if rid == "bad":
            raise RuntimeError("blastx failed on cox1.fas")
        return mc.Classification(rid, phylum="Chordata", genetic_code=2, support=1.0)
    monkeypatch.setattr(mc, "classify_seq", fake)
    monkeypatch.setattr(mc.mito_env, "find", lambda *a, **k: "/bin/true")
    import tempfile
    d = Path(tempfile.mkdtemp())
    (d / "featureProt").mkdir()
    got = mc.classify_records([("good", "ACGT"), ("bad", "ACGT")], str(d), threads=2)
    assert got[0].phylum == "Chordata" and "blastx failed" in got[1].note
