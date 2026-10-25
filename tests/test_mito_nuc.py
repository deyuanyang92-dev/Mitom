"""Tests for mito_nuc.py (rDNA and histone annotation of contigs by reference projection).

Run: python3 -m pytest -q tests/test_mito_nuc.py   (the end-to-end tests need blastn, tblastn, makeblastdb, exonerate)
"""
import csv
import random
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_env  # noqa: E402
import mito_nuc as nuc  # noqa: E402

rnd = random.Random(11)
TOOLS = all(mito_env.find(t) for t in ("blastn", "tblastn", "makeblastdb", "exonerate"))
needs_tools = pytest.mark.skipif(not TOOLS, reason="BLAST+ / exonerate not found")
CODONS = {aa: [c for c in (a + b + d for a in "ACGT" for b in "ACGT" for d in "ACGT")
               if str(Seq(c).translate()) == aa] for aa in "ACDEFGHIKLMNPQRSTVWY"}
H3 = ("MARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVG"
      "LFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA")


def rand(n):
    return "".join(rnd.choice("ACGT") for _ in range(n))


def mutate(s, rate):
    return "".join(rnd.choice("ACGT".replace(c, "")) if rnd.random() < rate else c for c in s)


def backtranslate(prot, seed):
    r = random.Random(seed)
    return "".join(r.choice(CODONS[a]) for a in prot) + "TAA"


UNIT = [("18S", 1800), ("ITS1", 350), ("5.8S", 160), ("ITS2", 300), ("28S", 3400)]


def unit(genes):
    """dict gene -> sequence for a reference rDNA unit"""
    return {g: rand(n) for g, n in UNIT} if genes is None else genes


def rdna_record(rid, parts, flank5, flank3, organism):
    seq, feats, pos = flank5, [], len(flank5)
    for g, _ in UNIT:
        s = parts[g]
        ftype = "misc_RNA" if g.startswith("ITS") else "rRNA"
        feats.append(SeqFeature(FeatureLocation(pos, pos + len(s), 1), type=ftype,
                                qualifiers={"product": [nuc.PRODUCT[g]]}))
        seq += s
        pos += len(s)
    seq += flank3
    r = SeqRecord(Seq(seq), id=rid, name=rid, description="")
    r.annotations.update(molecule_type="DNA", organism=organism, taxonomy=["Eukaryota", "Metazoa", "Annelida"])
    r.features = feats
    return r


def histone_record(rid, prot, organism, seed):
    cds = backtranslate(prot, seed)
    seq = rand(200) + cds + rand(200)
    r = SeqRecord(Seq(seq), id=rid, name=rid, description="")
    r.annotations.update(molecule_type="DNA", organism=organism, taxonomy=["Eukaryota", "Metazoa", "Annelida"])
    r.features = [SeqFeature(FeatureLocation(200, 200 + len(cds), 1), type="CDS",
                             qualifiers={"product": ["histone H3"], "codon_start": ["1"], "translation": [prot]})]
    return r


@pytest.fixture(scope="module")
def lib(tmp_path_factory):
    d = tmp_path_factory.mktemp("lib")
    base = unit(None)
    refs = []
    for i in range(4):          # four related references (~5 % divergent rRNA, unrelated spacers)
        parts = {g: (mutate(s, 0.05) if not g.startswith("ITS") else rand(len(s))) for g, s in base.items()}
        refs.append(rdna_record(f"REF{i}", parts, rand(300), rand(300), f"Genus{i} sp"))
    refs.append(histone_record("HREF1", H3, "Genus0 sp", 1))
    gb = d / "refs.gb"
    SeqIO.write(refs, str(gb), "genbank")
    n_r, n_h = nuc.build([gb], d / "lib")
    assert (n_r, n_h) == (4, 1)
    return {"dir": d / "lib", "base": base}


def read_rows(out):
    return list(csv.DictReader(open(out / "nuc_features.tsv"), delimiter="\t"))


def test_build_skips_non_histone_proteins(tmp_path):
    r = histone_record("X1", H3, "Genus0 sp", 2)
    r.features[0].qualifiers["product"] = ["histone H4 transcription factor"]
    SeqIO.write([r], str(tmp_path / "x.gb"), "genbank")
    assert nuc.build([tmp_path / "x.gb"], tmp_path / "lib") == (0, 0)


def test_orf_extends_to_stop_and_upstream_atg():
    cds = "ATG" + "GCT" * 20 + "TAA"
    seq = Seq("CCC" + cds + "GGG")
    # model covers codons 2..15 only, protein not from position 1
    exons, p5, p3 = nuc._orf(seq, [(7, 48)], 1, False)
    assert exons == [(4, 4 + len(cds) - 1)] and not p5 and not p3
    rc = seq.reverse_complement()
    n = len(rc)
    exons, p5, p3 = nuc._orf(rc, [(n - 48 + 1, n - 7 + 1)], -1, False)
    assert exons == [(n - (4 + len(cds) - 1) + 1, n - 4 + 1)] and not p5 and not p3


def test_orf_open_at_contig_end_is_partial():
    seq = Seq("GCT" * 10)                        # no ATG, no stop
    exons, p5, p3 = nuc._orf(seq, [(4, 27)], 1, False)
    assert exons == [(1, 30)] and p5 and p3


@needs_tools
def test_rdna_complete_unit_reverse_strand(tmp_path, lib):
    parts = {g: (mutate(s, 0.04) if not g.startswith("ITS") else rand(len(s))) for g, s in lib["base"].items()}
    q = rdna_record("Q1", parts, rand(5000), rand(4000), "Query sp")
    contig = SeqRecord(q.seq.reverse_complement(), id="contig1", description="")
    SeqIO.write([contig], str(tmp_path / "in.fa"), "fasta")
    nuc.annotate(tmp_path / "in.fa", tmp_path / "out", lib["dir"], ["rdna"], 2)
    rows = {r["gene"]: r for r in read_rows(tmp_path / "out")}
    n = len(q.seq)
    for f in q.features:
        g = next(k for k, v in nuc.PRODUCT.items() if v == f.qualifiers["product"][0])
        exp = (n - int(f.location.end) + 1, n - int(f.location.start))
        assert (int(rows[g]["start"]), int(rows[g]["end"]), rows[g]["strand"]) == (*exp, "-"), g
        assert rows[g]["partial_5"] == rows[g]["partial_3"] == "0"
    gb = list(SeqIO.parse(str(tmp_path / "out" / "nuc_annotation.gb"), "genbank"))
    assert len(gb) == 1 and len(gb[0].features) == 5
    assert (tmp_path / "out" / "genes" / "18S.fasta").exists()


@needs_tools
def test_rdna_contig_ending_in_18s_and_two_copies(tmp_path, lib):
    p1 = {g: (mutate(s, 0.04) if not g.startswith("ITS") else rand(len(s))) for g, s in lib["base"].items()}
    p2 = {g: (mutate(s, 0.04) if not g.startswith("ITS") else rand(len(s))) for g, s in lib["base"].items()}
    a = "".join(p1[g] for g, _ in UNIT)
    b = "".join(p2[g] for g, _ in UNIT)
    start = 700                                   # contig starts inside the first 18S
    seq = a[start:] + rand(3000) + b + rand(2000)
    SeqIO.write([SeqRecord(Seq(seq), id="arr", description="")], str(tmp_path / "in.fa"), "fasta")
    nuc.annotate(tmp_path / "in.fa", tmp_path / "out", lib["dir"], ["rdna"], 2)
    s18 = sorted((r for r in read_rows(tmp_path / "out") if r["gene"] == "18S"), key=lambda r: int(r["start"]))
    assert len(s18) == 2
    assert s18[0]["start"] == "1" and s18[0]["partial_5"] == "1" and int(s18[0]["end"]) == 1800 - start
    off = len(a) - start + 3000
    assert (int(s18[1]["start"]), int(s18[1]["end"]), s18[1]["partial_5"]) == (off + 1, off + 1800, "0")


@needs_tools
def test_leave_out_genus_removes_the_only_reference(tmp_path, lib):
    r = histone_record("Q2", H3, "Genus0 sp", 3)
    SeqIO.write([r], str(tmp_path / "in.fa"), "fasta")
    nuc.annotate(tmp_path / "in.fa", tmp_path / "out", lib["dir"], ["histone"], 2)
    rows = read_rows(tmp_path / "out")
    assert [(x["gene"], x["start"], x["end"], x["partial_5"], x["partial_3"]) for x in rows] == \
        [("H3", "201", str(200 + 3 * len(H3) + 3), "0", "0")]
    (tmp_path / "lo.tsv").write_text("Q2\tGenus0\n")
    nuc.annotate(tmp_path / "in.fa", tmp_path / "out2", lib["dir"], ["histone"], 2, tmp_path / "lo.tsv")
    assert read_rows(tmp_path / "out2") == []


@needs_tools
def test_histone_minus_strand_translation(tmp_path, lib):
    r = histone_record("Q3", H3, "Other sp", 4)
    SeqIO.write([SeqRecord(r.seq.reverse_complement(), id="c3", description="")], str(tmp_path / "in.fa"), "fasta")
    nuc.annotate(tmp_path / "in.fa", tmp_path / "out", lib["dir"], ["histone"], 2)
    gb = next(SeqIO.parse(str(tmp_path / "out" / "nuc_annotation.gb"), "genbank"))
    cds = [f for f in gb.features if f.type == "CDS"]
    assert len(cds) == 1 and cds[0].location.strand == -1
    assert cds[0].qualifiers["translation"][0] == H3 and cds[0].qualifiers["transl_table"] == ["1"]


@needs_tools
def test_consensus_from_fragmented_assembly(tmp_path, lib):
    """one rDNA unit split into overlapping SPAdes-like contigs (one with an error, low coverage), plus a
    contaminant 18S of another species: the consensus is the original gene, the contaminant is dropped"""
    parts = {g: (mutate(s, 0.04) if not g.startswith("ITS") else rand(len(s))) for g, s in lib["base"].items()}
    unit_seq = "".join(parts[g] for g, _ in UNIT)
    cuts = [(0, 1200), (900, 2600), (2300, 4300), (4000, 6010)]
    contigs = []
    for k, (a, b) in enumerate(cuts):
        s = unit_seq[a:b]
        contigs.append(SeqRecord(Seq(s if k % 2 else str(Seq(s).reverse_complement())),
                                 id=f"NODE_{k}_length_{b - a}_cov_{500 + k}.0", description=""))
    bad = list(unit_seq[300:1500])
    bad[100] = "A" if bad[100] != "A" else "C"     # one error, low coverage: outvoted
    contigs.append(SeqRecord(Seq("".join(bad)), id="NODE_9_length_1200_cov_3.0", description=""))
    other = mutate(parts["18S"], 0.12)
    contigs.append(SeqRecord(Seq(other), id="NODE_8_length_1800_cov_400.0", description=""))
    SeqIO.write(contigs, str(tmp_path / "in.fa"), "fasta")
    nuc.annotate(tmp_path / "in.fa", tmp_path / "out", lib["dir"], ["rdna"], 2)
    cons = {r.id: str(r.seq) for r in SeqIO.parse(str(tmp_path / "out" / "consensus" / "18S.fasta"), "fasta")}
    assert cons["18S"] == parts["18S"]
    for g in ("ITS1", "5.8S", "ITS2"):
        assert next(SeqIO.parse(str(tmp_path / "out" / "consensus" / f"{g}.fasta"), "fasta")).seq == parts[g], g
    rows = csv.DictReader(open(tmp_path / "out" / "consensus.tsv"), delimiter="\t")
    row = next(r for r in rows if r["gene"] == "18S")
    assert "NODE_8_" in row["dropped"] and row["partial_5"] == row["partial_3"] == "0"


@needs_tools
def test_assemble_joins_exact_kmer_overlap(tmp_path, lib):
    """SPAdes neighbours share a short exact overlap (k-mer); it must be joined, not repeated"""
    gene = lib["base"]["28S"]
    a, b = gene[:1500], gene[1475:]                 # 25 nt shared
    out = nuc._assemble([a, b], [100.0, 90.0], gene.lower(), tmp_path)
    assert out == gene



@needs_tools
def test_exonerate_protein2dna_gives_an_intronless_model(tmp_path):
    # protein2dna GFF has no exon lines, only one similarity line; the intronless model must still be returned
    rnd = random.Random(3)
    codon = {}
    for a in "ACGT":
        for b in "ACGT":
            for c in "ACGT":
                codon.setdefault(str(Seq(a + b + c).translate()), a + b + c)
    prot = "M" + "".join(rnd.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(79))
    flank = "".join(rnd.choice("ACGT") for _ in range(100))
    (tmp_path / "t.fa").write_text(f">t\n{flank}{''.join(codon[x] for x in prot)}TAA{flank}\n")
    m = nuc._exonerate(str(tmp_path / "t.fa"), prot, "protein2dna", tmp_path)
    assert m and [e[:2] for e in m[1]] == [(101, 100 + 3 * len(prot))]


@needs_tools
def test_curated_reference_decides_the_boundary(tmp_path):
    # four references end their 18S 12 nt early; one curated 18S record has the full gene and must win
    base = unit(None)
    refs = []
    for i in range(4):
        parts = {g: (mutate(s, 0.05) if not g.startswith("ITS") else rand(len(s))) for g, s in base.items()}
        r = rdna_record(f"REF{i}", parts, rand(300), rand(300), f"Genus{i} sp")
        f = r.features[0]
        f.location = FeatureLocation(f.location.start, f.location.end - 12, 1)
        refs.append(r)
    cur = SeqRecord(Seq(mutate(base["18S"], 0.05)), id="NR_1", name="NR_1", description="")
    cur.annotations.update(molecule_type="DNA", organism="Curata sp", taxonomy=["Eukaryota", "Metazoa"])
    cur.features = [SeqFeature(FeatureLocation(0, len(cur.seq), 1), type="rRNA",
                               qualifiers={"product": ["18S ribosomal RNA"]})]
    SeqIO.write(refs, str(tmp_path / "refs.gb"), "genbank")
    SeqIO.write([cur], str(tmp_path / "cur.gb"), "genbank")
    nuc.build([tmp_path / "refs.gb"], tmp_path / "lib", curated=[tmp_path / "cur.gb"])
    parts = {g: (mutate(s, 0.04) if not g.startswith("ITS") else rand(len(s))) for g, s in base.items()}
    q = rdna_record("Q9", parts, rand(500), rand(500), "Query sp")
    SeqIO.write([SeqRecord(q.seq, id="c9", description="")], str(tmp_path / "in.fa"), "fasta")
    nuc.annotate(tmp_path / "in.fa", tmp_path / "out", tmp_path / "lib", ["rdna"], 2)
    row = next(r for r in read_rows(tmp_path / "out") if r["gene"] == "18S")
    assert int(row["end"]) == 500 + len(base["18S"]) and row["evidence"].endswith("curated")
