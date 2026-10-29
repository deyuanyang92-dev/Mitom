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
    exons, p5, p3, _ = nuc._orf(seq, [(7, 48)], 1, False)
    assert exons == [(4, 4 + len(cds) - 1)] and not p5 and not p3
    rc = seq.reverse_complement()
    n = len(rc)
    exons, p5, p3, _ = nuc._orf(rc, [(n - 48 + 1, n - 7 + 1)], -1, False)
    assert exons == [(n - (4 + len(cds) - 1) + 1, n - 4 + 1)] and not p5 and not p3


def test_orf_open_at_contig_end_is_partial():
    seq = Seq("GCT" * 10)                        # no ATG, no stop
    exons, p5, p3, _ = nuc._orf(seq, [(4, 27)], 1, False)
    assert exons == [(1, 30)] and p5 and p3


def test_orf_minus_strand_partial_ends_and_codon_start():
    """AM921654.1 / AM921657.1 (table2asn StartCodon, NoStop): on the minus strand the 5' / 3' partial flags were
    swapped; AB429219.1: a partial CDS <2..>310 of a 312 nt sequence must run to the ends (<1..>312, codon_start 2)"""
    cds = "ATG" + "GCT" * 20 + "TAA"
    seq = Seq("CCC" + cds[:40])                  # gene ends at the sequence end: 3' partial
    rc = seq.reverse_complement()
    n = len(rc)
    exons, p5, p3, cs = nuc._orf(rc, [(1, n - 6)], -1, False)
    assert exons == [(1, n - 3)] and not p5 and p3 and cs == 1
    seq = Seq("C" + "GCT" * 20 + "TAAGG")        # no ATG, 1 nt before the first full codon
    exons, p5, p3, cs = nuc._orf(seq, [(5, 40)], 1, False)
    assert exons == [(1, 64)] and p5 and not p3 and cs == 2


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
    # the synthetic unit is random sequence, so the Rfam models cannot confirm the 18S start that sits at the first
    # contig's end: it is reported 5' partial (an unconfirmed end at a sequence end is never called complete)
    assert "NODE_8_" in row["dropped"] and (row["partial_5"], row["partial_3"]) == ("1", "0")


@needs_tools
def test_assemble_joins_exact_kmer_overlap(tmp_path, lib):
    """SPAdes neighbours share a short exact overlap (k-mer); it must be joined, not repeated"""
    gene = lib["base"]["28S"]
    a, b = gene[:1500], gene[1475:]                 # 25 nt shared
    out = nuc._assemble([a, b], [100.0, 90.0], gene.lower(), tmp_path)
    assert out == gene


@needs_tools
def test_assemble_joins_the_better_covered_branch(tmp_path, lib):
    """SRR27947587: copy variants split the 28S into contigs that overlap by k, with two branches over a variable
    stretch; the weaker branch was joined to one side, the better one to the other, and the gene stayed in two
    pieces (1097 nt consensus of a 3693 nt 28S)"""
    gene, k = lib["base"]["28S"], 77
    branch = gene[1500 - k:1520 + k]
    other = gene[1500 - k:1500] + rand(20) + gene[1520:1520 + k]
    out = nuc._assemble([gene[:1500], gene[1520:], branch, other], [100.0, 90.0, 60.0, 30.0], gene.lower(), tmp_path)
    assert out == gene


@needs_tools
def test_consensus_across_scaffold_gaps(tmp_path, lib):
    """SRR32909333 / SRR27947580: a SPAdes scaffold joins its contigs with runs of N, once over their k-mer overlap
    (the 77 nt then stand on both sides of the N), once over 100 nt that a short contig holds; the N were taken
    for sequence and the 28S consensus had 1653 N or lost 1444 nt"""
    gene, k = mutate(lib["base"]["28S"], 0.04), 77
    a, b = 1200, 2400
    recs = {"NODE_1_length_3577_cov_500.0": gene[:a] + "N" * 100 + gene[a - k:b] + "N" * 100 + gene[b + 100:],
            "NODE_2_length_254_cov_300.0": gene[b - k:b + 100 + k]}
    rows = [{"contig": c, "gene": "28S", "exons": f"1-{len(s)}", "strand": "+", "contig_cov": nuc._cov(c),
             "length": str(len(s)), "best_ref": "REF0", "partial_5": "1" if i else "0", "partial_3": "1" if i else "0"}
            for i, (c, s) in enumerate(recs.items())]
    recs = {c: SeqRecord(Seq(s), id=c) for c, s in recs.items()}
    cons = nuc.consensus(rows, recs, nuc.Library(lib["dir"], tmp_path / "lib"), tmp_path, 0.99)
    assert (cons[0]["seq"], cons[0]["n_fragments"]) == (gene, 2)


@needs_tools
def test_consensus_spacer_from_its_complete_copy(tmp_path, lib):
    """SRR27947585: the best-covered contig ends 27 nt into ITS1, a less covered one holds all of it; a spacer has no
    reference to order its pieces on and the 27 nt piece was the ITS1 consensus"""
    parts = {g: (mutate(s, 0.04) if not g.startswith("ITS") else rand(len(s))) for g, s in lib["base"].items()}
    u = "".join(parts[g] for g, _ in UNIT)
    i1 = len(parts["18S"])
    i2 = i1 + len(parts["ITS1"]) + len(parts["5.8S"])
    SeqIO.write([SeqRecord(Seq(u[:i1 + 27]), id=f"NODE_1_length_{i1 + 27}_cov_44.0", description=""),
                 SeqRecord(Seq(u[i1 - 100:i2 + 91]), id=f"NODE_2_length_{i2 - i1 + 191}_cov_31.0", description="")],
                str(tmp_path / "in.fa"), "fasta")
    nuc.annotate(tmp_path / "in.fa", tmp_path / "out", lib["dir"], ["rdna"], 2)
    assert next(SeqIO.parse(str(tmp_path / "out" / "consensus" / "ITS1.fasta"), "fasta")).seq == parts["ITS1"]


@needs_tools
def test_consensus_keeps_one_of_two_rdna_types(tmp_path, lib):
    """SRR27947581: two complete 18S 97 % alike at similar coverage (two rDNA types in one sample) and a short contig
    of one type; the votes of both types were mixed and the 18S consensus was a chimera of them"""
    parts = {g: (mutate(s, 0.04) if not g.startswith("ITS") else rand(len(s))) for g, s in lib["base"].items()}
    other = dict(parts, **{"18S": "".join(c if i % 40 else "ACGT"[("ACGT".index(c) + 1) % 4]
                                          for i, c in enumerate(parts["18S"]))})     # 2.5 % of sites differ
    ua, ub = ("".join(p[g] for g, _ in UNIT) for p in (parts, other))
    flank = rand(300)
    SeqIO.write([SeqRecord(Seq(flank + ua + flank), id="NODE_1_length_6610_cov_195.0", description=""),
                 SeqRecord(Seq(flank + ub + flank), id="NODE_2_length_6610_cov_164.0", description=""),
                 SeqRecord(Seq(ub[200:900]), id="NODE_3_length_700_cov_84.0", description="")],
                str(tmp_path / "in.fa"), "fasta")
    nuc.annotate(tmp_path / "in.fa", tmp_path / "out", lib["dir"], ["rdna"], 2)
    s18 = str(next(SeqIO.parse(str(tmp_path / "out" / "consensus" / "18S.fasta"), "fasta")).seq)
    assert s18 in (parts["18S"], other["18S"])



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


def _fake_lib(rdna):
    from types import SimpleNamespace
    return SimpleNamespace(rdna={r: {g: (a, b, True, True) for g, (a, b) in gs.items()} for r, gs in rdna.items()},
                           org={r: (f"{r} sp", r) for r in rdna}, curated=set())


def _hsp(qs, qe, ss, se, bits):
    n = qe - qs + 1
    return (qs, qe, ss, se, float(bits), "A" * n, "A" * n)


def test_rdna_gene_needs_alignment_inside_it():
    """23007-1 histone contig: 132 nt at the contig end hit the spacer after a sponge 28S (and the inside of a
    partial 28S annotated longer); 18S and 28S were projected kilobases away over the histone genes"""
    lib = _fake_lib({"A": {"18S": (1, 1800), "5.8S": (2150, 2310), "28S": (2610, 6010)},
                     "B": {"18S": (1, 1800), "5.8S": (2150, 2310), "28S": (2610, 6300)}})
    lib.rdna["B"]["28S"] = (2610, 6300, False, False)     # Pholoe KY753843: 28S ends marked unusable
    hits = {"A": [_hsp(6268, 6399, 6100, 6231, 176)], "B": [_hsp(6268, 6399, 6100, 6231, 150)]}
    assert nuc.rdna_loci(hits, lib, 6399) == []


def test_rdna_one_gene_per_place_from_two_loci():
    """23007-2: a reference ordered 28S - spacer - 18S chained the previous unit's 28S tail and the 18S into its
    own locus beside the unit's locus; the 18S was reported twice (4785 and 4789)"""
    lib = _fake_lib({"X": {"28S": (1, 3700), "18S": (4900, 6700), "5.8S": (7400, 7560)},
                     "P": {"18S": (1800, 3600), "5.8S": (4400, 4560), "28S": (4900, 8600)}})
    hits = {"X": [_hsp(1, 750, 2951, 3700, 3000), _hsp(4785, 6634, 4900, 6749, 400)],
            "P": [_hsp(4789, 10821, 1800, 7832, 3000)]}
    loci = nuc.rdna_loci(hits, lib, 10821)
    s18 = [(L["ref"], L["feats"]["18S"]["lo"]) for L in loci if "18S" in L["feats"]]
    assert s18 == [("P", 4789)]
    assert any(L["ref"] == "X" and L["feats"].get("28S", {}).get("hi") == 750 for L in loci)
    assert {"ITS1", "5.8S", "ITS2", "28S"} <= set(next(L for L in loci if L["ref"] == "P")["feats"])


def test_read_input_skips_records_without_sequence(tmp_path):
    """histone held-out set: NG_032568.1 (a CON record, no ORIGIN) crashed read_input with UndefinedSequenceError"""
    gb = tmp_path / "in.gb"
    gb.write_text("LOCUS       CON1                    50 bp    DNA     linear   CON 01-JAN-2020\n"
                  "DEFINITION  test.\nACCESSION   CON1\nVERSION     CON1.1\nFEATURES             Location/Qualifiers\n"
                  "     source          1..50\nCONTIG      join(gap(50))\n//\n")
    SeqIO.write([SeqRecord(Seq("ACGT" * 30), id="R2", name="R2", description="",
                           annotations={"molecule_type": "DNA"})], str(tmp_path / "r2.gb"), "genbank")
    gb.write_text(gb.read_text() + (tmp_path / "r2.gb").read_text())
    fa, org, n, nt, circ = nuc.read_input(gb, tmp_path)
    assert (n, nt, circ) == (1, 120, set()) and [r.id for r in SeqIO.parse(str(fa), "fasta")] == ["R2"]


def test_cut_circular_in_the_longest_stretch_without_genes():
    assert nuc._cut(100, [(10, 30), (50, 60)]) == 84       # free 61..100 + 1..9 (49 nt): its middle
    assert nuc._cut(100, [(90, 110)]) == 49                # a gene across the origin
    assert nuc._cut(100, [(40, 60), (140, 160)]) == 99     # one gene seen in two turns covers the same nt
    assert nuc._cut(100, []) is None and nuc._cut(100, [(1, 100)]) is None and nuc._cut(100, [(5, 150)]) is None


@needs_tools
def test_circular_input_turned_to_start_outside_the_genes(tmp_path, lib):
    """GetOrganelle '(circular)' rDNA unit: the path starts anywhere; SRR27947576 started inside the 28S, which came
    out 3' partial (3625 of 3694 nt) with its last 69 nt at the other end left unannotated"""
    parts = {g: (mutate(s, 0.04) if not g.startswith("ITS") else rand(len(s))) for g, s in lib["base"].items()}
    unit = "".join(parts[g] for g, _ in UNIT)
    circle = unit + rand(1500)
    seq = circle[len(unit) - 300:] + circle[:len(unit) - 300]       # starts 300 nt before the end of the 28S
    for name, d in (("EDGE_1+(circular)", "circ"), ("scaffold_1--5+", "lin")):
        (tmp_path / "in.fa").write_text(f">{name}\n{seq}\n")
        nuc.annotate(tmp_path / "in.fa", tmp_path / d, lib["dir"], ["rdna"], 2)
        rows = read_rows(tmp_path / d)
        s28 = [r for r in rows if r["gene"] == "28S"]
        if d == "circ":
            assert [(r["length"], r["partial_5"], r["partial_3"]) for r in s28] == [(str(len(parts["28S"])), "0", "0")]
            assert len(rows) == 5 and all(r["partial_5"] == r["partial_3"] == "0" for r in rows)
            assert next(SeqIO.parse(str(tmp_path / d / "consensus" / "28S.fasta"), "fasta")).seq == parts["28S"]
        else:                                                       # a linear sequence is never turned
            assert all("1" in (r["partial_5"], r["partial_3"]) for r in s28)


@needs_tools
def test_circular_genbank_record_turned(tmp_path, lib):
    parts = {g: (mutate(s, 0.04) if not g.startswith("ITS") else rand(len(s))) for g, s in lib["base"].items()}
    circle = "".join(parts[g] for g, _ in UNIT) + rand(1500)
    r = SeqRecord(Seq(circle[1000:] + circle[:1000]), id="C1", name="C1", description="rDNA unit",   # starts in the 18S
                  annotations={"molecule_type": "DNA", "topology": "circular", "organism": "Query sp"})
    SeqIO.write([r], str(tmp_path / "in.gb"), "genbank")
    nuc.annotate(tmp_path / "in.gb", tmp_path / "out", lib["dir"], ["rdna"], 2)
    s18 = [x for x in read_rows(tmp_path / "out") if x["gene"] == "18S"]
    assert [(x["length"], x["partial_5"], x["partial_3"]) for x in s18] == [(str(len(parts["18S"])), "0", "0")]


@needs_tools
def test_consensus_per_organism_for_genbank_records(tmp_path, lib):
    """a batch of GenBank records of several species: one consensus per organism, not one chimera of all"""
    recs = []
    for i in range(2):
        parts = {g: (mutate(s, 0.04) if not g.startswith("ITS") else rand(len(s))) for g, s in lib["base"].items()}
        recs.append(rdna_record(f"S{i}", parts, rand(500), rand(500), f"Species{i} sp"))
    SeqIO.write(recs, str(tmp_path / "in.gb"), "genbank")
    nuc.annotate(tmp_path / "in.gb", tmp_path / "out", lib["dir"], ["rdna"], 2)
    cons = list(csv.DictReader(open(tmp_path / "out" / "consensus.tsv"), delimiter="\t"))
    assert sorted({c["sample"] for c in cons}) == ["Species0 sp", "Species1 sp"]
    s18 = {r.id: str(r.seq) for r in SeqIO.parse(str(tmp_path / "out" / "consensus" / "18S.fasta"), "fasta")}
    assert s18 == {f"Species{i}_sp_18S": str(recs[i].features[0].extract(recs[i].seq)).upper() for i in range(2)}


def _feat(lo, hi, how, aln, sup):
    return {"lo": lo, "hi": hi, "p_lo": False, "p_hi": False, "how": how, "aln": set(aln), "sup": sup, "aligned": len(aln)}


def test_tidy_projected_end_stops_at_scaffold_gap():
    """AB272236.1: 100 N stand for a ~1.3 kb gap; the 28S 5' end was projected through it onto the 18S
    (NCBI: <2059..>4631). LC713089.1 / BK061357.1: a gap inside an aligned gene stays inside it"""
    seq = "A" * 1958 + "N" * 100 + "A" * 2573
    L = {"st": 1, "feats": {"28S": _feat(1249, 4631, "projected/projected", range(2100, 4100), 2000)}}
    nuc._tidy([L], seq, len(seq))
    f = L["feats"]["28S"]
    assert (f["lo"], f["hi"], f["p_lo"]) == (2059, 4631, True)
    L = {"st": 1, "feats": {"18S": _feat(1, 1646, "projected/aligned", list(range(1, 1400)) + list(range(1520, 1646)), 900)}}
    seq = "A" * 1482 + "N" * 26 + "A" * 4379
    nuc._tidy([L], seq, len(seq))
    assert (L["feats"]["18S"]["lo"], L["feats"]["18S"]["hi"]) == (1, 1646)


def test_tidy_drops_rrna_overlapping_another_on_either_strand():
    """MG590305.1: a reverse 28S projected over the end of the 18S into the 3.2 kb ITS1 (169 nt aligned);
    PP357355.1: a 28S projected from the IGS over the 18S start"""
    seq = "A" * 11000
    a = {"st": 1, "feats": {"18S": _feat(642, 2464, "curated/curated", range(642, 2465), 3000)}}
    b = {"st": -1, "feats": {"28S": _feat(2302, 4824, "projected/projected", range(3000, 3169), 150)}}
    c = {"st": 1, "feats": {"28S": _feat(1, 700, "projected/projected", range(600, 700), 90)}}
    nuc._tidy([a, b, c], seq, len(seq))
    assert "18S" in a["feats"] and not b["feats"] and not c["feats"]


def test_tidy_end_inside_gap_moves_to_its_edge():
    """PQ211133.1: the 28S start was projected into a 100 N gap (table2asn FeatureBeginsOrEndsInGap)"""
    seq = "A" * 2032 + "N" * 100 + "A" * 1035
    L = {"st": 1, "feats": {"28S": _feat(2058, 3167, "projected/aligned", range(2200, 3100), 900)}}
    nuc._tidy([L], seq, len(seq))
    assert (L["feats"]["28S"]["lo"], L["feats"]["28S"]["p_lo"]) == (2133, True)


def test_tidy_final_drops_rrna_placed_by_nothing():
    """MJ20260820041 CHY007 (GetOrganelle rDNA, minus strand): the contig ends 20 nt into ITS1; an 18S was projected
    there from the 5.8S across the reference's shorter ITS1 (0 nt aligned, no model hit) and cut the ITS1 to a
    'complete' 499 nt. In the final pass it goes and the ITS1 runs 5'-partial to the contig end (519 nt)."""
    seq = "A" * 4669
    def locus():
        f58 = _feat(3993, 4150, "hmm/hmm", range(3994, 4150), 300)
        s18 = dict(_feat(4650, 4669, "projected/projected", [], 0), p_hi=True)
        return {"st": -1, "feats": {"5.8S": f58, "18S": s18}}
    L = locus()
    nuc._tidy([L], seq, len(seq))                     # first pass: the Rfam models may still place it
    assert "18S" in L["feats"]
    L = locus()
    nuc._tidy([L], seq, len(seq), final=True)
    assert "18S" not in L["feats"]
    its1 = L["feats"]["ITS1"]
    assert (its1["lo"], its1["hi"], its1["p_hi"]) == (4151, 4669, True)   # minus strand: 5' end = high end
    L = locus()
    L["feats"]["18S"]["aligned"] = 15                 # some alignment inside: a real piece, kept
    nuc._tidy([L], seq, len(seq), final=True)
    assert "18S" in L["feats"]
    end18 = "GTAGGTGAACCTGCGGAAGGATCATTA"[-20:]       # the piece is the reference 18S 3' end: kept
    s2 = "A" * 4649 + str(Seq(end18).reverse_complement())
    L = locus()
    L["ref"] = "R1"
    nuc._tidy([L], s2, len(s2), final=True, ref_gene=lambda r, g: "C" * 1797 + end18 if g == "18S" else None)
    assert "18S" in L["feats"]
    L = locus()
    L["ref"] = "R1"                                   # CHY007's piece (ITS1 sequence) against the same reference
    s3 = "A" * 4649 + str(Seq("GAGAGCGATCGAGTCGGCGG").reverse_complement())
    nuc._tidy([L], s3, len(s3), final=True, ref_gene=lambda r, g: "C" * 1797 + end18 if g == "18S" else None)
    assert "18S" not in L["feats"]


def test_definition_line_follows_ncbi_grouping():
    """NCBI wording, e.g. KY753840.1 / MW412252.1 / PP239203.1 titles"""
    def f(g, a, b, p5=False, p3=False):
        loc = FeatureLocation(nuc.BeforePosition(a) if p5 else a, nuc.AfterPosition(b) if p3 else b, 1)
        return SeqFeature(loc, type="rRNA" if g in nuc.RRNA else "misc_RNA", qualifiers={"product": [nuc.PRODUCT[g]]})
    unit = [f("18S", 0, 1800, p5=True), f("ITS1", 1800, 2100), f("5.8S", 2100, 2260), f("ITS2", 2260, 2500),
            f("28S", 2500, 4000, p3=True)]
    assert nuc.definition("Perinereis wilsoni", {"specimen_voucher": "23007-1"}, unit) == (
        "Perinereis wilsoni voucher 23007-1 18S ribosomal RNA gene, partial sequence; internal transcribed spacer 1, "
        "5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence; and 28S ribosomal RNA gene, "
        "partial sequence")
    assert nuc.definition("X y", {}, [f("ITS1", 0, 300, p5=True), f("5.8S", 300, 460), f("ITS2", 460, 700)]) == (
        "X y internal transcribed spacer 1, partial sequence; and 5.8S ribosomal RNA gene and internal transcribed "
        "spacer 2, complete sequence")      # two groups: '; and' in 94 of 103 NCBI rDNA titles (2026-10-25)


def test_tidy_spacer_not_over_rrna_or_in_gap():
    """LC799786.1: ITS2 drawn inside a 100 N gap; FN995994.2: an ITS1 of one copy ran over the rRNAs of another"""
    seq = "A" * 1937 + "N" * 100 + "A" * 3000
    L = {"st": 1, "feats": {"5.8S": _feat(1780, 1938, "aligned/aligned", range(1780, 1939), 300),
                            "28S": _feat(2038, 5000, "projected/aligned", range(2100, 5000), 3000)}}
    nuc._tidy([L], seq, len(seq))
    assert "ITS2" not in L["feats"]
    seq = "A" * 3400
    a = {"st": 1, "feats": {"18S": _feat(1, 700, "projected/aligned", range(400, 700), 500),
                            "5.8S": _feat(971, 1127, "aligned/aligned", range(971, 1128), 300)}}
    b = {"st": 1, "feats": {"5.8S": _feat(2658, 2814, "aligned/aligned", range(2658, 2815), 290)}}
    nuc._tidy([a, b], seq, len(seq))
    assert (a["feats"]["ITS1"]["lo"], a["feats"]["ITS1"]["hi"]) == (701, 970)
    assert "ITS1" not in b["feats"]          # would run from 1 over the 18S and the other 5.8S


CM = all(mito_env.find(t) for t in ("cmsearch",)) and (nuc.LIB / "rfam" / "RF01960.cm").exists()


def _curated(acc, gene):
    import gzip
    row = next(r for r in csv.DictReader(open(nuc.LIB / "rdna_ref.tsv"), delimiter="\t")
               if r["ref"].startswith(acc) and r["gene"] == gene)
    with gzip.open(nuc.LIB / "rdna_ref.fa.gz", "rt") as fh:
        s = next(str(r.seq) for r in SeqIO.parse(fh, "fasta") if r.id == row["ref"])
    return s[int(row["start"]) - 1:int(row["end"])]


@pytest.mark.skipif(not CM, reason="cmsearch / Rfam models not found")
def test_cm_ends_place_the_curated_mature_ends(tmp_path):
    """human 18S (NR_003286) in random flanks; predicted ends 3 nt off (as the hesionid convention of polynoid
    records, ...GATCATTACCG) are moved to the mature ends, on both strands"""
    g = _curated("NR_003286", "18S")
    seq = Seq(rand(700) + g + rand(700))
    a, b = 701, 700 + len(g)
    seqs = {"c": seq, "r": seq.reverse_complement()}
    n = len(seq)
    jobs = [("c", 1, 5, a - 3), ("c", 1, 3, b + 3), ("r", -1, 5, n - a + 1 + 2), ("r", -1, 3, n - b + 1 - 2)]
    out = nuc.cm_ends(jobs, seqs, "18S", nuc.LIB / "rfam" / "RF01960.cm", tmp_path, 2)
    assert out == {0: a, 1: b, 2: n - a + 1, 3: n - b + 1}


@pytest.mark.skipif(not CM, reason="cmsearch / Rfam models not found")
def test_cm_rescue_58s_piece_at_sequence_start(tmp_path):
    """JN636100.1 / JN636101.1 (ITS2 amplicons): the sequence starts in the 5.8S; the 3' 60 nt of it were not
    annotated (too short for blastn)"""
    s58 = _curated("NR_003285", "5.8S")
    seq = Seq(s58[-60:] + rand(400) + _curated("NR_003287", "28S")[:600])
    L = {"st": 1, "feats": {"28S": {"lo": 461, "hi": len(seq), "p_lo": False, "p_hi": True, "how": "aligned/projected",
                                    "aln": set(range(461, 1061)), "sup": 900, "aligned": 600}}}
    assert nuc._cm_rescue_58s({"c": [L]}, {"c": seq}, nuc.LIB / "rfam", tmp_path, 2) == 1
    f = L["feats"]["5.8S"]
    assert (f["lo"], f["hi"], f["p_lo"]) == (1, 60, True) and L["feats"]["ITS2"]["lo"] == 61


@pytest.mark.skipif(not CM, reason="cmsearch / Rfam models not found")
def test_cm_refine_hmm_first_then_partial_when_unconfirmed(tmp_path):
    """the profile HMM places ends that are a few nt off (label hmm); an end at the sequence end that no model can
    confirm is partial (OK209792.1: the record starts inside the 18S, the projection said complete at nt 1)"""
    g = _curated("NR_003286", "18S")
    seq = Seq(rand(400) + g + rand(400))
    f = {"lo": 404, "hi": 400 + len(g) - 2, "p_lo": False, "p_hi": False, "how": "aligned/aligned"}
    L = {"st": 1, "feats": {"18S": f}}
    assert nuc.cm_refine({"c": [L]}, {"c": seq}, nuc.LIB, tmp_path, 2) == 2
    assert (f["lo"], f["hi"], f["how"]) == (401, 400 + len(g), "hmm/hmm")
    cut = Seq(g[25:] + rand(400))                  # starts 25 nt inside the 18S
    f = {"lo": 1, "hi": len(g) - 25, "p_lo": False, "p_hi": False, "how": "aligned/aligned"}
    L = {"st": 1, "feats": {"18S": f}}
    nuc.cm_refine({"c": [L]}, {"c": cut}, nuc.LIB, tmp_path, 2)
    assert f["p_lo"] and f["lo"] == 1 and not f["p_hi"]


def _h3(sites):
    """H3 test protein with residues 31/87/89/90 set (numbered without Met)"""
    p = list(H3)
    for pos, aa in zip(nuc.H3_SITES, sites):
        p[pos] = aa                               # index pos = residue pos (H3 starts with Met)
    return "".join(p)


def test_histone_class_h3_signatures():
    """MJ20260820041: the H3 consensus was H3.3 (S31 A87 I89 G90) in 3 annelid assemblies, canonical in 2"""
    assert nuc.histone_class("H3", _h3("ASVM"), 1, None)[0] == "replication-dependent"
    cls, why = nuc.histone_class("H3", _h3("SAIG"), 1, None)
    assert cls == "replacement" and "H3.3" in why
    assert nuc.histone_class("H3", _h3("SAIM"), 1, None)[0] == "replacement"       # 3 of 4 sites decide
    assert nuc.histone_class("H3", _h3("ASIG"), 1, None)[0] == ""                    # 2 vs 2: undecided
    assert nuc.histone_class("H3", H3[:60], 1, None)[0] == ""                       # fragment: sites not reached
    assert nuc.histone_class("H3", H3, 1, None)[0] == "replication-dependent"       # test protein is canonical


def test_histone_class_introns_copy_ratio_h2ax():
    assert nuc.histone_class("H4", "MSGRGKGG", 2, 300.0)[0] == "replacement"        # spliced: never a cluster gene
    assert nuc.histone_class("H4", "MSGRGKGG", 1, 300.0)[0] == "replication-dependent"
    assert nuc.histone_class("H4", "MSGRGKGG", 1, 1.2)[0] == "single-copy"
    assert nuc.histone_class("H4", "MSGRGKGG", 1, 3.0)[0] == ""                     # between: undecided
    assert nuc.histone_class("H4", "MSGRGKGG", 1, None)[0] == ""                    # no coverage: undecided
    assert nuc.histone_class("H3", _h3("SAIG"), 1, 300.0)[0] == "replacement"       # residues outrank coverage
    assert "H2A.X" in nuc.histone_class("H2A", "MSGRGKTGGKARAKAKSRSSRAGLQFPVGRVHRLLRKGNYSQEY", 1, None)[1]


def test_pick_histone_prefers_cluster_copy():
    """a complete single-copy H3.3 (14x) must not stand for H3 when cluster copies (4000x) are fragments"""
    rows = [{"contig": "h33", "contig_cov": "14", "partial_5": "0", "partial_3": "0", "length": "411",
             "class": "replacement"},
            {"contig": "frag", "contig_cov": "4000", "partial_5": "1", "partial_3": "0", "length": "300",
             "class": "replication-dependent"},
            {"contig": "unk", "contig_cov": "", "partial_5": "0", "partial_3": "0", "length": "411", "class": ""}]
    assert nuc.pick_histone(rows, 411)["contig"] == "frag"
    assert nuc.pick_histone(rows[::2], 411)["contig"] == "unk"     # undecided before replacement
    assert nuc.pick_histone(rows[:1], 411)["contig"] == "h33"      # only a variant: still reported
    h2ax = {"contig": "h2ax", "contig_cov": "15", "partial_5": "0", "partial_3": "0", "length": "402",
            "class": "single-copy"}
    assert nuc.pick_histone([h2ax, rows[0]], 375)["contig"] == "h2ax"   # single-copy before replacement
    assert nuc.pick_histone([h2ax, rows[1]], 375)["contig"] == "frag"   # cluster fragment before single copy


def test_single_copy_cov_from_spades_names():
    ids = ["NODE_1_length_50000_cov_12.0", "NODE_2_length_40000_cov_14.0", "NODE_3_length_2000_cov_4000.0",
           "NODE_4_length_500_cov_1.0"]
    assert nuc.single_copy_cov(ids) == 12.0
    assert nuc.single_copy_cov(["scaffold_1", "contig2"]) is None
