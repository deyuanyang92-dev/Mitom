"""Tests for mito_intron_ref.py (intron-site transfer).

Records (tests/data): Metridium senile NC_000933 (cox1 and nad5 group I introns), Exaiptasia diaphana NC_056771 (the
same two sites), Plakinastrella NC_010217 (cox1 with a 9-nt mini-exon; intron feature 6 nt off its CDS gap), Montipora
NC_006902 (11.5-kb nad5 intron hosting genes; intron feature 1 nt off its CDS gap), Platynereis NC_000931 (no
intron, negative control). Queries are the bare sequences rotated as the pipeline rotates them; predicted introns are
compared with NCBI's by sequence and flanks (rotation-independent).

Integration tests need tblastn and blastn (on PATH, or MITO_INTRON_BIN=/path/to/env/bin).
Run: python3 -m pytest -q tests/test_mito_intron_ref.py
"""
import os
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_env  # noqa: E402
import mito_intron_ref as M  # noqa: E402

DATA = HERE / "data"
SCRIPT = HERE.parent / "mito_intron_ref.py"
BIN = os.environ.get("MITO_INTRON_BIN") or mito_env.bin_with("tblastn", "blastn") or ""
needs_blast = pytest.mark.skipif(not BIN, reason="tblastn / blastn not available")


def ncbi(acc):
    return SeqIO.read(str(DATA / f"{acc}.gb"), "genbank")


def bare(rec, shift=0, revcomp=False):
    """Query as the pipeline step sees it: the record's sequence only, rotated by shift (and reverse-complemented)."""
    s = str(rec.seq)
    s = s[shift:] + s[:shift]
    if revcomp:
        s = str(Seq(s).reverse_complement())
    q = SeqRecord(Seq(s), id=rec.id, name=rec.name, description="query")
    q.annotations = {"molecule_type": "DNA", "topology": rec.annotations.get("topology", "circular")}
    return q


def intron_keys(rec, introns=None):
    """{(intron sequence, 10 nt before, 10 nt after)} on the transcript strand, from intron features (default) or
    from the gaps of joined CDS when the record leaves an intron to its join."""
    seq = str(rec.seq)
    L = len(seq)
    out = set()
    feats = [f for f in rec.features if f.type == "intron"] if introns is None else introns
    for f in feats:
        s, e, st = int(f.location.start), int(f.location.end), f.location.strand or 1
        if len(f.location.parts) > 1:                   # written across the origin: join(x..L,1..y)
            p = f.location.parts
            s, e = (int(p[0].start), int(p[-1].end) + L) if st == 1 else (int(p[-1].start), int(p[0].end) + L)
        ext = "".join(seq[i % L] for i in range(s - 10, e + 10))
        if st == -1:
            ext = str(Seq(ext).reverse_complement())
        out.add((ext[10:-10], ext[:10], ext[-10:]))
    return out


def join_gap_keys(rec, gene):
    """Intron keys of the gaps of the joined CDS of one gene (introns crossing the origin are given by the join)."""
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    L = len(rec.seq)
    gaps = []
    for f in rec.features:
        if f.type != "CDS" or M.gene_name(f.qualifiers, "CDS") != gene:
            continue
        for p, q in zip(f.location.parts, f.location.parts[1:]):
            if p.strand != q.strand:
                continue
            if (p.strand or 1) == 1:
                a, b = int(p.end), int(q.start) + (L if int(q.start) < int(p.end) else 0)
            else:
                a, b = int(q.end), int(p.start) + (L if int(p.start) < int(q.end) else 0)
            if b - a >= 100:
                if b <= L:
                    gaps.append(SeqFeature(FeatureLocation(a, b, strand=p.strand), type="intron"))
                else:
                    gaps.append(SeqFeature(FeatureLocation(a, L, strand=p.strand) + FeatureLocation(0, b - L, strand=p.strand),
                                           type="intron"))
    return intron_keys(rec, gaps)


def annotate(query, lib, code, own=None, lineage=None):
    tools = M.find_tools(BIN)
    models, rows = M.annotate_record(query, lib, code, lineage=lineage, tools=tools, own=own)
    written = M.apply_models(query, models, code)
    return models, rows, written


# ── no external programs ─────────────────────────────────────────────────────
def test_gene_names_and_intron_groups():
    assert M.gene_name({"gene": ["COX1a"]}, "CDS") == "cox1"
    assert M.gene_name({"gene": ["ND4L"]}, "CDS") == "nad4l"
    assert M.gene_name({"product": ["NADH dehydrogenase subunit 5"]}, "CDS") == "nad5"
    assert M.gene_name({"product": ["16S ribosomal RNA"]}, "rRNA") == "rrnL"
    assert M.gene_name({"gene": ["orf314"]}, "CDS") == ""
    assert M.intron_kind("mintron 1; group II intron") == "II"
    assert M.intron_kind("Group IB (predicted)") == "I"
    assert M.intron_kind("group 1 intron") == "I"
    assert M.intron_kind("Group II (domain V)") == "II"
    assert M.intron_kind("intron") == ""


def test_sites_of_metridium():
    rec = ncbi("NC_000933")
    genes = {g.gene: g for g in M.ref_genes(rec)}
    assert set(genes) == {"cox1", "nad5"}
    cox1 = genes["cox1"]
    j = cox1.junctions[0]
    assert (j.splice, j.length, j.offset, j.d5, j.d3) == ("cis", 853, 893, 0, 0)
    assert (j.offset // 3, j.offset % 3) == (297, 2)              # residue / codon phase of the insertion site
    assert len(j.e5) == len(j.i5) == len(j.i3) == len(j.e3) == M.PROBE
    s = str(rec.seq)
    assert j.e5 + j.i5 == s[4155 - 30:4155 + 30]                  # intron 4156..5008: probes centred on its ends
    assert j.i3 + j.e3 == s[5008 - 30:5008 + 30]
    feat = next(f for f in rec.features if f.type == "CDS" and f.qualifiers.get("gene") == ["COX1"])
    assert cox1.protein == feat.qualifiers["translation"][0]
    assert genes["nad5"].junctions[0].length == 1681


def test_annotation_quirks_are_recorded():
    """Intron features that differ from the CDS gap are kept as offsets (own-annotation mode reproduces them)."""
    g = M.ref_genes(ncbi("NC_006902"))[0]                          # Montipora nad5: intron 1 nt after its CDS gap
    assert (g.gene, g.junctions[0].length, g.junctions[0].d5, g.junctions[0].d3) == ("nad5", 11485, 1, 1)
    g = M.ref_genes(ncbi("NC_010217"))[0]                          # Plakinastrella cox1: 9-nt middle exon
    assert g.exon_lens == [717, 9, 843]
    assert [(j.length, j.d5, j.d3) for j in g.junctions] == [(388, 6, 0), (1118, 0, 0)]


def test_library_round_trip(tmp_path):
    genes = [g for a in ("NC_000933", "NC_010217", "NC_006902") for g in M.ref_genes(ncbi(a))]
    j = genes[0].junctions[0]
    j.half5, j.gap5, j.half3, j.gap3 = 159, 1, 7, 2                # trans-splicing fields survive the file
    path = tmp_path / "sites.tsv.gz"
    M.save_library(genes, path)
    back = M.load_library(path)
    assert [g.gid for g in back] == [g.gid for g in genes]
    for a, b in zip(genes, back):
        assert (a.protein, a.spliced, a.exon_lens, a.lineage) == (b.protein, b.spliced, b.exon_lens, b.lineage)
        for x, y in zip(a.junctions, b.junctions):
            assert (x.offset, x.length, x.e5, x.i5, x.i3, x.e3, x.d5, x.d3, x.half5, x.gap5, x.half3, x.gap3) == \
                (y.offset, y.length, y.e5, y.i5, y.i3, y.e3, y.d5, y.d3, y.half5, y.gap5, y.half3, y.gap3)


def test_strand_coordinates_and_origin():
    G = M.Genome("ACGTACGTAA", True)
    assert G.sub(-1, 0, 3) == "TTA"                               # reverse complement, transcript order
    assert G.sub(1, 8, 12) == "AAAC"                              # across the origin
    assert G.genomic(1, 8, 12) == [(8, 10), (0, 2)]
    assert G.genomic(-1, 0, 3) == [(7, 10)]
    assert G.genomic(-1, 8, 12) == [(0, 2), (8, 10)]              # minus strand, transcript order across the origin
    assert G.norm(1, 12, 15) == (1, 2, 5)
    assert M._signed_gap(M.Genome("A" * 20000, True), 100, 90, 3) < 0      # overlap is not a 20-kb intron
    assert M._signed_gap(M.Genome("A" * 20000, True), 19000, 500, 3) == 1499  # across the origin


def test_exclusions():
    g = M.RefGene("NC_1.1|cox1|1", "NC_1.1", ["KY742027"], "Decemunciger sp. AB-2017", [], 5, "cox1", "CDS", [], "", "")
    assert M.Exclude({"NC_1.2"})(g) and M.Exclude({"KY742027.1"})(g)          # any version, RefSeq copy
    assert M.Exclude(set(), species="decemunciger sp.")(g)
    assert M.Exclude(set(), genus="Decemunciger")(g)
    assert not M.Exclude({"KY774370.1"}, "decemunciger x", "other")(g)


def test_vote_prefers_close_majority():
    lib = {f"r{i}": M.RefGene(f"r{i}", f"A{i}", [], "x y", lin, 4, "cox1", "CDS", [], "", "") for i, lin in
           enumerate([list("abcdefg"), list("abcdefg"), list("abcdefg"), list("abcde")])}
    calls = [M.JCall("r0", 1, "present", 0, "cis", 1, 893, 1, 1747, ident=0.94),     # idiosyncratic close relative
             M.JCall("r1", 1, "present", 0, "cis", 1, 871, 1, 1743, ident=0.93),
             M.JCall("r2", 1, "present", 0, "cis", 1, 871, 1, 1743, ident=0.85),
             M.JCall("r3", 1, "present", 0, "cis", 1, 893, 1, 1747, ident=1.0)]
    rank = {k: i for i, k in enumerate(lib)}
    win = M.vote(calls, lib, list("abcdefg"), rank)
    assert (win.d, win.a) == (871, 1743)
    lib["r3"].own = True                                                             # the record's own annotation
    assert M.vote(calls, lib, list("abcdefg"), rank).gid == "r3"


def test_gold_placement_rules():
    """Exon 1 ends ..TT, intron ends ..GT: sliding 1 nt 5'-ward keeps the spliced sequence and gives exon-terminal T +
    omega-G (group I); a group II intron has no canonical placement there; a placement without a repeat cannot move;
    group typing: a covariance-model hit with E <= 1e-5 decides, a weaker one leaves both groups (or the reference's)."""
    exon1, intron, exon2 = "ACGACGAATT", "GCATGCATGC" * 10 + "AAGT", "CCATGCCATG"
    seq = exon1 + intron + exon2
    G = M.Genome(seq, False)
    ex = [(1, 0, len(exon1)), (1, len(exon1) + len(intron), len(seq))]
    assert M.gold_shift(G, ex, 0, {"I"}, 4, False)[0] == -1
    assert M.gold_shift(G, ex, 0, {"II"}, 4, False)[0] == 0
    moved = [(1, 0, len(exon1) - 1), (1, len(exon1) + len(intron) - 1, len(seq))]
    assert M.gold_shift(G, moved, 0, {"I"}, 4, False)[0] == 0                  # canonical already
    G2 = M.Genome("ACGACGAATC" + intron + exon2, False)                       # no repeat at the ends: no window
    assert M.gold_shift(G2, ex, 0, {"I", "II"}, 4, False)[0] == 0
    rc = M.Genome(str(Seq(seq).reverse_complement()), False)                  # minus strand, strand coordinates
    assert M.gold_shift(rc, [(-1, a, b) for _, a, b in ex], 0, {"I", "II"}, 4, False)[0] == -1
    assert M.allowed_groups("group II intron (Infernal/Rfam RF00029, E=1.6e-05)", "I") == {"I", "II"}
    assert M.allowed_groups("group I intron, subtype IB1 (model, E=4.7e-25)", "") == {"I"}
    assert M.allowed_groups("intron, type not identified", "II") == {"II"}
    assert M.allowed_groups("", "") == {"I", "II"}


# ── tblastn / blastn ─────────────────────────────────────────────────────────
@needs_blast
@pytest.mark.parametrize("shift,revcomp", [(4155 - 893, False), (12000, True)])
def test_own_annotation_exact(shift, revcomp):
    """Metridium rotated cox1-first (and reverse-complemented): its own sites give NCBI's introns exactly."""
    rec = ncbi("NC_000933")
    q = bare(rec, shift, revcomp)
    models, rows, written = annotate(q, [], 4, own=M.ref_genes(rec))
    assert sorted(written) == ["cox1", "nad5"]
    assert intron_keys(q) == intron_keys(rec)
    cds = {M.gene_name(f.qualifiers, "CDS"): f for f in q.features if f.type == "CDS"}
    assert "*" not in cds["cox1"].qualifiers["translation"][0]


@needs_blast
def test_own_annotation_quirks_and_mini_exon():
    """Plakinastrella: the 9-nt exon and the intron feature 6 nt off its CDS gap, as annotated."""
    rec = ncbi("NC_010217")
    q = bare(rec, 10830)
    annotate(q, [], 4, own=M.ref_genes(rec))
    assert intron_keys(q) == intron_keys(rec)


@needs_blast
def test_origin_crossing_intron_left_to_the_join():
    """Montipora nad5 intron (11.5 kb, hosts genes) put across the origin: no multi-interval intron feature (table2asn
    MultiIntervalIntron); the CDS join carries it, at NCBI's CDS gap."""
    rec = ncbi("NC_006902")
    q = bare(rec, 9000)
    models, rows, written = annotate(q, [], 4, own=M.ref_genes(rec))
    assert written == ["nad5"]
    assert not [f for f in q.features if f.type == "intron" and len(f.location.parts) > 1]
    cds_gap = M.ref_genes(rec)[0].junctions[0]
    assert len(join_gap_keys(q, "nad5")) == 1
    (seq, left, right), = join_gap_keys(q, "nad5")
    assert len(seq) == cds_gap.length and left.endswith(cds_gap.e5[-10:])
    assert any("crosses the origin" in n for f in q.features if f.type == "CDS" for n in f.qualifiers.get("note", []))


@needs_blast
def test_transfer_between_sea_anemones():
    """Exaiptasia introns from the Metridium sites (another family, same order): both found at NCBI's positions."""
    query, ref = ncbi("NC_056771"), ncbi("NC_000933")
    q = bare(query, 16279 - 890)
    lib = M.ref_genes(ref)
    models, rows, written = annotate(q, lib, 4, lineage=query.annotations["taxonomy"])
    assert sorted(written) == ["cox1", "nad5"]
    assert intron_keys(q) == intron_keys(query)


@needs_blast
def test_no_intron_in_negative_control():
    """Platynereis (annelid, no introns) with cnidarian and sponge sites: nothing written."""
    lib = [g for a in ("NC_000933", "NC_056771", "NC_010217", "NC_006902") for g in M.ref_genes(ncbi(a))]
    rec = ncbi("NC_000931")
    q = bare(rec, 1000)
    models, rows, written = annotate(q, lib, 5, lineage=rec.annotations["taxonomy"])
    assert written == [] and not [f for f in q.features if f.type == "intron"]
    assert not [r for r in rows if r["state"] == "present"]


@needs_blast
def test_command_line(tmp_path):
    lib = tmp_path / "lib" / "sites.tsv.gz"
    subprocess.run([sys.executable, "-W", "ignore", str(SCRIPT), "build", "-i", str(DATA / "NC_000933.gb"), "-o", str(lib)],
                   check=True, capture_output=True)
    assert lib.exists() and M.fasta_path(lib).exists()
    rec = ncbi("NC_056771")
    inp = tmp_path / "in"
    inp.mkdir()
    q = bare(rec, 16279 - 890)
    SeqIO.write(q, str(inp / "NC_056771.gbf"), "genbank")
    out = tmp_path / "out"
    subprocess.run([sys.executable, "-W", "ignore", str(SCRIPT), "annotate", "-i", str(inp), "-o", str(out), "--lib",
                    str(lib), "--genetic_code", "4", "--bin", BIN, "--lineage", ";".join(rec.annotations["taxonomy"])],
                   check=True, capture_output=True)
    assert intron_keys(SeqIO.read(str(out / "NC_056771.gbf"), "genbank")) == intron_keys(rec)
    assert (out / "intron_ref_report.tsv").read_text().count("\tpresent\t") == 2
    # own annotation through --ref_gb (same accession as the query), library ignored
    out2 = tmp_path / "out2"
    subprocess.run([sys.executable, "-W", "ignore", str(SCRIPT), "annotate", "-i", str(inp), "-o", str(out2),
                    "--ref_gb", str(DATA / "NC_056771.gb"), "--own_only", "--genetic_code", "4", "--bin", BIN],
                   check=True, capture_output=True)
    assert intron_keys(SeqIO.read(str(out2 / "NC_056771.gbf"), "genbank")) == intron_keys(rec)


@needs_blast
def test_merge_with_de_novo_model():
    """The record already holds a de novo model (Plakinastrella cox1, both introns, the second one 2 codons off); the
    transfer knows only the second site: that intron is re-placed, the first one (no library site) is kept."""
    from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
    rec = ncbi("NC_010217")
    q = bare(rec, 0)
    exons = [(10830, 11553), (11935, 11944), (13068, 13905)]          # NCBI: 10831..11553, 11936..11944, 13063..13905
    loc = CompoundLocation([FeatureLocation(a, b, strand=1) for a, b in exons])
    q.features = [SeqFeature(loc, type="CDS", qualifiers={"gene": ["COX1"]}),
                  SeqFeature(FeatureLocation(11553, 11935, strand=1), type="intron",
                             qualifiers={"gene": ["COX1"], "note": ["group I intron (de novo)"]}),
                  SeqFeature(FeatureLocation(11944, 13068, strand=1), type="intron", qualifiers={"gene": ["COX1"]})]
    lib = M.ref_genes(rec)
    lib[0].junctions = [j for j in lib[0].junctions if j.k == 2]       # the library has the second site only
    models, rows, written = annotate(q, lib, 4)
    assert written == ["COX1"]                                          # the record's own gene name is kept
    introns = sorted((int(f.location.start) + 1, int(f.location.end)) for f in q.features if f.type == "intron")
    assert introns == [(11554, 11935), (11945, 13062)]
    assert any("de novo model kept, 1 intron(s) re-placed" in n for f in q.features if f.type == "CDS"
               for n in f.qualifiers.get("note", []))


def test_leave_out_from_environment(monkeypatch):
    """mitom_lo: MITOM_EXCLUDE_REF leaves the record and its RefSeq / INSDC copy out, MITOM_EXCLUDE_LEVEL its species /
    genus; unset: the whole library."""
    lib = [g for a in ("NC_000933", "NC_056771", "NC_010217") for g in M.ref_genes(ncbi(a))]
    twin = M.RefGene("XY000001.1|cox1|1", "XY000001.1", [], "Metridium senile", lib[0].lineage, 4, "cox1", "CDS",
                     [], "", "")
    lib[0].same_as = ["XY000001"]                       # NC_000933.1 is a copy of XY000001.1
    lib.append(twin)
    rec = ncbi("NC_000933")
    monkeypatch.delenv("MITOM_EXCLUDE_REF", raising=False)
    assert M.env_exclude(lib, rec) is None
    monkeypatch.setenv("MITOM_EXCLUDE_REF", "NC_000933.1")
    ex = M.env_exclude(lib, rec)
    assert sorted({g.acc for g in lib if ex(g)}) == ["NC_000933.1", "XY000001.1"]
    monkeypatch.setenv("MITOM_EXCLUDE_REF", "XY000001.1")        # the copy named: the RefSeq record goes too
    assert sorted({g.acc for g in lib if M.env_exclude(lib, rec)(g)}) == ["NC_000933.1", "XY000001.1"]
    monkeypatch.setenv("MITOM_EXCLUDE_REF", "NC_056771.1")
    monkeypatch.setenv("MITOM_EXCLUDE_LEVEL", "genus")
    assert {g.acc for g in lib if M.env_exclude(lib, ncbi("NC_056771"))(g)} == {"NC_056771.1"}
    monkeypatch.setenv("MITOM_EXCLUDE_LEVEL", "species")
    rec_unknown = bare(ncbi("NC_000933"))
    rec_unknown.annotations["organism"] = "Metridium senile"
    monkeypatch.setenv("MITOM_EXCLUDE_REF", "ZZ999999.1")        # not in the library: the query's species
    assert {g.acc for g in lib if M.env_exclude(lib, rec_unknown)(g)} == {"NC_000933.1", "XY000001.1"}


@needs_blast
def test_lineage_option_is_a_fallback(tmp_path):
    """--lineage serves records without a taxonomy (FASTA input); a record's own taxonomy is never overridden."""
    lib = tmp_path / "lib" / "sites.tsv.gz"
    subprocess.run([sys.executable, "-W", "ignore", str(SCRIPT), "build", "-i", str(DATA / "NC_000933.gb"), "-o", str(lib)],
                   check=True, capture_output=True)
    rec = ncbi("NC_056771")
    inp = tmp_path / "in"
    inp.mkdir()
    no_tax = bare(rec, 16279 - 890)
    no_tax.id = no_tax.name = "QUERY1"
    with_tax = bare(rec, 16279 - 890)
    with_tax.id = with_tax.name = "QUERY2"
    with_tax.annotations["taxonomy"] = rec.annotations["taxonomy"]
    SeqIO.write(no_tax, str(inp / "q1.gbf"), "genbank")
    SeqIO.write(with_tax, str(inp / "q2.gbf"), "genbank")
    out = tmp_path / "out"
    subprocess.run([sys.executable, "-W", "ignore", str(SCRIPT), "annotate", "-i", str(inp), "-o", str(out), "--lib",
                    str(lib), "--genetic_code", "4", "--bin", BIN, "--lineage", "Eukaryota;Metazoa;Cnidaria"],
                   check=True, capture_output=True)
    import csv
    rows = list(csv.DictReader(open(out / "intron_ref_report.tsv"), delimiter="\t"))
    depth = {r["record"]: {int(r["lineage_shared"]) for r in rows if r["record"] == r_ and r["lineage_shared"]}
             for r_ in ("QUERY1", "QUERY2") for r in rows if r["record"] == r_}
    assert depth["QUERY1"] == {3}                                  # fallback: Eukaryota;Metazoa;Cnidaria
    assert min(depth["QUERY2"]) >= 5                               # its own lineage: Hexacorallia, Actiniaria ...
    for name in ("q1.gbf", "q2.gbf"):
        assert intron_keys(SeqIO.read(str(out / name), "genbank")) == intron_keys(rec)
