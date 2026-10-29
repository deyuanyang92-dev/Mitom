"""Tests for mito_intron.py against NCBI records with known group I introns.

Metridium senile NC_000933: COX1 and ND5 each contain a group I intron.
Nematostella NC_008164: ND5 is split (CDS join) but the record has no intron feature.
Platynereis dumerilii NC_000931: no introns (negative control).
Library: MFannot protein collections with the record's own species excluded.

Needs exonerate and cmsearch: on PATH, or MITO_INTRON_BIN=/path/to/env/bin.
Run: python3 -m pytest -q tests/test_mito_intron.py
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_env  # noqa: E402
DATA = HERE / "data"
SCRIPT = HERE.parent / "mito_intron.py"
NOTE = "exon-intron structure predicted by protein alignment (Exonerate, MFannot method)"
BIN = os.environ.get("MITO_INTRON_BIN") or mito_env.bin_with("exonerate", "cmsearch") or ""
pytestmark = pytest.mark.skipif(not BIN, reason="exonerate not available")


def run(inp: Path, out: Path, code: int) -> None:
    """MFannot library with the record's own species excluded (no circular validation)."""
    cmd = [sys.executable, "-W", "ignore", str(SCRIPT), "-i", str(inp), "-o", str(out),
           "--genetic_code", str(code), "--bin", BIN, "--exclude_own_species"]
    subprocess.run(cmd, check=True, capture_output=True)


def ncbi_cds(acc):
    rec = SeqIO.read(str(DATA / f"{acc}.gb"), "genbank")
    return {f.qualifiers["gene"][0]: (str(f.location), f.qualifiers["translation"][0])
            for f in rec.features if f.type == "CDS" and "gene" in f.qualifiers}


def intron_models(path: Path):
    rec = SeqIO.read(str(path), "genbank")
    cds = {f.qualifiers["gene"][0]: (str(f.location), f.qualifiers["translation"][0])
           for f in rec.features if f.type == "CDS" and NOTE in f.qualifiers.get("note", [])}
    introns = [(int(f.location.start) + 1, int(f.location.end), f.qualifiers["note"][0])
               for f in rec.features if f.type == "intron"]
    return cds, introns


def test_metridium_introns_recovered_from_mitoz_like_input(tmp_path):
    """Introns removed and split genes cut to exon 1 (as MitoZ leaves them) -> NCBI models restored."""
    rec = SeqIO.read(str(DATA / "NC_000933.gb"), "genbank")
    feats = []
    for f in rec.features:
        if f.type == "intron":
            continue
        if len(f.location.parts) > 1 or (f.type == "gene" and f.qualifiers.get("gene", [""])[0] in ("COX1", "ND5")):
            p = f.location.parts[0]
            f.location = FeatureLocation(p.start, p.end, strand=p.strand)
            f.qualifiers.pop("translation", None)
        feats.append(f)
    rec.features = feats
    (tmp_path / "in").mkdir()
    SeqIO.write(rec, str(tmp_path / "in" / "NC_000933.gbf"), "genbank")
    run(tmp_path / "in", tmp_path / "out", 4)
    cds, introns = intron_models(tmp_path / "out" / "NC_000933.gbf")
    ref = ncbi_cds("NC_000933")
    for gene in ("COX1", "ND5"):
        assert cds[gene] == ref[gene], gene               # same CDS coordinates and protein as NCBI
    assert [(s, e) for s, e, _ in introns] == [(4156, 5008), (7744, 9424)]
    assert all("group I intron" in n for _, _, n in introns)


def test_nematostella_missing_intron_added(tmp_path):
    (tmp_path / "in").mkdir()
    shutil.copy(DATA / "NC_008164.gb", tmp_path / "in" / "NC_008164.gbf")
    run(tmp_path / "in", tmp_path / "out", 4)
    cds, introns = intron_models(tmp_path / "out" / "NC_008164.gbf")
    assert cds["ND5"] == ncbi_cds("NC_008164")["ND5"]
    assert [(s, e) for s, e, _ in introns] == [(718, 2337)]


def test_intronless_record_unchanged(tmp_path):
    (tmp_path / "in").mkdir()
    shutil.copy(DATA / "NC_000931.gb", tmp_path / "in" / "NC_000931.gbf")
    run(tmp_path / "in", tmp_path / "out", 5)
    assert (tmp_path / "out" / "NC_000931.gbf").read_bytes() == (DATA / "NC_000931.gb").read_bytes()
    assert (tmp_path / "out" / "intron_report.tsv").read_text().count("\n") == 1   # header only


def mitoz_like(acc: str, tmp_path: Path) -> Path:
    """NCBI record as MitoZ would leave it: intron features removed, split CDS/genes cut to exon 1."""
    rec = SeqIO.read(str(DATA / f"{acc}.gb"), "genbank")
    split = {g for f in rec.features if f.type == "CDS" and len(f.location.parts) > 1 for g in f.qualifiers.get("gene", [])}
    feats = []
    for f in rec.features:
        if f.type == "intron":
            continue
        if len(f.location.parts) > 1 or (f.type == "gene" and set(f.qualifiers.get("gene", [])) & split):
            p = f.location.parts[0]
            f.location = FeatureLocation(p.start, p.end, strand=p.strand)
            f.qualifiers.pop("translation", None)
        feats.append(f)
    rec.features = feats
    d = tmp_path / "in"
    d.mkdir()
    SeqIO.write(rec, str(d / f"{acc}.gbf"), "genbank")
    return d


def report(out: Path):
    import csv
    return [r for r in csv.DictReader((out / "intron_report.tsv").open(), delimiter="\t") if r["intron_start"]]


def test_mini_exon_plakinastrella(tmp_path):
    """Sponge COX1: two introns separated by a mini-exon (NCBI join(10831..11547,11936..11944,13063..13905)) must
    not be merged. Since 2026-10-29 the introns around the inserted mini-exon get the boundary step again (exon 1 had
    run 9 nt into the intron, +3 aa): the gold placement, the canonical group I slide of NCBI's join with the same
    protein (exon 1 ends U, intron 1 ends omega-G)."""
    run(mitoz_like("NC_010217", tmp_path), tmp_path / "out", 4)
    cds, introns = intron_models(tmp_path / "out" / "NC_010217.gbf")
    assert len(introns) == 2
    rec = SeqIO.read(str(tmp_path / "out" / "NC_010217.gbf"), "genbank")
    cox1 = [f for f in rec.features if f.type == "CDS" and NOTE in f.qualifiers.get("note", [])][0]
    assert [(int(p.start) + 1, int(p.end)) for p in cox1.location.parts] == [(10831, 11546), (11935, 11944), (13063, 13905)]
    assert [x[:2] for x in introns] == [(11547, 11934), (11945, 13062)]
    ncbi = SeqIO.read(str(DATA / "NC_010217.gb"), "genbank")
    ref = [f for f in ncbi.features if f.type == "CDS" and len(f.location.parts) == 3][0]
    assert cox1.qualifiers["translation"][0] == ref.qualifiers["translation"][0]


def test_conflicting_model_not_written(tmp_path):
    """Montipora: a model fusing an ND4-region exon with nad5 exon 2 must not overwrite MitoZ's ND4
    (since 2026-10-21 that model is replaced by a re-aligned one, see test_next_chain_when_best_conflicts)."""
    run(mitoz_like("NC_006902", tmp_path), tmp_path / "out", 4)
    rows = report(tmp_path / "out")
    assert all(r["written"] == "no" for r in rows if "overlaps existing ND4" in r["evidence"])
    out_rec = SeqIO.read(str(tmp_path / "out" / "NC_006902.gbf"), "genbank") \
        if (tmp_path / "out" / "NC_006902.gbf").exists() else None
    if out_rec is not None:
        nd4 = [f for f in out_rec.features if f.type == "CDS" and f.qualifiers.get("gene") == ["ND4"]]
        assert nd4 and str(nd4[0].location) == ncbi_cds("NC_006902")["ND4"][0]   # MitoZ's ND4 untouched


def test_long_nd5_three_prime_end(tmp_path):
    """NC_056771 ND5 is 724 aa: the 3' end must reach the first in-frame stop (NCBI 3902)."""
    run(mitoz_like("NC_056771", tmp_path), tmp_path / "out", 4)
    cds, _ = intron_models(tmp_path / "out" / "NC_056771.gbf")
    assert cds["ND5"][0] == ncbi_cds("NC_056771")["ND5"][0]


def test_submission_qualifiers(tmp_path):
    """Every CDS written has /product; T/TA-terminated CDS carry /transl_except and the polyA note."""
    shutil.copy(DATA / "NC_008164.gb", tmp_path / "x.gbf")
    (tmp_path / "in").mkdir()
    shutil.move(str(tmp_path / "x.gbf"), tmp_path / "in" / "NC_008164.gbf")
    run(tmp_path / "in", tmp_path / "out", 4)
    rec = SeqIO.read(str(tmp_path / "out" / "NC_008164.gbf"), "genbank")
    ours = [f for f in rec.features if f.type == "CDS" and NOTE in f.qualifiers.get("note", [])]
    assert ours and all(f.qualifiers.get("product") for f in ours)
    for f in ours:
        tail = len(f.location.extract(rec.seq)) % 3
        assert bool(tail) == ("transl_except" in f.qualifiers)


def test_ref_gb_names_and_filters():
    """--ref_gb: CDS translations become homologs named like the MFannot library (nad5, cox1...)."""
    sys.path.insert(0, str(SCRIPT.parent))
    import mito_intron as mi
    lib = mi.load_ref_gb([str(DATA / "NC_000933.gb")], 4)
    genes = sorted({p.gene for p in lib.values()})
    assert "cox1" in genes and "nad5" in genes and "nad4L" in genes      # MFannot spelling
    assert all("*" not in p.seq for p in lib.values())
    assert {p.species for p in lib.values()} == {"metridium senile"}


def test_partial_cds_bare_fragment(tmp_path):
    """Travisia COX1 fragment, all annotation removed: the gene runs off the 5' end, so the CDS is
    5'-partial (no start codon) - the intron (NCBI 236..1100, group II) must still be written."""
    from Bio.SeqFeature import BeforePosition
    rec = SeqIO.read(str(DATA / "LC670763.gb"), "genbank")
    rec.features = [f for f in rec.features if f.type == "source"]
    (tmp_path / "in").mkdir()
    SeqIO.write(rec, str(tmp_path / "in" / "LC670763.gbf"), "genbank")
    run(tmp_path / "in", tmp_path / "out", 5)
    _, introns = intron_models(tmp_path / "out" / "LC670763.gbf")
    assert [(s, e) for s, e, _ in introns] == [(236, 1100)] and "group II" in introns[0][2]
    cds = [f for f in SeqIO.read(str(tmp_path / "out" / "LC670763.gbf"), "genbank").features if f.type == "CDS"][0]
    assert isinstance(cds.location.parts[0].start, BeforePosition)
    # the frame runs to the last base of the fragment without a stop: 3'-partial as well (NCBI 1101..>1523)
    from Bio.SeqFeature import AfterPosition
    assert isinstance(cds.location.parts[-1].end, AfterPosition) and int(cds.location.parts[-1].end) == 1523


def test_next_chain_when_best_conflicts(tmp_path):
    """2026-10-21 dev benchmark NC_006902.1: the best exonerate chain puts nad5 exon 1 on ND4 (nad4/nad5 are
    homologous) and is refused; the next chain that does not overlap another gene gives NCBI's ND5."""
    run(mitoz_like("NC_006902", tmp_path), tmp_path / "out", 4)
    cds, _ = intron_models(tmp_path / "out" / "NC_006902.gbf")
    assert cds["ND5"][0] == ncbi_cds("NC_006902")["ND5"][0]


def test_feature_notes_read_like_genbank_notes(tmp_path):
    """Review C finding 15: /note qualifiers reach the submission table, so they carry no script name or
    to-do: the CDS note is the neutral NOTE, intron notes name the 'MFannot method'."""
    rec = SeqIO.read(str(DATA / "LC670763.gb"), "genbank")
    rec.features = [f for f in rec.features if f.type == "source"]
    (tmp_path / "in").mkdir()
    SeqIO.write(rec, str(tmp_path / "in" / "LC670763.gbf"), "genbank")
    run(tmp_path / "in", tmp_path / "out", 5)
    out = SeqIO.read(str(tmp_path / "out" / "LC670763.gbf"), "genbank")
    cds = [f for f in out.features if f.type == "CDS"]
    introns = [f for f in out.features if f.type == "intron"]
    assert len(cds) == 1 and cds[0].qualifiers["note"] == [NOTE] and introns
    for f in cds + introns:
        for note in f.qualifiers.get("note", []):
            assert ".py" not in note and "mito_intron" not in note and "manually" not in note.lower()
    assert all("MFannot method" in f.qualifiers["note"][0] for f in introns)


def test_existing_stopless_cds_end_not_reused(tmp_path):
    """Native core review 2026-10: the core ends a CDS at its last aligned codon; core_ends() handed that stop-less
    end (NC_008164.1 ND5 ...GTA, the T stop 1 nt further) to the model, mito_curate P3 skips joins -> table2asn NoStop.
    An existing end is reused only when it is a stop / T-TA tail; otherwise the first in-frame complete stop (3468)."""
    rec = SeqIO.read(str(DATA / "NC_008164.gb"), "genbank")
    feats = []
    for f in rec.features:
        if f.type == "intron":
            continue
        if f.type in ("CDS", "gene") and f.qualifiers.get("gene") == ["ND5"]:
            f.location = FeatureLocation(2337, 3435, 1)                   # exon 2 only, ending on GTA (NCBI: ..3436 = T)
            f.qualifiers.pop("translation", None)
        feats.append(f)
    rec.features = feats
    (tmp_path / "in").mkdir()
    SeqIO.write(rec, str(tmp_path / "in" / "NC_008164.gbf"), "genbank")
    run(tmp_path / "in", tmp_path / "out", 4)
    cds, introns = intron_models(tmp_path / "out" / "NC_008164.gbf")
    assert cds["ND5"][0] == "join{[0:717](+), [2337:3468](+)}" and [(s, e) for s, e, _ in introns] == [(718, 2337)]
    import mito_intron as mi
    assert mi.ends_with_stop(str(rec.seq), [(1, 717), (2338, 3468)], 4)
    assert not mi.ends_with_stop(str(rec.seq), [(1, 717), (2338, 3435)], 4)
    assert mi.ends_with_stop(str(rec.seq), [(1, 717), (2338, 3436)], 4)    # NCBI's T tail


def test_chained_alignments_too_far_apart_are_not_an_intron():
    """Stress case 4a (linear record rotated by 8800 nt): the two pieces of nad5 at both sequence ends were chained
    into CDS join(<3..875,15777..16559) with a 14,901-nt 'intron'. A gene cut by the ends of a record has its
    N-terminal piece at the end and its C-terminal piece at the start: out of protein order, refused on a linear
    record (since 2026-10-29 by order, not by length); the 11,485-nt coral nad5 intron (Montipora NC_006902.1,
    test_next_chain_when_best_conflicts) and the 18,791-nt corallimorpharian one (NC_027100.1) pass."""
    sys.path.insert(0, str(SCRIPT.parent))
    import mito_intron as mi
    cut_pieces = [mi.Report(292, 552, 3, 875, 900, 1, [(3, 875)]), mi.Report(1, 291, 15777, 16559, 800, 1, [(15777, 16559)])]
    m, why = mi.build_model("p", "nad5", cut_pieces, 560)
    assert m is None and "out of protein order" in why
    m, why = mi.build_model("p", "nad5", cut_pieces, 560, circular=True)
    assert m is None and why == mi.ACROSS_ORIGIN               # modelled on a rotated copy (model_across_origin)
    far = [mi.Report(1, 291, 3, 875, 900, 1, [(3, 875)]), mi.Report(292, 552, 15777 + mi.MAX_CHAIN_INTRON,
                                                                      16559 + mi.MAX_CHAIN_INTRON, 800, 1,
                                                                      [(15777 + mi.MAX_CHAIN_INTRON, 16559 + mi.MAX_CHAIN_INTRON)])]
    m, why = mi.build_model("p", "nad5", far, 560)
    assert m is None and str(mi.MAX_CHAIN_INTRON) in why
    coral = [mi.Report(1, 240, 2440, 3159, 900, 1, [(2440, 3159)]), mi.Report(241, 612, 14645, 15760, 800, 1, [(14645, 15760)])]
    m, why = mi.build_model("p", "nad5", coral, 620)
    assert m is not None and m.exons == [(2440, 3159), (14645, 15760)]
    coralli = [mi.Report(0, 240, 0, 720, 900, 1, [(1, 720)]), mi.Report(240, 612, 19511, 20630, 800, 1, [(19512, 20630)])]
    m, why = mi.build_model("p", "nad5", coralli, 612, circular=True)
    assert m is not None and m.exons == [(1, 720), (19512, 20630)]
    # placozoan cox1 is trans-spliced: the + strand pieces chained across 17.7 kb leave the 32 residues of the exon on
    # the other strand unaligned (NC_041549.1, rotated record) - not one intron
    placo = [mi.Report(0, 240, 0, 720, 1218, 1, [(1, 720)]),
             mi.Report(272, 534, 18424, 19379, 1226, 1, [(18425, 18479), (18649, 19379)])]
    m, why = mi.build_model("p", "cox1", placo, 534, circular=True)
    assert m is None and "trans-spliced" in why and "241-272" in why
    # linear record, pieces at both ends: refused when chained across more than --maxintron; a fragment with a short
    # intron running off both ends (LC670763 cox1 join(<1..235,1101..>1523), test_partial_cds_bare_fragment) is kept;
    # on a circular record the (collinear) model is kept: a gene across the origin never gets here
    cut = mi.Model("p", "nad5", 1, [(3, 875), (10777, 11559)], 1700, 1.0)
    assert "cut by both ends" in mi.spans_record_ends(cut, 11560, False, 5000)
    assert mi.spans_record_ends(cut, 11560, True, 5000) == ""
    frag = mi.Model("p", "cox1", 1, [(1, 235), (1101, 1523)], 1700, 1.0)
    assert mi.spans_record_ends(frag, 1523, False, 5000) == ""


def test_chaining_is_collinear_and_crosses_the_origin():
    """group_reports chains alignments only in protein order along the transcript (both strands); on a circular record
    the step from the last piece before the origin to the first piece after it is the distance across the origin."""
    sys.path.insert(0, str(SCRIPT.parent))
    import mito_intron as mi
    L = 19000
    n_term = mi.Report(0, 160, 17200, 17690, 900, 1, [(17201, 17690)])       # exon 1 at the end of the record
    c_term = mi.Report(165, 511, 40, 1098, 1500, 1, [(41, 1098)])           # exon 2 at its start
    assert mi.chain_step([n_term], [c_term], L, circular=False) is None
    assert mi.chain_step([n_term], [c_term], L, circular=True) == 40 - 17690 + L
    assert [len(g) for g in mi.group_reports([c_term, n_term], 511, L, circular=True)] == [2]
    assert [len(g) for g in mi.group_reports([c_term, n_term], 511, L, circular=False)] == [1, 1]
    # minus strand: transcript runs to lower coordinates (exonerate tab > tae)
    a = mi.Report(0, 423, 12779, 10895, 1700, -1, [(10896, 11375), (11973, 12779)])
    b = mi.Report(471, 575, 10829, 10517, 450, -1, [(10518, 10829)])
    assert mi.chain_step([a], [b]) == 66 and mi.wrap_steps([a, b]) == []              # 10829 -> 10895
    assert mi.wrap_steps([b, mi.Report(600, 700, 15000, 14700, 300, -1, [(14701, 15000)])]) == [0]


def test_move_location_across_the_origin():
    """Locations moved on a circular sequence: a part crossing the origin becomes a join in transcript order (both
    strands), parts meeting at the old origin are merged, fuzzy ends are kept; moving back restores the location."""
    sys.path.insert(0, str(SCRIPT.parent))
    import mito_intron as mi
    from Bio.SeqFeature import BeforePosition, CompoundLocation
    L = 1000
    plus = CompoundLocation([FeatureLocation(BeforePosition(100), 300, 1), FeatureLocation(500, 700, 1)])
    moved = mi.move_location(plus, -250, L)                  # copy starting at position 251
    assert [(int(p.start), int(p.end)) for p in moved.parts] == [(850, 1000), (0, 50), (250, 450)]
    assert isinstance(moved.parts[0].start, BeforePosition)
    back = mi.move_location(moved, 250, L)
    assert [(int(p.start), int(p.end)) for p in back.parts] == [(100, 300), (500, 700)]
    assert isinstance(back.parts[0].start, BeforePosition)
    minus = CompoundLocation([FeatureLocation(0, 80, -1), FeatureLocation(900, 1000, -1)])   # transcript order
    m2 = mi.move_location(minus, 200, L)
    assert m2.strand == -1 and (int(m2.start), int(m2.end)) == (100, 280)
    m3 = mi.move_location(m2, -200, L)
    assert [(int(p.start), int(p.end)) for p in m3.parts] == [(0, 80), (900, 1000)]


def test_gene_across_the_origin_is_modelled(tmp_path):
    """MitoM rotates circular records to start at the core's cox1 piece; when that piece is exon 2, exon 1 and the
    intron lie at the end of the record (2026-10 annelid runs: cox1 of KY742027.1, LC677172.1, NC_028710.1, NC_087806.1
    refused as 'chained alignments 14-16 kb apart'). Metridium COX1 with the record rotated into its intron: the
    gene is modelled on a rotated copy and written across the origin, with NCBI's exons and protein."""
    sys.path.insert(0, str(SCRIPT.parent))
    import mito_intron as mi
    d = mitoz_like("NC_000933", tmp_path)
    rec = SeqIO.read(str(d / "NC_000933.gbf"), "genbank")
    rot = mi.rotated_copy(rec, 4200)                        # starts inside the COX1 intron (4156..5008)
    rot.annotations["molecule_type"] = "DNA"
    SeqIO.write(rot, str(d / "NC_000933.gbf"), "genbank")
    run(d, tmp_path / "out", 4)
    out = SeqIO.read(str(tmp_path / "out" / "NC_000933.gbf"), "genbank")
    cox1 = [f for f in out.features if f.type == "CDS" and f.qualifiers.get("gene") == ["COX1"]]
    assert len(cox1) == 1
    assert [(int(p.start) + 1, int(p.end)) for p in cox1[0].location.parts] == [(16507, 17399), (810, 1509)]
    assert cox1[0].qualifiers["translation"][0] == ncbi_cds("NC_000933")["COX1"][1]
    # the intron (17400..L + 1..809) has no feature: NCBI rejects an intron in two intervals (MultiIntervalIntron)
    assert not [f for f in out.features if f.type == "intron" and f.qualifiers.get("gene") == ["COX1"]]
    rows = report(tmp_path / "out")
    assert any("spans the origin" in r["evidence"] and "CDS join only" in r["evidence"] and r["written"] == "yes"
               and r["intron_start"] == "17400" and r["intron_end"] == "809" for r in rows)


def test_untyped_intron_with_group_ii_splice_sites():
    """Chaetopterus cox1 introns (NC_028710.1, PP895294-6) have no Rfam domain V/VI hit but the group II splice
    sites GUGCG...AC: an untyped intron whose best GUGYG...AY score reaches TYPE_II_MIN is placed by the group II
    rule (was scored as group I and held back as 'weak splice site')."""
    sys.path.insert(0, str(SCRIPT.parent))
    import mito_intron as mi
    ex1, ex2 = "ATG" + "GCT" * 30 + "GA", "T" + "GCT" * 30 + "TAA"
    intron = "GTGCG" + "ATTTAAATTA" * 20 + "TAC"
    s = ex1 + intron + ex2
    a = len(ex1)
    exons = [(1, a - 1), (a + len(intron), len(s))]             # exonerate's junction 1 nt off the splice sites
    adj, notes, scores = mi.adjust_introns(s, exons, [""], 4, lambda p: 0)
    assert adj == [(1, a), (a + len(intron) + 1, len(s))]
    assert scores[0] >= mi.TYPE_II_MIN and "group II splice sites" in notes[0]
    assert mi.evidence(s, adj, [""], scores, 4)[0] == "medium"


def test_core_reference_proteins(tmp_path):
    """The core's CDS note names its reference protein (refseq89m) and the residues it aligned: a gene whose CDS leaves
    >= STUB_MIN_AA residues unaligned at one end is an intron candidate even when MFannot's ORF pre-screen misses it
    (Decemunciger nad4, exon 2 of 52 aa), with that protein (its closest homolog) as the Exonerate query."""
    sys.path.insert(0, str(SCRIPT.parent))
    import mito_intron as mi
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature
    (tmp_path / "featureProt").mkdir()
    (tmp_path / "featureProt" / "nad4.fas").write_text(">NC_026855:nad4-1-2979-4322 Some worm\n" + "M" * 447 + "\n")
    (tmp_path / "featureProt" / "cox1.fas").write_text(">NC_011014:cox1-1-0-1544 Other worm\n" + "W" * 511 + "\n")
    rec = SeqRecord(Seq("A" * 3000), id="x")
    rec.features = [SeqFeature(FeatureLocation(0, 1191, 1), type="CDS", qualifiers={
                        "gene": ["ND4"], "note": ["ref NC_026855:nad4-1-2979-4322 score 900 residues 19-395/447"]}),
                    SeqFeature(FeatureLocation(1500, 3000, -1), type="CDS", qualifiers={
                        "gene": ["COX1"], "note": ["ref NC_011014:cox1-1-0-1544 score 2400 residues 3-500/511"]})]
    mi._REFSEQ.clear()
    refs = mi.core_references(rec, tmp_path)
    mi._REFSEQ.clear()
    assert refs["nad4"]["stub"] and refs["nad4"]["seq"] == "M" * 447 and refs["nad4"]["strand"] == 1
    assert not refs["cox1"]["stub"] and refs["cox1"]["strand"] == -1
    assert mi.core_references(rec, None) == {}


def test_mini_exon_search_capped(monkeypatch):
    """A homolog gap longer than MINI_EXON_MAX_AA, or an intron longer than MAX_INTRON, is not searched for a mini-exon
    (NC_015309.1: hours in find_mini_exon over a 19-kb intron with a 60-aa gap)."""
    sys.path.insert(0, str(SCRIPT.parent))
    import mito_intron as mi
    called = []
    monkeypatch.setattr(mi, "find_mini_exon", lambda *x: called.append(x))
    monkeypatch.setattr(mi, "homolog_gap", lambda *x: (60.0, 1.0))
    s = "ATG" * 4000
    mi.add_mini_exons(s, [(1, 300), (600, 900)], ["I"], ["n"], [20.0], 4, ["M" * 200])
    monkeypatch.setattr(mi, "homolog_gap", lambda *x: (3.0, 1.0))
    mi.add_mini_exons(s, [(1, 300), (6000, 6300)], ["I"], ["n"], [20.0], 4, ["M" * 200])
    assert not called
    mi.add_mini_exons(s, [(1, 300), (600, 900)], ["I"], ["n"], [20.0], 4, ["M" * 200])
    assert len(called) == 1


def test_core_only_gene_needs_typed_introns(monkeypatch, tmp_path):
    """A gene that is a candidate only through the core's partial CDS (no ORF pre-screen hit) is written only when
    Rfam types its introns: exonerate against a close homolog bridges divergent or frameshifted stretches with untyped
    'introns' >= 142 nt (sponge NC_006894.1 nad2 / nad5, 2026-10-29 benchmark). The same model from the pre-screen
    is written (medium: untyped intron with group II splice sites)."""
    import argparse
    sys.path.insert(0, str(SCRIPT.parent))
    import mito_intron as mi
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    ex1, ex2 = "ATG" + "GCT" * 30 + "GA", "T" + "GCT" * 30 + "TAA"
    intron = "GTGCG" + "ATTTAAATTA" * 20 + "TAC"
    s = ex1 + intron + ex2
    prot = str(Seq(ex1 + ex2).translate(table=4)).rstrip("*")
    rep = mi.Report(0, len(prot), 0, len(s) - 3, 300, 1, [(1, len(ex1)), (len(ex1) + len(intron) + 1, len(s) - 3)])
    monkeypatch.setattr(mi, "run_exonerate", lambda *x: [rep])
    monkeypatch.setattr(mi, "core_references", lambda *x: {"nad2": {"id": "X:nad2-1", "seq": prot, "strand": 1,
                                                                    "score": 500, "stub": True}})
    a = argparse.Namespace(genetic_code=4, exclude_own_species=False, refdir=None, threads=1, maxintron=5000,
                           cm_files=[], evalue=0.01, stop_rule="first")
    tools = {"exonerate": "exonerate", "cmsearch": None, "blastp": "blastp"}
    lib = {"L1": mi.LibProt("some species", "nad2", prot)}

    def rows_for(cands):
        monkeypatch.setattr(mi, "prescreen", lambda *x: dict(cands))
        rec = SeqRecord(Seq(s), id="x", annotations={"topology": "linear", "molecule_type": "DNA"})
        rows, across = mi.process_record(rec, lib, tmp_path / "db", a, tools, tmp_path)
        assert not across
        return [r for r in rows if r.get("intron_start")]
    core_only = rows_for({})
    assert core_only and core_only[0]["written"] == "no" and "core's partial CDS" in core_only[0]["evidence"]
    screened = rows_for({"nad2": {"prot": "L1", "strand": 1}})
    assert screened and screened[0]["written"] == "yes" and screened[0]["confidence"] == "medium"
    assert screened[0]["intron_start"] == len(ex1) + 1 and "group II splice sites" in screened[0]["note"]
    # the same untyped intron when the model aligns only part of its homolog (< UNTYPED_MIN_COVERAGE): not written
    lib["L2"] = mi.LibProt("some species", "nad2", prot + "W" * 100)
    partial = rows_for({"nad2": {"prot": "L2", "strand": 1}})
    assert partial and partial[0]["written"] == "no" and "of the homolog aligned" in partial[0]["evidence"]


def test_copies_of_a_duplicated_gene_are_not_exons():
    """Hexamermis NC_008828.1 has three identical atp6 copies; exonerate spliced the N-terminal part of one copy to the
    C-terminal part of another (CDS 6888..11167, a 3.7-kb 'intron'). An intron that holds a copy of the model's own
    exons is refused; a real intron (no copy) passes."""
    import random
    sys.path.insert(0, str(SCRIPT.parent))
    import mito_intron as mi
    rnd = random.Random(7)
    gene = "ATG" + "".join(rnd.choice("ACGT") for _ in range(297))
    spacer = "".join(rnd.choice("ACGT") for _ in range(2000))
    s = gene + spacer + gene                               # two identical copies
    second = len(gene) + len(spacer)
    chimera = [(1, 150), (second + 151, second + len(gene))]   # N-terminal half of copy 1 + C-terminal half of copy 2
    assert "copy of the exons" in mi.copied_exons(s, chimera)
    real = "".join(rnd.choice("ACGT") for _ in range(150)) + spacer[:600] + "".join(rnd.choice("ACGT") for _ in range(150))
    assert mi.copied_exons(real, [(1, 150), (751, 900)]) == ""


def test_three_prime_end_trimmed_before_next_gene():
    """Glycera KT989324.1: the first in-frame complete stop of cox1 lay 717 nt inside cox2, the model overlapped COX2
    and was refused. A complete stop more than STOP_OVERLAP_MAX nt inside the next feature gives way to a T / TA right
    before it (mito_curate P3); a stop that overlaps the next feature by less is kept."""
    sys.path.insert(0, str(SCRIPT.parent))
    import mito_intron as mi
    cds = "ATG" + "GCT" * 20 + "GCA" * 10                  # 93 nt, no stop
    s = cds + "T" + "CC" + "CCC" * 15 + "TAA" + "C" * 30     # T right before the next feature (position 95)
    assert mi.extend_to_stop(s, [(1, 63)], 4, 94, 5) == ([(1, 94)], False)
    s2 = cds + "G" + "CC" + "CCC" * 15 + "TAA" + "C" * 30    # no T / TA before it: the complete stop as before
    assert mi.extend_to_stop(s2, [(1, 63)], 4, 94, 5) == ([(1, len(cds) + 3 + 45 + 3)], False)
    s3 = cds + "T" + "CC" + "CCC" * 5 + "TAA" + "C" * 30     # stop 21 nt into the next feature: kept
    assert mi.extend_to_stop(s3, [(1, 63)], 4, 94, 5) == ([(1, len(cds) + 3 + 15 + 3)], False)


def test_trans_spliced_gene_pieces():
    """Placozoan cox1 is trans-spliced: a chain with exons on both strands, or two pieces chained across a long gap
    that leaves homolog residues unaligned, is not one gene; its runs of collinear same-strand alignments are modelled
    as pieces with their own (cis) introns, open where another piece continues the protein (NC_015309.1: 1 + 3 introns
    from the two + strand pieces; NC_008151.3: the + strand half). A piece without an intron is left as it is."""
    sys.path.insert(0, str(SCRIPT.parent))
    import mito_intron as mi
    a = mi.Report(0, 240, 1064, 2223, 1216, 1, [(1065, 1451), (1891, 2223)])
    b = mi.Report(275, 534, 21560, 26157, 1222, 1, [(21561, 21606), (22323, 22417), (24943, 25110), (25690, 26157)])
    m, why = mi.build_model("p", "cox1", [a, b], 534, circular=True)
    assert m is None and mi.TRANS_SPLICED in why
    pieces = mi.gene_modules("p", "cox1", [a, b], 534, circular=True)
    assert [p.exons for p in pieces] == [a.exons, b.exons]
    assert [(p.open5, p.open3) for p in pieces] == [(False, True), (True, False)]
    plus = mi.Report(0, 200, 0, 1500, 900, 1, [(1, 300), (1000, 1500)])
    minus = mi.Report(200, 535, 20000, 18000, 900, -1, [(18001, 19000)])
    m, why = mi.build_model("p", "cox1", [plus, minus], 535)
    assert m is None and mi.TRANS_SPLICED in why
    pieces = mi.gene_modules("p", "cox1", [plus, minus], 535)
    assert [p.exons for p in pieces] == [[(1, 300), (1000, 1500)]] and (pieces[0].open5, pieces[0].open3) == (False, True)


def test_short_repeat_of_exon_sequence_is_no_gene_copy():
    """True hexacoral nad5 introns can hold one 33-46-nt repeat of exon sequence (NC_047221.1 14 shared 20-mers,
    NC_082103.1 27); a duplicated gene copy shares hundreds (Hexamermis atp6 340): COPY_MIN_SHARED lies between."""
    import random
    sys.path.insert(0, str(SCRIPT.parent))
    import mito_intron as mi
    rnd = random.Random(11)
    ex1, ex2 = ("".join(rnd.choice("ACGT") for _ in range(n)) for n in (700, 1000))
    intron = "".join(rnd.choice("ACGT") for _ in range(1000)) + ex1[100:133] + "".join(rnd.choice("ACGT") for _ in range(1000))
    s = ex1 + intron + ex2
    assert mi.copied_exons(s, [(1, 700), (701 + len(intron), len(s))]) == ""


def test_median_homolog_gain_and_origin_overlap():
    """aa_score judges a shift by the median gain over the homologs (one distant homolog no longer vetoes it);
    chains across the origin may share WRAP_OVERLAP % of the protein (NC_040169.1 cox1: both pieces ran 37 residues
    into the intron, as exonerate cannot splice across the origin) - linear chains keep MFannot's 5 %."""
    sys.path.insert(0, str(SCRIPT.parent))
    import mito_intron as mi
    s = "ATG" + "GCT" * 9 + "AAA" * 10 + "TGG" * 12 + "TAA"
    old, new = [(1, 30), (61, len(s))], [(1, 33), (64, len(s))]           # ...A|W... vs ...AK|W...: another protein
    calls = iter([[10.0, 10.0, 10.0], [6.0, 13.0, 15.0]])                   # gains -4, +3, +5: median +3
    assert mi.aa_score(s, old, new, 4, lambda p: next(calls)) == 3.0
    calls = iter([[10.0, 10.0, 10.0], [6.0, 6.0, 15.0]])                    # gains -4, -4, +5: refused
    assert mi.aa_score(s, old, new, 4, lambda p: next(calls)) is None
    calls = iter([[10.0, 10.0, 10.0], [6.0, 6.0, 15.0]])                    # ... unless within a tolerance
    assert mi.aa_score(s, old, new, 4, lambda p: next(calls), tolerance=5) == -4.0
    L = 17697
    exon2 = mi.Report(222, 518, 24, 900, 1191, 1, [(25, 900)])           # rotated record: exon 2 at the start
    exon1 = mi.Report(2, 259, 15984, 16755, 1113, 1, [(15985, 16755)])   # exon 1 at the end, 37 residues shared
    assert [len(g) for g in mi.group_reports([exon2, exon1], 520, L, circular=True)] == [2]
    assert [len(g) for g in mi.group_reports([exon2, exon1], 520, L, circular=False)] == [1, 1]


def test_independent_junction_refinement():
    """refine_junction moves the two intron ends independently (frame kept) to the placement whose protein aligns
    best to the homologs, when the gain is at least REFINE_MIN_GAIN: exonerate put the 3' end of hexacoral nad5
    introns 11-12 nt upstream (NC_015642.1, NC_065774.1), which a same shift of both ends cannot correct."""
    import random
    from Bio.Seq import Seq
    sys.path.insert(0, str(SCRIPT.parent))
    import mito_intron as mi
    rnd = random.Random(5)
    codons = ["GCT", "GAA", "AAA", "CCA", "TGG", "TTT", "CAT", "ATG"]
    ex1 = "ATG" + "".join(rnd.choice(codons) for _ in range(60))
    ex2 = "".join(rnd.choice(codons) for _ in range(60)) + "TAA"
    junk = "GGCAATCCT"                                     # 3 codons the homolog lacks, at the start of exon 2
    intron = "ATTTAAATTA" * 30
    s = ex1 + intron + junk + ex2
    homolog = str(Seq(ex1 + ex2).translate(table=4)).rstrip("*")
    a = len(ex1)
    exons = [(1, a), (a + len(intron) + 1, len(s))]          # exonerate's model keeps the 3 extra codons
    fixed, gain = mi.refine_junction(s, exons, 0, 4, mi.homolog_scorer(homolog))
    assert fixed == [(1, a), (a + len(intron) + len(junk) + 1, len(s))] and gain >= mi.REFINE_MIN_GAIN
    same, gain = mi.refine_junction(s, fixed, 0, 4, mi.homolog_scorer(homolog))
    assert same == fixed and gain == 0.0


def test_group_ii_placement_with_small_homolog_loss_and_wide_window():
    """A typed group II intron with no GUGYG / AY placement within 6 nt that keeps the homolog score: one that loses at
    most II_AA_TOLERANCE is taken (Decemunciger nad1, GUGCG 6 nt upstream, -4), else within WIDE_II nt (Decemunciger
    nad4: GUGYG...AY 11 nt downstream of exonerate's junction)."""
    sys.path.insert(0, str(SCRIPT.parent))
    import mito_intron as mi
    ex1, ex2 = "ATG" + "GCT" * 40 + "CCAGCAGCTGCAGCC", "TGGTGGTGG" + "GCT" * 40 + "TAA"
    intron = "GTGCG" + "ATTTAAATTA" * 20 + "TAC"
    s = ex1 + intron + ex2
    a = len(ex1)
    exons = [(1, a - 11), (a + len(intron) - 10, len(s))]   # exonerate's junction 11 nt upstream of the splice sites
    adj, notes, scores = mi.adjust_introns(s, exons, ["II"], 4, lambda p: -1.0 if p.count("A") < 40 else 0.0)
    assert adj == [(1, a), (a + len(intron) + 1, len(s))] and "adjusted by +11 nt" in notes[0]
    assert scores[0] >= mi.TYPE_II_MIN


def test_weak_group_i_subtype_hit_yields_to_group_ii_sites(monkeypatch, tmp_path):
    """Endomyzostoma OR637230.1 cox1: the only covariance-model hit is a group I subtype (IA1, E=0.0013) but the
    intron has near-perfect group II splice sites (GUGUG...AC; NCBI's and the gold placement are group II). A subtype
    hit weaker than WEAK_SUBTYPE_E then gives no type and the group II rule places the intron."""
    import argparse
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    sys.path.insert(0, str(SCRIPT.parent))
    import mito_intron as mi
    ex1, ex2 = "ATG" + "GCT" * 30 + "GA", "T" + "GCT" * 30 + "TAA"
    intron = "GTGTG" + "ATTTAAATTA" * 20 + "TAC"
    s = ex1 + intron + ex2
    prot = str(Seq(ex1 + ex2).translate(table=5)).rstrip("*")
    rep = mi.Report(0, len(prot), 0, len(s) - 3, 300, 1, [(1, len(ex1) - 1), (len(ex1) + len(intron), len(s) - 3)])
    monkeypatch.setattr(mi, "run_exonerate", lambda *x: [rep])
    monkeypatch.setattr(mi, "core_references", lambda *x: {})
    monkeypatch.setattr(mi, "prescreen", lambda *x: {"cox1": {"prot": "L1", "strand": 1}})
    monkeypatch.setattr(mi, "intron_types", lambda *x: {"m0i0": ("I", 0.0013, "IA1")})
    a = argparse.Namespace(genetic_code=5, exclude_own_species=False, refdir=None, threads=1, maxintron=5000,
                           cm_files=[], evalue=0.01, stop_rule="first")
    rec = SeqRecord(Seq(s), id="x", annotations={"topology": "linear", "molecule_type": "DNA"})
    rows, _ = mi.process_record(rec, {"L1": mi.LibProt("sp", "cox1", prot)}, tmp_path / "db", a,
                                {"exonerate": "e", "cmsearch": "c", "blastp": "b"}, tmp_path)
    row = [r for r in rows if r.get("intron_start")][0]
    assert row["intron_start"] == len(ex1) + 1 and row["intron_end"] == len(ex1) + len(intron)
    assert "type not identified" in row["note"] and "group II splice sites" in row["note"]
