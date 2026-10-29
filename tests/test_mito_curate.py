"""Tests for mito_curate.py on real MitoZ output (Aphrodita KY753833; Lepidonotus clava OW387151).

Rules: Fourdrilis et al. 2018 criteria as summarised by Ghiselli et al. 2021 (doi:10.1098/rstb.2020.0159).
Run: python3 -m pytest -q tests/test_mito_curate.py
"""
import csv
import subprocess
import sys
from pathlib import Path

from Bio import SeqIO

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
SCRIPT = HERE.parent / "mito_curate.py"


def curate(tmp_path, fixture, *extra):
    (tmp_path / "in").mkdir(exist_ok=True)
    (tmp_path / "in" / fixture).write_bytes((DATA / fixture).read_bytes())
    subprocess.run([sys.executable, "-W", "ignore", str(SCRIPT), "-i", str(tmp_path / "in"), "-o",
                    str(tmp_path / "out"), "--genetic_code", "5", *extra], check=True, capture_output=True)
    rec = SeqIO.read(str(tmp_path / "out" / fixture), "genbank")
    rows = list(csv.DictReader((tmp_path / "out" / "curation_report.tsv").open(), delimiter="\t"))
    return rec, rows


def feat(rec, ftype, start1):
    return next(f for f in rec.features if f.type == ftype and int(f.location.start) + 1 == start1)


def test_rrna_trimmed_not_overlapping_trnas(tmp_path):
    rec, _ = curate(tmp_path, "KY753833_mitoz.gbf")
    rrna = sorted((int(f.location.start) + 1, int(f.location.end), f.qualifiers["gene"][0], f.qualifiers["product"][0])
                  for f in rec.features if f.type == "rRNA")
    assert rrna == [(9833, 10643, "rrnS", "12S ribosomal RNA"), (10708, 11997, "rrnL", "16S ribosomal RNA")]
    trnas = [(int(f.location.start) + 1, int(f.location.end)) for f in rec.features if f.type == "tRNA"]
    for s, e, *_ in rrna:
        assert not any(ts <= e and te >= s for ts, te in trnas)          # criterion vi: no overlap


def test_rrna_fill_extends_into_free_space_only(tmp_path):
    rec, _ = curate(tmp_path, "KY753833_mitoz.gbf", "--rrna_bounds", "fill")
    rrna = {f.qualifiers["gene"][0]: (int(f.location.start) + 1, int(f.location.end)) for f in rec.features if f.type == "rRNA"}
    assert rrna["rrnS"] == (9833, 10643) and rrna["rrnL"] == (10708, 11997)  # already flank-to-flank here


def test_cox1_abbreviated_stop_before_trna(tmp_path):
    """MitoZ cox1 1..1543 has no stop; the first complete stop lies inside trnN (1557..): the gene
    ends with TA right before trnN (NCBI KY753833 likewise ends COX1 with TA before tRNA-Asn)."""
    rec, rows = curate(tmp_path, "KY753833_mitoz.gbf")
    cox1 = feat(rec, "CDS", 1)
    assert int(cox1.location.end) == 1556
    assert cox1.qualifiers["transl_except"] == ["(pos:1555..1556,aa:TERM)"]
    assert any(r["rule"] == "P3" and r["gene"].upper() == "COX1" for r in rows)


def test_different_frame_overlap_keeps_start(tmp_path):
    """7994..8281 and 8275..9642 overlap in different frames (criterion ii allows it): start kept."""
    rec, _ = curate(tmp_path, "KY753833_mitoz.gbf")
    assert int(feat(rec, "CDS", 8275).location.start) + 1 == 8275


def test_internal_stops_reported_not_changed(tmp_path):
    before = SeqIO.read(str(DATA / "OW387151_mitoz.gbf"), "genbank")
    rec, rows = curate(tmp_path, "OW387151_mitoz.gbf")
    bad = [r for r in rows if r["rule"] == "P4" and "internal stop" in r["reason"]]
    assert bad and bad[0]["gene"].upper() == "COX1"
    loc = {str(f.location) for f in before.features if f.type == "CDS"}
    assert any(str(f.location) in loc for f in rec.features if f.type == "CDS" and f.qualifiers["gene"][0].upper() == "COX1")


def test_ncbi_validation_before_after(tmp_path):
    """Built-in NCBI table2asn check before vs after curation: errors fixed, none introduced."""
    import shutil
    import mito_env
    t2a = mito_env.find("table2asn")
    if not t2a:
        import pytest
        pytest.skip("table2asn not available")
    curate(tmp_path, "KY753833_mitoz.gbf", "--table2asn", t2a)
    rows = list(csv.DictReader((tmp_path / "out" / "validation_compare.tsv").open(), delimiter="\t"))
    err = [r for r in rows if r["severity"] == "ERROR"]
    assert any(r["status"] == "fixed" for r in err)
    assert not any(r["status"] == "new" for r in err)


def _spacer_record(tmp_path, spacer):
    """+ strand gene ATG(GCC)x30 + T, then `spacer`, then a 66-nt tRNA whose frame holds the first TAA."""
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    gene = "ATG" + "GCC" * 30 + "T"                                  # 94 nt, abbreviated stop T
    k = -(len(gene) + len(spacer) + 30) % 3                          # put the TAA in the gene's frame
    trna = ("GGCCGC" * 5 + "G" * k + "TAA" + "GCGGCC" * 5 + "CCC")[:66]
    pad = "CCCGGG" * 20
    seq = gene + spacer + trna + pad
    t0 = len(gene) + len(spacer)
    frame_stop = next(i for i in range(0, len(seq) - 2, 3) if seq[i:i + 3] in ("TAA", "TAG"))
    rec = SeqRecord(Seq(seq), id="X", description="")
    rec.annotations.update(molecule_type="DNA", topology="linear")
    cds = FeatureLocation(0, frame_stop + 3, 1)                      # MitoZ: greedy stop inside the tRNA
    rec.features = [SeqFeature(FeatureLocation(0, len(seq), 1), type="source"),
                    SeqFeature(cds, type="gene", qualifiers={"gene": ["COX3"]}),
                    SeqFeature(cds, type="CDS", qualifiers={"gene": ["COX3"], "codon_start": ["1"], "transl_table": ["5"]}),
                    SeqFeature(FeatureLocation(t0, t0 + len(trna), 1), type="tRNA",
                               qualifiers={"gene": ["trnS(uga)"], "product": ["tRNA-Ser"]})]
    d = tmp_path / "in"
    d.mkdir()
    SeqIO.write(rec, str(d / "x.gbf"), "genbank")
    subprocess.run([sys.executable, "-W", "ignore", str(SCRIPT), "-i", str(d), "-o", str(tmp_path / "out"),
                    "--genetic_code", "5", "--no_validate"], check=True, capture_output=True)
    out = SeqIO.read(str(tmp_path / "out" / "x.gbf"), "genbank")
    return next(f for f in out.features if f.type == "CDS")


def test_abbreviated_stop_with_short_spacer(tmp_path):
    """2026-10-20 dev benchmark: NCBI ends incomplete-stop genes 0 (143x), -1 (9x), 1-3 nt (10x) before the next
    same-strand tRNA; P3 only accepted 0, so e.g. NC_018362.1 nad2 (T, 1 nt, trnE) kept MitoZ's stop inside the tRNA."""
    cds = _spacer_record(tmp_path, "CC")
    assert int(cds.location.end) == 94 and cds.qualifiers["transl_except"] == ["(pos:94,aa:TERM)"]


def test_abbreviated_stop_no_spacer_unchanged(tmp_path):
    assert int(_spacer_record(tmp_path, "").location.end) == 94


def test_rrna_never_contains_a_trna(tmp_path):
    """2026-10-21 table2asn on dev outputs: 16x SEQ_FEAT.BadRRNAcomponentOverlapTRNA, e.g. NC_028229.1 rrnS 3368..4148
    containing trnE 3373..3429 - R1 only trimmed neighbours overlapping an rRNA end. The rRNA keeps the longer side."""
    import random
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    rnd = random.Random(2)
    rec = SeqRecord(Seq("".join(rnd.choice("ACGT") for _ in range(2000))), id="X", description="")
    rec.annotations.update(molecule_type="DNA", topology="linear")
    def add(kind, s, e, gene, product):
        loc = FeatureLocation(s - 1, e, 1)
        rec.features += [SeqFeature(loc, type="gene", qualifiers={"gene": [gene]}),
                         SeqFeature(loc, type=kind, qualifiers={"gene": [gene], "product": [product]})]
    rec.features = [SeqFeature(FeatureLocation(0, 2000, 1), type="source")]
    add("tRNA", 100, 160, "trnA(ugc)", "tRNA-Ala")
    add("rRNA", 150, 1000, "s-rRNA", "12S ribosomal RNA")          # overlaps trnA end, contains trnE
    add("tRNA", 170, 230, "trnE(uuc)", "tRNA-Glu")
    add("tRNA", 1001, 1060, "trnS(uga)", "tRNA-Ser")
    d = tmp_path / "in"
    d.mkdir()
    SeqIO.write(rec, str(d / "x.gbf"), "genbank")
    subprocess.run([sys.executable, "-W", "ignore", str(SCRIPT), "-i", str(d), "-o", str(tmp_path / "out"),
                    "--genetic_code", "5", "--no_validate"], check=True, capture_output=True)
    out = SeqIO.read(str(tmp_path / "out" / "x.gbf"), "genbank")
    r = next(f for f in out.features if f.type == "rRNA")
    assert (int(r.location.start) + 1, int(r.location.end)) == (231, 1000)
    assert r.qualifiers["gene"] == ["rrnS"]


def test_ncbi_gene_names(tmp_path):
    """User 2026-10-21: names consistent and NCBI-conform. /gene is the plain symbol (NCBI genome annotation guide:
    tRNA labels go to /note); Leu/Ser by anticodon (L1 CUN, L2 UUR, S1 AGN, S2 UCN); rRNA rrnS/rrnL also when the
    feature is a join (3 dev records kept MitoZ's l-rRNA / s-rRNA); CDS synonyms -> COX1..CYTB, ND1..ND6, ND4L."""
    import random
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    rnd = random.Random(3)
    rec = SeqRecord(Seq("".join(rnd.choice("ACGT") for _ in range(3000))), id="X", description="")
    rec.annotations.update(molecule_type="DNA", topology="circular")
    rec.features = [SeqFeature(FeatureLocation(0, 3000, 1), type="source")]
    def add(kind, loc, gene, **q):
        rec.features += [SeqFeature(loc, type="gene", qualifiers={"gene": [gene]}),
                         SeqFeature(loc, type=kind, qualifiers={"gene": [gene], **{k: [v] for k, v in q.items()}})]
    add("tRNA", FeatureLocation(100, 165, 1), "trnL(uaa)", product="tRNA-Leu")
    add("tRNA", FeatureLocation(200, 265, -1), "trnS(gcu)", product="tRNA-Ser")
    add("tRNA", FeatureLocation(300, 365, 1), "trnM(cau)", product="tRNA-Met")
    add("rRNA", CompoundLocation([FeatureLocation(2500, 3000, 1), FeatureLocation(0, 90, 1)]), "l-rRNA",
        product="16S ribosomal RNA")
    add("CDS", FeatureLocation(400, 1000, 1), "cob", product="Cytb", codon_start="1", transl_table="5")
    d = tmp_path / "in"
    d.mkdir()
    SeqIO.write(rec, str(d / "x.gbf"), "genbank")
    subprocess.run([sys.executable, "-W", "ignore", str(SCRIPT), "-i", str(d), "-o", str(tmp_path / "out"),
                    "--genetic_code", "5", "--no_validate"], check=True, capture_output=True)
    out = SeqIO.read(str(tmp_path / "out" / "x.gbf"), "genbank")
    got = [(f.type, f.qualifiers["gene"][0], f.qualifiers.get("product", [""])[0], f.qualifiers.get("note", [""])[0])
           for f in out.features if f.type in ("tRNA", "rRNA", "CDS")]
    assert ("tRNA", "trnL2", "tRNA-Leu", "anticodon: UAA") in got
    assert ("tRNA", "trnS1", "tRNA-Ser", "anticodon: GCU") in got
    assert ("tRNA", "trnM", "tRNA-Met", "anticodon: CAU") in got
    assert ("rRNA", "rrnL", "16S ribosomal RNA", "") in got
    assert any(x[:3] == ("CDS", "CYTB", "cytochrome b") for x in got)
    genes = sorted(f.qualifiers["gene"][0] for f in out.features if f.type == "gene")
    assert genes == ["CYTB", "rrnL", "trnL2", "trnM", "trnS1"]


def _one_cds(tmp_path, nt, code, gene_loc=None):
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    seq = "CCCGGG" * 10 + nt + "CCCGGG" * 10
    rec = SeqRecord(Seq(seq), id="X", description="")
    rec.annotations.update(molecule_type="DNA", topology="linear")
    loc = FeatureLocation(60, 60 + len(nt), 1)
    rec.features = [SeqFeature(FeatureLocation(0, len(seq), 1), type="source"),
                    SeqFeature(gene_loc or loc, type="gene", qualifiers={"gene": ["ND4L"]}),
                    SeqFeature(loc, type="CDS", qualifiers={"gene": ["ND4L"], "codon_start": ["1"], "transl_table": [str(code)]})]
    d = tmp_path / "in"
    d.mkdir()
    SeqIO.write(rec, str(d / "x.gbf"), "genbank")
    subprocess.run([sys.executable, "-W", "ignore", str(SCRIPT), "-i", str(d), "-o", str(tmp_path / "out"),
                    "--genetic_code", str(code), "--no_validate"], check=True, capture_output=True)
    return SeqIO.read(str(tmp_path / "out" / "x.gbf"), "genbank")


def test_illegal_start_marked_partial(tmp_path):
    """2026-10-21 table2asn: SEQ_FEAT.StartCodon ERROR on echinoderm nad4l (ATT/ATC, code 9 starts are ATG/GTG only);
    NCBI writes the same CDS as <start with note 'start codon not determined' (NC_001627.1, NC_061414.1)."""
    out = _one_cds(tmp_path, "ATT" + "GCC" * 40 + "TAA", 9)
    cds = next(f for f in out.features if f.type == "CDS")
    assert str(cds.location).startswith("[<60") and "start codon not determined" in cds.qualifiers["note"]
    gene = next(f for f in out.features if f.type == "gene")
    assert str(gene.location) == str(cds.location)


def test_gene_feature_matches_cds(tmp_path):
    """table2asn SEQ_FEAT.CDSgeneRange: gene features must span exactly their CDS / RNA."""
    from Bio.SeqFeature import FeatureLocation
    out = _one_cds(tmp_path, "ATG" + "GCC" * 40 + "TAA", 5, gene_loc=FeatureLocation(55, 150, 1))
    cds, gene = (next(f for f in out.features if f.type == t) for t in ("CDS", "gene"))
    assert (int(gene.location.start), int(gene.location.end)) == (int(cds.location.start), int(cds.location.end))


def test_circular_defline_complete():
    """table2asn SEQ_INST.CompleteCircleProblem on every circular record: the .fsa needs [completeness=complete]."""
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    sys.path.insert(0, str(HERE.parent))
    import ncbi_validate as nv
    r = SeqRecord(Seq("ACGT"), id="X")
    r.annotations.update(organism="Aus bus", topology="circular")
    assert "[completeness=complete]" in nv.fasta_defline(r, 5)
    r.annotations["topology"] = "linear"
    assert "completeness" not in nv.fasta_defline(r, 5)


def test_missing_stop_marked_partial(tmp_path):
    """2026-10-21 dev benchmark: SEQ_FEAT.NoStop ERROR (e.g. NC_044095.1 ND6, NC_085563.1 ND6) when a CDS has neither
    a complete nor an abbreviated (T/TA) stop; written 3'-partial with 'stop codon not determined' like P5."""
    out = _one_cds(tmp_path, "ATG" + "GCC" * 40, 5)
    cds = next(f for f in out.features if f.type == "CDS")
    assert str(cds.location).endswith(">183](+)") and "stop codon not determined" in cds.qualifiers["note"]


def _origin_record(tmp_path, strand):
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    body = "CCCGGG" * 20
    gene = "ATG" + "GCC" * 21 + "TAA"                      # 69 nt; the last 4 nt lie past the origin
    plus = gene[-4:] + body + gene[:-4]
    seq = plus if strand == 1 else str(Seq(plus).reverse_complement())
    L = len(seq)
    loc = FeatureLocation(L - 65, L - 2, 1) if strand == 1 else FeatureLocation(2, 65, -1)
    rec = SeqRecord(Seq(seq), id="X", description="")
    rec.annotations.update(molecule_type="DNA", topology="circular")
    rec.features = [SeqFeature(FeatureLocation(0, L, 1), type="source"),
                    SeqFeature(loc, type="gene", qualifiers={"gene": ["ATP8"]}),
                    SeqFeature(loc, type="CDS", qualifiers={"gene": ["ATP8"], "codon_start": ["1"], "transl_table": ["5"]})]
    d = tmp_path / "in"
    d.mkdir(parents=True)
    SeqIO.write(rec, str(d / "x.gbf"), "genbank")
    subprocess.run([sys.executable, "-W", "ignore", str(SCRIPT), "-i", str(d), "-o", str(tmp_path / "out"),
                    "--genetic_code", "5", "--no_validate"], check=True, capture_output=True)
    return SeqIO.read(str(tmp_path / "out" / "x.gbf"), "genbank"), gene


def test_cds_extended_across_origin(tmp_path):
    """2026-10-21 dev benchmark NC_006902.1 / NC_022829.1 (ATP8 overlaps COX1, cut by the cox1-first rotation),
    NC_025335.1 (ND4L): CDS stops 2 nt before the end of a circular record; its stop lies past the origin."""
    for strand in (1, -1):
        out, gene = _origin_record(tmp_path / str(strand), strand)
        cds = next(f for f in out.features if f.type == "CDS")
        assert len(cds.location.parts) == 2 and str(cds.extract(out.seq)) == gene, strand
        g = next(f for f in out.features if f.type == "gene")
        assert str(g.location) == str(cds.location)


def test_tbl_minus_strand_join_in_biological_order():
    """2026-10-21: NC_035061.1 trnY complement(join(16132..16186,1..8)) was written to the .tbl as 16186-16132
    then 8-1; a feature table lists intervals in transcript order (8-1, then 16186-16132). Biopython already
    keeps minus-strand parts in transcript order."""
    import io
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    sys.path.insert(0, str(HERE.parent))
    import ncbi_validate as nv
    r = SeqRecord(Seq("A" * 100), id="X")
    r.annotations["molecule_type"] = "DNA"
    r.features = [SeqFeature(CompoundLocation([FeatureLocation(0, 10, -1), FeatureLocation(90, 100, -1)]), type="tRNA")]
    h = io.StringIO()
    SeqIO.write(r, h, "genbank")
    assert "complement(join(91..100,1..10))" in h.getvalue()
    assert nv.feature_lines(r)[1:3] == ["10\t1\ttRNA", "100\t91"]


# ── in-process helpers (review C fixes) ──
def _mc():
    sys.path.insert(0, str(HERE.parent))
    import mito_curate
    return mito_curate


def _args(**kw):
    import argparse
    d = dict(genetic_code=5, rrna_bounds="fill", max_fill=300, max_start_shift=30, aggressive_start=False,
             aggressive_start_genes="")
    d.update(kw)
    return argparse.Namespace(**d)


def _mkrec(seq, feats, topology="linear"):
    """feats: (type, location, gene[, extra qualifiers]); a matching gene feature is added for CDS / RNA."""
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    rec = SeqRecord(Seq(seq), id="X", description="")
    rec.annotations.update(molecule_type="DNA", topology=topology)
    rec.features = [SeqFeature(FeatureLocation(0, len(seq), 1), type="source")]
    for ftype, loc, gene, *extra in feats:
        q = {"gene": [gene], **(extra[0] if extra else {})}
        if ftype == "CDS":
            q.setdefault("codon_start", ["1"])
            q.setdefault("transl_table", ["5"])
        rec.features += [SeqFeature(loc, type="gene", qualifiers={"gene": [gene]}), SeqFeature(loc, type=ftype, qualifiers=q)]
    return rec


def _curate_rec(rec, **kw):
    a = _args(**kw)
    rows = _mc().Curator(rec, a.genetic_code, a).run()
    return rows


def _by(rec, ftype, gene):
    return next(f for f in rec.features if f.type == ftype and f.qualifiers.get("gene", [""])[0] == gene)


def test_wrapped_trna_does_not_block_p1():
    """Review C finding 6: span() is the hull of a compound location, so an origin-spanning tRNA
    join(L-29..L,1..35) counted as a preceding feature ending at L and disabled P1 for every same-strand CDS."""
    from Bio.SeqFeature import CompoundLocation, FeatureLocation
    L = 1200
    seq = "GCC" * 100 + "ATG" + "GCC" + "GCC" * 40 + "TAA" + "GCC" * 257        # ATG 301..303, annotated start 307 = GCC
    assert len(seq) == L
    cds = FeatureLocation(306, 429, 1)
    trna = CompoundLocation([FeatureLocation(L - 30, L, 1), FeatureLocation(0, 35, 1)])
    rec = _mkrec(seq, [("tRNA", trna, "trnY", {"product": ["tRNA-Tyr"]}), ("CDS", cds, "ND1")], "circular")
    rows = _curate_rec(rec)
    assert any(r["rule"] == "P1" for r in rows), rows
    assert int(_by(rec, "CDS", "ND1").location.start) + 1 == 301


def _rrna_rec(feats, n=3000):
    import random
    rnd = random.Random(11)
    seq = "".join(rnd.choice("ACGT") for _ in range(n))
    from Bio.SeqFeature import FeatureLocation
    out = []
    for kind, s, e, gene in feats:
        out.append((kind, FeatureLocation(s - 1, e, 1), gene, {"product": [{"tRNA": "tRNA-Gln", "rRNA": "16S ribosomal RNA"}[kind]]}))
    return _mkrec(seq, out)


def test_rrna_not_trimmed_for_inside_gene_when_loss_large():
    """Review C finding 7: a spurious tRNA call inside a correctly annotated 1.3-kb l-rRNA cut the rRNA to 830 nt
    (36 % lost). Trim only when the discarded side is <= 10 % of the rRNA; else keep it and log R3."""
    rec = _rrna_rec([("rRNA", 1001, 2300, "l-rRNA"), ("tRNA", 1401, 1470, "trnQ(uug)")])
    rows = _curate_rec(rec)
    r = _by(rec, "rRNA", "rrnL")
    assert (int(r.location.start) + 1, int(r.location.end)) == (1001, 2300)
    r3 = [x for x in rows if x["rule"] == "R3" and "inside rRNA" in x["reason"]]
    assert len(r3) == 1 and "trnQ" in r3[0]["reason"] and "not trimmed" in r3[0]["reason"], rows
    assert not any(x["rule"] == "R1" for x in rows)


def test_rrna_trimmed_for_inside_gene_near_the_end():
    rec = _rrna_rec([("rRNA", 1001, 2300, "l-rRNA"), ("tRNA", 2201, 2270, "trnQ(uug)")])
    rows = _curate_rec(rec)
    r = _by(rec, "rRNA", "rrnL")
    assert (int(r.location.start) + 1, int(r.location.end)) == (1001, 2200)
    assert any(x["rule"] == "R1" for x in rows)


def test_translation_kept_for_every_cds(tmp_path):
    """Review C finding 9: pcg() popped /translation from every simple CDS and main() now always rewrites the file, so
    KY753833 went from 13 /translation qualifiers to 0 (changed or not). Recomputed from the final location."""
    before = SeqIO.read(str(DATA / "KY753833_mitoz.gbf"), "genbank")
    rec, _ = curate(tmp_path, "KY753833_mitoz.gbf")
    n_in = sum("translation" in f.qualifiers for f in before.features if f.type == "CDS")
    cds = [f for f in rec.features if f.type == "CDS"]
    assert n_in == 13 and sum("translation" in f.qualifiers for f in cds) == 13
    for f in cds:
        prot = f.qualifiers["translation"][0]
        assert prot.startswith("M") and "*" not in prot, f.qualifiers["gene"]
        assert prot == str(f.extract(rec.seq)[: len(f.extract(rec.seq)) // 3 * 3].translate(table=5)).rstrip("*"), \
            f.qualifiers["gene"]


def test_translation_alternative_start_is_methionine():
    """NCBI convention: a (non-partial) CDS beginning with a start codon of the table (ATT under code 5) is
    translated with M; checked on a + CDS, a - CDS and a join."""
    from Bio.Seq import Seq
    from Bio.SeqFeature import CompoundLocation, FeatureLocation
    gene = "ATT" + "GCC" * 40 + "TAA"                                       # 126 nt
    pad = "CCCGGG" * 10
    ex1, ex2 = "ATT" + "GCC" * 20, "GCC" * 20 + "TAA"
    seq = pad + gene + pad + str(Seq(gene).reverse_complement()) + pad + ex1 + "A" * 20 + ex2 + pad
    j = 60 + 126 + 60 + 126 + 60                                             # exon 1 starts here (0-based)
    rec = _mkrec(seq, [("CDS", FeatureLocation(60, 186, 1), "ND1"),
                       ("CDS", FeatureLocation(246, 372, -1), "ND2"),
                       ("CDS", CompoundLocation([FeatureLocation(j, j + 63, 1), FeatureLocation(j + 83, j + 146, 1)]), "ND3")])
    _curate_rec(rec)
    for gene_name in ("ND1", "ND2", "ND3"):
        assert _by(rec, "CDS", gene_name).qualifiers["translation"] == ["M" + "A" * 40], gene_name


def test_translation_not_written_with_internal_stop_and_partial_start_plain():
    """An internal stop remains -> no /translation (P4 reports it); a 5'-partial CDS keeps the plain translation (I)."""
    from Bio.SeqFeature import BeforePosition, FeatureLocation
    nt = "ATG" + "GCC" * 10 + "TAA" + "GCC" * 100 + "TAA"
    rec = _mkrec("CCCGGG" * 10 + nt + "CCCGGG" * 10, [("CDS", FeatureLocation(60, 60 + len(nt), 1), "ND1")])
    rows = _curate_rec(rec)
    assert "translation" not in _by(rec, "CDS", "ND1").qualifiers
    assert any(r["rule"] == "P4" and "internal stop" in r["reason"] for r in rows)
    nt = "ATT" + "GCC" * 20 + "TAA"
    rec = _mkrec("CCCGGG" * 10 + nt + "CCCGGG" * 10, [("CDS", FeatureLocation(BeforePosition(60), 60 + len(nt), 1), "ND1")])
    _curate_rec(rec)
    assert _by(rec, "CDS", "ND1").qualifiers["translation"] == ["I" + "A" * 20]


def test_p6_stop_tested_in_frame_only():
    """Review C finding 10: an out-of-frame TAA at the end of a CDS whose length is not a multiple of 3 suppressed P6
    (table2asn then gives SEQ_FEAT.NoStop). Stop test only when length % 3 == 0; otherwise T / TA tail or no stop."""
    from Bio.SeqFeature import FeatureLocation
    nt = "ATG" + "GCC" * 40 + "ATAA"                                       # 127 nt, no in-frame stop
    rec = _mkrec("CCCGGG" * 10 + nt + "CCCGGG" * 10, [("CDS", FeatureLocation(60, 60 + len(nt), 1), "ND1")])
    rows = _curate_rec(rec)
    assert any(r["rule"] == "P6" for r in rows), rows
    assert any(r["rule"] == "P4" and "not a multiple of 3" in r["reason"] for r in rows), rows
    assert "stop codon not determined" in _by(rec, "CDS", "ND1").qualifiers["note"]


def _cut_cds_record(stop_at, first_feature=True):
    """Circular record; + strand ND4L = ATG + GCC*40 ends at L-1 (no stop); its first in-frame stop lies past the origin,
    TAA at 0-based s[stop_at..] (stop_at % 3 == 2). First feature COX1 covers position 1 (cox1-first rotation)."""
    from Bio.SeqFeature import FeatureLocation
    first = list(("GCC" * 134)[:400])
    first[stop_at:stop_at + 3] = "TAA"
    nd4l = "ATG" + "GCC" * 40
    seq = "".join(first) + "CCCGGG" * 30 + nd4l + "G"
    feats = [("CDS", FeatureLocation(580, 703, 1), "ND4L")]
    if first_feature:
        feats.append(("tRNA", FeatureLocation(0, 400, 1), "trnY", {"product": ["tRNA-Tyr"]}))
    rec = _mkrec(seq, feats, "circular")
    return rec, _curate_rec(rec)


def test_origin_extension_capped_by_feature_at_position_1():
    """Review C finding 11: O1 extended a CDS up to 300 nt past the origin, i.e. up to 300 nt into the feature at
    position 1 (cox1-first rotation); genuine mitochondrial overlaps are <= ~50 nt (dev cases 13-19 nt)."""
    rec, rows = _cut_cds_record(14)                               # stop 17 nt into the first feature: genuine overlap
    assert any(r["rule"] == "O1" for r in rows)
    assert len(_by(rec, "CDS", "ND4L").location.parts) == 2
    rec, rows = _cut_cds_record(149)                              # stop 152 nt into the first feature: not O1
    assert not any(r["rule"] == "O1" for r in rows), rows
    nd4l = _by(rec, "CDS", "ND4L")
    assert len(nd4l.location.parts) == 1 and str(nd4l.location).endswith(">704](+)")   # P6: cut codon's nt included
    assert any(r["rule"] == "P6" for r in rows)
    rec, rows = _cut_cds_record(149, first_feature=False)         # nothing at position 1: the 300 nt limit applies
    assert any(r["rule"] == "O1" for r in rows)


def test_cds_extension_does_not_overlap_following_rrna():
    """Review C finding 12: R1 ran before the PCG loop, so P3 later extended a CDS (no stop at 1453, first stop at 1531)
    into the rRNA that R1 had just cleared (78-nt CDS/rRNA overlap). The PCG loop now runs first."""
    from Bio.SeqFeature import FeatureLocation
    seq = ("CCCGGG" * 167)[:1000] + "ATG" + "GCC" * 150 + "GCC" * 25 + "TAA" + ("CCCGGG" * 300)[:1469]
    assert len(seq) == 3000 and seq[1528:1531] == "TAA"
    rec = _mkrec(seq, [("CDS", FeatureLocation(1000, 1453, 1), "COX1"),
                       ("rRNA", FeatureLocation(1453, 2500, 1), "s-rRNA", {"product": ["12S ribosomal RNA"]})])
    rows = _curate_rec(rec)
    cds, rrna = _by(rec, "CDS", "COX1"), _by(rec, "rRNA", "rrnS")
    assert int(cds.location.end) == 1531 and any(r["rule"] == "P3" for r in rows)
    assert int(rrna.location.start) >= int(cds.location.end)


def test_failed_curation_leaves_record_unchanged(tmp_path, monkeypatch):
    """Review C finding 14: Curator.run() mutates the record in place, so after an exception the row said 'left
    unchanged' but renames / trims done before the failure were still written."""
    mc = _mc()
    (tmp_path / "in").mkdir()
    (tmp_path / "in" / "k.gbf").write_bytes((DATA / "KY753833_mitoz.gbf").read_bytes())

    def boom(self):
        self.names(next(f for f in self.rec.features if f.type == "CDS"))
        self.rec.features[1].qualifiers["gene"] = ["MUTATED"]
        self.rec.features.pop()
        raise RuntimeError("boom")

    monkeypatch.setattr(mc.Curator, "run", boom)
    assert mc.main(["-i", str(tmp_path / "in"), "-o", str(tmp_path / "out"), "--genetic_code", "5", "--no_validate"]) == 0
    before, after = (SeqIO.read(str(tmp_path / d / "k.gbf"), "genbank") for d in ("in", "out"))

    def sig(r):
        return [(f.type, str(f.location), f.qualifiers.get("gene")) for f in r.features]
    assert sig(after) == sig(before)
    rows = list(csv.DictReader((tmp_path / "out" / "curation_report.tsv").open(), delimiter="\t"))
    assert any(r["rule"] == "error" and "left unchanged" in r["reason"] for r in rows)


def test_origin_spanning_feature_sorts_by_its_larger_coordinate_part():
    """Review C finding 20: sort key int(location.start) put join(L-29..L,1..35) (start 0) before every feature."""
    from Bio.SeqFeature import CompoundLocation, FeatureLocation
    seq = "CCCGGG" * 100 + "ATG" + "GCC" * 20 + "TAA" + "CCCGGG" * 100
    L = len(seq)
    trna = CompoundLocation([FeatureLocation(L - 30, L, 1), FeatureLocation(0, 35, 1)])
    nt = 3 + 60 + 3
    rec = _mkrec(seq, [("tRNA", trna, "trnY", {"product": ["tRNA-Tyr"]}), ("CDS", FeatureLocation(600, 600 + nt, 1), "ND1")],
                 "circular")
    _curate_rec(rec)
    order = [(f.type, f.qualifiers.get("gene", [""])[0]) for f in rec.features]
    assert order == [("source", ""), ("gene", "ND1"), ("CDS", "ND1"), ("gene", "trnY"), ("tRNA", "trnY")], order


def test_rrna_trimmed_for_inside_trna_near_the_start_up_to_250_nt():
    """Dev regression after review fix 7 (NC_023929.1 / NC_023931.1): MitoZ's l-rRNA began 199 nt before trnV, which
    therefore counted as 'inside' (loss 199 nt > 10 % of 1370). A tRNA-sized overhang at either end is a normal
    core error and is trimmed: the discarded side may be up to max(10 %, 250 nt)."""
    rec = _rrna_rec([("rRNA", 1001, 2370, "l-rRNA"), ("tRNA", 1011, 1200, "trnV(uac)")])
    rows = _curate_rec(rec)
    r = _by(rec, "rRNA", "rrnL")
    assert (int(r.location.start) + 1, int(r.location.end)) == (1201, 2370)
    assert any(x["rule"] == "R1" for x in rows) and not any(x["rule"] == "R3" for x in rows)


def test_rrna_overhang_is_trimmed_before_p1_pushes_a_cds_start():
    """Dev regression after review fix 12 (NC_060482.1 nad6 exact -> +9 nt, NC_045069.1 nad1 exact -> +30 nt): with the
    CDS loop before R1, the raw rRNA overhang still overlapped the CDS start and P1 moved the start downstream.
    Order is rRNA (R1) -> CDS -> rRNA again (so P3/O1 extensions are trimmed as well)."""
    from Bio.SeqFeature import FeatureLocation
    seq = ("CCCGGG" * 167)[:1000] + ("CCCGGG" * 84)[:500] + "ATG" + "GCC" * 32 + "ATG" + "GCC" * 60 + "TAA" + ("CCCGGG" * 300)[:1215]
    assert len(seq) == 3000 and seq[1500:1503] == "ATG" and seq[1599:1602] == "ATG"   # in-frame ATG at 1501 and 1600
    rec = _mkrec(seq, [("rRNA", FeatureLocation(1000, 1560, 1), "l-rRNA", {"product": ["16S ribosomal RNA"]}),
                       ("CDS", FeatureLocation(1500, 1785, 1), "ND6")])
    rows = _curate_rec(rec)
    cds, rrna = _by(rec, "CDS", "ND6"), _by(rec, "rRNA", "rrnL")
    assert int(cds.location.start) + 1 == 1501 and not any(r["rule"] == "P1" for r in rows), rows
    assert int(rrna.location.end) == 1500 and any(r["rule"] == "R1" for r in rows)


def test_minus_strand_frame_anchored_at_the_genome_end():
    """Dev regression NC_015829.1 (tardigrade) after review fix 6: nad4L/nad4 on the minus strand overlap by 7 nt in
    DIFFERENT frames (NCBI keeps the overlap; the frame of a minus-strand CDS is set by its genome END), but the
    frame test used the genome starts, called them 'same frame', and P1 moved the nad4 start 18 nt downstream.
    No P1 and no P2 for such a pair; a true same-frame overlap (multiple of 3) is still reported."""
    from Bio.SeqFeature import FeatureLocation
    from Bio.Seq import Seq
    nad4 = "ATG" + "GCC" * 5 + "ATG" + "GCC" * 94 + "T"      # 304 nt (T-tail stop, as NC_015829.1 nad4: length % 3 == 1);
    #                                                          in-frame ATG 18 nt downstream = where P1 would move the start
    nad4l = "ATG" + "GCC" * 60 + "TAA"                       # 186 nt
    # genome (minus strand genes): ... [nad4l rc] overlapping the last 7 nt of [nad4 rc]
    tx = nad4l[:-7] + nad4                                   # nad4l stop TAA + 4 nt overlap nad4 start -> 7 nt overlap
    rc = str(Seq(tx).reverse_complement())
    seq = ("CCCGGG" * 200)[:1000] + rc + ("CCCGGG" * 200)[:1000]
    L = len(seq)
    g_nad4_s, g_nad4_e = 1000, 1000 + len(nad4)              # nad4 occupies the 3' part of tx -> low genome coords
    g_nad4l_s, g_nad4l_e = 1000 + len(nad4) - 7, 1000 + len(tx)
    assert g_nad4l_s < g_nad4_e and (g_nad4l_s - g_nad4_s) % 3 == 0        # the OLD (start-based) test says 'same frame'
    assert (g_nad4l_e - g_nad4_e) % 3 != 0                                  # the end-based frame differs
    rec = _mkrec(seq, [("CDS", FeatureLocation(g_nad4_s, g_nad4_e, -1), "ND4"),
                       ("CDS", FeatureLocation(g_nad4l_s, g_nad4l_e, -1), "ND4L")])
    rows = _curate_rec(rec)
    assert int(_by(rec, "CDS", "ND4").location.end) == g_nad4_e, rows
    assert not any(r["rule"] in ("P1", "P2") for r in rows), rows
    assert L > 0


def _with_atp8_copy(diverge: int):
    """NC_000931 with a copy of atp8 at the record end; diverge = codons whose third base is changed in the copy."""
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    rec = SeqIO.read(str(DATA / "NC_000931.gb"), "genbank")
    atp8 = next(f for f in rec.features if f.type == "CDS" and f.qualifiers["gene"][0].upper() == "ATP8")
    s, e = int(atp8.location.start), int(atp8.location.end)
    seq = str(rec.seq)
    copy = list(seq[s:e])
    for k in range(3, 3 + 3 * diverge, 3):                                 # after the start codon, before the stop
        copy[k + 2] = {"A": "C", "C": "A", "G": "T", "T": "G"}[copy[k + 2]]
    rec.seq = type(rec.seq)(seq[:len(seq) - (e - s)] + "".join(copy))
    copy_loc = FeatureLocation(len(rec.seq) - (e - s), len(rec.seq), 1)
    rec.features = [f for f in rec.features if int(f.location.end) <= copy_loc.start or f.type == "source"]
    rec.features.append(SeqFeature(copy_loc, type="gene", qualifiers={"gene": ["ATP8"]}))
    rec.features.append(SeqFeature(copy_loc, type="CDS", qualifiers=dict(atp8.qualifiers)))
    return rec, s


def test_duplicate_pcg_copy_removed():
    """D1: a second, diverged (< 95 % identical) copy of atp8 on a repeated region is removed; one CDS + gene stays."""
    rec, s = _with_atp8_copy(diverge=12)                                   # 12 of ~160 nt: 92.5 % identical
    rows = _curate_rec(rec)
    cds = [f for f in rec.features if f.type == "CDS" and f.qualifiers["gene"][0].upper() == "ATP8"]
    genes = [f for f in rec.features if f.type == "gene" and f.qualifiers["gene"][0].upper() == "ATP8"]
    assert len(cds) == 1 and len(genes) == 1 and int(cds[0].location.start) == s
    assert any(r["rule"] == "D1" for r in rows)


def test_identical_complete_second_copy_kept():
    """D1 (2026-10-29): a complete second copy >= 95 % identical to the gene is a real gene copy and stays, as NCBI
    annotates it (heron / cormorant / turaco nad6, Hexamermis atp6: 30 of 38 identical copies in 1,050 records)."""
    for diverge in (0, 2):
        rec, _ = _with_atp8_copy(diverge)
        rows = _curate_rec(rec)
        assert len([f for f in rec.features if f.type == "CDS" and f.qualifiers["gene"][0].upper() == "ATP8"]) == 2
        assert len([f for f in rec.features if f.type == "gene" and f.qualifiers["gene"][0].upper() == "ATP8"]) == 2
        assert not any(r["rule"] == "D1" for r in rows)


def test_identical_copy_without_stop_yet_kept():
    """D1: the core ends a gene at its last aligned codon (the stop comes from P3 later); two identical copies whose
    stop is not yet annotated are a real gene copy too (NC_044480.1 cob, NC_010636.1 cox2 were removed)."""
    from Bio.SeqFeature import FeatureLocation
    rec, _ = _with_atp8_copy(0)
    for f in [x for x in rec.features if x.type in ("CDS", "gene") and x.qualifiers["gene"][0].upper() == "ATP8"]:
        s, e = int(f.location.start), int(f.location.end)
        f.location = FeatureLocation(s, e - 3, 1)                         # without the stop codon
    rows = _curate_rec(rec)
    assert len([f for f in rec.features if f.type == "CDS" and f.qualifiers["gene"][0].upper() == "ATP8"]) == 2
    assert not any(r["rule"] == "D1" for r in rows)


def test_duplicated_block_copies_kept():
    """D1: a block of two genes repeated (NC_009093.1 mite: 9 genes twice, both copies in NCBI) keeps both copies."""
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    rec = SeqIO.read(str(DATA / "NC_000931.gb"), "genbank")
    cds = {f.qualifiers["gene"][0].upper(): f for f in rec.features if f.type == "CDS"}
    a8, a6 = cds["ATP8"], cds["ATP6"]
    s, e = int(a8.location.start), int(a6.location.end)
    seq = str(rec.seq)
    rec.seq = type(rec.seq)(seq[:len(seq) - (e - s)] + seq[s:e])          # copy of atp8-atp6 at the record end
    off = len(rec.seq) - (e - s) - s
    rec.features = [f for f in rec.features if int(f.location.end) <= len(rec.seq) - (e - s) or f.type == "source"]
    for f in (a8, a6):
        loc = FeatureLocation(int(f.location.start) + off, int(f.location.end) + off, 1)
        rec.features.append(SeqFeature(loc, type="gene", qualifiers={"gene": f.qualifiers["gene"]}))
        rec.features.append(SeqFeature(loc, type="CDS", qualifiers=dict(f.qualifiers)))
    rows = _curate_rec(rec)
    for name in ("ATP8", "ATP6"):
        assert len([f for f in rec.features if f.type == "CDS" and f.qualifiers["gene"][0].upper() == name]) == 2
    assert not any(r["rule"] == "D1" for r in rows)


def test_aggressive_start_genes_only_for_listed_genes():
    """P1 upstream search for nad5 (listed, as `mitom run` sets for annelids) but not cob: both start with a valid
    ATT 9 codons downstream of an in-frame ATA that follows the preceding tRNA. Without the list nothing moves (an
    upstream ATG would be taken for any gene: test_upstream_atg_for_non_atg_start)."""
    from Bio.SeqFeature import FeatureLocation
    body = "ATA" + "GCC" * 8 + "ATT" + "GCC" * 100 + "TAA"
    seq = "C" * 1000 + body + "C" * 1000
    for gene, moved, genes in (("ND5", True, "ND5,ND4L"), ("CYTB", False, "ND5,ND4L"), ("ND5", False, "")):
        rec = _mkrec(seq, [("tRNA", FeatureLocation(900, 965, 1), "trnF", {"product": ["tRNA-Phe"]}),
                           ("CDS", FeatureLocation(1027, 1000 + len(body), 1), gene)])
        _curate_rec(rec, aggressive_start_genes=genes)
        cds = _by(rec, "CDS", gene)
        assert (int(cds.location.start) == 1000) == moved, (gene, cds.location)


def test_p3_short_end_extended_to_first_stop_across_empty_gap():
    """Metazoa benchmark 2026-10: the cores end nad5 / nad4 / nad2 / cox1 100-500 nt early in divergent phyla (a T
    inside the gene taken as an abbreviated stop); NCBI ends them at the first complete stop. Nothing annotated in
    between -> extended; a tRNA in between -> not."""
    from Bio.SeqFeature import FeatureLocation
    seq = "ATG" + "GCC" * 49 + "T" + "GC" + "GCC" * 49 + "TAA" + "GCC" * 100      # T at 151, TAA ends at 303
    rec = _mkrec(seq, [("CDS", FeatureLocation(0, 151, 1), "ND5")])
    rows = _curate_rec(rec)
    assert int(_by(rec, "CDS", "ND5").location.end) == 303, rows
    assert any(r["rule"] == "P3" and "stopped short" in r["reason"] for r in rows)
    rec = _mkrec(seq, [("CDS", FeatureLocation(0, 151, 1), "ND5"),
                       ("tRNA", FeatureLocation(200, 260, -1), "trnF", {"product": ["tRNA-Phe"]})])
    _curate_rec(rec)
    assert int(_by(rec, "CDS", "ND5").location.end) == 151


def test_p7_frame_repair_keeps_stop_and_moves_start_into_its_frame():
    """Metazoa benchmark 2026-10 (NC_010432.1 cox3, NC_013831.1 nad4): the start is 1-2 nt off the frame of the
    stop -> internal stops (NCBI ERROR). The stop is kept and the in-frame start nearest the old start is used."""
    from Bio.SeqFeature import FeatureLocation
    seq = "GCC" * 10 + "ATG" + "CTAAGC" * 30 + "TAA" + "GCC" * 20             # ATG 31..33, TAA ends at 216;
    # read from 32 the body gives TAA (internal stops)
    rec = _mkrec(seq, [("CDS", FeatureLocation(31, 216, 1), "COX3")])          # starts at 32: wrong frame
    rows = _curate_rec(rec)
    f = _by(rec, "CDS", "COX3")
    assert (int(f.location.start) + 1, int(f.location.end)) == (31, 216), rows
    assert any(r["rule"] == "P7" for r in rows)


def test_d1_antisense_copy_with_internal_stops_removed():
    from Bio.SeqFeature import FeatureLocation
    gene = "ATG" + "GCC" * 60 + "TAA"
    seq = "GCC" * 10 + gene + "GCC" * 20
    rec = _mkrec(seq, [("CDS", FeatureLocation(30, 216, 1), "ND3"), ("CDS", FeatureLocation(60, 180, -1), "ND3")])
    rows = _curate_rec(rec)
    cds = [f for f in rec.features if f.type == "CDS"]
    assert len(cds) == 1 and cds[0].location.strand == 1, rows
    assert any(r["rule"] == "D1" and "antisense" in r["reason"] for r in rows)


def test_p7_in_frame_internal_stop_left_to_p4():
    """Review 2026-10: P7 must not cut a CDS whose start and end are in frame (a nonsense codon or frameshift the
    user has to look at); only ends in different frames are repaired."""
    from Bio.SeqFeature import FeatureLocation
    seq = "GCC" * 10 + "ATG" + "GCC" * 30 + "TAA" + "GCC" * 40 + "TAA" + "GCC" * 20   # internal stop 123 nt before the end
    rec = _mkrec(seq, [("CDS", FeatureLocation(30, 249, 1), "ND2")])
    rows = _curate_rec(rec)
    f = _by(rec, "CDS", "ND2")
    assert (int(f.location.start) + 1, int(f.location.end)) == (31, 249), rows
    assert not any(r["rule"] == "P7" for r in rows)
    assert any(r["rule"] == "P4" and "internal stop" in r["reason"] for r in rows), rows


def test_p7_start_search_stops_at_the_preceding_gene():
    """The stop-kept ORF may not reach into the preceding tRNA for its start (P1's bound)."""
    from Bio.SeqFeature import FeatureLocation
    seq = "GCC" * 20 + "TTG" + "GCC" * 60 + "CTAAGC" * 40 + "TAA" + "GCC" * 10
    rec = _mkrec(seq, [("tRNA", FeatureLocation(0, 150, 1), "trnF", {"product": ["tRNA-Phe"]}),
                       ("CDS", FeatureLocation(244, 486, 1), "ND1")])
    rows = _curate_rec(rec)
    assert int(_by(rec, "CDS", "ND1").location.start) + 1 == 245, rows     # the only start (TTG at 61) is in the tRNA
    assert not any(r["rule"] == "P7" for r in rows)


def test_p3_extension_capped():
    from Bio.SeqFeature import FeatureLocation
    seq = "GCC" * 10 + "ATG" + "GCC" * 100 + "TA" + "CGG" * 220 + "C" + "TAA" + "GCC" * 20
    rec = _mkrec(seq, [("CDS", FeatureLocation(30, 335, 1), "COX1")])
    _curate_rec(rec)
    assert int(_by(rec, "CDS", "COX1").location.end) == 335          # the first stop lies 664 nt further: TA stop kept


def test_p3_annotated_complete_stop_overlapping_next_gene_kept():
    """NCBI records hold both patterns (complete stops overlapping the next gene, T/TA ends before it); a complete
    stop the core annotated inside the next gene is kept, not retracted to a chance T."""
    from Bio.SeqFeature import FeatureLocation
    seq = ("GCC" * 10 + "ATG" + "GCC" * 40 + "TCC" + "G" + "ATG" + "GCC" + "GCC" + "GCT" + "AAG" + "GCC" * 20 + "TAA"
           + "GCC" * 10)
    rec = _mkrec(seq, [("CDS", FeatureLocation(30, 171, 1), "COX3"), ("CDS", FeatureLocation(157, 235, 1), "ND3")])
    rows = _curate_rec(rec)
    assert int(_by(rec, "CDS", "COX3").location.end) == 171, rows       # TAA 169..171, 14 nt inside ND3
    assert int(_by(rec, "CDS", "ND3").location.start) + 1 == 158


def test_p3_gene_branch_respects_antisense_trna_between():
    """A tRNA on the other strand between two genes bounds the end like a same-strand one (no extension across it)."""
    from Bio.SeqFeature import FeatureLocation
    seq = ("GCC" * 10 + "ATG" + "GCC" * 41 + "GCC" * 23 + "T" + "ATG" + "GCC" * 5 + "GCT" + "AAG" + "GCC" * 13 + "TAA"
           + "GCC" * 10)
    rec = _mkrec(seq, [("CDS", FeatureLocation(30, 156, 1), "COX1"),
                       ("tRNA", FeatureLocation(156, 225, -1), "trnK", {"product": ["tRNA-Lys"]}),
                       ("CDS", FeatureLocation(226, 292, 1), "COX2")])
    rows = _curate_rec(rec)
    assert int(_by(rec, "CDS", "COX1").location.end) == 156, rows       # not moved to the chance T at 226


def test_p1_overlapping_start_kept_when_next_start_codon_is_far():
    """173-record benchmark: NCBI keeps a start overlapping the preceding tRNA when the next in-frame start codon
    is far (moves > 15 nt wrong 28 / right 14); a near one (<= 15 nt) is taken (right 31 / wrong 12)."""
    from Bio.SeqFeature import FeatureLocation
    far = "GCC" * 20 + "ATG" + "GCC" * 19 + "ATG" + "GCC" * 40 + "TAA" + "GCC" * 10   # tRNA 1..62 overlaps ATG 61..63 by 2
    rec = _mkrec(far, [("tRNA", FeatureLocation(0, 62, 1), "trnF", {"product": ["tRNA-Phe"]}),
                       ("CDS", FeatureLocation(60, 246, 1), "ND1")])
    rows = _curate_rec(rec)
    assert int(_by(rec, "CDS", "ND1").location.start) + 1 == 61, rows          # next ATG is 60 nt away: kept
    near = "GCC" * 20 + "ATG" + "GCC" * 2 + "ATG" + "GCC" * 57 + "TAA" + "GCC" * 10   # next ATG 9 nt away
    rec = _mkrec(near, [("tRNA", FeatureLocation(0, 62, 1), "trnF", {"product": ["tRNA-Phe"]}),
                        ("CDS", FeatureLocation(60, 246, 1), "ND1")])
    rows = _curate_rec(rec)
    assert int(_by(rec, "CDS", "ND1").location.start) + 1 == 70, rows
    assert any(r["rule"] == "P1" for r in rows)


def test_p3_first_stop_deep_inside_following_rrna_not_taken():
    """173-record truth: a complete stop never reaches more than 2 nt into a following same-strand rRNA; the 'first
    complete stop within 90 nt' rule must not extend a T/TA end 50 nt into the rRNA (NC_067900.1 nad1)."""
    from Bio.SeqFeature import FeatureLocation
    seq = "GCC" * 10 + "ATG" + "GCC" * 50 + "TA" + "C" + "GCC" * 16 + "TAA" + "GCC" * 100   # TAA 235..237, in frame
    rec = _mkrec(seq, [("CDS", FeatureLocation(30, 185, 1), "ND1"),
                       ("rRNA", FeatureLocation(185, 537, 1), "rrnL", {"product": ["16S ribosomal RNA"]})])
    rows = _curate_rec(rec)
    assert int(_by(rec, "CDS", "ND1").location.end) == 185, rows
    assert int(_by(rec, "rRNA", "rrnL").location.start) + 1 == 186


def test_p3_complete_stop_up_to_30_nt_inside_next_gene_taken():
    """173-record truth: a complete stop 11-30 nt inside the next same-strand protein gene is what NCBI annotates
    (17 vs 6 T/TA ends); deeper than 30 nt the T/TA end wins."""
    from Bio.SeqFeature import FeatureLocation
    seq = ("GCC" * 10 + "ATG" + "GCC" * 40 + "TAC" + "G" + "ATG" + "GCC" * 4 + "GCT" + "AAG" + "GCC" * 20 + "TAA"
           + "GCC" * 10)                                          # COX1-frame TAA at 175..177, 20 nt inside ND3 (158..241)
    rec = _mkrec(seq, [("CDS", FeatureLocation(30, 155, 1), "COX1"), ("CDS", FeatureLocation(157, 241, 1), "ND3")])
    rows = _curate_rec(rec)
    assert int(_by(rec, "CDS", "COX1").location.end) == 177, rows
    assert (int(_by(rec, "CDS", "ND3").location.start) + 1, int(_by(rec, "CDS", "ND3").location.end)) == (158, 241)


def test_p3_end_without_stop_signal_takes_the_tail_before_a_deep_stop():
    """Native core ends at the aligned codon; the first in-frame stop 51 nt inside the antisense tRNA must not be
    taken (173-record truth: never deeper than 10 nt) - the T right before the tRNA is the abbreviated stop."""
    from Bio.SeqFeature import FeatureLocation
    seq = "GCC" * 10 + "ATG" + "GCC" * 40 + "TGC" + "CGC" * 16 + "TAA" + "GCC" * 20      # T at 154; TAA 205..207
    rec = _mkrec(seq, [("CDS", FeatureLocation(30, 153, 1), "ND1"),
                       ("tRNA", FeatureLocation(154, 224, -1), "trnS", {"product": ["tRNA-Ser"]})])
    rows = _curate_rec(rec)
    assert int(_by(rec, "CDS", "ND1").location.end) == 154, rows
    assert any(r["rule"] == "P3" and "before the next feature" in r["reason"] for r in rows)


def test_o1_abbreviated_stop_at_the_sequence_end_beats_a_deep_continuation():
    """NC_065213.1 / NC_065650.1 nad2: the gene reaches the sequence end without a stop codon; the first stop lies
    23-29 nt inside cox1 at the origin, NCBI ends the gene with the TA at the sequence end."""
    from Bio.SeqFeature import FeatureLocation
    cox1 = "ATG" + "GCC" * 5 + "GCT" + "AAG" + "GCC" * 20 + "TAA"     # at the origin; a ND2-frame TAA at 21..23
    nad2 = "ATG" + "GCC" * 30 + "T"                                     # ends with T at the sequence end (T + AT: no stop)
    seq = cox1 + "GCC" * 10 + nad2
    L = len(seq)
    rec = _mkrec(seq, [("CDS", FeatureLocation(0, len(cox1), 1), "COX1"),
                       ("CDS", FeatureLocation(L - len(nad2), L - 1, 1), "ND2")], "circular")   # core end: aligned codon
    rows = _curate_rec(rec)
    nd2 = _by(rec, "CDS", "ND2")
    assert len(nd2.location.parts) == 1 and int(nd2.location.end) == L, rows
    assert nd2.qualifiers.get("transl_except") == [f"(pos:{L},aa:TERM)"]


# ── native core review 2026-10 (stop-less CDS: every path must end in P3 / O1 / P6 with a stop or a sequence-end partial) ──
def test_partial_marked_single_interval_cds_is_curated():
    """Review record A: the core writes a CDS cut at a sequence end 5'-partial (`<100..700`, TAA right after); simple()
    refused fuzzy ends, so only P4 'not curated' was logged and the stop was never added (table2asn NoStop)."""
    from Bio.SeqFeature import BeforePosition, FeatureLocation
    seq = "C" * 100 + "GCT" * 200 + "TAA" + "C" * 497
    rec = _mkrec(seq, [("CDS", FeatureLocation(BeforePosition(100), 700, 1), "COX2")], "circular")
    rows = _curate_rec(rec)
    cds = _by(rec, "CDS", "COX2")
    assert int(cds.location.end) == 703 and str(cds.location).startswith("[<100:")
    assert {r["rule"] for r in rows if r["gene"] == "COX2"} == {"P3", "P5"}, rows
    assert "stop codon not determined" not in cds.qualifiers.get("note", [])


def test_p3_short_end_extended_past_an_rrna_filled_to_it():
    """Review record B: core end 150 nt short, rRNA 200 nt downstream. R1 (fill) first moves the rRNA start to the
    CDS end, so the branch-4 bound (rRNA start + 10) blocked the extension and P6 marked the gene 3'-partial. The
    rRNA start is soft (R1 trims it back afterwards), as in the branch-3 rRNA exception."""
    from Bio.SeqFeature import FeatureLocation
    seq = "C" * 100 + "ATG" + "GCT" * 300 + "C" * 150 + "TAA" + "C" * 47 + "G" * 1000 + "C" * 397
    assert seq[1153:1156] == "TAA" and len(seq) == 2600
    rec = _mkrec(seq, [("CDS", FeatureLocation(100, 1003, 1), "ND1"),
                       ("rRNA", FeatureLocation(1203, 2203, 1), "l-rRNA", {"product": ["16S ribosomal RNA"]})], "circular")
    rows = _curate_rec(rec)
    assert int(_by(rec, "CDS", "ND1").location.end) == 1156, rows
    assert int(_by(rec, "rRNA", "rrnL").location.start) + 1 == 1157
    assert not any(r["rule"] == "P6" for r in rows)
    # a tRNA in between still bounds the extension (review record C: designed conservative behaviour)
    rec = _mkrec(seq, [("CDS", FeatureLocation(100, 1003, 1), "ND1"),
                       ("tRNA", FeatureLocation(1050, 1115, -1), "trnL1", {"product": ["tRNA-Leu"]})], "circular")
    rows = _curate_rec(rec)
    assert int(_by(rec, "CDS", "ND1").location.end) == 1003 and any(r["rule"] == "P6" for r in rows)


def test_o1_continues_a_cds_ending_further_from_the_sequence_end():
    """Review record F: the core ends nad2 at its last aligned codon 12 nt before the end of the rotated circular
    record, the stop lies 6 nt past the origin; O1 required an end <= 2 nt before L, so P6 marked it 3'-partial."""
    from Bio.SeqFeature import FeatureLocation
    seq = "C" * 6 + "TAA" + "ATG" + "GCT" * 199 + "C" * 491 + "ATG" + "GCT" * 295 + "C" * 12
    assert len(seq) == 2000
    rec = _mkrec(seq, [("CDS", FeatureLocation(9, 609, 1), "COX1"), ("CDS", FeatureLocation(1100, 1988, 1), "ND2")], "circular")
    rows = _curate_rec(rec)
    nd2 = _by(rec, "CDS", "ND2")
    assert [(int(p.start) + 1, int(p.end)) for p in nd2.location.parts] == [(1101, 2000), (1, 9)], rows
    assert any(r["rule"] == "O1" and r["gene"] == "ND2" for r in rows)
    # an in-frame stop before the sequence end (P3's business), or a feature in between: not continued
    seq2 = seq[:1991] + "TAA" + seq[1994:]                                  # in frame: (1991 - 1100) % 3 == 0
    rec = _mkrec(seq2, [("CDS", FeatureLocation(9, 609, 1), "COX1"), ("CDS", FeatureLocation(1100, 1988, 1), "ND2")], "circular")
    rows = _curate_rec(rec)
    assert len(_by(rec, "CDS", "ND2").location.parts) == 1 and int(_by(rec, "CDS", "ND2").location.end) == 1994


def test_p6_cds_cut_at_the_sequence_end_extended_to_it():
    """Stress case 4b (linear NC_065213.1): nad2 15558..>16559 on a 16560-nt record left the T at 16560 out; a 3'-partial
    end short of the sequence end is table2asn WARNING SEQ_FEAT.PartialProblemNotSpliceConsensus3Prime."""
    from Bio.SeqFeature import FeatureLocation
    seq = "C" * 60 + "ATG" + "GCT" * 40 + "T"                               # 184 nt; CDS 61..183, cut codon's T at 184
    rec = _mkrec(seq, [("CDS", FeatureLocation(60, 183, 1), "ND2")], "linear")
    rows = _curate_rec(rec)
    nd2 = _by(rec, "CDS", "ND2")
    assert str(nd2.location) == "[60:>184](+)" and str(_by(rec, "gene", "ND2").location) == "[60:>184](+)"
    assert any(r["rule"] == "P6" and "extended 1 nt to the sequence end" in r["reason"] for r in rows), rows
    assert nd2.qualifiers["translation"] == ["M" + "A" * 40]
    # minus strand: cut at position 1
    from Bio.Seq import Seq
    rc = str(Seq(seq).reverse_complement())
    rec = _mkrec(rc, [("CDS", FeatureLocation(1, 124, -1), "ND2")], "linear")
    _curate_rec(rec)
    assert str(_by(rec, "CDS", "ND2").location) == "[<0:124](-)"


def test_trna_number_kept_when_the_anticodon_is_unknown():
    """The core numbers trnL1/L2, trnS1/S2 from the MiTFi model even when MiTFi finds no anticodon; names() recomputed
    the name from the (empty) anticodon and wrote trnL."""
    from Bio.SeqFeature import FeatureLocation
    seq = "C" * 300
    rec = _mkrec(seq, [("tRNA", FeatureLocation(10, 75, 1), "trnL1", {"product": ["tRNA-Leu"]}),
                       ("tRNA", FeatureLocation(100, 165, 1), "trnS2", {"product": ["tRNA-Ser"]}),
                       ("tRNA", FeatureLocation(200, 265, 1), "trnL1", {"product": ["tRNA-Leu"], "note": ["anticodon: UAA"]})])
    _curate_rec(rec)
    names = [(f.qualifiers["gene"][0], int(f.location.start)) for f in rec.features if f.type == "tRNA"]
    assert names == [("trnL1", 10), ("trnS2", 100), ("trnL2", 200)]       # a known anticodon still decides
    assert _by(rec, "gene", "trnS2").qualifiers["gene"] == ["trnS2"]


def test_abbreviated_overlap_vertebrate_convention():
    """--abbreviated_overlap (Chordata): a complete stop running 1-2 nt into the next same-strand tRNA becomes the
    T / TA before it; without the option the complete stop is kept (other phyla annotate it so)."""
    from Bio.SeqFeature import FeatureLocation
    seq = "GCC" * 10 + "ATG" + "GCC" * 40 + "TAA" + "GCCGCC" + "G" * 59 + "GCC" * 10
    # CDS 31..156 (ATG .. TAA at 154..156); tRNA starts at 156: the TAA overlaps it by 1 nt, TA ends at 155
    for flag, end in ((True, 155), (False, 156)):
        rec = _mkrec(seq, [("CDS", FeatureLocation(30, 156, 1), "COX3"),
                           ("tRNA", FeatureLocation(155, 220, 1), "trnG", {"product": ["tRNA-Gly"]})])
        _curate_rec(rec, abbreviated_overlap=flag)
        assert int(_by(rec, "CDS", "COX3").location.end) == end, flag


def test_upstream_atg_for_non_atg_start():
    """P1: a non-ATG start (ATT) with an in-frame ATG upstream, open frame, after the preceding tRNA -> the ATG (leave-out
    benchmark +37 starts); --no_upstream_atg (nematodes, arthropods) keeps ATT; an ATG start is never moved."""
    from Bio.SeqFeature import FeatureLocation
    for first, flag, moved in (("ATG", False, True), ("ATG", True, False), ("ATA", False, False)):
        body = first + "GCC" * 8 + "ATT" + "GCC" * 100 + "TAA"
        seq = "C" * 1000 + body + "C" * 1000
        rec = _mkrec(seq, [("tRNA", FeatureLocation(900, 965, 1), "trnF", {"product": ["tRNA-Phe"]}),
                           ("CDS", FeatureLocation(1027, 1000 + len(body), 1), "CYTB")])
        _curate_rec(rec, no_upstream_atg=flag)
        assert (int(_by(rec, "CDS", "CYTB").location.start) == 1000) == moved, (first, flag)


def test_p3_complete_stop_a_few_nt_into_the_trna_when_no_abbreviated_stop():
    """No T / TA ends right before the downstream tRNA (its predicted start is approximate): the first complete stop is
    taken when it runs <= OVERLAP_FALLBACK nt into the tRNA, instead of a 3'-partial gene."""
    from Bio.SeqFeature import FeatureLocation
    seq = "GCC" * 10 + "ATG" + "GCC" * 40 + "GCC" + "TAA" + "GCC" * 30
    # CDS 31..153 (no stop), GCC 154..156, TAA 157..159; tRNA predicted from 155: no T / TA before it, TAA 5 nt inside
    rec = _mkrec(seq, [("CDS", FeatureLocation(30, 153, 1), "COX2"),
                       ("tRNA", FeatureLocation(154, 220, 1), "trnK", {"product": ["tRNA-Lys"]})])
    rows = _curate_rec(rec)
    assert int(_by(rec, "CDS", "COX2").location.end) == 159, rows


def test_o1_stop_completed_one_nt_past_the_origin_is_an_abbreviated_stop():
    """NC_023772.1 cob ends TA at the sequence end, the A at position 1 completes TAA: NCBI writes TA + transl_except
    (11 of 12 such genes in 255 RefSeq records), not join(...,1..1)."""
    from Bio.SeqFeature import FeatureLocation
    cds = "ATG" + "GCC" * 30 + "TA"
    seq = "A" + "GCC" * 30 + cds
    L = len(seq)
    rec = _mkrec(seq, [("CDS", FeatureLocation(L - len(cds), L - 2, 1), "CYTB")], "circular")
    rows = _curate_rec(rec)
    cob = _by(rec, "CDS", "CYTB")
    assert len(cob.location.parts) == 1 and int(cob.location.end) == L, rows
    assert cob.qualifiers.get("transl_except") == [f"(pos:{L - 1}..{L},aa:TERM)"]


def test_p3_core_complete_stop_a_few_nt_into_the_trna_kept():
    """NC_023772.1 nad4: the core's alignment reaches a complete TAG that runs 10 nt into the predicted tRNA-His;
    NCBI keeps it (P3 used to cut back to a T / TA before the tRNA)."""
    from Bio.SeqFeature import FeatureLocation
    seq = "GCC" * 10 + "ATG" + "GCC" * 40 + "GCT" + "TAG" + "GCC" * 30
    # CDS 31..159 with TAG 157..159; tRNA predicted from 150 (TAG 10 nt inside); a T ends at 148 (GCT)
    rec = _mkrec(seq, [("CDS", FeatureLocation(30, 159, 1), "ND4"),
                       ("tRNA", FeatureLocation(149, 215, 1), "trnH", {"product": ["tRNA-His"]})])
    rows = _curate_rec(rec)
    assert int(_by(rec, "CDS", "ND4").location.end) == 159, rows


def test_rrna_across_origin_trimmed_and_filled():
    """R1 skipped rRNAs written as a join across the origin. MitoM's frame starts at cox1, so the echinoderm rrnL
    before cox1 crosses it (NC_013881.1: 19 nt into cox1, 46 nt to nad2 left unfilled): one interval on the circle."""
    import random
    from Bio.SeqFeature import CompoundLocation, FeatureLocation
    n = 3000
    rnd = random.Random(5)
    seq = "".join(rnd.choice("ACGT") for _ in range(n))
    rrnl = CompoundLocation([FeatureLocation(1650, n, 1), FeatureLocation(0, 19, 1)])     # 1651..3000 + 1..19
    rec = _mkrec(seq, [("CDS", FeatureLocation(0, 600, 1), "COX1"),
                       ("CDS", FeatureLocation(1000, 1604, 1), "ND2"),
                       ("rRNA", rrnl, "l-rRNA", {"product": ["16S ribosomal RNA"]})], "circular")
    rows = _curate_rec(rec)
    r = _by(rec, "rRNA", "rrnL")
    assert len(r.location.parts) == 1, r.location                   # trimmed to the cox1 start: no longer across
    assert (int(r.location.start) + 1, int(r.location.end)) == (1605, n)
    assert any(x["rule"] == "R1" for x in rows)
    assert str(_by(rec, "gene", "rrnL").location) == str(r.location)


def test_rrna_filled_across_the_origin_is_a_join():
    """A gap that lies across the origin is filled like any other: the rRNA becomes a join over the origin."""
    import random
    from Bio.SeqFeature import FeatureLocation
    n = 3000
    rnd = random.Random(6)
    seq = "".join(rnd.choice("ACGT") for _ in range(n))
    rec = _mkrec(seq, [("CDS", FeatureLocation(20, 620, 1), "COX1"),
                       ("tRNA", FeatureLocation(1500, 1570, 1), "trnV", {"product": ["tRNA-Val"]}),
                       ("rRNA", FeatureLocation(1570, n - 10, 1), "l-rRNA", {"product": ["16S ribosomal RNA"]})],
                 "circular")
    _curate_rec(rec)
    r = _by(rec, "rRNA", "rrnL")
    assert [(int(p.start), int(p.end)) for p in r.location.parts] == [(1570, n), (0, 20)], r.location


def test_rrna_end_left_in_place_keeps_its_partial_mark():
    """DQ517436.1: the rRNA HMM step adds rrnL <1..481 (the partial sequence starts inside the gene); R1 trims the 3'
    end to the next tRNA and must keep the < of the 5' end it does not move. An end R1 fills to a gene is determined;
    a < inside the sequence (model columns missing, not a cut gene) goes as before when R1 changes the rRNA."""
    import random
    from Bio.SeqFeature import BeforePosition, FeatureLocation
    rnd = random.Random(12)
    seq = "".join(rnd.choice("ACGT") for _ in range(3000))
    rec = _mkrec(seq, [("rRNA", FeatureLocation(BeforePosition(0), 481, 1), "l-rRNA", {"product": ["16S ribosomal RNA"]}),
                       ("tRNA", FeatureLocation(479, 549, 1), "trnA", {"product": ["tRNA-Ala"]})])
    rows = _curate_rec(rec)
    r = _by(rec, "rRNA", "rrnL")
    assert isinstance(r.location.start, BeforePosition) and (int(r.location.start) + 1, int(r.location.end)) == (1, 479)
    assert str(_by(rec, "gene", "rrnL").location) == str(r.location) and any(x["rule"] == "R1" for x in rows)
    rec = _mkrec(seq, [("tRNA", FeatureLocation(0, 70, 1), "trnV", {"product": ["tRNA-Val"]}),
                       ("rRNA", FeatureLocation(BeforePosition(100), 900, 1), "l-rRNA", {"product": ["16S ribosomal RNA"]})])
    _curate_rec(rec)
    r = _by(rec, "rRNA", "rrnL")
    assert not isinstance(r.location.start, BeforePosition) and int(r.location.start) == 70     # filled to trnV
    rec = _mkrec(seq, [("tRNA", FeatureLocation(0, 70, 1), "trnV", {"product": ["tRNA-Val"]}),
                       ("rRNA", FeatureLocation(BeforePosition(500), 981, 1), "l-rRNA", {"product": ["16S ribosomal RNA"]}),
                       ("tRNA", FeatureLocation(979, 1049, 1), "trnA", {"product": ["tRNA-Ala"]})])
    _curate_rec(rec)
    r = _by(rec, "rRNA", "rrnL")
    assert not isinstance(r.location.start, BeforePosition) and (int(r.location.start), int(r.location.end)) == (500, 979)


def test_core_notes_for_itself_are_not_in_the_final_genbank():
    """N1: the core's own notes (reference protein and score, MitoM-ML / MiTFi scores, rRNA model columns) are not for
    NCBI; biological notes (anticodon, polyA-completed stop) stay."""
    import random
    from Bio.SeqFeature import FeatureLocation
    rnd = random.Random(13)
    seq = "".join(rnd.choice("ACGT") for _ in range(3000))
    rec = _mkrec(seq, [
        ("CDS", FeatureLocation(0, 300, 1), "COX1", {"note": ["ref NC_023943:cox1-1-0-1533 score 2461 residues 1-511/511"]}),
        ("tRNA", FeatureLocation(400, 465, 1), "trnN", {"product": ["tRNA-Asn"],
                                                        "note": ["anticodon: GUU", "MitoM-ML p 0.969 ends p 0.852"]}),
        ("tRNA", FeatureLocation(470, 535, 1), "trnD", {"product": ["tRNA-Asp"], "note": ["MiTFi trnD.cm bits 45.1 E 1.2e-09"]}),
        ("rRNA", FeatureLocation(600, 1500, 1), "rrnL", {"product": ["16S ribosomal RNA"],
                                                         "note": ["nhmmer Annelida.rrnL.hmm bits 404 model 116-1192/1307"]})])
    _curate_rec(rec)
    notes = {f.qualifiers["gene"][0]: f.qualifiers.get("note", []) for f in rec.features if f.type in ("CDS", "tRNA", "rRNA")}
    assert notes["trnN"] == ["anticodon: GUU"] and notes["trnD"] == [] and notes["rrnL"] == []
    assert not any(n.startswith("ref ") for n in notes["COX1"])


def test_p9_nad3_extra_nucleotide_is_a_join():
    """Bird / turtle nad3: one untranslated nucleotide at CDS position 174 (Mindell et al. 1998). NCBI writes
    join(start..173,175..end) with the frameshift note; a single-frame CDS had internal stops (table2asn InternalStop in
    109 nad3 of the 2026-10 runs). The gene spans the whole join; the protein is translated across it."""
    from Bio.SeqFeature import FeatureLocation
    coding = "ATG" + "GTA" + "AGC" + "GCC" * 54 + "GCT" + "AAC" + "GCC" * 56 + "TAA"   # 351 nt, frame +1 closed early
    gene = coding[:173] + "T" + coding[173:]                                           # extra T at CDS position 174
    seq = "C" * 60 + gene + "C" * 60
    rec = _mkrec(seq, [("CDS", FeatureLocation(60, 60 + len(gene), 1), "ND3", {"transl_table": ["2"]})])
    rows = _curate_rec(rec, genetic_code=2)
    cds = _by(rec, "CDS", "ND3")
    assert [(int(p.start), int(p.end)) for p in cds.location.parts] == [(60, 233), (234, 60 + len(gene))], rows
    assert _mc().ND3_NOTE in cds.qualifiers["note"] and any(r["rule"] == "P9" for r in rows)
    assert (int(_by(rec, "gene", "ND3").location.start), int(_by(rec, "gene", "ND3").location.end)) == (60, 60 + len(gene))
    from Bio.Seq import Seq
    assert cds.qualifiers["translation"] == [str(Seq(coding[:-3]).translate(table=2))]


def test_p8_start_across_a_stop_moves_past_it():
    """NC_027667.1 nad1 / insect nad4L: the start placed upstream across an in-frame stop; the first start codon after the
    last internal stop is taken (>= half the gene kept)."""
    from Bio.SeqFeature import FeatureLocation
    body = "ATG" + "GCC" * 5 + "TAA" + "GCC" * 3 + "ATG" + "GCC" * 100 + "TAA"
    seq = "C" * 60 + body + "C" * 60
    rec = _mkrec(seq, [("CDS", FeatureLocation(60, 60 + len(body), 1), "ND1")])
    rows = _curate_rec(rec)
    cds = _by(rec, "CDS", "ND1")
    assert int(cds.location.start) == 60 + 30 and int(cds.location.end) == 60 + len(body), rows
    assert any(r["rule"] == "P8" for r in rows) and "pseudo" not in cds.qualifiers


def test_p10_untranslatable_cds_is_pseudo():
    """A stop in the middle, no frameshift join, neither open part >= half: /pseudo on CDS and gene, no translation
    (table2asn InternalStop otherwise)."""
    from Bio.SeqFeature import FeatureLocation
    body = "ATG" + "GCC" * 40 + "TAA" + "GCC" * 60 + "TAA"
    seq = "C" * 60 + body + "C" * 60
    rec = _mkrec(seq, [("CDS", FeatureLocation(60, 60 + len(body), 1), "ND2")])
    rows = _curate_rec(rec)
    cds = _by(rec, "CDS", "ND2")
    assert "pseudo" in cds.qualifiers and "pseudo" in _by(rec, "gene", "ND2").qualifiers, rows
    assert "translation" not in cds.qualifiers and (int(cds.location.start), int(cds.location.end)) == (60, 60 + len(body))


def test_o2_core_join_across_the_origin_cut_back():
    """Annelid nad2 before cox1 in the cox1-first frame: the core joined it 14 nt past the origin without a stop codon
    (table2asn NoStop); it ends with the T at the sequence end (abbreviated stop), as NCBI writes it."""
    from Bio.SeqFeature import CompoundLocation, FeatureLocation
    cox1 = "ATG" + "GCC" * 200 + "TAA"
    nad2 = "ATG" + "GCC" * 99 + "T"
    seq = cox1 + "C" * 100 + nad2
    L = len(seq)
    loc = CompoundLocation([FeatureLocation(L - len(nad2), L, 1), FeatureLocation(0, 14, 1)])
    rec = _mkrec(seq, [("CDS", FeatureLocation(0, len(cox1), 1), "COX1"), ("CDS", loc, "ND2")], "circular")
    rows = _curate_rec(rec)
    cds = _by(rec, "CDS", "ND2")
    assert len(cds.location.parts) == 1 and (int(cds.location.start), int(cds.location.end)) == (L - len(nad2), L), rows
    assert cds.qualifiers.get("transl_except") == [f"(pos:{L},aa:TERM)"] and any(r["rule"] == "O2" for r in rows)


def test_p1_start_moved_past_a_whole_trna():
    """Unionid nad1 (NC_085464.1): the start lay upstream across trnL2 (table2asn BadCDScomponentOverlapTRNA); it moves to
    the first in-frame start codon after the tRNA."""
    from Bio.SeqFeature import FeatureLocation
    body = "ATG" + "GCC" * 36 + "ATG" + "GCC" * 270 + "TAA"
    seq = "C" * 59 + body + "C" * 60
    rec = _mkrec(seq, [("CDS", FeatureLocation(59, 59 + len(body), 1), "ND1"),
                       ("tRNA", FeatureLocation(99, 165, 1), "trnL2", {"product": ["tRNA-Leu"]})])
    rows = _curate_rec(rec)
    assert int(_by(rec, "CDS", "ND1").location.start) + 1 == 171, rows


def test_p3_end_before_a_whole_trna_without_stop():
    """A gene whose first stop lies past a whole same-strand tRNA and no T / TA right before it: it ends before the tRNA,
    3'-partial (the gene must not contain the tRNA: table2asn BadCDScomponentOverlapTRNA)."""
    from Bio.SeqFeature import FeatureLocation
    body = "ATG" + "GCC" * 212 + "TAA"
    seq = "C" * 60 + body + "C" * 60
    rec = _mkrec(seq, [("CDS", FeatureLocation(60, 60 + len(body), 1), "COX3"),
                       ("tRNA", FeatureLocation(609, 675, 1), "trnG", {"product": ["tRNA-Gly"]})])
    rows = _curate_rec(rec)
    cds = _by(rec, "CDS", "COX3")
    assert int(cds.location.end) == 609 and _mc().NO_STOP_NOTE in cds.qualifiers["note"], rows


def test_t1_weak_trna_inside_a_gene_removed_strong_kept():
    """A weak MiTFi call (E 9.3e-02) inside a protein gene on its strand is dropped; a strong one (E 4.8e-08) is kept."""
    from Bio.SeqFeature import FeatureLocation
    body = "ATG" + "GCC" * 300 + "TAA"
    seq = "C" * 60 + body + "C" * 60
    for note, kept in (("MiTFi trnC.cm bits 17.8 E 9.3e-02", False), ("MiTFi trnP.cm bits 49.3 E 4.8e-08", True)):
        rec = _mkrec(seq, [("CDS", FeatureLocation(60, 60 + len(body), 1), "ND2"),
                           ("tRNA", FeatureLocation(500, 566, 1), "trnC", {"product": ["tRNA-Cys"], "note": [note]})])
        rows = _curate_rec(rec)
        assert any(f.type == "tRNA" for f in rec.features) == kept, rows


def test_c1_incomplete_gene_around_a_trna_removed():
    """NC_061193.1: a 120-nt 'atp8' without start or stop codon around tRNA-Glu (flatworms lack atp8): removed."""
    from Bio.SeqFeature import FeatureLocation
    body = "TTT" + "GCC" * 39
    seq = "C" * 60 + body + "C" * 60
    rec = _mkrec(seq, [("CDS", FeatureLocation(60, 180, 1), "ATP8"),
                       ("tRNA", FeatureLocation(95, 160, 1), "trnE", {"product": ["tRNA-Glu"]})])
    rows = _curate_rec(rec)
    assert not any(f.type == "CDS" for f in rec.features) and any(r["rule"] == "C1" for r in rows), rows


def test_partial_end_keeps_no_abbreviated_stop():
    """A CDS cut by the end of a linear sequence (3'-partial) whose last nucleotide is T: no T stop / transl_except (the
    codon may continue past the sequenced part; table2asn PartialProblemHasStop)."""
    from Bio.SeqFeature import AfterPosition, ExactPosition, FeatureLocation
    body = "ATG" + "GCC" * 100 + "T"
    seq = "C" * 60 + body
    rec = _mkrec(seq, [("CDS", FeatureLocation(ExactPosition(60), AfterPosition(len(seq)), 1), "ND2",
                        {"note": ["stop codon not determined"]})])
    _curate_rec(rec)
    cds = _by(rec, "CDS", "ND2")
    assert "transl_except" not in cds.qualifiers and isinstance(cds.location.end, AfterPosition)


def test_o1_no_stop_past_the_origin_ends_with_t_at_the_sequence_end():
    """KT429007.1 / NC_012618.1 / EF494507.1: nad2 before cox1 (cox1-first frame) has no stop within 60 nt past the origin
    (cox1's frame differs); NCBI ends it with the in-frame T at the sequence end + transl_except - it was written
    3'-partial ('Partial CDS on complete sequence')."""
    from Bio.SeqFeature import FeatureLocation
    cox1 = "ATG" + "GCC" * 200 + "TAA"
    nad2 = "ATG" + "GCC" * 99 + "T"                    # 301 nt: in-frame T at the sequence end
    seq = cox1 + "C" * 100 + nad2
    L = len(seq)
    rec = _mkrec(seq, [("CDS", FeatureLocation(0, len(cox1), 1), "COX1"),
                       ("CDS", FeatureLocation(L - len(nad2), L - 1, 1), "ND2")], "circular")
    rows = _curate_rec(rec)
    cds = _by(rec, "CDS", "ND2")
    assert (int(cds.location.start), int(cds.location.end)) == (L - len(nad2), L), rows
    assert cds.qualifiers.get("transl_except") == [f"(pos:{L},aa:TERM)"] and _mc().NO_STOP_NOTE not in cds.qualifiers.get("note", [])


def test_t1_weak_trna_at_a_gene_end_kept_and_the_gene_end_moved():
    """NC_060482.1 / NC_060708.1: the weak trnD (MiTFi E 0.096, as MITOS2) lies in the last third of nad6, which P3 had
    extended to the stop inside trnD; NCBI ends nad6 right before trnD. The tRNA is kept, the gene ends before it."""
    from Bio.SeqFeature import FeatureLocation
    body = "ATG" + "GCC" * 146 + "T" + "AAGCC" * 12 + "GCCGCCTAA"         # T at 502 then tRNA 503..562, TAA at 563..571
    seq = "C" * 60 + body + "C" * 60
    rec = _mkrec(seq, [("CDS", FeatureLocation(60, 60 + len(body), 1), "ND6"),
                       ("tRNA", FeatureLocation(502, 562, 1), "trnD", {"product": ["tRNA-Asp"],
                                                                      "note": ["MiTFi trnD.cm bits 19.8 E 9.6e-02"]})])
    rows = _curate_rec(rec)
    assert any(f.type == "tRNA" for f in rec.features) and not any(r["rule"] == "T1" for r in rows), rows
    assert int(_by(rec, "CDS", "ND6").location.end) <= 504, rows


def test_p9_join_ending_with_ta_before_the_next_trna_vertebrate_convention():
    """Bird nad3 join (skip at CDS position 174) whose TAA runs 1 nt into tRNA-Arg: with --abbreviated_overlap it ends
    with TA right before the tRNA (NCBI, 9 of 107 nad3). The tail is counted on the join (one nucleotide skipped)."""
    from Bio.SeqFeature import FeatureLocation
    from Bio.Seq import Seq
    coding = "ATG" + "GTA" + "AGC" + "GCC" * 54 + "GCT" + "AAC" + "GCC" * 56 + "TAA"
    gene = coding[:173] + "T" + coding[173:]
    seq = "C" * 60 + gene + "C" * 60
    end = 60 + len(gene)                                   # 1-based end of the TAA
    rec = _mkrec(seq, [("CDS", FeatureLocation(60, end, 1), "ND3", {"transl_table": ["2"]}),
                       ("tRNA", FeatureLocation(end - 1, end + 64, 1), "trnR", {"product": ["tRNA-Arg"]})])
    rows = _curate_rec(rec, genetic_code=2, abbreviated_overlap=True)
    cds = _by(rec, "CDS", "ND3")
    assert [(int(p.start), int(p.end)) for p in cds.location.parts] == [(60, 233), (234, end - 1)], rows
    assert cds.qualifiers["transl_except"] == [f"(pos:{end - 2}..{end - 1},aa:TERM)"]
    assert cds.qualifiers["translation"] == [str(Seq(coding[:-3]).translate(table=2))]


def test_cds_across_the_origin_without_stop_marked_3prime_partial():
    """Leave-out runs: the core joined nad4L across the origin on - with only 38 nt before the origin (O2 cannot cut it
    back) and no stop codon (table2asn NoStop, NC_006386.1). The joined CDS is checked like one interval: 3'-partial
    with 'stop codon not determined'."""
    from Bio.Seq import Seq
    from Bio.SeqFeature import BeforePosition, CompoundLocation, FeatureLocation
    plus = "ATG" + "GCC" * 89                                  # 270 nt, no stop
    seq_plus = plus[38:] + "C" * 600 + plus[:38]               # + strand: 38 nt before the origin, 232 after
    L = len(seq_plus)
    rc = str(Seq(seq_plus).reverse_complement())               # the gene on - : 1..38 then L-231..L
    loc = CompoundLocation([FeatureLocation(0, 38, -1), FeatureLocation(L - 232, L, -1)])
    rec = _mkrec(rc, [("CDS", loc, "ND4L")], "circular")
    rows = _curate_rec(rec)
    cds = _by(rec, "CDS", "ND4L")
    assert str(cds.extract(rec.seq)) == plus
    assert isinstance(cds.location.parts[-1].start, BeforePosition) and _mc().NO_STOP_NOTE in cds.qualifiers["note"], rows
