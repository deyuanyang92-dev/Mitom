"""--topology auto (2026-10-30): the core's evidence at the sequence ends (mito_core.end_evidence, exonerate replaced by
a string search), the decision (mitom_topology.decide / apply_topology), where a label comes from and --topology
circular / linear (mitom_inputs), and the LOCUS line of the restored header (mitom_metadata.keep_topology).

Run: python3 -m pytest -q tests/test_topology.py
"""
import argparse
import logging
import random
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_core as mc  # noqa: E402
import mitom_topology as mt  # noqa: E402
from mito_common import CoreFeature  # noqa: E402
from mito_taxa import linear_mtdna  # noqa: E402
from mitom_inputs import SampleInfo, split_or_sanitize_fastas  # noqa: E402
from mitom_metadata import keep_topology  # noqa: E402

LOG = logging.getLogger("test")
rng = random.Random(7)
CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in ("TAA", "TAG")]
GENE = "".join(rng.choice(CODONS) for _ in range(300))           # a 300-codon "cox1", open frame (code 5)
FILLER = "".join(rng.choice("ACGT") for _ in range(11000))


def fake_exonerate(gene: str, trim: int = 0):
    """exonerate_best stand-in: the whole gene where the window holds it (ends joined), else its longest 5' part at
    the window's end or 3' part at its start (the gene cut by the junction); trim: nt of a 5' part not aligned at its
    end (exonerate drops the last codons of a gene cut by a sequence end)."""
    def best(exonerate, window, prots, code, tmp, log, model=None):
        if model == "protein2genome":                     # spliced: the gene's two parts with an insertion between
            for k in range(90, len(gene) - 89, 3):
                i = window.find(gene[:k])
                j = window.find(gene[k:], i + k + 1) if i >= 0 else -1
                if i >= 0 and j > i + k:
                    return "ref1", 0, 300, i, j + len(gene) - k, 900.0
            return None
        i = window.find(gene)
        if i >= 0:
            return "ref1", 0, 300, i, i + len(gene), 900.0
        for k in range(len(gene) - 3, 89, -3):
            if window.find(gene[:k]) >= 0:
                i = window.find(gene[:k])
                return "ref1", 0, (k - trim) // 3, i, i + k - trim, 3.0 * k
            if window.find(gene[-k:]) >= 0:
                i = window.find(gene[-k:])
                return "ref1", 300 - k // 3, 300, i, i + k, 3.0 * k
        return None
    return best


def evidence(monkeypatch, genome, feats, gene=GENE, trim=0):
    monkeypatch.setattr(mc, "exonerate_best", fake_exonerate(gene, trim))
    monkeypatch.setattr(mc, "reference_proteins", lambda db, ids: {i: "M" * 300 for i in ids})
    rows = mc.end_evidence(genome, feats, Path("."), 5, {"exonerate": "exonerate"}, Path("."), LOG)
    return [dict(zip(mc.ENDS_HEADER.rstrip("\n").split("\t")[1:], r.rstrip("\n").split("\t"))) for r in rows]


def cds(start, end, strand, residues):
    return CoreFeature("CDS", "cox1", start, end, strand, f"ref ref1 score 900 residues {residues}/300")


def test_gene_across_the_ends_of_a_linear_record_joins_them(monkeypatch):
    """An assembly started inside cox1 (labelled linear: several sequences in one FASTA): the 3' part at the start,
    the 5' part at the end; aligned across the junction from either piece -> joined."""
    genome = GENE[600:] + FILLER + GENE[:600]
    n = len(genome)
    rows = evidence(monkeypatch, genome, [cds(n - 600, n, 1, "1-200"), cds(0, 300, 1, "201-300")])
    assert [(r["end"], r["before_nt"], r["after_nt"], r["verdict"]) for r in rows] == \
        [("end", "600", "300", "joined"), ("start", "600", "300", "joined")]


def test_minus_strand_gene_across_the_ends(monkeypatch):
    """The same on the minus strand (NC_035740.1 cox1: its 3' end lies at the end of the record)."""
    rc = str(Seq(GENE[600:] + FILLER + GENE[:600]).reverse_complement())
    n = len(rc)
    rows = evidence(monkeypatch, rc, [cds(0, 600, -1, "1-200")])
    assert (rows[0]["end"], rows[0]["before_nt"], rows[0]["after_nt"], rows[0]["verdict"]) == ("start", "600", "300", "joined")
    circ = evidence(monkeypatch, rc, [cds(n - 300, n + 600, -1, "1-300")])         # the core's join across the origin
    assert (circ[0]["end"], circ[0]["verdict"]) == ("across", "joined")


def test_gene_cut_at_the_end_of_a_partial_sequence(monkeypatch):
    """A partial genome: cox1's 5' 600 nt at the end, its 3' part missing -> cut (200 of 300 residues)."""
    genome = FILLER + GENE[:600]
    rows = evidence(monkeypatch, genome, [cds(len(genome) - 600, len(genome), 1, "1-200")])
    assert (rows[0]["residues"], rows[0]["after_nt"], rows[0]["verdict"]) == ("1-200/300", "0", "cut")


def test_gene_split_by_an_intron_across_the_ends(monkeypatch):
    """A circular record started inside a gene with an intron (Travisia cox1: exon 2 at the start, exon 1 and the
    intron at the end; every MJ20260820041 CHY record): protein2dna stops at the junction, protein2genome aligns the
    gene across it with the intron -> joined_intron (circular), not cut (linear)."""
    intron = "".join(random.Random(3).choice("ACGT") for _ in range(1100))
    genome = GENE[300:] + FILLER + GENE[:300] + intron                # exon 2 | ... | exon 1, intron | (junction)
    rows = evidence(monkeypatch, genome, [cds(0, 600, 1, "101-300")])
    assert (rows[0]["end"], rows[0]["verdict"]) == ("start", "joined_intron")
    assert mt.decide("circular", rows)[0] == "circular" and "intron" in mt.decide("circular", rows)[1]


def test_intron_longer_than_half_the_rest_of_the_circle(monkeypatch):
    """The intron can hold other genes and be most of the record (hexacoral nad5: 11-12 kb of an 18-kb record), wholly
    on one side of the junction: the spliced re-check gets the rest of the circle on either side, on both strands."""
    intron = "".join(random.Random(4).choice("ACGT") for _ in range(5000))
    genome = GENE[300:] + "A" * 200 + GENE[:300] + intron                # rest of the circle 5,500 nt, intron 5,000
    rows = evidence(monkeypatch, genome, [cds(0, 600, 1, "101-300")])
    assert (rows[0]["end"], rows[0]["verdict"]) == ("start", "joined_intron")
    rc = str(Seq(genome).reverse_complement())                          # the same record on the minus strand
    n = len(rc)
    rows = evidence(monkeypatch, rc, [cds(n - 600, n, -1, "101-300")])
    assert rows[0]["verdict"] == "joined_intron"


def test_partial_sequence_is_still_cut_with_the_spliced_model(monkeypatch):
    """The rest of the gene is nowhere in the record: protein2genome finds no spliced alignment either -> cut."""
    genome = GENE[300:] + FILLER
    rows = evidence(monkeypatch, genome, [cds(0, 600, 1, "101-300")])
    assert rows[0]["verdict"] == "cut"


def test_alignment_trimmed_short_of_the_end(monkeypatch):
    """exonerate stops ~20 nt before the end of a gene cut by it (synthetic tests 2026-10-30: NC_006990.1 cox1,
    NC_045069.1 nad2 21 nt short): still checked, as the missing residues would not fit before the end."""
    genome = GENE[600:] + FILLER + GENE[:600]
    n = len(genome)
    rows = evidence(monkeypatch, genome, [cds(n - 600, n - 21, 1, "1-193")], trim=21)
    assert (rows[0]["end"], rows[0]["verdict"]) == ("end", "joined")
    part = FILLER + GENE[:600]
    rows = evidence(monkeypatch, part, [cds(len(part) - 600, len(part) - 21, 1, "1-193")], trim=21)
    assert (rows[0]["after_nt"], rows[0]["verdict"]) == ("-21", "cut")


def test_stop_codon_before_the_junction_is_a_complete_gene(monkeypatch):
    """A gene whose divergent 3' end does not align (120 residues missing) but whose stop codon lies before the
    junction is complete, not cut by it (a circular genome stays circular)."""
    for stop, verdict in (("TAA", "complete"), ("GCC", "cut")):
        genome = FILLER + GENE[:540] + stop + GENE[543:600]
        n = len(genome)
        rows = evidence(monkeypatch, genome, [cds(n - 600, n - 60, 1, "1-180")], gene=GENE[:540] + stop + GENE[543:],
                        trim=60)
        assert rows[0]["verdict"] == verdict, stop


def test_alignment_across_needs_both_halves_or_a_long_open_span(monkeypatch):
    """An alignment across the junction is joined when each half aligns alone and the pieces meet, or when both sides
    are long (>= 90 nt) in an open frame - a divergent half may align only with the rest (synthetic NC_045069.1 nad2,
    510 / 372 nt). A short extension whose halves do not meet is no evidence either way (5 synthetic partial genomes:
    35-42 nt across by chance)."""
    genome = GENE[600:] + FILLER + GENE[:600]
    n = len(genome)
    whole = fake_exonerate(GENE)

    def one_half(exonerate, window, prots, code, tmp, log):        # the 3' half too divergent to align alone
        hit = whole(exonerate, window, prots, code, tmp, log)
        return hit if hit and hit[1] < 200 else None
    monkeypatch.setattr(mc, "exonerate_best", one_half)
    monkeypatch.setattr(mc, "reference_proteins", lambda db, ids: {i: "M" * 300 for i in ids})
    rows = mc.end_evidence(genome, [cds(n - 600, n, 1, "1-200")], Path("."), 5, {"exonerate": "x"}, Path("."), LOG)
    assert rows[0].split("\t")[-1].strip() == "joined"
    short = GENE[870:] + FILLER + GENE[:600]                          # 30 nt across: halves needed, the 3' one is missing
    n = len(short)

    def short_half(exonerate, window, prots, code, tmp, log):
        i, k = window.find(GENE[:600] + GENE[870:]), 630
        if i >= 0:
            return "ref1", 0, 210, i, i + k, 630.0
        i = window.find(GENE[:600])
        return ("ref1", 0, 200, i, i + 600, 600.0) if i >= 0 else None
    monkeypatch.setattr(mc, "exonerate_best", short_half)
    rows = mc.end_evidence(short, [cds(n - 600, n, 1, "1-200")], Path("."), 5, {"exonerate": "x"}, Path("."), LOG)
    assert rows[0].split("\t")[-1].strip() == "complete"


def test_complete_gene_at_an_end_is_no_evidence(monkeypatch):
    """cox1 starting at position 1 (most NCBI records): complete, whatever lies before it."""
    genome = GENE + FILLER
    rows = evidence(monkeypatch, genome, [cds(0, 900, 1, "1-300")])
    assert (rows[0]["end"], rows[0]["before_nt"], rows[0]["verdict"]) == ("start", "0", "complete")
    assert evidence(monkeypatch, genome, [cds(5000, 5900, 1, "1-300")]) == []        # not at an end: not checked


def test_decide():
    joined = {"gene": "cox1", "end": "end", "residues": "1-300/300", "before_nt": "600", "after_nt": "300",
              "verdict": "joined"}
    cut = {**joined, "residues": "1-200/300", "after_nt": "0", "verdict": "cut"}
    assert mt.decide("linear", [joined])[0] == "circular"
    assert mt.decide("circular", [cut]) == ("linear", "cox1 stops at the junction of the sequence ends with residues "
                                            "1-200/300 of its protein, the rest is not on the other side: the ends do "
                                            "not join")
    assert mt.decide("linear", [cut]) == ("linear", "as labelled")
    assert mt.decide("circular", []) == ("circular", "as labelled")
    hydrozoa = "Eukaryota;Metazoa;Cnidaria;Hydrozoa;Hydroidolina"
    assert mt.decide("circular", [], hydrozoa) == ("linear", "Hydrozoa: linear mitochondrial genome, no gene across "
                                                             "the ends")
    assert mt.decide("circular", [joined], hydrozoa)[0] == "circular"              # NC_035740.1 (Scyphozoa)
    assert linear_mtdna("Metazoa;Porifera;Calcarea") == "Calcarea" and linear_mtdna("Metazoa;Cnidaria;Anthozoa") == ""


def test_apply_topology_annotates_changed_records_again(tmp_path, monkeypatch):
    split = tmp_path / "split"
    split.mkdir()
    samples, results = [], {}
    for i, (name, verdict) in enumerate((("asm1", "joined"), ("asm2", "complete")), 1):
        iid = f"S{i:05d}"
        fa = split / f"{iid}.fasta"
        fa.write_text(f">{iid} topology=linear\n{FILLER}\n")
        samples.append(SampleInfo(name, name, iid, "linear", fa, fa, None, "FASTA without topology=: 2 sequences in the file"))
        wd = tmp_path / "01.anno_round1" / iid
        wd.mkdir(parents=True)
        (wd / f"{iid}.ends.tsv").write_text(mc.ENDS_HEADER + f"{iid}\tcox1\t+\tend\tref1\t1-300/300\t600\t300\t{verdict}\n")
        results[iid] = {"ok": True, "workdir": str(wd), "best_gb": str(wd / f"{iid}.gbf")}
    seen = []

    def fake_stage(*, samples, phase_dirname, out_root, args, fasta_selector):
        seen.extend((s.sample_name, s.topology, fasta_selector(s).read_text().split("\n")[0]) for s in samples)
        return {s.internal_id: {"ok": True, "workdir": "new", "best_gb": "new.gbf"} for s in samples}
    monkeypatch.setattr(mt, "run_annotation_stage", fake_stage)
    args = argparse.Namespace(topology="auto", engine="native")
    new, res = mt.apply_topology(samples, results, args, tmp_path, tmp_path / "reports")
    assert seen == [("asm1", "circular", ">S00001 topology=circular")]
    assert [(s.sample_name, s.topology) for s in new] == [("asm1", "circular"), ("asm2", "linear")]
    assert res["S00001"]["best_gb"] == "new.gbf" and res["S00002"]["best_gb"].endswith("S00002.gbf")
    rows = (tmp_path / "reports" / "topology_decision.tsv").read_text().splitlines()
    assert rows[0] == mt.DECISION_HEADER.rstrip("\n")
    assert rows[1].split("\t")[:6] == ["asm1", "11000", "linear", "FASTA without topology=: 2 sequences in the file",
                                       "circular", "yes"]
    assert rows[2].split("\t")[4:] == ["linear", "no", "as labelled"]
    seen.clear()
    mt.apply_topology(samples, results, argparse.Namespace(topology="label", engine="native"), tmp_path, tmp_path / "r2")
    assert seen == []                                                               # --topology label: nothing changes
    monkeypatch.setattr(mt, "run_annotation_stage", lambda **kw: {s.internal_id: {"ok": False} for s in kw["samples"]})
    new, res = mt.apply_topology(samples, results, args, tmp_path, tmp_path / "r3")    # annotation as circular failed
    assert [s.topology for s in new] == ["linear", "linear"] and res["S00001"] is results["S00001"]
    row = (tmp_path / "r3" / "topology_decision.tsv").read_text().splitlines()[1].split("\t")
    assert row[4:6] == ["linear", "no"] and row[6].endswith("annotation as circular failed, kept linear")


def test_label_source_and_forced_topology(tmp_path):
    fa = tmp_path / "asm.fa"
    fa.write_text(f">a1\n{FILLER}\n>a2 topology=circular\n{FILLER}\n")
    one = tmp_path / "one.fa"
    one.write_text(f">g1\n{FILLER}\n")
    samples, _ = split_or_sanitize_fastas([fa, one], tmp_path / "s1", True, 15)
    assert [(s.topology, s.topology_source) for s in samples] == [
        ("circular", "FASTA without topology=: decided from the sequence"), ("circular", "FASTA header"),
        ("circular", "FASTA without topology=: decided from the sequence")]
    samples, _ = split_or_sanitize_fastas([fa, one], tmp_path / "s2", True, 15, topology="linear")
    assert [(s.topology, s.label) for s in samples] == [("linear", "circular"), ("linear", "circular"),
                                                        ("linear", "circular")]
    assert "topology=linear" in samples[1].internal_fasta.read_text().split("\n")[0]
    args = argparse.Namespace(topology="linear", engine="native")
    results = {s.internal_id: {"ok": True, "workdir": str(tmp_path)} for s in samples}
    mt.apply_topology(samples, results, args, tmp_path, tmp_path / "rep")
    rows = [r.split("\t") for r in (tmp_path / "rep" / "topology_decision.tsv").read_text().splitlines()[1:]]
    assert [r[2:] for r in rows[1:]] == [["circular", "FASTA header", "linear", "yes", "set by --topology linear"],
                                         ["circular", "FASTA without topology=: decided from the sequence", "linear",
                                          "yes", "set by --topology linear"]]


def test_linear_lineage_is_annotated_linear_from_round_1(tmp_path, monkeypatch):
    """A circular-labelled hydrozoan: linear in round 1 (no second annotation); annotated again as circular only when
    a gene runs across its ends (NC_035740.1, a scyphozoan whose cox1 does)."""
    fa = tmp_path / "S00001.fasta"
    fa.write_text(f">S00001 topology=circular\n{FILLER}\n")
    si = SampleInfo("NC_1", "NC_1", "S00001", "circular", fa, fa, None, "GenBank LOCUS line", "circular")
    monkeypatch.setattr(mt, "sample_lineage", lambda s: "Eukaryota;Metazoa;Cnidaria;Hydrozoa")
    args = argparse.Namespace(topology="auto", engine="native")
    [lin] = mt.linear_lineages_first([si], args)
    assert lin.topology == "linear" and lin.label == "circular" and "topology=linear" in lin.internal_fasta.read_text()
    assert mt.linear_lineages_first([si], argparse.Namespace(topology="label", engine="native")) == [si]
    wd = tmp_path / "w"
    wd.mkdir()
    redo = []
    monkeypatch.setattr(mt, "run_annotation_stage", lambda **kw: redo.extend(kw["samples"]) or
                        {s.internal_id: {"ok": True, "workdir": "x"} for s in kw["samples"]})
    new, _ = mt.apply_topology([lin], {"S00001": {"ok": True, "workdir": str(wd)}}, args, tmp_path, tmp_path / "r")
    assert redo == [] and new[0].topology == "linear"
    row = (tmp_path / "r" / "topology_decision.tsv").read_text().splitlines()[1].split("\t")
    assert row[2:] == ["circular", "GenBank LOCUS line", "linear", "yes", "Hydrozoa: linear mitochondrial genome, no "
                       "gene across the ends"]
    (wd / "S00001.ends.tsv").write_text(mc.ENDS_HEADER + "S00001\tcox1\t-\tstart\tr\t1-530/530\t1548\t42\tjoined\n")
    new, _ = mt.apply_topology([lin], {"S00001": {"ok": True, "workdir": str(wd)}}, args, tmp_path, tmp_path / "r")
    assert [s.topology for s in redo] == ["circular"] and new[0].topology == "circular"


def test_restored_header_keeps_the_annotated_topology():
    meta = "LOCUS       NC_035740              17024 bp    DNA     circular INV 03-APR-2023\nDEFINITION  x.\n"
    new = "LOCUS       S00001                 17024 bp    DNA     linear   INV 01-JAN-1980\n"
    assert keep_topology(meta, new).splitlines()[0] == \
        "LOCUS       NC_035740              17024 bp    DNA     linear   INV 03-APR-2023"
    back = keep_topology(keep_topology(meta, new), meta)
    assert back == meta


def test_intron_features_restored_after_rotation():
    """An intron across the origin is given by the CDS join only (no multi-interval intron); once the record starts
    at the gene's first exon the gap is inside the record and gets its feature back, with the de novo description;
    a gap that still crosses the origin gets none."""
    from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    import mitom_orient
    from mito_intron import INTRON_NOTE
    note = f"{INTRON_NOTE} 1: group II intron (Rfam RF00029)"
    loc = CompoundLocation([FeatureLocation(0, 273, 1), FeatureLocation(1415, 2687, 1)])
    rec = SeqRecord(Seq("A" * 17000), id="r", annotations={"topology": "circular"},
                    features=[SeqFeature(loc, type="CDS", qualifiers={"gene": ["COX1"], "note": [note]})])
    assert mitom_orient.add_exon_gap_introns(rec) == 1
    intron = next(f for f in rec.features if f.type == "intron")
    assert (int(intron.location.start), int(intron.location.end)) == (273, 1415)
    assert intron.qualifiers["note"] == ["group II intron (Rfam RF00029)"]
    assert not any(x.startswith(INTRON_NOTE) for f in rec.features for x in f.qualifiers.get("note", []))
    across = CompoundLocation([FeatureLocation(16000, 17000, 1), FeatureLocation(0, 300, 1)])   # still across
    rec2 = SeqRecord(Seq("A" * 17000), id="r2", annotations={"topology": "circular"},
                     features=[SeqFeature(across, type="CDS", qualifiers={"gene": ["COX1"]})])
    assert mitom_orient.add_exon_gap_introns(rec2) == 0


ALL = {"cox1", "cox2", "cox3", "cob", "atp6", "atp8", "nad1", "nad2", "nad3", "nad4", "nad4l", "nad5", "nad6",
       "rrnS", "rrnL"}


def test_unlabelled_record_decided_by_its_gene_set():
    """No topology label and no gene at the ends (a complete genome started at a gene boundary, or a contig ending in
    non-coding sequence): the complete gene set (>= 12 of 13 protein genes - atp8 is missing in several lineages -
    and both rRNAs) makes it circular, an incomplete one linear; ends that decide win; a label is kept."""
    assert mt.decide("circular", [], unlabelled=True, genes=ALL)[0] == "circular"
    assert mt.decide("circular", [], unlabelled=True, genes=ALL - {"atp8"})[0] == "circular"
    top, why = mt.decide("circular", [], unlabelled=True, genes=ALL - {"atp8", "nad6"})
    assert top == "linear" and "11 of 13 protein genes" in why
    # phyla that keep all tRNAs: fewer than 21 = a contig that lost some; cnidarians (1-2 tRNAs) not checked
    assert mt.decide("circular", [], unlabelled=True, genes=ALL | {"tRNA:19"}, phylum="Annelida")[0] == "linear"
    assert mt.decide("circular", [], unlabelled=True, genes=ALL | {"tRNA:22"}, phylum="Annelida")[0] == "circular"
    assert mt.decide("circular", [], unlabelled=True, genes=ALL | {"tRNA:2"}, phylum="Cnidaria")[0] == "circular"
    assert mt.decide("circular", [], unlabelled=True, genes=ALL - {"rrnS"})[0] == "linear"
    cut = [{"gene": "cox1", "verdict": "cut", "residues": "100-513/513"}]
    assert mt.decide("circular", cut, unlabelled=True, genes=ALL)[0] == "linear"           # a gene cut at an end
    joined = [{"gene": "nad5", "verdict": "joined", "before_nt": "300", "after_nt": "900"}]
    assert mt.decide("linear", joined, unlabelled=True, genes=set())[0] == "circular"
    assert mt.decide("circular", [], genes=set()) == ("circular", "as labelled")           # labelled: kept
    assert mt.decide("linear", [], unlabelled=True, genes=ALL) == ("linear", "as labelled")  # < 10 kb: a contig


def test_gene_set_of_an_annotation(tmp_path):
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    r = SeqRecord(Seq("A" * 300), id="x", annotations={"molecule_type": "DNA"})
    r.features = [SeqFeature(FeatureLocation(0, 90, 1), type="CDS", qualifiers={"gene": ["COX1"]}),
                  SeqFeature(FeatureLocation(90, 200, 1), type="rRNA", qualifiers={"product": ["16S ribosomal RNA"]}),
                  SeqFeature(FeatureLocation(200, 260, 1), type="tRNA", qualifiers={"gene": ["trnL1"]})]
    SeqIO.write(r, str(tmp_path / "x.gb"), "genbank")
    assert mt.gene_set(tmp_path / "x.gb") == {"cox1", "rrnL", "tRNA:1"}
    assert mt.gene_set(tmp_path / "missing.gb") == set()


def test_threads_fit_the_cores():
    """--threads x records at once > cores made Snakemake run one record at a time: threads per record lowered."""
    from mitom_stage import fit_threads
    a = argparse.Namespace(threads=48, max_tasks=16, smk_cores=64)
    fit_threads(a, "snakemake", 148)
    assert a.threads == 4
    a = argparse.Namespace(threads=48, max_tasks=16, smk_cores=64)
    fit_threads(a, "snakemake", 1)                      # one record: all cores for it
    assert a.threads == 48
    a = argparse.Namespace(threads=2, max_tasks=24, smk_cores=64)
    fit_threads(a, "snakemake", 148)                    # fits: unchanged
    assert a.threads == 2


def test_threads_fit_the_cpus_with_the_threads_backend(monkeypatch):
    import mitom_stage
    monkeypatch.setattr(mitom_stage.os, "cpu_count", lambda: 32)
    a = argparse.Namespace(threads=8, max_tasks=16)
    mitom_stage.fit_threads(a, "threads", 100)
    assert a.threads == 2


def test_header_words_give_the_topology():
    from mitom_inputs import header_hint
    assert header_hint("NC_064058.1 Sirsoe methanicola mitochondrion, complete genome")[0] == "circular"
    assert header_hint("1+,2-,3+(circular)")[0] == "circular"                       # GetOrganelle path
    assert header_hint("KU672397.1 Hirudo verbana mitochondrion, partial genome")[0] == "linear"
    assert header_hint("LC726566.1 Metaphire soulensis mitochondrial DNA, nearly complete genome")[0] == "linear"
    assert header_hint("contig_7 non-circular")[0] is None and header_hint("contig_7 len=15000")[0] is None


def test_terminal_overlap_found_and_removed(tmp_path):
    """An assembler's circular overlap (the end repeats the start): circular, the repeated end removed; a chance
    match shorter than MIN_OVERLAP, or none, is no overlap."""
    from mitom_inputs import terminal_overlap
    rnd = random.Random(5)
    body = "".join(rnd.choice("ACGT") for _ in range(12000))
    assert terminal_overlap(body + body[:300]) == 300
    assert terminal_overlap(body + body[:20]) == 0 and terminal_overlap(body) == 0
    fa = tmp_path / "asm.fa"
    fa.write_text(f">ctg1\n{body + body[:300]}\n>ctg2\n{body[:11000]}\n")
    samples, _ = split_or_sanitize_fastas([fa], tmp_path / "s", True, 15)
    first = SeqIO.read(str(samples[0].internal_fasta), "fasta")
    assert len(first.seq) == 12000 and samples[0].topology == "circular"
    assert samples[0].topology_source == "terminal overlap of 300 nt (removed)"
    assert samples[1].topology_source == "FASTA without topology=: decided from the sequence"
    assert terminal_overlap("ACGT" * 4000) == 0                    # a tandem repeat at both ends: no overlap
