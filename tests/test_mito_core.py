"""mito_core.py: MitoM's own annotation core (pure functions; the programs themselves are not run here)."""
import logging
import sys
from pathlib import Path

from Bio.Data import CodonTable

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_core as mc  # noqa: E402

T5 = CodonTable.unambiguous_dna_by_id[5]
LOG = logging.getLogger("test")


def test_place_cds_nearest_upstream_start_and_adjacent_stop():
    """MitoZ's rule: the nearest start codon upstream of the alignment (within missing residues + 10 codons); the end
    is the last aligned codon, plus the stop codon when it follows directly (other stops: mito_curate P3)."""
    window = "GCC" * 5 + "ATG" + "GCC" * 3 + "ATT" + "GCC" * 20 + "TAA" + "GCC" * 5
    tab = 15 + 3 + 9 + 3                    # alignment starts one codon after the ATT (missing 1 residue)
    tae = tab + 60                          # ... and ends right before the TAA
    s, e, partial = mc.place_cds(window, tab, tae, 1, T5)
    assert (s, e, partial) == (27, tae + 3, False)                     # ATT at 27 (nearest), adjacent TAA included
    s2, e2, partial2 = mc.place_cds(window, tab, tae - 6, 0, T5)       # alignment ends 2 codons early: no stop added
    assert s2 == 27 and e2 == tae - 6 and not partial2


def test_place_cds_stop_blocks_upstream_search_and_downstream_fallback():
    window = "TAA" + "GCC" * 4 + "GCC" + "ATG" + "GCC" * 20 + "TAG" + "GCC" * 3
    tab = 3 + 12                             # aligned start at the GCC before the ATG, no start upstream before TAA
    s, e, partial = mc.place_cds(window, tab, tab + 66, 5, T5)       # alignment ends right before the TAG
    assert s == 18 and not partial and window[e - 3:e] == "TAG"       # downstream ATG taken, adjacent TAG included


def test_place_cds_partial_when_no_start():
    window = "GCC" * 30 + "TAA"
    s, e, partial = mc.place_cds(window, 30, 90, 0, T5)
    assert partial and s == 30 and e == 93


def test_trna_name_by_anticodon_then_model():
    assert mc.trna_name("L", "tag", "trnL2.cm") == "trnL1"             # anticodon decides, not the model
    assert mc.trna_name("L", "taa", "trnL1.cm") == "trnL2"
    assert mc.trna_name("S", "gct", "trnS2.cm") == "trnS1" and mc.trna_name("S", "tga", "x") == "trnS2"
    assert mc.trna_name("S", "", "trnS1.cm") == "trnS1"                # no anticodon: the model's number
    assert mc.trna_name("N", "gtt", "trnA.cm") == "trnN"               # MiTFi renames by anticodon (A model, Asn)
    assert mc.trna_name("", "", "sequence.fa_trnW.cmout") == "trnW"


def test_parse_mitfi_selects_tophits_and_copies_and_resolves_overlaps():
    text = "\n".join([
        ">X|5|1541|1601|1|1|72|42.2|2.4E-6|0.0|41|28|GUU|N|trnN.cm|1|1|?", "SEQ", "STRUCT",
        ">X|5|1545|1605|1|1|67|30.0|1E-4|0.0|41|27|UGC|A|trnA.cm|1|1|?", "SEQ", "STRUCT",     # overlaps trnN: dropped
        ">X|5|9750|9812|-1|1|67|43.3|4.2E-6|0.0|25|31|UUC|E|trnA.cm|0|1|?", "SEQ", "STRUCT",  # copy: kept
        ">X|5|9900|9950|1|1|67|20.0|5E-2|0.0|25|0||?|trnQ.cm|0|0|?", "SEQ", "STRUCT",         # neither: dropped
    ])
    feats = mc.parse_mitfi(text, LOG)
    assert [(f.name, f.start, f.end, f.strand, f.anticodon, f.acpos, f.model) for f in feats] == \
        [("trnN", 1540, 1601, 1, "gtt", 27, "trnN"), ("trnE", 9749, 9812, -1, "ttc", 30, "trnA")]


def test_loci_groups_same_strand_hits():
    h = lambda s, e, strand, bits: mc.Hit("cox1", "r", s, e, strand, bits, 1, 10, 500)  # noqa: E731
    L = mc.loci([h(0, 300, 1, 50), h(350, 900, 1, 400), h(5000, 5300, -1, 90), h(1200, 1500, 1, 30)])
    assert [[x.start for x in locus] for locus in L] == [[0, 350], [5000], [1200]]    # best locus first


def test_write_genbank_names(tmp_path):
    feats = [mc.Feature("CDS", "cox1", 0, 90, 1, "ref r", score=10), mc.Feature("rRNA", "rrnL", 100, 200, -1),
             mc.Feature("tRNA", "trnL1", 300, 360, 1, anticodon="tag", acpos=30),
             mc.Feature("tRNA", "trnW", 380, 440, -1, anticodon="tca", acpos=30),
             mc.Feature("CDS", "nad5", 400, 490, 1, partial5=True)]
    gb = mc.write_genbank("A" * 500, feats, 5, "X1", "circular", tmp_path / "x.gbf")
    from Bio import SeqIO
    rec = SeqIO.read(str(gb), "genbank")
    got = [(f.type, f.qualifiers.get("gene", [""])[0], str(f.location)) for f in rec.features if f.type != "source"]
    assert got[0] == ("gene", "COX1", "[0:90](+)") and got[1] == ("CDS", "COX1", "[0:90](+)")
    assert ("rRNA", "l-rRNA", "[100:200](-)") in got and ("tRNA", "trnL1", "[300:360](+)") in got
    assert any(f.type == "CDS" and f.qualifiers["gene"] == ["ND5"] and str(f.location).startswith("[<400")
               for f in rec.features)
    trnas = [f for f in rec.features if f.type == "tRNA"]
    assert trnas[0].qualifiers["product"] == ["tRNA-Leu"]
    assert trnas[0].qualifiers["note"][0] == "anticodon: UAG"           # the MITOS2-style note mito_curate reads
    assert "anticodon" not in trnas[1].qualifiers                        # no positional qualifier: it would not follow re-orientation


def test_window_to_genome_minus_strand():
    """Benchmark 2026-10-22: every minus-strand protein gene was placed at a mirrored position (gs = n - (ws + e))."""
    from Bio.Seq import Seq
    genome = "".join("ACGT"[(i * 7 + i // 3) % 4] for i in range(400))
    ws, we = 100, 300
    rc_window = str(Seq(genome[ws:we]).reverse_complement())
    s, e = 40, 70                                                      # a stretch of the minus-strand window
    gs, ge = mc.window_to_genome(ws, we, -1, s, e)
    assert (gs, ge) == (230, 260)
    assert str(Seq(genome[gs:ge]).reverse_complement()) == rc_window[s:e]
    assert mc.window_to_genome(ws, we, 1, s, e) == (140, 170)


def test_coding_frame_picks_the_frame_without_stops():
    body = "ATG" + ("TTA" + "ATA" + "GCA" + "CTT" + "AAA") * 8 + "TAA"       # the shifted frames hold TAA / TAG codons
    window = "GG" + body + "GCC" * 5                                        # the gene starts at 2
    stops = set(T5.stop_codons)
    end = 2 + len(body) - 3
    n = [sum(1 for c in range(2 + d, end - 2, 3) if window[c:c + 3] in stops) for d in range(3)]
    assert n[0] == 0 and n[1] > 0 and n[2] > 0, n                         # the test sequence is as intended
    assert mc.coding_frame(window, 2, end, T5) == 2
    assert mc.coding_frame(window, 1, end, T5) == 2                        # off by one: corrected
    assert mc.coding_frame(window, 3, end, T5) == 5                        # the other way: next codon of the true frame


def test_rrna_hmm_lookup_by_phylum(tmp_path):
    (tmp_path / "Annelida.rrnS.hmm").write_text("")
    (tmp_path / "Annelida.rrnL.hmm").write_text("")
    assert mc.rrna_hmm(tmp_path, "annelida", "rrnL") == tmp_path / "Annelida.rrnL.hmm"
    assert mc.rrna_hmm(tmp_path, "Mollusca", "rrnL") is None and mc.rrna_hmm(tmp_path, None, "rrnS") is None


def test_clean_protein_for_exonerate():
    """RefSeq proteins carry J / B / Z / U / O and a trailing *: exonerate stops on them (NC_037188.1, nad5 reference
    with J at residue 324)."""
    assert mc.clean_protein("MAJLBZUO*") == "MAXLXXXX"
    assert mc.clean_protein("mkvl") == "MKVL"


def test_location_across_origin_and_minus_partial():
    """Review 2026-10-22: genes across the origin of a circular record are found on seq + seq[:ORIGIN_EXT] and written
    as joins; a 5'-partial minus-strand CDS gets > at its genome end (was a complete location)."""
    n = 1000
    plus = mc.location(mc.Feature("CDS", "cox1", 900, 1150, 1), n)
    assert str(plus) == "join{[900:1000](+), [0:150](+)}"
    minus = mc.location(mc.Feature("tRNA", "trnF", 980, 1030, -1), n)
    assert str(minus) == "join{[0:30](-), [980:1000](-)}"
    assert str(mc.location(mc.Feature("CDS", "nad5", 400, 490, -1, partial5=True), n)) == "[400:>490](-)"
    assert str(mc.location(mc.Feature("CDS", "nad5", 400, 490, 1, partial5=True), n)) == "[<400:490](+)"


def test_parse_mitfi_origin_copy_and_model_name_without_anticodon():
    """A tRNA in the first ORIGIN_EXT nt is reported twice (i and i + n): one feature. Without an anticodon MiTFi
    writes the model name in the amino-acid column: trnD, not Thr (NC_058285.1)."""
    text = "\n".join([
        ">X|5|101|160|1|1|60|40.0|1E-6|0.0|25|-1|---|trnD|trnD.cm|1|0|?", "S", "T",
        ">X|5|1101|1160|1|1|60|40.0|1E-6|0.0|25|-1|---|trnD|trnD.cm|0|1|?", "S", "T",     # same tRNA, appended copy
        ">X|5|501|560|1|1|60|30.0|5E-3|0.0|25|30|UAG|L1|trnL2.cm|0|1|?", "S", "T",        # copy, E > 0.001: dropped
        ">X|5|601|660|1|1|60|35.0|1E-5|0.0|25|30|UAG|L1|trnL2.cm|1|0|?", "S", "T",        # MiTFi's own L1
    ])
    feats = mc.parse_mitfi(text, LOG, 1000)
    assert [(f.name, f.start, f.end) for f in feats] == [("trnD", 100, 160), ("trnL1", 600, 660)]


def test_drop_paralogs_keeps_real_overlaps():
    """Review: with atp8 absent, its reference hit inside atp6 was annotated as atp8 (score 107 vs 991)."""
    f = lambda name, s, e, score: mc.Feature("CDS", name, s, e, 1, score=score)  # noqa: E731
    kept = mc.drop_paralogs([f("atp6", 1000, 1700, 991), f("atp8", 1050, 1212, 107), f("nad4l", 3000, 3300, 400),
                             f("nad4", 3293, 4600, 1500)], 10000, LOG)
    assert [x.name for x in kept] == ["atp6", "nad4l", "nad4"]              # 7-nt nad4l / nad4 overlap is real


def test_exonerate_code_table_string_for_codes_exonerate_lacks():
    assert mc.exonerate_code(5) == "5"
    t33 = mc.exonerate_code(33)                                             # exonerate 2.4.0: no built-in code 33
    assert len(t33) == 64 and t33.startswith("FFLLSSSSYYY*CCWW") and t33 == \
        "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG"


def test_coding_frame_keeps_frame_with_single_stop():
    """exonerate aligns through one premature stop; a short gene's other frames may hold none (review: atp8-size gene
    moved into a nonsense frame)."""
    body = "ATG" + "GCC" * 28 + "TAA" + "GCC" * 22
    assert mc.coding_frame("GG" + body + "GCC" * 3, 2, 2 + len(body), T5) == 2


def test_tbl_hits_model_columns(tmp_path):
    tbl = tmp_path / "x.tbl"
    tbl.write_text("#\nX - rrnL - cm 52 1230 15964 17183 + no 1 0.30 0.0 908.0 1e-200 ! -\n"
                   "X - rrnL - cm 1 40 900 861 - no 1 0.30 0.0 20.0 1e-2 ! -\n")
    h = mc.tbl_hits(tbl)
    assert (h[0].start, h[0].end, h[0].strand, h[0].hmm_from, h[0].hmm_to, h[0].score) == (15963, 17183, 1, 52, 1230, 908.0)
    assert (h[1].start, h[1].end, h[1].strand) == (860, 900, -1)


def test_place_cds_without_downstream_fallback():
    """Echinoderm nad4l (NO_START_DOWN): no start codon at / before the alignment start -> 5'-partial there, not the ATG
    downstream (NC_023772.1: NCBI starts it right after the preceding tRNA)."""
    window = "GCC" * 30 + "ATG" + "GCC" * 20 + "TAA"
    s, e, partial = mc.place_cds(window, 60, 150, 0, T5, downstream=False)
    assert (s, partial) == (60, True)
    s2, _, partial2 = mc.place_cds(window, 60, 150, 0, T5)
    assert (s2, partial2) == (90, False)
    assert "nad4l" in mc.NO_START_DOWN["echinodermata"]


def test_read_fasta_matches_seqio(tmp_path):
    """read_fasta (no Bio.SeqIO, peak memory) gives SeqIO's ids, descriptions and wrapped sequences"""
    from Bio import SeqIO
    fa = tmp_path / "x.fa"
    fa.write_text(">NC_1.1 Foo bar topology=linear\nACGT\nAC GT\n\n>NC_2.1\nTTTT\n")
    got = [(r.id, r.description, str(r.seq)) for r in mc.read_fasta(fa)]
    assert got == [(r.id, r.description, str(r.seq)) for r in SeqIO.parse(str(fa), "fasta")]
    assert got[0] == ("NC_1.1", "NC_1.1 Foo bar topology=linear", "ACGTACGT")


def test_orf_before_neighbour_both_strands():
    """ADJACENT: atp8 = the ORF ending in the start of atp6 (stop 7 nt into atp6), start codon nearest the typical length"""
    from Bio.Data import CodonTable
    from Bio.Seq import Seq
    table = CodonTable.unambiguous_dna_by_id[5]
    atp8 = "ATG" + "GCT" * 50 + "TAA"                           # 156 nt
    genome = "TAG" * 10 + atp8[:-7] + atp8[-7:] + "GCT" * 60 + "TAA" * 5
    a6 = 30 + len(atp8) - 7                                     # atp6 starts 7 nt before the atp8 stop ends
    nb = mc.Feature("CDS", "atp6", a6, a6 + 180, 1)
    rc = str(Seq(genome).reverse_complement())
    f = mc.orf_before(genome, rc, nb, "atp8", table, 162, 90, 240)
    assert (f.start, f.end, f.strand) == (30, 30 + len(atp8), 1)
    N = len(genome)                                             # the same on the minus strand
    nb_m = mc.Feature("CDS", "atp6", N - nb.end, N - nb.start, -1)
    g = mc.orf_before(rc, genome, nb_m, "atp8", table, 162, 90, 240)
    assert (g.start, g.end, g.strand) == (N - f.end, N - f.start, -1)


def test_rrna_weak_hmm_span_falls_back_to_the_covariance_model(tmp_path, monkeypatch):
    """A divergent genome (tunicates against the vertebrate-dominated Chordata model: NC_009833.1, NC_013561.1) gives
    the HMM only a partial span; the covariance model, which finds these rRNAs, is tried instead of dropping the gene."""
    import logging
    hmm = tmp_path / "Chordata.rrnS.hmm"
    hmm.write_text("")
    (tmp_path / "Chordata.rrnL.hmm").write_text("")
    weak = mc.mito_rrna_hmm.Hit(start=100, end=300, strand=1, hmm_from=10, hmm_to=200, score=50.0)   # 191 of 900
    good = mc.mito_rrna_hmm.Hit(start=1000, end=1800, strand=1, hmm_from=1, hmm_to=880, score=400.0)
    monkeypatch.setattr(mc.mito_rrna_hmm, "nhmmer", lambda h, fa: [weak])
    monkeypatch.setattr(mc.mito_rrna_hmm, "model_len", lambda h: 900)
    monkeypatch.setattr(mc.mito_rrna_hmm, "best_spans", lambda hits, n: (hits[0] if hits else None, None))
    monkeypatch.setattr(mc, "cmsearch", lambda *a, **k: tmp_path / "x.out")
    monkeypatch.setattr(mc, "tbl_hits", lambda tbl: [good])
    monkeypatch.setattr(mc, "cm_len", lambda cm: 900)
    feats = mc.annotate_rrnas(tmp_path / "g.fa", tmp_path, {"cmsearch": "cmsearch"}, tmp_path,
                              logging.getLogger("t"), phylum="Chordata", rrna_lib=tmp_path)
    assert [(f.name, f.start, f.end) for f in feats] == [("rrnS", 1000, 1800), ("rrnL", 1000, 1800)]
    assert all("cmsearch" in f.note for f in feats)


def test_rrna_second_copy_added_when_it_scores_like_the_gene(tmp_path, monkeypatch):
    """Duplicated rRNAs (NC_044480.1 rrnS + rrnL, NC_008046.1 rrnS: both copies in NCBI) were annotated once: the
    runner-up span over the same model columns is a second copy when it scores >= COPY_FRAC of the gene."""
    import logging
    (tmp_path / "Chordata.rrnS.hmm").write_text("")
    (tmp_path / "Chordata.rrnL.hmm").write_text("")
    gene = mc.mito_rrna_hmm.Hit(start=100, end=1000, strand=1, hmm_from=1, hmm_to=880, score=700.0)
    copy = mc.mito_rrna_hmm.Hit(start=9000, end=9900, strand=1, hmm_from=1, hmm_to=880, score=690.0)
    weak = mc.mito_rrna_hmm.Hit(start=5000, end=5900, strand=-1, hmm_from=1, hmm_to=880, score=300.0)
    spans = {"Chordata.rrnS.hmm": (gene, copy), "Chordata.rrnL.hmm": (gene, weak)}
    monkeypatch.setattr(mc.mito_rrna_hmm, "nhmmer", lambda h, fa: [h.name])
    monkeypatch.setattr(mc.mito_rrna_hmm, "model_len", lambda h: 900)
    monkeypatch.setattr(mc.mito_rrna_hmm, "best_spans", lambda hits, n: spans[hits[0]])
    feats = mc.annotate_rrnas(tmp_path / "g.fa", tmp_path, {}, tmp_path, logging.getLogger("t"), phylum="Chordata",
                              rrna_lib=tmp_path)
    assert [(f.name, f.start, f.end) for f in feats] == [("rrnS", 100, 1000), ("rrnS", 9000, 9900), ("rrnL", 100, 1000)]


def test_trna_rescue_adds_mitfi_calls_of_families_ml_did_not_call():
    """MitoM-ML misses weak tRNAs of divergent genomes (acanthocephalans: 17 % recall); MiTFi's calls of the families it
    did not call are added, lowest E first, one per family, when they overlap no call by more than 10 nt."""
    t = lambda name, s, e, st=1, ev="1.0e-05": mc.Feature("tRNA", name, s, e, st, f"MiTFi {name}.cm bits 30.0 E {ev}")  # noqa: E731
    ml = [mc.Feature("tRNA", "trnA", 100, 170, 1, "MitoM-ML p 0.9"), mc.Feature("tRNA", "trnC", 300, 370, -1, "MitoM-ML p 0.8"),
          mc.Feature("tRNA", "trnW", 990, 1050, 1, "MitoM-ML p 0.9")]          # across the origin of a 1,000-nt record
    mitfi = [t("trnA", 105, 175),                         # family already called
             t("trnD", 800, 870, ev="5.0e-01"), t("trnD", 500, 570, ev="5.0e-02"),   # one per family, lowest E first
             t("trnE", 160, 230, -1),                     # 10 nt on trnA (other strand): allowed
             t("trnF", 150, 220),                         # 20 nt on trnA: rejected
             t("trnY", 0, 40)]                            # 40 nt on trnW across the origin: rejected
    out = mc.rescue_trnas(ml, mitfi, 1000)
    assert [(f.name, f.start) for f in out[3:]] == [("trnE", 160), ("trnD", 500)]
    assert out[:3] == ml


def test_vertebrate_trnas_take_mitfi_ends_where_mitfi_calls_the_same_trna():
    """The learned ends never saw a vertebrate; MiTFi's are NCBI's more often there (held out: 91.3 -> 92.7 % exact)."""
    ml = [mc.Feature("tRNA", "trnF", 100, 170, 1, "MitoM-ML p 0.9", anticodon="gaa"),
          mc.Feature("tRNA", "trnV", 300, 370, 1, "MitoM-ML p 0.9"), mc.Feature("tRNA", "trnL2", 500, 575, -1, "MitoM-ML p 0.9")]
    mitfi = [mc.Feature("tRNA", "trnF", 99, 169, 1, "MiTFi trnF.cm bits 50.0 E 1.0e-10"),
             mc.Feature("tRNA", "trnV", 300, 370, -1, "MiTFi trnV.cm bits 50.0 E 1.0e-10"),     # other strand: not the same
             mc.Feature("tRNA", "trnL1", 501, 575, -1, "MiTFi trnL1.cm bits 50.0 E 1.0e-10")]   # other name: not the same
    out = mc.mitfi_ends(ml, mitfi, 16000)
    assert [(f.name, f.start, f.end) for f in out] == [("trnF", 99, 169), ("trnV", 300, 370), ("trnL2", 500, 575)]
    assert out[0].anticodon == "gaa" and out[0].note.endswith("ends by MiTFi") and out[1].note == "MitoM-ML p 0.9"
    assert mc.vertebrate("Eukaryota;Metazoa;Chordata;Craniata;Vertebrata;Euteleostomi") and mc.vertebrate("Metazoa;Chordata")
    assert not mc.vertebrate("Eukaryota;Metazoa;Chordata;Tunicata;Ascidiacea") and not mc.vertebrate("Metazoa;Annelida")
