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


def test_trna_second_opinion_keeps_sure_calls_and_weak_ones_mitfi_also_calls():
    """p < 0.9 needs MiTFi's call of the same tRNA (37 % of such calls were wrong without it on the tuning halves);
    a weaker second call of a family needs p >= 0.9 even with MiTFi; without MiTFi only the second-copy rule applies."""
    t = lambda name, s, e, p, st=1: mc.Feature("tRNA", name, s, e, st, f"MitoM-ML p {p}", score=p)   # noqa: E731
    ml = [t("trnA", 100, 170, 0.95), t("trnC", 300, 370, 0.6), t("trnD", 500, 570, 0.7),
          t("trnM", 700, 770, 0.97), t("trnM", 900, 970, 0.85), t("trnS1", 1100, 1170, 0.99), t("trnS1", 1300, 1370, 0.92)]
    mitfi = [mc.Feature("tRNA", "trnC", 302, 371, 1, "MiTFi"), mc.Feature("tRNA", "trnM", 900, 970, 1, "MiTFi"),
             mc.Feature("tRNA", "trnD", 500, 570, -1, "MiTFi")]                     # other strand: no support
    keep = lambda fs: [(f.name, f.start) for f in fs]                              # noqa: E731
    assert keep(mc.trna_second_opinion(ml, mitfi, 16000)) == [("trnA", 100), ("trnC", 300), ("trnM", 700),
                                                              ("trnS1", 1100), ("trnS1", 1300)]
    assert keep(mc.trna_second_opinion(ml, None, 16000)) == [("trnA", 100), ("trnC", 300), ("trnD", 500), ("trnM", 700),
                                                             ("trnS1", 1100), ("trnS1", 1300)]


def test_mitfi_path_runs_cmsearch_and_mitfi_as_mitos2(tmp_path, monkeypatch):
    """E-values at the record length (cmsearch -Z: the record is searched with ORIGIN_EXT nt appended) and MiTFi
    -overlap 50 -method l, as runmitos runs them: NC_060482.1 trnD at E 0.096 (0.12 with the appended nt) and
    NC_060484.1 trnL2 overlapping trnP by 12 nt were lost; with both, MiTFi's output equals MITOS2's own tRNAout.nc."""
    calls = []

    class P:
        stdout = ""

    def fake_run(cmd, **kw):
        calls.append(cmd)
        return P()
    monkeypatch.setattr(mc.mito_trna, "run_or_raise", fake_run)      # the tRNA step is in mito_trna.py
    (tmp_path / "ncRNA").mkdir()
    for m in ("trnA", "trnD"):
        (tmp_path / "ncRNA" / f"{m}.cm").write_text("")
    fa = tmp_path / "g.fa"
    fa.write_text(">g\nACGT\n")
    tools = {"cmsearch": "cmsearch", "java": "java", "mitfi": (Path("mitfi_import.jar"), Path("codes"))}
    mc.annotate_trnas("A" * 16845, fa, tmp_path, 5, 1, tools, tmp_path, LOG, n=13845)
    cms = [c for c in calls if c[0] == "cmsearch"]
    assert len(cms) == 4 and all(c[c.index("-Z") + 1] == f"{2 * 13845 / 1e6:.6f}" for c in cms)   # fast + sensitive
    java = [c for c in calls if c[0] == "java"]
    assert java and all(c[c.index("-overlap") + 1] == "50" and c[c.index("-method") + 1] == "l" for c in java)


def test_gene_free_leaves_out_gene_interiors_only():
    """The sensitive pass searches all but the gene interiors: the core's gene ends may cover a tRNA (arthropod trnH
    inside the core's nad4 end, abutting nad5: NC_069264.1, NC_073524.1)."""
    genes = [mc.Feature("CDS", "cox1", 100, 1000, 1), mc.Feature("rRNA", "rrnL", 1200, 2500, 1)]
    # 3,000-nt circular record searched with 500 nt appended: cox1 covers 3100..3500 again
    assert mc.gene_free(genes, 3000, 3500, 60) == [(0, 160), (940, 1260), (2440, 3160)]
    abutting = [mc.Feature("CDS", "nad5", 100, 1000, -1), mc.Feature("CDS", "nad4", 1000, 2000, -1)]
    assert mc.gene_free(abutting, 3000, 3000, 60) == [(0, 160), (940, 1060), (1940, 3000)]
    short = [mc.Feature("CDS", "nad6", 100, 200, 1), mc.Feature("CDS", "cob", 250, 1000, -1)]   # nad6 < 2 x 60 nt
    assert mc.gene_free(short, 3000, 3000, 60) == [(0, 310), (940, 3000)]
    assert mc.gene_free([], 3000, 3500) == [(0, 3500)]


CMOUT = """Hit scores:
 rank     E-value  score  bias  sequence  start    end   mdl trunc   gc  description
 ----   --------- ------ -----  -------- ------ ------   --- ----- ----  -----------
  (1) !   3.5e-05   38.9   0.0  s15364      124     72 -  cm    no 0.43  -
 ------ inclusion threshold ------
  (2) ?     0.033   22.1   0.0  s1473        72    124 +  cm    no 0.43  -


Hit alignments:
>> s15364
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       cyksc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ------ ----- ----
  (1) !   3.5e-05   38.9   0.0  cm        1       68 []         124          72 - ..   36.4    no 0.43

             (((((((,,<<<<_____>>>>,<<<<<_______>>>>>,,,~~~~~~))))))): CS
    trnG   1 auuuuuuUAGUAuAAuuAgUACAuuUGaCUUCCAAuCAaaaAG*[17]*aaaaaauA 68
             A: U::  AGUAU  +  GUACA: ::A+UUCCA U:: :A+G       ::A :UA
  s15364 124 AGUUAGGCAGUAU-UAGGGUACAGGCAAUUUCCAUUUGUCAGG*[ 3]*ACUAGCUA 72

>> s1473
 rank     E-value  score  bias mdl mdl from   mdl to       seq from      seq to       cyksc trunc   gc
 ----   --------- ------ ----- --- -------- --------    ----------- -----------      ------ ----- ----
  (2) ?     0.033   22.1   0.0  cm        1       68 []          72         124 + ..   19.4    no 0.43

            (((((((,,<<<<_____>>>>,<<<<<_______>>>>>,,,~~~~~~))))))): CS
   trnG   1 auuuuuuUAGUAuAAuuAgUACAuuUGaCUUCCAAuCAaaaAG*[17]*aaaaaauA 68
            A: U::  AGUAU  +  GUACA: ::A+UUCCA U:: :A+G       ::A :UA
  s1473  72 AGUUAGGCAGUAU-UAGGGUACAGGCAAUUUCCAUUUGUCAGG*[ 3]*ACUAGCUA 124
"""


def test_shift_cmout_moves_stretch_hits_to_genome_coordinates():
    """MiTFi reads the hit table, the alignment headers and the alignment rows; after the shift a stretch search reads
    like a search of the whole sequence (NC_060482.1: trnG 1545..1597 as NCBI, found in stretch s1473 at 72..124)."""
    out = mc.shift_cmout(CMOUT, "S00001").splitlines()
    rows = [line.split() for line in out if line.startswith("  (")]
    assert [r[5:9] for r in rows[:2]] == [["S00001", "15488", "15436", "-"], ["S00001", "1545", "1597", "+"]]  # table
    assert [r[9:12] for r in rows[2:]] == [["15488", "15436", "-"], ["1545", "1597", "+"]]          # alignment headers
    assert [line.split()[:2] for line in out if line.startswith(">>")] == [[">>", "S00001"], [">>", "S00001"]]
    aln = [line.split() for line in out if line.lstrip().startswith("S00001")]
    assert [(a[1], a[-1]) for a in aln] == [("15488", "15436"), ("1545", "1597")]
    model = [line for line in CMOUT.splitlines() if "trnG   1" in line or line.endswith(": CS")]
    assert [line for line in out if "trnG   1" in line or line.endswith(": CS")] == model    # the rest untouched


def test_trna_sensitive_pass_searches_the_gene_free_stretches_with_the_cm_filter(tmp_path, monkeypatch):
    """--max on the whole record took minutes (NC_014885.1: 17 models, 297 s with 2 threads); the sensitive pass now
    runs cmsearch --nohmm outside the protein gene / rRNA interiors, E-values still for the whole record."""
    calls, fastas = [], {}

    class P:
        stdout = ""

    def fake_run(cmd, **kw):
        calls.append(cmd)
        if cmd[0] == "cmsearch":
            fastas[cmd[-1]] = Path(cmd[-1]).read_text()
        return P()
    monkeypatch.setattr(mc.mito_trna, "run_or_raise", fake_run)      # the tRNA step is in mito_trna.py
    (tmp_path / "ncRNA").mkdir()
    (tmp_path / "ncRNA" / "trnD.cm").write_text("")
    genome = "A" * 3500                                   # 3,000-nt circular record + 500 nt appended
    fa = tmp_path / "S00001.fa"
    fa.write_text(f">S00001\n{genome}\n")
    tools = {"cmsearch": "cmsearch", "java": "java", "mitfi": (Path("mitfi_import.jar"), Path("codes"))}
    genes = [mc.Feature("CDS", "cox1", 100, 1000, 1), mc.Feature("rRNA", "rrnL", 1200, 2500, 1)]
    mc.annotate_trnas(genome, fa, tmp_path, 5, 1, tools, tmp_path, LOG, n=3000, genes=genes)
    fast, sens = [c for c in calls if c[0] == "cmsearch"]
    assert "--noF4b" in fast and fast[-1] == str(fa)
    assert "--nohmm" in sens and "--notrunc" in sens and "--max" not in sens
    assert sens[sens.index("-Z") + 1] == f"{2 * 3000 / 1e6:.6f}" == fast[fast.index("-Z") + 1]
    stretches = fastas[sens[-1]].split()                  # all but the gene interiors (100 nt of each end searched)
    assert stretches[0::2] == [">s0", ">s900", ">s2400"] and [len(s) for s in stretches[1::2]] == [200, 400, 800]


def test_parse_mitfi_keeps_trnas_overlapping_as_mitos2_allows():
    """NC_060484.1: trnL2 8233..8305 and trnP 8294..8357 overlap by 12 nt, both top hits; MITOS2 keeps both."""
    text = "\n".join([
        ">X|5|8294|8357|1|1|67|40.8|1.7E-6|0.0|36|29|UGG|P|trnP.cm|1|1|l", "SEQ", "STRUCT",
        ">X|5|8233|8305|1|1|71|35.7|1.3E-5|0.0|37|29|UAA|L2|trnL2.cm|1|1|l", "SEQ", "STRUCT",
    ])
    assert [f.name for f in mc.parse_mitfi(text, LOG)] == ["trnL2", "trnP"]


def test_rrna_method_cm_uses_only_the_mitos2_covariance_models(tmp_path, monkeypatch):
    """--rrna_method cm: MITOS2's rrnS / rrnL covariance models only (no nhmmer); hmm (default): the per-phylum HMM
    first. Same accuracy on 450 leave-out records (2026-10-28); hmm is faster."""
    calls = []
    monkeypatch.setattr(mc, "rrna_hmm", lambda lib, phylum, model: tmp_path / f"{model}.hmm")
    monkeypatch.setattr(mc.mito_rrna_hmm, "nhmmer", lambda h, fa: calls.append(("nhmmer", h.name)) or [])
    monkeypatch.setattr(mc.mito_rrna_hmm, "best_spans", lambda hits, n: (None, None))
    monkeypatch.setattr(mc.mito_rrna_hmm, "model_len", lambda h: 900)
    monkeypatch.setattr(mc, "cmsearch", lambda cm, model, fa, out, sensitive, log, rrna=False, z=None:
                        calls.append(("cmsearch", model.name)) or out)
    monkeypatch.setattr(mc, "tbl_hits", lambda tbl: [])
    monkeypatch.setattr(mc, "cm_len", lambda cm: 900)
    tools = {"cmsearch": "cmsearch"}
    mc.annotate_rrnas(tmp_path / "g.fa", tmp_path, tools, tmp_path, LOG, "Annelida", None, None, "cm")
    assert calls == [("cmsearch", "rrnS.cm"), ("cmsearch", "rrnL.cm")]
    calls.clear()
    mc.annotate_rrnas(tmp_path / "g.fa", tmp_path, tools, tmp_path, LOG, "Annelida", None, None)
    assert calls[0] == ("nhmmer", "rrnS.hmm") and ("nhmmer", "rrnL.hmm") in calls


def test_pipeline_passes_rrna_method_to_the_core(tmp_path):
    import argparse
    import mitom_pipeline as bm
    a = argparse.Namespace(genetic_code=5, threads=2, mitos_refdir=None, mitos_bin=None, phylum=None, rrna_hmm_lib=None,
                           rrna_method="cm")
    cmd = bm.native_core_cmd(tmp_path / "x.fa", tmp_path, "S00001", a)
    assert cmd[cmd.index("--rrna_method") + 1] == "cm"
    cmd = bm.native_core_cmd(tmp_path / "x.fa", tmp_path, "S00001", {"genetic_code": 5, "threads": 2, "rrna_method": "hmm"})
    assert cmd[cmd.index("--rrna_method") + 1] == "hmm"


def test_single_second_copy_kept_only_when_identical():
    """2026-10-29: one gene with a second locus is kept when the copy is >= COPY_IDENTITY identical (length within 3 nt),
    as NCBI annotates such copies (heron nad6, Hexamermis atp6); a diverged or longer second locus is still dropped."""
    g = "ATG" + "GCT" * 30 + "TAA"
    genome = "C" * 50 + g + "C" * 50 + g + "C" * 50
    a = mc.Feature("CDS", "atp8", 50, 50 + len(g), 1)
    b = mc.Feature("CDS", "atp8", 100 + len(g), 100 + 2 * len(g), 1)
    assert mc.copy_identity(genome, a, b) == 1.0
    minus = mc.Feature("CDS", "atp8", 100 + len(g), 100 + 2 * len(g), -1)        # other strand: not the same sequence
    assert mc.copy_identity(genome, a, minus) < mc.COPY_IDENTITY
    longer = mc.Feature("CDS", "atp8", 100 + len(g), 100 + 2 * len(g) + 6, 1)
    assert mc.copy_identity(genome, a, longer) == 0.0
    diverged = genome[:100 + len(g)] + "ATG" + "GCA" * 30 + "TAA" + genome[100 + 2 * len(g):]
    assert 0 < mc.copy_identity(diverged, a, b) < mc.COPY_IDENTITY


def test_trna_step_lives_in_mito_trna(monkeypatch):
    """The core's tRNA step is mito_trna.core_trnas (its names stay importable from mito_core): MitoM-ML calls with
    p < 0.9 need MiTFi's call of the same tRNA, MiTFi's calls of the families MitoM-ML did not call are added; without
    MiTFi every MitoM-ML call stays; without MitoM-ML the tRNAs are MiTFi's."""
    import mito_common
    import mito_trna
    assert mc.annotate_trnas is mito_trna.annotate_trnas and mc.trna_name is mito_trna.mitfi_trna_name
    assert mc.Feature is mito_common.CoreFeature
    ml = [mc.Feature("tRNA", "trnA", 100, 170, 1, "MitoM-ML", score=0.95),
          mc.Feature("tRNA", "trnC", 300, 370, 1, "MitoM-ML", score=0.5)]
    mitfi = [mc.Feature("tRNA", "trnD", 500, 570, 1, "MiTFi trnD.cm bits 30.0 E 1.0e-05", score=30.0)]
    monkeypatch.setattr(mito_trna, "annotate_trnas_ml", lambda *a: list(ml))
    monkeypatch.setattr(mito_trna, "annotate_trnas", lambda *a, **k: list(mitfi))

    def names(tools, use_ml):
        return [f.name for f in mito_trna.core_trnas("A" * 16000, "A" * 16000, Path("g.fa"), Path("."), 5, 1, tools,
                                                     Path("."), LOG, 16000, None, "Metazoa;Annelida", True, use_ml, [])]
    assert names({"mitfi": ("jar", "codes")}, True) == ["trnA", "trnD"]
    assert names({}, True) == ["trnA", "trnC"]
    assert names({"mitfi": ("jar", "codes")}, False) == ["trnD"]
