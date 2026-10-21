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
    assert [(f.name, f.start, f.end, f.strand, f.anticodon, f.acpos) for f in feats] == \
        [("trnN", 1540, 1601, 1, "gtt", 27), ("trnE", 9749, 9812, -1, "ttc", 30)]


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
    assert trnas[0].qualifiers["anticodon"] == ["(pos:331..333,aa:Leu,seq:tag)"]            # NCBI form, genome coordinates
    assert trnas[1].qualifiers["anticodon"] == ["(pos:complement(408..410),aa:Trp,seq:tca)"]
