"""mito_rna.py: relative depth, library strand and poly(A) 3' ends from SAM lines."""
import random

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import mito_rna as R

GENES = [("COX1", 0, 600), ("ND2", 700, 1300), ("orf99", 1400, 1700)]


def record():
    rnd = random.Random(1)
    rec = SeqRecord(Seq("".join(rnd.choice("ACGT") for _ in range(2000))), id="M.1")
    rec.features = [SeqFeature(FeatureLocation(s, e, 1), type="CDS", qualifiers={"gene": [g]}) for g, s, e in GENES]
    rec.features.append(SeqFeature(FeatureLocation(1700, 1771, 1), type="tRNA", qualifiers={"gene": ["trnK"]}))
    return rec


def read(name, pos, length, flag, seq=None, cigar=None):
    seq = seq or "A" * length
    return "\t".join([name, str(flag), "M.1", str(pos + 1), "60", cigar or f"{length}M", "*", "0", "0", seq,
                      "I" * len(seq)])


def test_reverse_stranded_library_relative_depth_and_polya():
    rec = record()
    sam = []
    k = 0
    for g, s, e, per in (("COX1", 0, 600, 200), ("ND2", 700, 1300, 100), ("orf99", 1400, 1700, 50)):
        for i in range(per):                               # read 1 antisense (flag 16 + 64): reverse-stranded
            pos = s + (i * 7) % (e - s - 60)
            sam.append(read(f"r{k}", pos, 60, 16 | 64))
            k += 1
    for i in range(5):                                     # poly(A) at the end of ND2 (sense read 2, 3' A clip)
        seq = "C" * 50 + "AAAAAAAAAA"
        sam.append(read(f"p{i}", 1250, 60, 128, seq, "50M10S"))
    rows, sites, info = R.analyse(rec, ["@HD\tVN:1.6"] + sam)
    assert info["stranded"] and info["reverse"]
    by = {r["feature"]: r for r in rows}
    assert by["ND2"]["sense_share"] == "1.00"              # strands swapped: reads count on the gene's strand
    assert float(by["orf99"]["relative_depth"]) < float(by["ND2"]["relative_depth"])
    assert by["orf99"]["transcribed_like_genes"] == "yes"
    assert any(s["feature_ending_here"] == "ND2" and s["reads"] >= 5 for s in sites)


def test_low_mapq_reads_are_not_used():
    rec = record()
    line = read("x", 10, 60, 0).split("\t")
    line[4] = "3"
    rows, sites, info = R.analyse(rec, ["\t".join(line)])
    assert info["reads_used"] == 0 and info["reads_low_mapq"] == 1


def test_too_few_reads_give_no_expression_calls():
    rec = record()
    rows, sites, info = R.analyse(rec, [read("a", 10, 60, 16 | 64), read("b", 720, 60, 16 | 64)])
    assert info["pcg_median_depth"] < R.MIN_MEDIAN_DEPTH and all(r["transcribed_like_genes"] == "" for r in rows)


def test_internal_priming_rule():
    import mito_rna
    g = "CCGTCGTCGG" + "AAAAAAAAAG" + "CCGTCGTCGGCGTCGT"     # site at the G after an A9 run (plus strand)
    assert mito_rna.internal_priming(g, 19, 1)
    assert not mito_rna.internal_priming(g, 9, 1)            # ends in CCGTCGTCGG: a real end
    rc = g[::-1].translate(str.maketrans("ACGT", "TGCA"))     # the same site on the minus strand
    assert mito_rna.internal_priming(rc, len(g) - 1 - 19, -1)


def test_internal_priming_site_not_counted(tmp_path):
    """A site on a genomic A-run is listed (internal_priming yes) but not counted for the feature."""
    import mito_rna
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    body = "ACGTGCTAGCTAGGCTACGATCGATCGGCTAGCTAG" * 5 + "AAAAAAAAA"     # 189 nt ending on A9
    rec = SeqRecord(Seq(body + "GCGCGCGCTTCCGG" * 10), id="t")
    rec.features.append(SeqFeature(FeatureLocation(0, len(body)), type="CDS",
                                   qualifiers={"gene": ["cox1"], "product": ["x"]}))
    end = len(body)                         # the read aligns up to the A-run's end, then 20 clipped A's
    read = body[end - 50:end] + "A" * 20
    sam = ["@SQ\tSN:t\tLN:%d" % len(rec.seq)] + [
        f"r{i}\t0\tt\t{end - 49}\t60\t50M20S\t*\t0\t0\t{read}\t*" for i in range(5)]
    rows, sites, info = mito_rna.analyse(rec, sam)
    assert info["polya_sites"] == 0 and info["internal_priming_sites"] == 1
    assert sites[0]["internal_priming"] == "yes"
    assert next(r for r in rows if r["feature"] == "cox1")["polya_reads_at_3prime_end"] == 0


def test_refine_moves_end_to_polya_site_with_incomplete_stop():
    """A strong poly(A) site inside the gene's end, leaving an in-frame TA, moves the 3' end; a site at the
    annotated end of another gene confirms it."""
    import mito_rna
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    cds = "ATG" + "GCT" * 30 + "TA" + "GCCGCCGCCTAA"              # site after ...TA (95 nt); annotated end at TAA
    seq = cds + "C" * 40 + "ATG" + "GCC" * 20 + "TAA" + "C" * 30
    rec = SeqRecord(Seq(seq), id="t")
    for name, s, e in (("ND6", 0, len(cds)), ("COX1", len(cds) + 40, len(cds) + 40 + 66)):
        q = {"gene": [name], "product": ["x"], "transl_table": ["5"], "codon_start": ["1"]}
        rec.features.append(SeqFeature(FeatureLocation(s, e, 1), type="gene", qualifiers={"gene": [name]}))
        rec.features.append(SeqFeature(FeatureLocation(s, e, 1), type="CDS", qualifiers=q))
    sites = [{"position": "95", "strand": "+", "reads": "20", "internal_priming": ""},
             {"position": str(len(cds) + 40 + 66), "strand": "+", "reads": "5", "internal_priming": ""}]
    corr = mito_rna.refine(rec, [], sites, stranded=False, run="SRR1")
    assert [(c["gene"], c["new_end"], c["incomplete_stop"]) for c in corr] == [("ND6", 95, "TA")]
    nd6 = next(f for f in rec.features if f.type == "CDS" and f.qualifiers["gene"] == ["ND6"])
    assert int(nd6.location.end) == 95 and nd6.qualifiers["transl_except"] == ["(pos:94..95,aa:TERM)"]
    cox1 = next(f for f in rec.features if f.type == "CDS" and f.qualifiers["gene"] == ["COX1"])
    assert any("confirmed by 5 poly(A)" in n for n in cox1.qualifiers["note"])


def test_refine_utr_extension_and_duplication():
    """Complete stop + site 20 nt beyond: confirmed (3' UTR). No stop + site beyond leaving an in-frame T: extended.
    A gene present twice: duplicated, depth not measurable."""
    import mito_rna
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord
    a = "ATG" + "GCTCGA" * 15 + "TAA"                       # 96 nt, complete stop
    b = "ATG" + "GCCAGC" * 15 + "GTC"                       # 96 nt, no stop; then CCG CAT + T = incomplete stop
    dupl = "ATGACCTTGGAGCAAGTCCTGAACGGTTCAGCGATCCATTGA"
    seq = a + "C" * 40 + b + "CCGCATT" + "G" * 40 + dupl + "C" * 30 + dupl + "A" * 5
    rec = SeqRecord(Seq(seq), id="t")
    loc = {"A": (0, 96), "B": (136, 232), "D": (279, 279 + len(dupl))}
    for name, (s, e) in loc.items():
        q = {"gene": [name], "product": ["x"], "transl_table": ["5"], "codon_start": ["1"]}
        rec.features.append(SeqFeature(FeatureLocation(s, e, 1), type="CDS", qualifiers=q))
    sites = [{"position": "116", "strand": "+", "reads": "8", "internal_priming": ""},      # 20 nt after A's stop
             {"position": "239", "strand": "+", "reads": "15", "internal_priming": ""}]     # B + CCGCATT
    corr = mito_rna.refine(rec, [], sites, stranded=False)
    assert [(c["gene"], c["kind"], c["new_end"], c["incomplete_stop"]) for c in corr] == [("B", "extended", 239, "T")]
    f = {x.qualifiers["gene"][0]: x for x in rec.features}
    assert any("ends 20 nt after the stop codon" in n for n in f["A"].qualifiers["note"])
    assert any("occurs elsewhere" in n for n in f["D"].qualifiers["note"])


def test_minus_strand_polya_never_shortens_a_plus_strand_cds_even_when_unstranded():
    """A poly(A) site's strand is the transcript's whatever the library is: 20 minus-strand sites at 95 leave the CDS."""
    seq = "ATG" + "GCT" * 30 + "TAT" + "GCT" * 3 + "TAA"      # 108 nt, complete; a cut at 95 would end ...TA
    rec = SeqRecord(Seq(seq), id="M.1")
    rec.features = [SeqFeature(FeatureLocation(0, 108, 1), type="CDS", qualifiers={"gene": ["X"], "transl_table": ["5"]})]
    sites = [{"position": "95", "strand": "-", "reads": "20", "internal_priming": ""}]
    assert R.refine(rec, [], sites, stranded=False) == []
    assert (rec.features[0].location.start, rec.features[0].location.end) == (0, 108)


def test_joined_feature_measured_over_its_parts_and_other_contigs_ignored():
    from Bio.SeqFeature import CompoundLocation
    rnd = random.Random(2)
    rec = SeqRecord(Seq("".join(rnd.choice("ACGT") for _ in range(1000))), id="M.1")
    loc = CompoundLocation([FeatureLocation(900, 1000, 1), FeatureLocation(0, 100, 1)])   # across the origin
    rec.features = [SeqFeature(loc, type="CDS", qualifiers={"gene": ["COX1"]})]
    sam = ["@SQ\tSN:M.1\tLN:1000", "@SQ\tSN:X\tLN:900"]
    for i in range(50):
        sam += [read(f"a{i}", 900, 100, 0), read(f"b{i}", 0, 100, 0)]
    sam += [read(f"f{i}", 0, 50, 0, None, "50M10S").replace("M.1", "X") for i in range(10)]   # foreign contig
    rows, sites, info = R.analyse(rec, sam, "M.1")
    r = rows[0]
    assert r["mean_depth"] == "50.0" and r["covered_share"] == "1.00"
    assert info["reads_used"] == 100
