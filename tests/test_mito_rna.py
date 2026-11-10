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
