import random

import mito_intron_exon as ie


def test_locate_finds_planted_exon_and_rejects_noise():
    rng = random.Random(3)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    prot = "".join(rng.choice(aa) for _ in range(60))
    homs = []
    for _ in range(8):
        h = list(prot)
        for i in rng.sample(range(60), 30):
            h[i] = rng.choice(aa)
        homs.append("".join(h))
    codon = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA",
             "L": "CTT", "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT",
             "W": "TGG", "Y": "TAT"}
    flank = "".join(rng.choice("ACGT") for _ in range(300))
    exon = "".join(codon[x] for x in prot[:40])
    use = flank + exon + flank
    got = ie.locate(use, prot, homs, 0, 40, 0, len(use), 5, [], 142)
    assert got and any(abs(s - 300) <= 9 for s, _ in got)
    noise = "".join(rng.choice("ACGT") for _ in range(900))
    assert ie.locate(noise, prot, homs, 0, 40, 0, len(noise), 5, [], 142) == []


def test_complete_chain_without_holes_is_unchanged():
    ex = [(1, 90)]
    assert ie.complete_chain("A" * 100, ex, [(0, 30, 1, 90)], "M" * 30, [], 5, [], 5000, 142) == (ex, "")
