import mito_nuc_lib as L


def hsp(q0, q1, s0, s1, bits):
    return (q0, q1, s0, s1, bits, "A" * (q1 - q0 + 1), "A" * (q1 - q0 + 1))


def test_tandem_copies_are_not_hybridised():
    # reference: 18S 1-1800, 5.8S 2000-2150, 28S 2400-5000. Two tandem copies on one contig (unit 6000 nt).
    c1 = [hsp(1, 1800, 1, 1800, 900), hsp(2001, 2150, 2000, 2149, 200)]
    c2 = [hsp(6401, 9000, 2400, 5000, 1200)]  # only 28S strong in copy 2
    c1b = [hsp(6001, 7800, 1, 1800, 400), hsp(8001, 8150, 2000, 2149, 100)]
    out = L.chains(c1 + c1b + c2)
    spans = [(c[0][0], c[-1][1]) for c in out]
    assert all(not (a <= 2001 and b >= 6401) for a, b in spans), spans
    assert len(out) >= 2
