"""Score rRNA ends of tools against experimentally mapped mature-rRNA ends (rrna_end_anchors.tsv, literature survey).

ANCHORS: (accession, gene, end, accepted 1-based positions on the accession, provenance class); an interval
(a, b) accepts every position from a to b. 'paper' = number from the paper's text, figure or the authors'
RNA/cDNA clone mapped onto the genome (codes T, F, M); 'derived' = an abstract statement converted with the
record's own GenBank features, or the GenBank value only (codes A, G): those equal NCBI by construction.
5' end = first nucleotide of the mature rRNA, 3' end = last templated nucleotide (transcript direction).
usage: anchor_score.py ROOT TOOL [TOOL ...]   (ROOT/runs/<tool>/<acc>/final.gbf, ROOT/gb/<acc>.gb = NCBI)
"""
import os
import sys

sys.path.insert(0, "/media/deyuan/18T_B/yzz/MitoM/benchmark")
from mito_bench_eval import norm_features   # noqa: E402
from Bio import SeqIO                        # noqa: E402

ANCHORS = [
    ("NC_001453.1", "rrnS", "5'", [74], "derived"),
    ("NC_001453.1", "rrnS", "3'", [(945, 951)], "derived"),
    ("NC_001453.1", "rrnL", "5'", [4260], "derived"),
    ("NC_001453.1", "rrnL", "3'", [(5787, 5789)], "derived"),
    ("NC_001572.1", "rrnS", "3'", [956], "derived"),
    ("NC_002177.1", "rrnL", "3'", [13099], "paper"),                 # 118 of 151 ESTs; minor sites 13089-13100
    ("NC_002081.1", "rrnS", "5'", [70, 69], "paper"),
    ("NC_002081.1", "rrnS", "3'", [1018], "paper"),
    ("NC_002081.1", "rrnL", "5'", [1090, 1091], "paper"),
    ("NC_002081.1", "rrnL", "3'", [2759, 2764], "paper"),
    ("JN248424.2", "rrnL", "3'", [(12032, 12037)], "paper"),          # 3' RACE, most clones
    ("NC_006817.1", "rrnL", "5'", [14017], "paper"),
    ("NC_006817.1", "rrnL", "3'", [12685], "paper"),
    ("NC_006817.1", "rrnS", "3'", [14090], "paper"),
    ("NC_001322.1", "rrnS", "5'", [14942], "derived"),
    ("NC_024511.2", "rrnL", "5'", [14058], "paper"),
    ("NC_024511.2", "rrnL", "3'", [12735, 12734], "paper"),
    ("NC_024511.2", "rrnS", "5'", [14914], "paper"),
    ("MK374364.1", "rrnL", "5'", [13742], "paper"),
    ("MK374364.1", "rrnL", "3'", [12466], "paper"),
    ("MK374364.1", "rrnS", "5'", [14624], "paper"),
    ("MK374364.1", "rrnS", "3'", [13811], "paper"),
    ("NC_001328.1", "rrnS", "5'", [898], "derived"),
    ("NC_001328.1", "rrnL", "5'", [10403], "derived"),
    ("NC_001327.1", "rrnS", "5'", [901], "derived"),
    ("NC_001327.1", "rrnL", "5'", [11346], "derived"),
    ("NC_006886.2", "rrnS", "5'", [14070], "derived"),
    ("NC_006886.2", "rrnL", "5'", [15501], "derived"),
    ("NC_005089.1", "rrnS", "5'", [70], "derived"),
    ("NC_005089.1", "rrnS", "3'", [1024, 1023], "derived"),
    ("NC_005089.1", "rrnL", "5'", [1094], "derived"),
    ("NC_005089.1", "rrnL", "3'", [(2675, 2682)], "derived"),
]


def dist(pos, accepted):
    best = None
    for a in accepted:
        lo, hi = a if isinstance(a, tuple) else (a, a)
        d = 0 if lo <= pos <= hi else min(abs(pos - lo), abs(pos - hi))
        best = d if best is None else min(best, d)
    return best


def ends(f):
    """1-based (5' end, 3' end) of a Feat (0-based start, exclusive end, may pass the origin)."""
    n = f.n or 10 ** 9
    first, last = f.start % n + 1, (f.end - 1) % n + 1
    return (first, last) if f.strand > 0 else (last, first)


def main():
    root, tools = sys.argv[1], sys.argv[2:]
    cols = ["NCBI"] + tools
    res = {c: [] for c in cols}
    for acc, gene, end, accepted, prov in ANCHORS:
        line = [f"{acc} {gene} {end} {prov:7s} {accepted}"]
        for c in cols:
            path = os.path.join(root, "gb", acc + ".gb") if c == "NCBI" else os.path.join(root, "runs", c, acc, "final.gbf")
            d = None
            if os.path.exists(path):
                fs = [f for f in norm_features(SeqIO.read(path, "genbank")) if f.kind == "rRNA" and f.name == gene]
                if fs:
                    d = min(dist(ends(f)[0 if end == "5'" else 1], accepted) for f in fs)
            res[c].append((prov, end, d))
            line.append(f"{c}={'missing' if d is None else d}")
        print("  ".join(line))
    print()
    for prov in ("paper", "derived", "all"):
        for end in ("5'", "3'", "both"):
            sel = [(p, e) for p, e, _ in res["NCBI"] if (prov == "all" or p == prov) and (end == "both" or e == end)]
            if not sel:
                continue
            out = [f"{prov:7s} {end:4s} n={len(sel):2d}"]
            for c in cols:
                ds = [d for p, e, d in res[c] if (prov == "all" or p == prov) and (end == "both" or e == end)]
                ok0 = sum(d == 0 for d in ds if d is not None)
                ok1 = sum(d <= 1 for d in ds if d is not None)
                ok3 = sum(d <= 3 for d in ds if d is not None)
                out.append(f"{c}: exact {ok0}, <=1 {ok1}, <=3 {ok3}")
            print(" | ".join(out))


if __name__ == "__main__":
    main()
