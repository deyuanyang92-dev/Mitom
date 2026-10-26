"""Detection if the 13 CM-verified NCBI rRNA errors are scored against the covariance model's strand / name /
position instead of NCBI's (secondary number; the main tables keep NCBI as truth).
usage: ncbi_fixed_det.py ROOT TOOL"""
import os, sys
from Bio import SeqIO
sys.path.insert(0, "/media/deyuan/18T_B/yzz/MitoM/benchmark")
from mito_bench_eval import norm_features

root, tool = sys.argv[1:3]
# accession, CM-supported name, strand, 1-based start, end
CM = [("NC_028407.1", "rrnS", -1, 12360, 13182), ("NC_028085.1", "rrnL", -1, 537, 1204),
      ("NC_028085.1", "rrnS", -1, 3146, 3780), ("NC_032061.1", "rrnL", -1, 11833, 12976),
      ("NC_071876.1", "rrnL", -1, 6885, 8495), ("NC_006162.1", "rrnS", -1, 7428, 8156),
      ("NC_005840.1", "rrnS", -1, 8180, 8868), ("NC_067052.1", "rrnS", -1, 13985, 14719),
      ("NC_028708.1", "rrnS", -1, 12329, 13506), ("NC_044703.1", "rrnS", 1, 9610, 10339),
      ("NC_027943.1", "rrnS", 1, 12838, 13745), ("NC_027943.1", "rrnL", 1, 14914, 16240),
      ("NC_068907.1", "rrnL", 1, 12273, 13590), ("PP998669.1", "rrnL", -1, 8182, 9304),
      ("PP998669.1", "rrnS", -1, 9386, 9984), ("OQ807204.1", "rrnS", 1, 15689, 16503)]
ok = 0
for acc, name, st, a, b in CM:
    p = os.path.join(root, "runs", tool, acc, "final.gbf")
    if not os.path.exists(p):
        continue
    rec = SeqIO.read(p, "genbank")
    want = {i % len(rec) for i in range(a - 1, b)}
    hit = False
    for f in norm_features(rec):
        if f.kind == "rRNA" and f.name == name and f.strand == st:
            n = f.n or 10 ** 9
            got = {i % n for i in range(f.start, f.end)}
            hit |= len(got & want) >= 0.5 * len(want)
    ok += hit
    print(acc, name, "+" if st > 0 else "-", f"{a}..{b}", "MitoM agrees with the CM" if hit else "MitoM does not")
print(f"{tool}: {ok} of the CM-verified genes in {root} found as the CM says")
