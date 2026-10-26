"""NCBI rRNAs that MitoM puts on the other strand / calls the other gene / places elsewhere: which strand and gene
does the full MITOS2 covariance model (cmsearch default pipeline, not --hmmonly) support?
usage: ncbi_errors_cm.py OUTDIR   (-> result_ncbi_errors_cm.txt)
"""
import subprocess, sys
from Bio import SeqIO
W = "/media/deyuan/18T_B/yzz/mitom_work/{}/gb/"
OUT = sys.argv[1]
CS = "/home/deyuan/miniforge3/envs/mitotools/bin/cmsearch"
CM = "/home/deyuan/mitos_refs/refseq89m/ncRNA/{}.cm"
CASES = [("NC_028407.1", "rrnS", 12360, 13182), ("NC_028085.1", "rrnL", 537, 1204), ("NC_028085.1", "rrnS", 3146, 3780),
         ("NC_032061.1", "rrnL", 11833, 12976), ("NC_071876.1", "rrnL", 6885, 8495), ("NC_006162.1", "rrnS", 7428, 8156),
         ("NC_005840.1", "rrnS", 8180, 8868), ("NC_067052.1", "rrnL", 13985, 14719), ("NC_028708.1", "rrnL", 12329, 13506),
         ("NC_044703.1", "rrnL", 9610, 10339), ("NC_027943.1", "genome", None, None),
         ("NC_068907.1", "rrnS", 12273, 13590, "bench_annelida"), ("PP998669.1", "rrnL", 8182, 9304, "bench_annelida"),
         ("PP998669.1", "rrnS", 9386, 9984, "bench_annelida"), ("OQ807204.1", "rrnS", 15689, 16503, "bench_annelida")]


def hits(model, fa, tag):
    tbl = f"{OUT}/{tag}_{model}.tbl"
    subprocess.run([CS, "--cpu", "2", "--notextw", "-E", "10", "--tblout", tbl, CM.format(model), fa],
                   capture_output=True, check=True)
    out = []
    for line in open(tbl):
        if not line.startswith("#"):
            c = line.split()
            out.append((float(c[14]), c[9], int(c[7]), int(c[8]), int(c[5]), int(c[6])))
    return sorted(out, reverse=True)


for acc, g, a, b, *root in CASES:
    rec = SeqIO.read(W.format(root[0] if root else "bench_large") + acc + ".gb", "genbank")
    seq = str(rec.seq if a is None else (rec.seq + rec.seq)[a - 1:b])          # b > length: across the origin
    tag = f"{acc}_{g}"
    fa = f"{OUT}/{tag}.fa"
    open(fa, "w").write(f">{acc}\n{seq}\n")
    res = []
    for model in ("rrnS", "rrnL"):
        hs = hits(model, fa, tag)
        best = {"+": 0.0, "-": 0.0}
        for bits, st, *_ in hs:
            best[st] = max(best[st], bits)
        top = hs[0] if hs else None
        res.append(f"{model}: + {best['+']:.0f} / - {best['-']:.0f} bits"
                   + (f" (top {top[2]}..{top[3]} {top[1]} model {top[4]}-{top[5]})" if top and a is None else ""))
    ann = [f for f in rec.features if f.type == "rRNA" and a is not None and any(int(p.start) + 1 == a for p in f.location.parts)]
    st = ("+" if ann[0].location.strand == 1 else "-") if ann else "?"
    print(acc, g, f"NCBI {a}..{b} strand {st}" if a else "whole genome", "|", " | ".join(res), flush=True)
