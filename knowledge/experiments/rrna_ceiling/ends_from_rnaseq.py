"""Mature rRNA 5' / 3' ends from paired-end RNA-seq fragments mapped to a mitogenome.

Fragments = properly paired, end-to-end alignments (no soft clipping) to the genome + its first EXT nt (circular).
For every rRNA of the GenBank record (NCBI annotation, used only to say where to look), counts of fragment left
ends and right ends are taken in +-WIN nt around each annotated end; for a + strand rRNA the 5' end is the
left-end peak and the 3' end the right-end peak (minus strand: mirrored).  Library strandedness is inferred from
the reads, and when stranded only fragments from the rRNA's strand are counted.

usage: ends_from_rnaseq.py RECORD.gb FRAGS.tsv OUT.tsv
FRAGS.tsv = 'left right strand' per fragment (0-based left, exclusive right, strand of the transcript if known
or '.'), produced by frags_from_bam.py.
"""
import csv, sys
from collections import Counter
sys.path.insert(0, "/media/deyuan/18T_B/yzz/MitoM/benchmark")
from mito_bench_eval import norm_features   # noqa: E402
from Bio import SeqIO                        # noqa: E402

WIN = 60


def main():
    gb, frags, out = sys.argv[1:4]
    rec = SeqIO.read(gb, "genbank")
    n = len(rec)
    L = {"+": Counter(), "-": Counter(), ".": Counter()}
    R = {"+": Counter(), "-": Counter(), ".": Counter()}
    nfr = 0
    for line in open(frags):
        a, b, s = line.split()
        a, b = int(a) % n, (int(b) - 1) % n          # last base, folded
        L[s][a] += 1
        R[s][b] += 1
        nfr += 1
    stranded = (sum(L["+"].values()) + sum(L["-"].values())) > 0
    rows = []
    for f in norm_features(rec):
        if f.kind != "rRNA":
            continue
        st = "+" if f.strand > 0 else "-"
        keyL = L[st] if stranded else L["."]
        keyR = R[st] if stranded else R["."]
        s0, e0 = f.start % n, (f.end - 1) % n           # first and last base (0-based)
        # 5' end of a + strand rRNA = left end; of a - strand rRNA = right end
        for end, cnt, centre in (("5'", keyL if st == "+" else keyR, s0 if st == "+" else e0),
                                 ("3'", keyR if st == "+" else keyL, e0 if st == "+" else s0)):
            win = [((centre + d) % n, cnt[(centre + d) % n]) for d in range(-WIN, WIN + 1)]
            ranked = sorted(win, key=lambda x: -x[1])
            peak, top = ranked[0]
            second = ranked[1][1]
            vals = sorted(c for _, c in win)
            bg = vals[len(vals) // 2]
            # offset of the peak from the NCBI end in transcript direction
            off = ((peak - centre + n // 2) % n) - n // 2
            if st == "-":
                off = -off
            rows.append(dict(name=f.name, strand=st, end=end, ncbi=centre + 1, peak=peak + 1, offset=off, top=top,
                             second=second, background=bg, fragments=nfr, stranded=int(stranded)))
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0].keys()) if rows else ["name"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
