"""Reference library of NCBI-annotated mitochondrial rRNAs with flanks.

For every rrnS / rrnL feature of the RefSeq pool: region = rRNA +- FLANK nt (circular), in the rRNA's strand,
with the rRNA occupying region positions FLANK .. FLANK + len - 1.  Also the gap to the nearest annotated feature
on each side (5' / 3', in transcript direction; 0 = contiguous, negative = overlap) to describe the submitter's
convention.

usage: build_lib.py Metazoa.gb OUTDIR [workers]
"""
import sys, os, re, io
from multiprocessing import Pool
sys.path.insert(0, "/media/deyuan/18T_B/yzz/MitoM/benchmark")
from mito_names import gene_name          # noqa: E402
from Bio import SeqIO                     # noqa: E402

FLANK = 150


def span(loc, n):
    parts = list(loc.parts)
    strand = loc.strand or 1
    wrapped = any((int(b.start) < int(a.start)) if strand != -1 else (int(b.start) > int(a.start))
                  for a, b in zip(parts, parts[1:]))
    if not wrapped:
        return int(loc.start), int(loc.end)
    if strand != -1:
        return int(parts[0].start), int(parts[-1].end) + n
    return int(parts[-1].start), int(parts[0].end) + n


def circ(seq, s, e):
    n = len(seq)
    t = seq * 3
    return t[s + n:e + n]


COMP = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")


def rc(s):
    return s.translate(COMP)[::-1]


def gaps(feats, i, n):
    """(gap5, gap3) of feature i to the nearest other feature, in transcript direction."""
    s, e, strand = feats[i][1], feats[i][2], feats[i][3]
    best_left, best_right = None, None
    for j, (_, s2, e2, _) in enumerate(feats):
        if j == i:
            continue
        for k in (-n, 0, n):
            a, b = s2 + k, e2 + k
            if a >= e or b <= s:          # no overlap: distance
                if b <= s:
                    d = s - b
                    best_left = d if best_left is None or d < best_left else best_left
                if a >= e:
                    d = a - e
                    best_right = d if best_right is None or d < best_right else best_right
            else:                          # overlap at an end: negative gap
                if a < s < b <= e:
                    d = s - b
                    best_left = d if best_left is None or d < best_left else best_left
                if s <= a < e < b:
                    d = a - e
                    best_right = d if best_right is None or d < best_right else best_right
    if strand == -1:
        return best_right, best_left
    return best_left, best_right


def records(chunk_text):
    return SeqIO.parse(io.StringIO(chunk_text), "genbank")


def work(chunk_text):
    out = []
    for rec in records(chunk_text):
        try:
            seq = str(rec.seq).upper()
        except Exception:          # CONTIG records without sequence
            continue
        n = len(seq)
        if n < 5000:
            continue
        org = rec.annotations.get("organism", "")
        tax = rec.annotations.get("taxonomy", [])
        feats = []
        for f in rec.features:
            if f.type not in ("CDS", "rRNA", "tRNA"):
                continue
            try:
                s, e = span(f.location, n)
            except Exception:
                continue
            g = gene_name(f.qualifiers, f.type) if f.type in ("CDS", "rRNA") else "tRNA"
            feats.append((g or f.type, s, e, f.location.strand or 1, f.type))
        plain = [(g, s, e, st) for g, s, e, st, _ in feats]
        types = [t for *_, t in feats]
        for i, (g, s, e, strand) in enumerate(plain):
            if g not in ("rrnS", "rrnL") or types[i] != "rRNA":
                continue
            L = e - s
            if L < 200 or L > 4000:
                continue
            region = circ(seq, s - FLANK, e + FLANK)
            if strand == -1:
                region = rc(region)
            g5, g3 = gaps(plain, i, n)
            # neighbour type on each side (transcript direction)
            out.append((f"{rec.id}|{g}|{s}", rec.id, g, org, ";".join(tax), L, g5, g3, region))
    return out


def chunks(path, size=40_000_000):
    buf, cur = [], 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            buf.append(line)
            cur += len(line)
            if line.startswith("//") and cur >= size:
                yield "".join(buf)
                buf, cur = [], 0
    if buf:
        yield "".join(buf)


def main():
    gb, outdir = sys.argv[1], sys.argv[2]
    workers = int(sys.argv[3]) if len(sys.argv) > 3 else 24
    os.makedirs(outdir, exist_ok=True)
    rows = []
    with Pool(workers) as pool:
        for res in pool.imap_unordered(work, chunks(gb, 8_000_000)):
            rows.extend(res)
    rows.sort()
    with open(os.path.join(outdir, "lib.fa"), "w") as fa, open(os.path.join(outdir, "lib.tsv"), "w") as tsv:
        tsv.write("ref_id\taccession\tgene\torganism\ttaxonomy\tlength\tgap5\tgap3\n")
        for rid, acc, g, org, tax, L, g5, g3, region in rows:
            fa.write(f">{rid}\n{region}\n")
            tsv.write(f"{rid}\t{acc}\t{g}\t{org}\t{tax}\t{L}\t{'' if g5 is None else g5}\t{'' if g3 is None else g3}\n")
    print("rRNAs:", len(rows))


if __name__ == "__main__":
    main()
