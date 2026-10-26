"""rRNA end conventions in RefSeq (training data: the pool minus every species of the benchmark sets).

Per rRNA end (transcript direction): phylum, gene, end, type of the nearest annotated feature on that side
(tRNA / CDS / rRNA, any strand, within 400 nt), gap to it (0 = contiguous, < 0 = overlap).
Writes conv_ends.tsv and prints the contiguity rate per (phylum, gene, end, neighbour type).
usage: conventions.py Metazoa.gb SET1.tsv [SET2.tsv ...]
"""
import csv, io, re, sys
from collections import Counter, defaultdict
from multiprocessing import Pool
sys.path.insert(0, "/media/deyuan/18T_B/yzz/MitoM/benchmark")
from mito_names import gene_name   # noqa: E402
from Bio import SeqIO              # noqa: E402

PHYLA = set()
EXCL = set()


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


def binom(org):
    return " ".join(re.sub(r"[^A-Za-z0-9. ]", " ", org).split()[:2]).lower()


def work(args):
    text, phyla, excl = args
    out = []
    for rec in SeqIO.parse(io.StringIO(text), "genbank"):
        if binom(rec.annotations.get("organism", "")) in excl:
            continue
        tax = rec.annotations.get("taxonomy", [])
        phylum = next((t for t in tax if t in phyla), "other")
        try:
            n = len(rec.seq)
        except Exception:
            continue
        feats = []
        for f in rec.features:
            if f.type not in ("CDS", "rRNA", "tRNA"):
                continue
            try:
                s, e = span(f.location, n)
            except Exception:
                continue
            g = gene_name(f.qualifiers, f.type) if f.type != "tRNA" else "tRNA"
            feats.append((f.type, g, s, e, f.location.strand or 1))
        for i, (t, g, s, e, st) in enumerate(feats):
            if t != "rRNA" or g not in ("rrnS", "rrnL"):
                continue
            left = right = None
            for j, (t2, g2, s2, e2, st2) in enumerate(feats):
                if j == i:
                    continue
                for k in (-n, 0, n):
                    a, b = s2 + k, e2 + k
                    if a < s and b <= e and b > s - 400:           # neighbour on the left (may overlap the start)
                        d = s - b
                        if left is None or abs(d) < abs(left[0]):
                            left = (d, t2)
                    if b > e and a >= s and a < e + 400:
                        d = a - e
                        if right is None or abs(d) < abs(right[0]):
                            right = (d, t2)
            five, three = (left, right) if st == 1 else (right, left)
            for end, nb in (("5'", five), ("3'", three)):
                out.append((rec.id, phylum, g, end, nb[1] if nb else "none", "" if nb is None else nb[0]))
    return out


def chunks(path, size=8_000_000):
    buf, cur = [], 0
    for line in open(path, encoding="utf-8", errors="replace"):
        buf.append(line)
        cur += len(line)
        if line.startswith("//") and cur >= size:
            yield "".join(buf)
            buf, cur = [], 0
    if buf:
        yield "".join(buf)


def main():
    gb = sys.argv[1]
    phyla, excl = set(), set()
    for tsv in sys.argv[2:]:
        for r in csv.DictReader(open(tsv), delimiter="\t"):
            phyla.add(r["phylum"])
            excl.add(binom(r["species"]))
    rows = []
    with Pool(48) as pool:
        for res in pool.imap_unordered(work, ((c, phyla, excl) for c in chunks(gb))):
            rows.extend(res)
    with open("conv_ends.tsv", "w") as fh:
        fh.write("accession\tphylum\tgene\tend\tnb_type\tgap\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    tab = defaultdict(Counter)
    for acc, ph, g, end, nbt, gap in rows:
        c = tab[(ph, g, end, nbt)]
        c["n"] += 1
        if gap != "":
            c["c0"] += gap == 0
            c["le3"] += abs(gap) <= 3
            c["gap4_300"] += 3 < gap <= 300
            c["ovl4"] += gap < -3
    print("rRNA ends", len(rows))
    print("phylum\tgene\tend\tnb\tn\tcontig%\twithin3%\tgap4-300%\toverlap>3%")
    for k in sorted(tab, key=lambda k: (k[0], k[1], k[2], -tab[k]["n"])):
        c = tab[k]
        if c["n"] < 8:
            continue
        print("\t".join(k) + f"\t{c['n']}\t{100*c['c0']/c['n']:.0f}\t{100*c['le3']/c['n']:.0f}\t{100*c['gap4_300']/c['n']:.0f}"
              f"\t{100*c['ovl4']/c['n']:.0f}")


if __name__ == "__main__":
    main()
