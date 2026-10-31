"""Topology test set: one FASTA without topology labels, truth known per sequence.
  circ_orig   complete genome as in NCBI (truth: NCBI topology)
  circ_rot    the same rotated to a random start (truth: NCBI topology)
  frag        a random piece of 3-12 kb (truth: linear)
  gap         the genome with a random 500-2000 nt piece removed, opened at the gap (truth: linear)
Genomes: every record of the given GenBank files that NCBI annotated (>= 10 CDS), at most --max per file, seeded.
    make_set.py OUTDIR --max 40 a.gb b.gb ..."""
import argparse
import random
from pathlib import Path

from Bio import SeqIO

p = argparse.ArgumentParser()
p.add_argument("out")
p.add_argument("gb", nargs="+")
p.add_argument("--max", type=int, default=40)
p.add_argument("--seed", type=int, default=1)
a = p.parse_args()
rnd = random.Random(a.seed)
out = Path(a.out)
out.mkdir(parents=True, exist_ok=True)
fa = open(out / "unlabelled.fasta", "w")
truth = open(out / "truth.tsv", "w")
truth.write("seq_id\tsource\tkind\ttruth\tlength\tstart\tend\torganism\n")
seen = set()
for g in a.gb:
    recs = [r for r in SeqIO.parse(g, "genbank") if sum(f.type == "CDS" for f in r.features) >= 10]
    recs = [r for r in recs if str(r.seq) not in seen]
    rnd.shuffle(recs)
    for r in recs[:a.max]:
        seen.add(str(r.seq))
        s, n, topo = str(r.seq), len(r), r.annotations.get("topology", "circular")
        org = r.annotations.get("organism", "")
        acc = r.id.split(".")[0]
        items = [(f"{acc}_orig", "circ_orig", topo, s, 0, n)]
        k = rnd.randrange(n)
        items.append((f"{acc}_rot", "circ_rot", topo, s[k:] + s[:k], k, k))
        L = rnd.randint(3000, min(12000, n - 2000))
        k = rnd.randrange(n - L)
        items.append((f"{acc}_frag", "frag", "linear", s[k:k + L], k, k + L))
        d = rnd.randint(500, 2000)
        k = rnd.randrange(n - d)
        items.append((f"{acc}_gap", "gap", "linear", s[k + d:] + s[:k], k + d, k))
        for sid, kind, t, seq, st, en in items:
            fa.write(f">{sid}\n{seq}\n")
            truth.write(f"{sid}\t{r.id}\t{kind}\t{t}\t{len(seq)}\t{st}\t{en}\t{org}\n")
fa.close(); truth.close()
