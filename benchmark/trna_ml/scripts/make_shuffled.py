#!/usr/bin/env python3
"""Negative-control genomes: dinucleotide-shuffled copies (esl-shuffle -d) of records from one split.

Usage: make_shuffled.py <dataset> <out_dataset> --split test --n 150 [--seed 7]
Writes <out>/records.tsv (acc prefixed SHUF_, split test), <out>/fasta/, an empty trnas_labeled.tsv and
<out>/accs.txt.
"""
import argparse
import random
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import INF, read_tsv, write_tsv  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset")
    ap.add_argument("out")
    ap.add_argument("--split", default="test")
    ap.add_argument("--n", type=int, default=150)
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    recs = [r for r in read_tsv(f"{a.dataset}/records.tsv") if r["split"] == a.split]
    random.Random(a.seed).shuffle(recs)
    recs = recs[:a.n]
    (Path(a.out) / "fasta").mkdir(parents=True, exist_ok=True)
    out = []
    for r in recs:
        acc = "SHUF_" + r["acc"]
        s = subprocess.run([f"{INF}/esl-shuffle", "-d", "--seed", str(a.seed), f"{a.dataset}/fasta/{r['acc']}.fa"],
                           capture_output=True, text=True, check=True).stdout
        seq = "".join(s.split("\n")[1:])
        Path(f"{a.out}/fasta/{acc}.fa").write_text(f">{acc}\n{seq}\n")
        out.append(dict(r, acc=acc, n_trna=0, split="test"))
    write_tsv(f"{a.out}/records.tsv", out, list(recs[0].keys()))
    write_tsv(f"{a.out}/trnas_labeled.tsv", [], list(read_tsv(f"{a.dataset}/trnas_labeled.tsv")[0].keys()))
    Path(f"{a.out}/accs.txt").write_text("".join(r["acc"] + "\n" for r in out))


if __name__ == "__main__":
    main()
