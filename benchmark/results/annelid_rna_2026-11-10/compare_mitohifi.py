#!/usr/bin/env python3
"""A fresh MitoHiFi assembly against the deposited mitogenome at the sites HiFi reads showed to be wrong
(hifi_check/<acc>.tsv): which version does the new assembly carry - the deposited one or the reads' one? Plus every
difference between the two assemblies (minimap2 -cx asm5 --cs).

Usage: compare_mitohifi.py --deposited dep.gb --new new.fasta --sites hifi_check/<acc>.tsv [--minimap2 PATH]
"""
import argparse
import csv
import re
import subprocess
import tempfile

from Bio import SeqIO

FLANK = 12
RC = str.maketrans("ACGT", "TGCA")


def windows(g, sites):
    """(position, kind, bases, deposited window, reads' window) per site; near-by edits applied together."""
    def act(e):
        return e[0] - 1 if e[1] == "insertion" else e[0]
    edits = {act(e): e for e in sites}
    out = []
    for e in sites:
        ppos, kind, seq = e
        lo, hi = ppos - FLANK, ppos + FLANK + (len(seq) if kind == "deletion" else 0)
        alt, i = "", lo
        while i < hi:
            x = edits.get(i)
            if x and x[1] == "insertion":
                alt += g[i] + x[2]
                i += 1
            elif x and x[1] == "deletion":
                i += len(x[2])
            else:
                alt += g[i]
                i += 1
        out.append((ppos, kind, seq, g[lo:hi], alt))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--deposited", required=True)
    p.add_argument("--new", required=True)
    p.add_argument("--sites", required=True)
    p.add_argument("--minimap2", default="minimap2")
    a = p.parse_args()
    dep = SeqIO.read(a.deposited, "genbank")
    g = str(dep.seq).upper()
    new = [str(r.seq).upper() for r in SeqIO.parse(a.new, "fasta")]
    nn = [s + s[:200] for s in new] + [(s + s[:200])[::-1].translate(RC) for s in new]    # circular, both strands
    sites = [(int(r["position"]), r["kind"], r["bases"]) for r in csv.DictReader(open(a.sites), delimiter="\t")
             if r["verdict"].startswith("assembly")]
    print("position\tkind\tbases\tnew assembly carries")
    for ppos, kind, seq, dw, rw in windows(g, sites):
        d, r = any(dw in s for s in nn), any(rw in s for s in nn)
        print(f"{ppos}\t{kind}\t{seq}\t" + ("the reads' version" if r and not d else "the deposited (wrong) version"
                                            if d and not r else "both" if d and r else "neither"))
    with tempfile.NamedTemporaryFile("w", suffix=".fa") as t:
        t.write(f">dep\n{g}\n")
        t.flush()
        paf = subprocess.run([a.minimap2, "-cx", "asm5", "--cs", t.name, a.new], capture_output=True,
                             text=True).stdout
    for ln in paf.splitlines():
        c = ln.split("\t")
        cs = next((x[5:] for x in c[12:] if x.startswith("cs:Z:")), "")
        ops = re.findall(r"[*+-][a-z]+", cs)
        print(f"alignment {c[0]} {c[2]}-{c[3]} to deposited {c[7]}-{c[8]} ({c[4]}): "
              f"{sum(o[0] == '*' for o in ops)} mismatches, {sum(o[0] == '+' for o in ops)} insertions, "
              f"{sum(o[0] == '-' for o in ops)} deletions relative to the deposited sequence")


if __name__ == "__main__":
    main()
