#!/usr/bin/env python3
"""Sites where the RNA-seq differs from the mitogenome (mito_rnaseq_annot.py step 4), checked in HiFi DNA reads of the
same specimen: an assembly base error (DNA = RNA) or an RNA-level change (DNA = genome).

For every flagged insertion / deletion the RNA consensus (bwa mem RNA BAM, samtools mpileup: the most frequent indel
string) is applied; a window of FLANK nt either side is built twice - the genome's and the RNA's (all RNA-consensus
edits inside the window applied, so close sites are judged together) - and exact matches are counted in the HiFi reads
and the RNA reads (both strands). Verdict: assembly error when >= 80 % of the HiFi reads that match either version
carry the RNA's; RNA-level when >= 80 % carry the genome's; else mixed (heteroplasmy?).

Usage: hifi_site_check.py -g annotated.gb --evidence rnaseq_evidence.tsv --rna-bam rna.bam --hifi hifi_mito.fa
       [--min-share 0.2] -o sites.tsv
"""
import argparse
import collections
import csv
import re
import subprocess

from Bio import SeqIO

FLANK = 12
RC = str.maketrans("ACGT", "TGCA")
SITE = re.compile(r"RNA (insertion|deletion) at (\d+) in (\d+) of (\d+) reads")


def consensus(bam, contig, kind, pos, samtools, ref=None):
    """Most frequent indel string at a site (insertion after 1-based pos; deletion of 1-based pos ...)."""
    ppos = pos if kind == "insertion" else pos - 1
    out = subprocess.run([samtools, "mpileup", "-B", "-Q", "20", "-d", "0", "-r", f"{contig}:{ppos}-{ppos}"]
                         + (["-f", ref] if ref else []) + [bam],
                         capture_output=True, text=True).stdout.split("\t")
    if len(out) < 5:
        return None
    b = re.sub(r"\^.", "", out[4])
    sign = "+" if kind == "insertion" else "-"
    c = collections.Counter(m.group(2)[: int(m.group(1))].upper()
                            for m in re.finditer(re.escape(sign) + r"(\d+)([ACGTNacgtn]+)", b))
    return (ppos, kind, c.most_common(1)[0][0]) if c else None


def load(path):
    out = []
    for r in SeqIO.parse(path, "fasta"):
        s = str(r.seq).upper()
        out += [s, s[::-1].translate(RC)]
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-g", "--genbank", required=True)
    p.add_argument("--evidence", required=True)
    p.add_argument("--rna-bam", required=True)
    p.add_argument("--rna-fa", required=True, help="RNA reads mapped to the mitogenome (FASTA)")
    p.add_argument("--hifi", required=True, help="HiFi reads mapped to the mitogenome (FASTA)")
    p.add_argument("--min-share", type=float, default=0.2)
    p.add_argument("--samtools", default="samtools")
    p.add_argument("--ref", default=None, help="FASTA of the mitogenome (shows deleted bases instead of N)")
    p.add_argument("-o", "--out", required=True)
    a = p.parse_args()
    rec = SeqIO.read(a.genbank, "genbank")
    g = str(rec.seq).upper()
    sites = []
    for e in csv.DictReader(open(a.evidence), delimiter="\t"):
        if e["step"] != "4":
            continue
        m = SITE.search(e["detail"])
        if not m:
            continue
        kind, pos, c, d = m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4))
        if c >= a.min_share * d:
            cons = consensus(a.rna_bam, rec.id, kind, pos, a.samtools, a.ref)
            if cons:
                sites.append((*cons, e["feature"], c / d))
    def act(e):                          # 0-based index where the edit acts: after it (insertion) / from it (deletion)
        return e[0] - 1 if e[1] == "insertion" else e[0]
    edits = {act(s): s for s in sites if s[4] >= 0.8}           # near-all sites: applied as context to neighbours
    hifi, rna = load(a.hifi), load(a.rna_fa)
    rows = []
    for ppos, kind, seq, feat, share in sites:
        lo, hi = ppos - FLANK, ppos + FLANK + (len(seq) if kind == "deletion" else 0)
        gen = g[lo:hi]
        alt = ""
        i = lo
        here = act((ppos, kind, seq))
        while i < hi:                       # rebuild the window with this edit and near-all edits inside it
            e = (ppos, kind, seq) if i == here else edits.get(i)
            if e and e[1] == "insertion":
                alt += g[i] + e[2]
                i += 1
            elif e and e[1] == "deletion":
                i += len(e[2])
            else:
                alt += g[i]
                i += 1
        hg, ha = sum(gen in s for s in hifi), sum(alt in s for s in hifi)
        rg, ra = sum(gen in s for s in rna), sum(alt in s for s in rna)
        tot = hg + ha
        verdict = ("too few HiFi reads" if tot < 5 else "assembly error (DNA = RNA)" if ha >= 0.8 * tot else
                   "RNA-level change (DNA = genome)" if hg >= 0.8 * tot else "mixed in DNA (heteroplasmy?)")
        rows.append({"position": ppos, "kind": kind, "bases": seq, "feature": feat, "rna_share": f"{share:.2f}",
                     "hifi_genome": hg, "hifi_rna_version": ha, "rna_genome": rg, "rna_rna_version": ra,
                     "verdict": verdict})
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]) if rows else ["position"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print("\t".join(str(v) for v in r.values()))


if __name__ == "__main__":
    main()
