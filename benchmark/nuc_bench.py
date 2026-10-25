#!/usr/bin/env python3
"""nuc_bench.py - accuracy of mitom nuc (rDNA boundaries, histone CDS models) against complete NCBI annotations.

Truth: the GenBank records of `mitom nuc fetch` (rdna.gb, histone.gb).
  rDNA: records with exactly one complete 18S, ITS1, 5.8S, ITS2 and 28S feature on one strand, the first record
  of each species; every boundary of 18S / 5.8S / 28S is compared (exact, <= 1, 5, 20 nt).
  histones: complete histone H1/H2A/H2B/H3/H4 CDS (the library filter of mito_nuc.build); a model is exact when
  its exons and gene name equal the truth.
References of the query's own genus are left out (mito_nuc annotate --leave_out). barrnap (0.9, --kingdom euk) is
run on the same rDNA queries when --barrnap gives its bin directory.

Usage:
  python3 benchmark/nuc_bench.py --gb nuc_refs/rdna.gb nuc_refs/histone.gb -o nuc_bench [--barrnap ENV/bin]
  held-out records with the shipped library: ... --gb new_rdna.gb new_histone.gb --lib models/nuc
Writes in -o: library/ (built from the same records), rdna_queries.fa, histone_queries.fa, mitom_rdna/,
mitom_histone/, barrnap.gff, report.txt.
"""
import argparse
import csv
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_nuc as nuc  # noqa: E402

TOL = (0, 1, 5, 20)


def rdna_truth(gbs):
    """acc -> {gene: (start, end, strand)}, organism genus; one record per species"""
    truth, genus, species = {}, {}, set()
    for gb in gbs:
        for rec in SeqIO.parse(str(gb), "genbank"):
            d = defaultdict(list)
            for f in rec.features:
                g = nuc.rrna_name(f)
                if not g and f.type == "misc_RNA":
                    lab = nuc._label(f)
                    g = "ITS1" if re.search(r"internal transcribed spacer 1\b|ITS1\b", lab) else \
                        "ITS2" if re.search(r"internal transcribed spacer 2\b|ITS2\b", lab) else None
                if g and not any(nuc._partial(f)):
                    d[g].append(f)
            if not all(len(d.get(g, [])) == 1 for g in nuc.UNIT) or \
                    len({d[g][0].location.strand for g in nuc.UNIT}) != 1:
                continue
            org = rec.annotations.get("organism", "")
            if org in species:
                continue
            species.add(org)
            truth[rec.id] = (rec, {g: (int(d[g][0].location.start) + 1, int(d[g][0].location.end),
                                       d[g][0].location.strand) for g in nuc.RRNA})
            genus[rec.id] = org.split()[0]
    return truth, genus


def histone_truth(gbs):
    truth, genus = {}, {}
    for gb in gbs:
        for rec in SeqIO.parse(str(gb), "genbank"):
            cds = []
            for f in rec.features:
                h = nuc.histone_name(f)
                if not h or "translation" not in f.qualifiers or any(nuc._partial(f)) or \
                        int(f.qualifiers.get("codon_start", ["1"])[0]) != 1 or \
                        len(f.qualifiers["translation"][0]) > 400:
                    continue
                cds.append((h, f.location.strand, sorted((int(p.start) + 1, int(p.end)) for p in f.location.parts)))
            if cds:
                truth[rec.id] = (rec, cds)
                genus[rec.id] = rec.annotations.get("organism", "").split()[0]
    return truth, genus


def write_queries(truth, genus, fa, lo):
    SeqIO.write([r for r, _ in truth.values()], str(fa), "fasta")
    lo.write_text("".join(f"{a}\t{g}\n" for a, g in genus.items()))


def rdna_scores(truth, pred):
    """pred: acc -> {gene: (start, end, strand)}; returns per-gene counts of boundaries within each tolerance"""
    hit = defaultdict(lambda: [0] * len(TOL))
    n = Counter()
    for acc, (_, t) in truth.items():
        for g in nuc.RRNA:
            n[g] += 2
            p = pred.get(acc, {}).get(g)
            if not p or p[2] != t[g][2]:
                continue
            for i in (0, 1):
                for k, tol in enumerate(TOL):
                    hit[g][k] += abs(p[i] - t[g][i]) <= tol
    return hit, n


def mitom_rdna(path):
    """longest feature per (contig, gene)"""
    best = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        k = (r["contig"], r["gene"])
        if k not in best or int(r["length"]) > int(best[k]["length"]):
            best[k] = r
    pred = defaultdict(dict)
    for (c, g), r in best.items():
        pred[c][g] = (int(r["start"]), int(r["end"]), 1 if r["strand"] == "+" else -1)
    return pred


def barrnap_rdna(path):
    """best (lowest e-value) barrnap hit per (sequence, gene)"""
    pred, ev = defaultdict(dict), {}
    for line in open(path):
        if line.startswith("#"):
            continue
        x = line.rstrip("\n").split("\t")
        g = re.search(r"Name=([^;]+)", x[8]).group(1).replace("_rRNA", "").replace("5_8S", "5.8S")
        e = float(x[5]) if x[5] != "." else 1.0
        if (x[0], g) not in ev or e < ev[(x[0], g)]:
            ev[(x[0], g)] = e
            pred[x[0]][g] = (int(x[3]), int(x[4]), 1 if x[6] == "+" else -1)
    return pred


def histone_scores(truth, path):
    pred = defaultdict(list)
    for r in csv.DictReader(open(path), delimiter="\t"):
        pred[r["contig"]].append((r["gene"], 1 if r["strand"] == "+" else -1,
                                  [tuple(map(int, e.split("-"))) for e in r["exons"].split(";")]))
    c, used = Counter(), set()
    for acc, (_, cds) in truth.items():
        for g, st, ex in cds:
            c["truth"] += 1
            ps = [p for p in pred[acc] if p[1] == st and p[2][0][0] <= ex[-1][1] and p[2][-1][1] >= ex[0][0]]
            if not ps:
                continue
            p = max(ps, key=lambda p: min(ex[-1][1], p[2][-1][1]) - max(ex[0][0], p[2][0][0]))
            used.add((acc, p[2][0][0]))
            c["found"] += 1
            c["name"] += p[0] == g
            c["exact"] += p[0] == g and p[2] == ex
    c["extra"] = sum(1 for a, ps in pred.items() for p in ps if (a, p[2][0][0]) not in used)
    return c


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gb", nargs="+", type=Path, required=True)
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--barrnap", help="bin directory of barrnap 0.9 (and its nhmmer)")
    ap.add_argument("--lib", type=Path, help="existing library (held-out test); default: built from --gb")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    lib = a.lib or a.out / "library"
    if not a.lib:
        nuc.build(a.gb, lib)
    me = [sys.executable, str(HERE.parent / "mito_nuc.py"), "annotate", "--lib", str(lib),
          "--threads", str(a.threads), "--consensus_identity", "0"]
    rep = []

    rt, rg = rdna_truth(a.gb)
    write_queries(rt, rg, a.out / "rdna_queries.fa", a.out / "rdna_leave_out.tsv")
    subprocess.run(me + ["-i", str(a.out / "rdna_queries.fa"), "-o", str(a.out / "mitom_rdna"), "--genes", "rdna",
                         "--leave_out", str(a.out / "rdna_leave_out.tsv")], check=True, capture_output=True)
    tools = {"mitom nuc": mitom_rdna(a.out / "mitom_rdna" / "nuc_features.tsv")}
    if a.barrnap:
        env = dict(os.environ, PATH=a.barrnap + os.pathsep + os.environ["PATH"])
        with open(a.out / "barrnap.gff", "w") as o:
            subprocess.run(["barrnap", "--kingdom", "euk", "--threads", str(a.threads), str(a.out / "rdna_queries.fa")],
                           stdout=o, stderr=subprocess.DEVNULL, env=env, check=True)
        tools["barrnap"] = barrnap_rdna(a.out / "barrnap.gff")
    rep.append(f"rDNA: {len(rt)} records (one per species), references of the same genus left out")
    rep.append("tool        gene   " + "  ".join(f"<={t:<3}" for t in TOL) + "   (% of boundaries)")
    for name, pred in tools.items():
        hit, n = rdna_scores(rt, pred)
        for g in (*nuc.RRNA, "all"):
            h = [sum(hit[x][k] for x in nuc.RRNA) for k in range(len(TOL))] if g == "all" else hit[g]
            m = sum(n.values()) if g == "all" else n[g]
            rep.append(f"{name:11} {g:5}  " + "  ".join(f"{100 * x / m:5.1f}" for x in h))

    ht, hg = histone_truth(a.gb)
    write_queries(ht, hg, a.out / "histone_queries.fa", a.out / "histone_leave_out.tsv")
    subprocess.run(me + ["-i", str(a.out / "histone_queries.fa"), "-o", str(a.out / "mitom_histone"), "--genes",
                         "histone", "--leave_out", str(a.out / "histone_leave_out.tsv")], check=True,
                   capture_output=True)
    c = histone_scores(ht, a.out / "mitom_histone" / "nuc_features.tsv")
    n = c["truth"]
    rep.append(f"histones: {n} complete CDS in {len(ht)} records, references of the same genus left out")
    rep.append(f"found {100 * c['found'] / n:.1f} %  name {100 * c['name'] / n:.1f} %  exact model and name "
               f"{100 * c['exact'] / n:.1f} %  predictions not in the truth {c['extra']}")
    (a.out / "report.txt").write_text("\n".join(rep) + "\n")
    print("\n".join(rep))


if __name__ == "__main__":
    main()
