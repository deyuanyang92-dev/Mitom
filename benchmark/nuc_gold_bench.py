#!/usr/bin/env python3
"""nuc_gold_bench.py - rDNA boundaries against curated RefSeq rRNA of the same species (independent of submitters).

Gold: the curated NR_ records of `mitom nuc fetch` (curated.gb; human, rat, macaque, rabbit, cow, hamster 5.8S,
Drosophila 18S, C. elegans 18S / 5.8S) mapped onto GenBank records of the same species (blastn >= 98 % identity, the
gold gene end aligned); every gold boundary inside a record is scored for
  NCBI submitters - the record's own rRNA features (ends marked partial count as not annotated),
  mitom nuc       - with the record's genus left out (default) or the full library (--no_leave_out),
  barrnap 0.9     - when --barrnap gives its bin directory.
Test records (fetched 2026-10-23): for each gold species
  "<species>"[Organism] AND (18S[Title] OR 28S[Title] OR 5.8S[Title] OR "ribosomal DNA"[Title] OR rDNA[Title] OR
  "ribosomal RNA gene"[Title]) AND 140:60000[SLEN] NOT srcdb_refseq[PROP] NOT biomol_mrna[PROP]   (first 400 each)

Usage:
  python3 benchmark/nuc_gold_bench.py --curated nuc_refs/curated.gb --gb gold_species.gb -o gold_bench \
      [--lib models/nuc] [--no_leave_out] [--barrnap ENV/bin]
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

TOL = (0, 1, 5)


def gold_boundaries(curated, recs, out):
    """(record, gene, side, strand) -> gold positions (copies of one species may differ by a few nt; any counts)"""
    org = {r.id: r.annotations["organism"] for r in recs}
    with open(out / "gold.fa", "w") as o:
        for r in SeqIO.parse(str(curated), "genbank"):
            g = next((nuc.rrna_name(f) for f in r.features if nuc.rrna_name(f)), None)
            if g:
                o.write(f">{r.id}|{g}|{r.annotations['organism'].replace(' ', '_')}\n{r.seq}\n")
    nuc.run([nuc.tool("makeblastdb"), "-in", str(out / "q.fa"), "-dbtype", "nucl"])
    truth = defaultdict(set)
    for line in nuc.run([nuc.tool("blastn"), "-query", str(out / "gold.fa"), "-db", str(out / "q.fa"), "-evalue",
                         "1e-20", "-max_target_seqs", "5000", "-num_threads", "8",
                         "-outfmt", "6 qseqid sseqid pident qstart qend sstart send qlen"]).splitlines():
        q, s, pid, qs, qe, ss, se, ql = line.split("\t")
        qs, qe, ss, se, ql = map(int, (qs, qe, ss, se, ql))
        _, g, o = q.split("|")
        if o.replace("_", " ") != org[s] or float(pid) < 98:
            continue
        st = 1 if se >= ss else -1
        if qs == 1:
            truth[(s, g, "start" if st == 1 else "end", st)].add(ss)
        if qe == ql:
            truth[(s, g, "end" if st == 1 else "start", st)].add(se)
    return truth


def submitters(recs):
    p = defaultdict(list)
    for r in recs:
        for f in r.features:
            g = nuc.rrna_name(f)
            if g:
                p5, p3 = nuc._partial(f)
                st = f.location.strand or 1
                lo_ok, hi_ok = (not p5, not p3) if st == 1 else (not p3, not p5)
                p[(r.id, g, st)].append((int(f.location.start) + 1 if lo_ok else None,
                                         int(f.location.end) if hi_ok else None))
    return p


def mitom(path):
    p = defaultdict(list)
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["type"] == "rRNA":
            st = 1 if r["strand"] == "+" else -1
            lo_p, hi_p = (r["partial_5"], r["partial_3"]) if st == 1 else (r["partial_3"], r["partial_5"])
            p[(r["contig"], r["gene"], st)].append((None if lo_p == "1" else int(r["start"]),
                                                    None if hi_p == "1" else int(r["end"])))
    return p


def barrnap(path):
    p = defaultdict(list)
    for line in open(path):
        if not line.startswith("#"):
            x = line.split("\t")
            g = re.search(r"Name=([^;]+)", x[8]).group(1).replace("_rRNA", "").replace("5_8S", "5.8S")
            p[(x[0], g, 1 if x[6] == "+" else -1)].append((int(x[3]), int(x[4])))
    return p


def score(truth, pred):
    """per boundary: |distance| to the nearest predicted end of that gene and strand, None = not annotated"""
    d = {}
    for (s, g, side, st), pos in truth.items():
        ends = [(a if side == "start" else b) if st == 1 else (b if side == "start" else a)
                for a, b in pred.get((s, g, st), [])]
        ends = [e for e in ends if e is not None]
        d[(s, g, side, st)] = min(abs(e - t) for e in ends for t in pos) if ends else None
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--curated", type=Path, required=True)
    ap.add_argument("--gb", nargs="+", type=Path, required=True)
    ap.add_argument("-o", "--out", type=Path, required=True)
    ap.add_argument("--lib", type=Path, default=nuc.LIB)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--no_leave_out", action="store_true", help="full library (only the record itself is skipped)")
    ap.add_argument("--barrnap", help="bin directory of barrnap 0.9")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    species = {r.annotations["organism"] for r in SeqIO.parse(str(a.curated), "genbank")}
    recs = [r for gb in a.gb for r in SeqIO.parse(str(gb), "genbank")
            if r.annotations.get("organism") in species and not r.id.startswith("NR_")]
    SeqIO.write(recs, str(a.out / "q.fa"), "fasta")
    truth = gold_boundaries(a.curated, recs, a.out)
    cmd = [sys.executable, str(HERE.parent / "mito_nuc.py"), "annotate", "-i", str(a.out / "q.fa"), "-o",
           str(a.out / "mitom"), "--lib", str(a.lib), "--genes", "rdna", "--threads", str(a.threads),
           "--consensus_identity", "0"]
    if not a.no_leave_out:
        (a.out / "leave_out.tsv").write_text("".join(f"{r.id}\t{r.annotations['organism'].split()[0]}\n" for r in recs))
        cmd += ["--leave_out", str(a.out / "leave_out.tsv")]
    subprocess.run(cmd, check=True, capture_output=True)
    tools = {"NCBI submitters": submitters(recs), "mitom nuc": mitom(a.out / "mitom" / "nuc_features.tsv")}
    if a.barrnap:
        env = dict(os.environ, PATH=a.barrnap + os.pathsep + os.environ["PATH"])
        with open(a.out / "barrnap.gff", "w") as o:
            subprocess.run(["barrnap", "--kingdom", "euk", "--threads", str(a.threads), str(a.out / "q.fa")],
                           stdout=o, stderr=subprocess.DEVNULL, env=env, check=True)
        tools["barrnap"] = barrnap(a.out / "barrnap.gff")
    d = {name: score(truth, p) for name, p in tools.items()}
    annotated = [k for k, v in d["NCBI submitters"].items() if v is not None]
    rep = [f"gold: {len(truth)} boundaries in {len({k[0] for k in truth})} records "
           f"({dict(Counter(recs_org for recs_org in (next(r.annotations['organism'] for r in recs if r.id == k[0]) for k in truth)))}); "
           f"mitom {'full library' if a.no_leave_out else 'genus left out'}",
           "set                    n     tool              exact  <=1    <=5    not annotated (%)"]
    for label, keys in (("where NCBI annotates", annotated), ("all gold boundaries", list(truth))):
        for name, dd in d.items():
            v = [dd[k] for k in keys]
            rep.append(f"{label:21} {len(v):4}  {name:16} " + "  ".join(
                f"{100 * sum(x is not None and x <= t for x in v) / len(v):5.1f}" for t in TOL) +
                f"  {100 * sum(x is None for x in v) / len(v):5.1f}")
    with open(a.out / "per_boundary.tsv", "w") as o:
        o.write("record\tgene\tside\t" + "\t".join(d) + "\n")
        for k in truth:
            o.write(f"{k[0]}\t{k[1]}\t{k[2]}\t" + "\t".join("" if dd[k] is None else str(dd[k]) for dd in d.values()) + "\n")
    (a.out / "report.txt").write_text("\n".join(rep) + "\n")
    print("\n".join(rep))


if __name__ == "__main__":
    main()
