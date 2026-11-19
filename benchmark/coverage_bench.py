#!/usr/bin/env python3
"""coverage_bench.py - which way of computing mitogenome read depth is fastest and most accurate?

Reads are simulated with wgsim from a copy of the reference ROTATED by half its length, so reads that span the
reference's origin are continuous in the simulation and the true depth of every base is known from the read names
(wgsim names: <chr>_<start>_<end>_...; read 1 covers start.., read 2 ends at end). Compared (all end in a per-base
depth table):
  samtools-bwa    bwa mem | samtools view -F12 | sort | index | depth      (using_bwa_draw_coveragev1.2.py; MitoZ's
                  depth file is made the same way)
  stream-<m>      mito_coverage.py: <m> SAM on a pipe, reference span added to a difference array, circular
                  reference extended and folded; m = bwa, bwa-mem2, minimap2, bowtie2
  samtools-bwa-circ  as samtools-bwa but on the extended reference, folded afterwards (what the first one lacks)
bowtie2 with --sensitive --dovetail --no-discordant is GetOrganelle's mapper setting.
Metrics: wall seconds (index + map + count), Pearson r with the true depth, sum|depth - truth| / sum(truth),
fraction of bases with depth >= 0.9 x truth, depth at the origin +-50 bp relative to the truth.

Usage: python3 benchmark/coverage_bench.py REF.gb OUT [--pairs 200000 --div 0.0 --threads 8 --reps 2 --tag name]
Needs wgsim, samtools, bwa, bwa-mem2, minimap2, bowtie2 on PATH (conda env mapbench)."""
import argparse
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import mito_coverage as mc  # noqa: E402

READ = 150


def simulate(rec, work: Path, pairs: int, div: float, seed: int):
    s = str(rec.seq)
    n = len(s)
    h = n // 2
    rot = work / "rot.fa"
    rot.write_text(">rot\n" + s[h:] + s[:h] + "\n")
    r1, r2 = work / "r1.fq", work / "r2.fq"
    subprocess.run(["wgsim", "-N", str(pairs), "-1", str(READ), "-2", str(READ), "-r", str(div), "-R", "0", "-e",
                    "0.002", "-S", str(seed), str(rot), str(r1), str(r2)], check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)
    diff = np.zeros(n + 1, dtype=np.int64)          # truth in rotated coordinates, then rotate back
    with open(r1) as fh:
        for i, line in enumerate(fh):
            if i % 4:
                continue
            m = re.match(r"@rot_(\d+)_(\d+)_", line)
            a, b = int(m.group(1)) - 1, int(m.group(2))
            diff[a] += 1
            diff[a + READ] -= 1
            diff[b - READ] += 1
            diff[b] -= 1
    t = np.cumsum(diff[:n])
    return np.roll(t, h), r1, r2          # rot index i is original (i + h) % n


def samtools_depth(rec, r1, r2, threads, circular_fold: bool) -> np.ndarray:
    n = len(rec.seq)
    ext = min(mc.ORIGIN_EXT, n - 1) if circular_fold else 0
    with tempfile.TemporaryDirectory() as tmp:
        fa = Path(tmp) / "ref.fa"
        s = str(rec.seq)
        fa.write_text(f">{rec.id}\n{s + s[:ext]}\n")
        null = subprocess.DEVNULL
        subprocess.run(["bwa", "index", str(fa)], check=True, stdout=null, stderr=null)
        bam = Path(tmp) / "a.bam"
        mem = subprocess.Popen(["bwa", "mem", "-t", str(threads), str(fa), str(r1), str(r2)], stdout=subprocess.PIPE,
                               stderr=null)
        view = subprocess.Popen(["samtools", "view", "-bF", "12", "-"], stdin=mem.stdout, stdout=subprocess.PIPE,
                                stderr=null)
        mem.stdout.close()
        srt = subprocess.Popen(["samtools", "sort", "-@", str(threads), "-o", str(bam)], stdin=view.stdout,
                               stderr=null)
        view.stdout.close()
        srt.wait()
        subprocess.run(["samtools", "index", str(bam)], check=True)
        out = subprocess.run(["samtools", "depth", str(bam)], check=True, capture_output=True, text=True).stdout
    d = np.zeros(n + ext, dtype=np.int64)
    for line in out.splitlines():
        f = line.split("\t")
        d[int(f[1]) - 1] = int(f[2])
    res = d[:n].copy()
    if ext:
        res[:ext] += d[n:]
    return res


def score(d, truth, n):
    m = truth > 0
    r = float(np.corrcoef(d, truth)[0, 1])
    rel = float(np.abs(d - truth).sum() / truth.sum())
    ok = float((d[m] >= 0.9 * truth[m]).mean())
    idx = np.r_[np.arange(n - 50, n), np.arange(0, 50)]
    org = float(d[idx].mean() / max(truth[idx].mean(), 1e-9))
    return r, rel, ok, org


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ref")
    ap.add_argument("out")
    ap.add_argument("--pairs", type=int, default=200000)
    ap.add_argument("--div", type=float, default=0.0, help="wgsim SNP rate between the simulated individual and ref")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--reps", type=int, default=2)
    ap.add_argument("--seed", type=int, default=11)
    ap.add_argument("--tag", default="run")
    a = ap.parse_args()
    rec = next(SeqIO.parse(a.ref, "genbank"))
    rec.annotations["topology"] = "circular"
    n = len(rec.seq)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        truth, r1, r2 = simulate(rec, Path(tmp), a.pairs, a.div, a.seed)
        methods = {"samtools-bwa": lambda: samtools_depth(rec, r1, r2, a.threads, False),
                   "samtools-bwa-circ": lambda: samtools_depth(rec, r1, r2, a.threads, True)}
        for m in mc.MAPPERS:
            methods[f"stream-{m}"] = (lambda m=m: mc.depth_from_reads([rec], [str(r1), str(r2)], m, a.threads)[rec.id])
        rows = []
        for name, fn in methods.items():
            times, d = [], None
            for _ in range(a.reps):
                t0 = time.time()
                d = fn()
                times.append(time.time() - t0)
            r, rel, ok, org = score(d, truth, n)
            rows.append((a.tag, name, a.pairs, a.div, a.threads, min(times), r, rel, ok, org))
            print(f"{a.tag}\t{name}\t{min(times):.1f}s\tr={r:.4f}\trel_err={rel:.4f}\t>=0.9truth={ok:.3f}\torigin={org:.3f}",
                  flush=True)
    path = out / "results.tsv"
    new = not path.exists()
    with open(path, "a") as h:
        if new:
            h.write("scenario\tmethod\tpairs\tdivergence\tthreads\twall_s\tpearson_r\trel_abs_err\tfrac_ge_0.9truth\t"
                    "origin_depth_ratio\n")
        for r in rows:
            h.write("\t".join(f"{x:.4f}" if isinstance(x, float) else str(x) for x in r) + "\n")


if __name__ == "__main__":
    main()
