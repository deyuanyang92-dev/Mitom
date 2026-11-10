#!/usr/bin/env python3
"""Map paired reads to one mitogenome with each mapper; SAM + time / memory per mapper.

Usage: run_mappers.py -r ref.fa -1 R1.fq.gz -2 R2.fq.gz -o outdir [-t 8] [--mappers a,b] [--bin envbin]
Writes outdir/<mapper>.sam and outdir/<mapper>.time (GNU time -v), outdir/<mapper>.cmd.
Settings: each mapper's defaults for RNA-seq / short reads, plus the variant that matters for mitochondrial
transcripts (soft clipping for poly(A) tails; no splicing - animal mitochondrial genes have no spliceosomal introns).
"""
from __future__ import annotations

import argparse
import os
import shlex
import subprocess
from pathlib import Path


def mappers(B: str, ref: str, idx: Path, r1: str, r2: str, t: int, sam: str):
    s = shlex.quote
    hx, bx, sx, ux, gx = idx / "hisat2", idx / "bt2", idx / "star", idx / "subread", idx / "segemehl.idx"
    return {
        # name: (index command or None, map command writing SAM to `sam`)
        "minimap2_sr": (None, f"{B}/minimap2 -ax sr -Y -t {t} {s(ref)} {s(r1)} {s(r2)} > {sam}"),
        "minimap2_splice_sr": (None, f"{B}/minimap2 -ax splice:sr -Y -t {t} {s(ref)} {s(r1)} {s(r2)} > {sam}"),
        "bwa_mem": (f"{B}/bwa index -p {idx}/bwa {s(ref)}", f"{B}/bwa mem -Y -t {t} {idx}/bwa {s(r1)} {s(r2)} > {sam}"),
        "bwa_mem2": (f"{B}/bwa-mem2 index -p {idx}/bwa2 {s(ref)}",
                     f"{B}/bwa-mem2 mem -Y -t {t} {idx}/bwa2 {s(r1)} {s(r2)} > {sam}"),
        "bowtie2_e2e": (f"{B}/bowtie2-build -q {s(ref)} {bx}", f"{B}/bowtie2 -p {t} -x {bx} -1 {s(r1)} -2 {s(r2)} -S {sam}"),
        "bowtie2_local": (f"{B}/bowtie2-build -q {s(ref)} {bx}",
                          f"{B}/bowtie2 --very-sensitive-local -p {t} -x {bx} -1 {s(r1)} -2 {s(r2)} -S {sam}"),
        "hisat2": (f"{B}/hisat2-build -q {s(ref)} {hx}", f"{B}/hisat2 -p {t} -x {hx} -1 {s(r1)} -2 {s(r2)} -S {sam}"),
        "hisat2_nosplice": (f"{B}/hisat2-build -q {s(ref)} {hx}",
                            f"{B}/hisat2 --no-spliced-alignment -p {t} -x {hx} -1 {s(r1)} -2 {s(r2)} -S {sam}"),
        "star": (f"mkdir -p {sx} && {B}/STAR --runMode genomeGenerate --genomeDir {sx} --genomeFastaFiles {s(ref)} "
                 f"--genomeSAindexNbases 6 --outFileNamePrefix {sx}/ > /dev/null",
                 f"{B}/STAR --runThreadN {t} --genomeDir {sx} --readFilesIn {s(r1)} {s(r2)} --readFilesCommand zcat "
                 f"--outFileNamePrefix {sam}. > /dev/null && mv {sam}.Aligned.out.sam {sam}"),
        "star_nosplice": (f"mkdir -p {sx} && {B}/STAR --runMode genomeGenerate --genomeDir {sx} "
                          f"--genomeFastaFiles {s(ref)} --genomeSAindexNbases 6 --outFileNamePrefix {sx}/ > /dev/null",
                          f"{B}/STAR --runThreadN {t} --genomeDir {sx} --readFilesIn {s(r1)} {s(r2)} "
                          f"--readFilesCommand zcat --alignIntronMax 1 --alignMatesGapMax 1000 "
                          f"--outFileNamePrefix {sam}. > /dev/null && mv {sam}.Aligned.out.sam {sam}"),
        "bbmap_local": (None, f"{B}/bbmap.sh -Xmx4g ref={s(ref)} path={idx}/bbmap in={s(r1)} in2={s(r2)} out={sam} "
                              f"local=t threads={t} overwrite=t nodisk=f 2> {sam}.log"),
        "subread": (f"{B}/subread-buildindex -o {ux} {s(ref)} > /dev/null",
                    f"{B}/subread-align -t 0 -T {t} -i {ux} -r {s(r1)} -R {s(r2)} -o {sam} --SAMoutput > /dev/null"),
        "segemehl": (f"{B}/segemehl.x -x {gx} -d {s(ref)} 2> /dev/null",
                     f"{B}/segemehl.x -t {t} -i {gx} -d {s(ref)} -q {s(r1)} -p {s(r2)} > {sam} 2> /dev/null"),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-r", "--ref", required=True)
    p.add_argument("-1", dest="r1", required=True)
    p.add_argument("-2", dest="r2", required=True)
    p.add_argument("-o", "--out", required=True)
    p.add_argument("-t", "--threads", type=int, default=8)
    p.add_argument("--mappers", default="")
    p.add_argument("--drop_unmapped", action="store_true",
                   help="after the timed mapping, keep only aligned records (real data: most reads are nuclear)")
    p.add_argument("--bin", default=os.path.expanduser("~/miniforge3/envs/mapbench/bin"))
    a = p.parse_args()
    out = Path(a.out)
    idx = out / "index"
    idx.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, PATH=f"{a.bin}:{os.environ['PATH']}")
    built = set()
    todo = mappers(a.bin, a.ref, idx, a.r1, a.r2, a.threads, "SAM")
    for name in (a.mappers.split(",") if a.mappers else todo):
        sam = out / f"{name}.sam"
        if sam.exists() and sam.stat().st_size:
            continue
        icmd, mcmd = mappers(a.bin, a.ref, idx, a.r1, a.r2, a.threads, shlex.quote(str(sam)))[name]
        if icmd and icmd not in built:
            subprocess.run(icmd, shell=True, check=True, env=env)
            built.add(icmd)
        (out / f"{name}.cmd").write_text(mcmd + "\n")
        r = subprocess.run(["/usr/bin/time", "-v", "-o", str(out / f"{name}.time"), "bash", "-c", mcmd], env=env)
        if r.returncode:
            print(f"[run_mappers] {name} failed ({r.returncode})")
            sam.unlink(missing_ok=True)
        else:
            if a.drop_unmapped:
                tmp = sam.with_suffix(".mapped.sam")
                subprocess.run([f"{a.bin}/samtools", "view", "-h", "-F", "4", "-o", str(tmp), str(sam)], check=True)
                tmp.replace(sam)
            print(f"[run_mappers] {name} done")


if __name__ == "__main__":
    main()
