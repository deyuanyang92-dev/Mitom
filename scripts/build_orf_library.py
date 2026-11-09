"""Reference library of mitochondrial ORFs for mito_orf.py (models/orf_library.fna.gz).

Every NCBI RefSeq animal (Metazoa) mitochondrial genome of 8-100 kb: the CDS that are not one of the 13 standard
protein genes (annotated ORFs, e.g. bivalve FORF / MORF, octocoral mutS), and the open reading frames of >= 40 codons
in the sequence no gene covers (mito_orf.candidates: the longest unannotated stretch of each frame between two stops,
from its first start codon; the record's own genetic code). Records with fewer than 10 annotated CDS are left out
(unannotated records: their ORFs are the 13 genes). These are the possible homologs
of an ORF found in a new genome: mito_orf.py writes an ORF only when it has a homolog that is itself open after
enough substitutions that a non-coding stretch would have gained a stop (see mito_orf.py).

Nucleotides in reading order; header >accession|start|end|strand|kind|code|phylum|organism (kind: cds = annotated,
ncr = ORF in unannotated sequence; 1-based coordinates on the record; code = genetic code).

Usage: python3 scripts/build_orf_library.py <work> [--out models/orf_library.faa.gz]
NCBI esearch / efetch (tool=mitom, no e-mail), 200 records per request; <work>/chunk_*.gb kept for a rerun.
"""
import argparse
import gzip
import re
import subprocess
import sys
import time
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from mito_orf import candidates, code_of, gene_key  # noqa: E402

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TERM = "Metazoa[Organism] AND mitochondrion[filter] AND refseq[filter] AND 8000:100000[SLEN]"
MIN_AA = 40


def curl(*args):
    for _ in range(5):
        p = subprocess.run(["curl", "-s", "--max-time", "900", *args], capture_output=True, text=True)
        if p.returncode == 0 and p.stdout.strip():
            return p.stdout
        time.sleep(10)
    raise SystemExit(f"NCBI request failed: {args[:3]}")


def fetch(work: Path) -> None:
    x = curl("-G", f"{EUTILS}/esearch.fcgi", "--data-urlencode", "db=nuccore", "--data-urlencode", "tool=mitom",
             "--data-urlencode", "usehistory=y", "--data-urlencode", "retmax=0", "--data-urlencode", f"term={TERM}")
    n = int(re.search(r"<Count>(\d+)</Count>", x).group(1))
    web, key = re.search(r"<WebEnv>(\S+?)</WebEnv>", x).group(1), re.search(r"<QueryKey>(\d+)</QueryKey>", x).group(1)
    print(f"{n} records", file=sys.stderr)
    for i in range(0, n, 200):
        f = work / f"chunk_{i // 200:04d}.gb"
        if f.exists() and f.stat().st_size and f.read_text()[-200:].rstrip().endswith("//"):
            continue
        f.write_text(curl(f"{EUTILS}/efetch.fcgi?db=nuccore&tool=mitom&rettype=gbwithparts&retmode=text"
                          f"&WebEnv={web}&query_key={key}&retstart={i}&retmax=200"))
        time.sleep(0.5)


def phylum_of(rec) -> str:
    tax = rec.annotations.get("taxonomy", [])
    return tax[tax.index("Metazoa") + 1] if "Metazoa" in tax and tax.index("Metazoa") + 1 < len(tax) else ""


def main():
    p = argparse.ArgumentParser()
    p.add_argument("work")
    p.add_argument("--out", default=str(Path(__file__).resolve().parents[1] / "models" / "orf_library.fna.gz"))
    p.add_argument("--no_fetch", action="store_true")
    a = p.parse_args()
    work = Path(a.work)
    work.mkdir(parents=True, exist_ok=True)
    if not a.no_fetch:
        fetch(work)
    n_rec = n_cds = n_ncr = 0
    seen = set()
    with gzip.open(a.out, "wt") as out:
        for f in sorted(work.glob("chunk_*.gb")):
            for rec in SeqIO.parse(str(f), "genbank"):
                if rec.id in seen:
                    continue
                seen.add(rec.id)
                n_rec += 1
                if sum(c.type == "CDS" for c in rec.features) < 10:
                    continue
                code = code_of(rec, 0)
                if not code:
                    continue
                org = rec.annotations.get("organism", "").replace("|", " ").replace(" ", "_")
                tag = f"{code}|{phylum_of(rec)}|{org}"
                for c in rec.features:
                    if c.type != "CDS" or gene_key(c):
                        continue
                    nt = str(c.extract(rec.seq)).upper()
                    if len(nt) >= 3 * MIN_AA:
                        st = "-" if c.location.strand == -1 else "+"
                        out.write(f">{rec.id}|{int(c.location.start) + 1}|{int(c.location.end)}|{st}|cds|{tag}\n"
                                  f"{nt}\n")
                        n_cds += 1
                for o in candidates(rec, code, MIN_AA):
                    out.write(f">{rec.id}|{o.start + 1}|{o.end}|{'+' if o.strand == 1 else '-'}|ncr|{tag}\n"
                              f"{o.piece_nt}\n")
                    n_ncr += 1
    print(f"{n_rec} records: {n_cds} non-standard CDS, {n_ncr} ORFs in unannotated sequence -> {a.out}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
