"""Recall and false positives of mito_pcg_hmm.py.

Truth: annotated NCBI mitogenomes (the RefSeq download of scripts/build_orf_library.py and the NCBI annelid download
of ../intron_annelid_2026-11-09/), up to --per_phylum records per phylum with all 13 genes, plus every Sabellida
record with atp8. For each record and gene: the gene's CDS and gene features are removed, the gene is searched with
the references of other genera only (leave-genus-out), and the result is compared with the NCBI gene: correct =
same strand, >= 50 % overlap (start / end differences reported); wrong place; not found.
False positives: atp8 searched in Platyhelminthes and Nematoda records without an annotated atp8 (most have none).
Usage: python3 eval_pcg_hmm.py <refseq dir> <annelid dir> <out.tsv> [--per_phylum 15] [--jobs 32]
"""
import argparse
import csv
import random
import sys
from collections import defaultdict
from multiprocessing import Pool
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import mito_pcg_hmm as P  # noqa: E402
from mito_common import iter_gb_files  # noqa: E402
from mito_names import PCG, gene_name  # noqa: E402

REFS = None


def init():
    global REFS
    REFS = P.load_refs()


def job(args):
    rec, kind = args
    genus = rec.annotations.get("organism", "").split(" ")[0]
    out = []
    if kind == "negative":
        rows = P.search_record(rec, REFS, 5, ["atp8"], exclude_genus={genus}, write=False)
        for r in rows:
            out.append([rec.id, kind, phylum(rec), "atp8", "", r["status"], r.get("start", ""), r.get("end", ""),
                        r.get("evalue", ""), "", ""])
        return out
    for g in sorted(PCG):
        feats = [f for f in rec.features if f.type in ("CDS", "gene") and gene_name(f.qualifiers, f.type) == g]
        cds = [f for f in feats if f.type == "CDS"]
        if not cds:
            continue
        t = cds[0]
        r2 = rec[:]
        r2.annotations = dict(rec.annotations)
        r2.features = [f for f in rec.features if f not in feats]
        rows = P.search_record(r2, REFS, 5, [g], exclude_genus={genus}, write=False)
        r = rows[0] if rows else {"status": "not found"}
        ts, te, tst = int(t.location.start) + 1, int(t.location.end), t.location.strand
        cls, ds, de = r["status"], "", ""
        if r["status"] == "found":
            ps, pe, pst = int(r["start"]), int(r["end"]), 1 if r["strand"] == "+" else -1
            ov = max(0, min(pe, te) - max(ps, ts) + 1)
            cls = "correct" if pst == tst and ov >= 0.5 * (te - ts + 1) else "wrong_place"
            ds, de = (ps - ts, pe - te) if tst == 1 else (te - pe, ts - ps)
        out.append([rec.id, kind, phylum(rec), g, len(t) // 3, cls, r.get("start", ""), r.get("end", ""),
                    r.get("evalue", ""), ds, de])
    return out


def phylum(rec):
    tax = rec.annotations.get("taxonomy", [])
    return tax[tax.index("Metazoa") + 1] if "Metazoa" in tax and tax.index("Metazoa") + 1 < len(tax) else ""


def main():
    p = argparse.ArgumentParser()
    p.add_argument("refseq")
    p.add_argument("annelid")
    p.add_argument("out")
    p.add_argument("--per_phylum", type=int, default=15)
    p.add_argument("--jobs", type=int, default=32)
    a = p.parse_args()
    by = defaultdict(list)
    sab, neg = [], []
    seen = set()
    for d in (a.refseq, a.annelid):
        for f in iter_gb_files(Path(d)):
            for rec in SeqIO.parse(str(f), "genbank"):
                if rec.id in seen:
                    continue
                seen.add(rec.id)
                genes = {gene_name(x.qualifiers, "CDS") for x in rec.features if x.type == "CDS"}
                tax = rec.annotations.get("taxonomy", [])
                if "Sabellida" in tax and "atp8" in genes:
                    sab.append(rec)
                elif genes >= set(PCG) and d == a.refseq:
                    by[phylum(rec)].append(rec)
                elif d == a.refseq and ("Platyhelminthes" in tax or "Nematoda" in tax) and "atp8" not in genes \
                        and len(genes) >= 11:
                    neg.append(rec)
    rnd = random.Random(1)
    truth = [r for ph, rs in sorted(by.items()) for r in rnd.sample(rs, min(a.per_phylum, len(rs)))] + sab
    neg = rnd.sample(neg, min(200, len(neg)))
    print(f"{len(truth)} truth records ({len(sab)} Sabellida), {len(neg)} negatives", file=sys.stderr)
    with Pool(a.jobs, initializer=init) as pool, open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["record", "kind", "phylum", "gene", "codons", "class", "start", "end", "evalue", "d_start", "d_end"])
        for rows in pool.imap_unordered(job, [(r, "truth") for r in truth] + [(r, "negative") for r in neg]):
            w.writerows(rows)
            fh.flush()


if __name__ == "__main__":
    main()
