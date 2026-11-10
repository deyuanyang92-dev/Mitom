"""atp8 / atp6 recall of mito_pcg_hmm.py by real phylum (leave-genus-out).

Up to --per_phylum NCBI records per phylum (names in PHYLA, found in the record's taxonomy; Sabellida of the annelids
apart, the tube worms' order) with the gene annotated; the gene is removed and searched again with references of
other genera only; correct = same strand, >= 50 % overlap. Records: the RefSeq download of scripts/build_orf_library.py
and the NCBI annelid download of ../intron_annelid_2026-11-09/.
Usage: python3 eval_atp.py <refseq dir> <annelid dir> <out.tsv> [--per_phylum 15] [--jobs 40]
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
from mito_names import gene_name  # noqa: E402

PHYLA = ["Porifera", "Cnidaria", "Ctenophora", "Placozoa", "Xenacoelomorpha", "Chordata", "Echinodermata",
         "Hemichordata", "Arthropoda", "Nematoda", "Tardigrada", "Onychophora", "Priapulida", "Nematomorpha",
         "Kinorhyncha", "Annelida", "Mollusca", "Brachiopoda", "Bryozoa", "Nemertea", "Platyhelminthes", "Rotifera",
         "Phoronida", "Entoprocta", "Chaetognatha", "Gastrotricha", "Cycliophora"]
GENES = ("atp8", "atp6")
REFS = None


def group(rec):
    t = rec.annotations.get("taxonomy", [])
    if "Sabellida" in t:
        return "Annelida: Sabellida"
    return next((p for p in PHYLA if p in t), "")


def init():
    global REFS
    REFS = P.load_refs()


def job(args):
    rec, g = args
    genus = rec.annotations.get("organism", "").split(" ")[0]
    feats = [f for f in rec.features if f.type in ("CDS", "gene") and gene_name(f.qualifiers, f.type) == g]
    t = [f for f in feats if f.type == "CDS"][0]
    r2 = rec[:]
    r2.annotations = dict(rec.annotations)
    r2.features = [f for f in rec.features if f not in feats]
    r = (P.search_record(r2, REFS, 5, [g], exclude_genus={genus}, write=False) or [{"status": "not found"}])[0]
    ts, te, tst = int(t.location.start) + 1, int(t.location.end), t.location.strand
    cls, ds, de = r["status"], "", ""
    if cls == "found":
        ps, pe, pst = int(r["start"]), int(r["end"]), 1 if r["strand"] == "+" else -1
        ov = max(0, min(pe, te) - max(ps, ts) + 1)
        cls = "correct" if pst == tst and ov >= 0.5 * (te - ts + 1) else "wrong_place"
        ds, de = (ps - ts, pe - te) if tst == 1 else (te - pe, ts - ps)
    return [rec.id, group(rec), rec.annotations.get("organism", ""), g, len(t) // 3, cls, r.get("evalue", ""),
            r.get("prior", ""), ds, de]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("refseq")
    p.add_argument("annelid")
    p.add_argument("out")
    p.add_argument("--per_phylum", type=int, default=15)
    p.add_argument("--jobs", type=int, default=40)
    a = p.parse_args()
    pool = defaultdict(list)
    seen = set()
    for d in (a.refseq, a.annelid):
        for f in iter_gb_files(Path(d)):
            for rec in SeqIO.parse(str(f), "genbank"):
                if rec.id in seen:
                    continue
                seen.add(rec.id)
                genes = {gene_name(x.qualifiers, "CDS") for x in rec.features if x.type == "CDS"}
                for g in GENES:
                    if g in genes and group(rec):
                        pool[(group(rec), g)].append(rec)
    rnd = random.Random(1)
    jobs = [(r, g) for (ph, g), rs in sorted(pool.items(), key=lambda x: x[0])
            for r in rnd.sample(rs, min(a.per_phylum, len(rs)))]
    print(f"{len(jobs)} searches in {len({k[0] for k in pool})} groups", file=sys.stderr)
    with Pool(a.jobs, initializer=init) as pl, open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["record", "group", "organism", "gene", "codons", "class", "evalue", "prior", "d_start", "d_end"])
        for row in pl.imap_unordered(job, jobs):
            w.writerow(row)


if __name__ == "__main__":
    main()
