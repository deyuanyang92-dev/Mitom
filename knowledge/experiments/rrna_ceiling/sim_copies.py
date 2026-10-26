"""Would adding the second copy of an rRNA (mito_rrna_hmm.best_spans' runner-up) find NCBI's duplicated rRNAs?

The core keeps only the best span per gene; best_spans already returns a runner-up that covers mostly the same model
columns elsewhere in the genome (a copy). Per record and gene: nhmmer with the leave-out library (as mitom_lo),
the copy's mito_rrna_hmm._verdict (CDS and rRNA positions taken = MitoM's final CDS + the best span), and whether
the copy hits an NCBI rRNA of that gene (>= 50 % of it) that the best span does not.
usage: sim_copies.py SET.tsv ROOT LIB|split=LIB_PATTERN TOOL OUT.tsv   (LIB_PATTERN with {split}, e.g. rrna_lib_{split}_v2)
"""
import csv
import os
import sys
import tempfile
from multiprocessing import Pool
from pathlib import Path

sys.path.insert(0, "/media/deyuan/18T_B/yzz/MitoM")
import mito_rrna_hmm as R   # noqa: E402
from Bio import SeqIO       # noqa: E402

SET, ROOT, LIBARG, TOOL, OUT = sys.argv[1:6]


def lib_for(split):
    return Path(ROOT) / (LIBARG.format(split=split) if "{split}" in LIBARG else LIBARG)


def ncbi_rrnas(rec, gene):
    out = []
    for f in rec.features:
        if f.type == "rRNA" and R.gene_of(f) == gene:
            out.append({i % len(rec) for p in f.location.parts for i in range(int(p.start), int(p.end))})
    return out


def one(row):
    acc, phy = row["accession"], row["phylum"]
    lib = lib_for(row["split"])
    man = R.read_manifest(lib)
    rec = SeqIO.read(os.path.join(ROOT, "gb", acc + ".gb"), "genbank")
    n = len(rec)
    circ = rec.annotations.get("topology") == "circular"
    seq = str(rec.seq)
    fin = os.path.join(ROOT, "runs", TOOL, acc, "final.gbf")
    cds = []
    if os.path.exists(fin):
        for f in SeqIO.read(fin, "genbank").features:
            if f.type == "CDS":
                cds.append({i % n for p in f.location.parts for i in range(int(p.start), int(p.end))})
    rows = []
    with tempfile.TemporaryDirectory() as td:
        fa = Path(td) / "q.fa"
        fa.write_text(f">{acc}\n{seq + seq[:min(R.ORIGIN_EXT, n)] if circ else seq}\n")
        for g in ("rrnS", "rrnL"):
            m = man.get((phy, g))
            hmm = lib / f"{phy}.{g}.hmm"
            if not m or not hmm.exists():
                continue
            best, copy = R.best_spans(R.nhmmer(hmm, fa), n if circ else None)
            truth = ncbi_rrnas(rec, g)
            out = {"accession": acc, "phylum": phy, "gene": g, "n_ncbi": len(truth), "best": "", "copy": "",
                   "copy_score": "", "best_score": "", "copy_cover": "", "verdict": "", "copy_hits_ncbi": "",
                   "best_hits_ncbi": ""}
            if best is None:
                rows.append(out)
                continue
            leng = R.model_len(hmm)
            bpos = {i % n for i in range(best.start, best.end)}
            hit_b = [i for i, t in enumerate(truth) if len(t & bpos) >= 0.5 * len(t)]
            out.update(best=R._loc(best, n), best_score=f"{best.score:.0f}", best_hits_ncbi=",".join(map(str, hit_b)))
            if copy is not None:
                cpos = {i % n for i in range(copy.start, copy.end)}
                hit_c = [i for i, t in enumerate(truth) if len(t & cpos) >= 0.5 * len(t) and i not in hit_b]
                out.update(copy=R._loc(copy, n), copy_score=f"{copy.score:.0f}",
                           copy_cover=f"{(copy.hmm_to - copy.hmm_from + 1) / leng:.2f}",
                           verdict=R._verdict(copy, m, leng, cds + [bpos], n), copy_hits_ncbi=",".join(map(str, hit_c)))
            rows.append(out)
    return rows


def main():
    meta = list(csv.DictReader(open(SET), delimiter="\t"))
    with Pool(36) as p:
        res = p.map(one, meta, chunksize=4)
    rows = [r for rs in res for r in rs]
    with open(OUT, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    dup = [r for r in rows if r["n_ncbi"] >= 2]
    cop = [r for r in rows if r["copy"]]
    add = [r for r in cop if r["verdict"].startswith("added")]
    print(f"gene searches {len(rows)}; NCBI has >= 2 copies: {len(dup)}; runner-up copy found: {len(cop)}; "
          f"copy passes the verdict: {len(add)}; of these hit an NCBI rRNA the best span does not: "
          f"{sum(bool(r['copy_hits_ncbi']) for r in add)}")


if __name__ == "__main__":
    main()
