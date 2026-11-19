"""Rotation invariance of trna_ml.predict: calls on the record vs calls on rotations that put each tRNA at the
origin (start at position 1, and straddling the origin by 30 nt). Usage: rot.py REPO out.tsv acc.gb [...]"""
import sys, json, time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
REPO = Path(sys.argv[1]); sys.path.insert(0, str(REPO))
import os; os.chdir(REPO)
import trna_ml, mito_env
from Bio import SeqIO
REF = mito_env.refseq_dir(mito_env.mitos_refdir(None)) / "ncRNA"
CMS, CMA = mito_env.find("cmsearch"), mito_env.find("cmalign")
FIX = os.environ.get("MITOM_TRNA_CIRC", "") not in ("", "0")

def run(args):
    seq, lin, code, off = args
    L = len(seq)
    rot = seq[off:] + seq[:off]
    t = time.time()
    calls = trna_ml.predict(rot, lin, code, REF, CMS, CMA, cpu=1, circ_fix=FIX)
    out = set()
    for c in calls:
        s = (c["start"] - 1 + off) % L + 1
        out.add((c["family"], c["strand"], s, c["end"] - c["start"] + 1))
    return off, out, time.time() - t

def main():
    rows = []
    for gb in sys.argv[3:]:
        rec = SeqIO.read(gb, "genbank")
        seq, lin = str(rec.seq).upper(), ";".join(rec.annotations.get("taxonomy", []))
        code = 5
        for f in rec.features:
            if f.type == "CDS" and "transl_table" in f.qualifiers:
                code = int(f.qualifiers["transl_table"][0]); break
        L = len(seq)
        _, base, t0 = run((seq, lin, code, 0))
        offs = sorted({(s - 1) % L for _, _, s, _ in base} | {(s - 1 + 30) % L for _, _, s, _ in base}
                      | {(s - 1 + n - 20) % L for _, _, s, n in base})
        with ProcessPoolExecutor(24) as ex:
            res = list(ex.map(run, [(seq, lin, code, o) for o in offs]))
        for off, calls, t in res:
            lost, gained = base - calls, calls - base
            rows.append((rec.id, L, off, len(base), len(calls), sorted(lost), sorted(gained)))
            if lost or gained:
                print(rec.id, off, "lost", sorted(lost), "gained", sorted(gained), flush=True)
        print(rec.id, "L", L, "base", len(base), "rotations", len(offs), "differing",
              sum(1 for r in rows if r[0] == rec.id and (r[5] or r[6])), f"t0 {t0:.0f}s", flush=True)
    with open(sys.argv[2], "w") as fh:
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
main()
