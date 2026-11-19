import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
REPO = Path(sys.argv[1]); sys.path.insert(0, str(REPO))
import trna_ml, mito_env
from Bio import SeqIO
REF = mito_env.refseq_dir(mito_env.mitos_refdir(None)) / "ncRNA"
CMS, CMA = mito_env.find("cmsearch"), mito_env.find("cmalign")
rec = SeqIO.read(sys.argv[2], "genbank"); seq = str(rec.seq).upper(); L = len(seq)
lin = ";".join(rec.annotations.get("taxonomy", []))
base = trna_ml.predict(seq, lin, 5, REF, CMS, CMA, cpu=1)
def one(off):
    rot = seq[off:] + seq[:off]; r = {}
    for fix in (False, True):
        cs = trna_ml.predict(rot, lin, 5, REF, CMS, CMA, cpu=1, circular=False, circ_fix=fix)
        r[fix] = {(c["family"], c["strand"], (c["start"] - 1 + off) % L + 1, c["end"] - c["start"] + 1) for c in cs if c["end"] <= L}
    return off, r
offs = sorted({(c["start"] - 1 + d) % L for c in base for d in (0, 30)})
with ProcessPoolExecutor(24) as ex:
    for off, r in ex.map(one, offs):
        print(off, len(r[False]), len(r[True]), "only old", sorted(r[False] - r[True]), "only new", sorted(r[True] - r[False]), flush=True)
