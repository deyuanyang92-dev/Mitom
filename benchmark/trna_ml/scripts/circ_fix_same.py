import sys
from pathlib import Path
REPO = Path(sys.argv[1]); sys.path.insert(0, str(REPO))
import trna_ml, mito_env
from Bio import SeqIO
REF = mito_env.refseq_dir(mito_env.mitos_refdir(None)) / "ncRNA"
CMS, CMA = mito_env.find("cmsearch"), mito_env.find("cmalign")
for gb in sys.argv[2:]:
    rec = SeqIO.read(gb, "genbank"); seq = str(rec.seq); lin = ";".join(rec.annotations.get("taxonomy", []))
    r = [{(c["family"], c["strand"], c["start"], c["end"], c["p"], c["bound_p"]) for c in trna_ml.predict(seq, lin, 5, REF, CMS, CMA, cpu=4, circ_fix=f)} for f in (False, True)]
    print(rec.id, len(r[0]), "identical" if r[0] == r[1] else f"DIFF old-only {sorted(r[0]-r[1])} new-only {sorted(r[1]-r[0])}", flush=True)
