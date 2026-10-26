"""Shared helpers for the invertebrate mt-tRNA project."""
import csv
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
INF = Path("/opt/micromamba/envs/invertrnax/bin")
PAN = Path("/opt/mitos_refs/refseq89m/ncRNA")
FAMS = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L1", "L2", "M", "N", "P", "Q", "R", "S1", "S2", "T", "V", "W", "Y"]
COMP = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")

# codon families of the two Leu / Ser tRNAs: anticodon (5'->3', DNA) -> family
L_S = {"tag": "L1", "taa": "L2", "tct": "S1", "tga": "S2", "gct": "S1", "tcc": "S1"}


def rc(s):
    return s.translate(COMP)[::-1]


def fam_of(aa, anticodon):
    if aa in ("L", "S"):
        return L_S.get(anticodon, "")
    return aa


def read_tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(p, rows, fields=None):
    rows = list(rows)
    with open(p, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields or list(rows[0]), delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def sh(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, check=False, **kw)


def hairpin_loops(ss):
    """Loop column spans (start, end exclusive) of the hairpins in a WUSS SS_cons line, 5'->3'."""
    out, i = [], 0
    n = len(ss)
    while i < n:
        if ss[i] in "<([{":
            j = i
            while j < n and ss[j] in "<([{":
                j += 1
            k = j
            while k < n and ss[k] not in "<([{>)]}":
                k += 1
            if k < n and ss[k] in ">)]}":
                out.append((j, k))
            i = k
        else:
            i += 1
    return out


def parse_sto(text):
    """Stockholm -> (seqs {name: aligned}, SS_cons, RF)."""
    seqs, ss, rf = {}, "", ""
    for line in text.splitlines():
        if line.startswith("#=GC SS_cons"):
            ss += line.split()[-1]
        elif line.startswith("#=GC RF"):
            rf += line.split()[-1]
        elif line and not line.startswith("#") and not line.startswith("//"):
            k, v = line.split()[:2]
            seqs[k] = seqs.get(k, "") + v
    return seqs, ss, rf


def anticodon_from_alignment(aln, ss):
    """Anticodon = the middle 3 nt of the 7-nt anticodon loop (2nd of the 3 cloverleaf hairpins)."""
    loops = hairpin_loops(ss)
    if len(loops) < 2:
        return ""
    lo = loops[1] if len(loops) >= 3 else loops[0]
    cols = [c for c in range(lo[0], lo[1])]
    nts = [aln[c] for c in cols if aln[c] not in ".-~"]
    loop_cons = [c for c in cols if ss[c] not in "~"]
    if len(nts) != 7 or len(loop_cons) < 5:
        return ""
    return "".join(nts[2:5]).lower().replace("u", "t")
