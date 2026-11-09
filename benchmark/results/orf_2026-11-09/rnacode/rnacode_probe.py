"""Prototype: RNAcode p-value of a candidate ORF + its homologs (protein-guided codon MSA, mito codons recoded)."""
import sys, subprocess, tempfile, os, re
sys.path.insert(0, "/media/cody/18T_R1/yzz/Mitom")
from Bio.Seq import Seq
import mito_orf as m
RNACODE = os.path.expanduser("~/miniforge3/envs/rnacode_test/bin/RNAcode")
MAFFT = os.path.expanduser("~/miniforge3/envs/mitotools/bin/mafft")

def recode(nt, code):
    """In frame: codons whose meaning differs from the standard code -> a standard synonym of the mito meaning."""
    from Bio.Data import CodonTable
    mt = CodonTable.unambiguous_dna_by_id[code].forward_table
    st = CodonTable.unambiguous_dna_by_id[1].forward_table
    syn = {}
    for c, aa in mt.items():
        if st.get(c) != aa:
            syn[c] = next(k for k, v in st.items() if v == aa)
    stops_mt = set(CodonTable.unambiguous_dna_by_id[code].stop_codons)
    out = []
    for i in range(0, len(nt) - len(nt) % 3, 3):
        c = nt[i:i+3]
        out.append(syn.get(c, c))
    return "".join(out)

def codon_msa(seqs, code):
    """seqs: {name: nt in frame}; protein MSA by MAFFT, back-translated."""
    prots = {k: str(Seq(v[: len(v)//3*3]).translate(table=code)).replace("*", "X") for k, v in seqs.items()}
    with tempfile.TemporaryDirectory() as td:
        f = os.path.join(td, "p.fa")
        open(f, "w").write("".join(f">{k}\n{v}\n" for k, v in prots.items()))
        aln = subprocess.run([MAFFT, "--quiet", "--auto", f], capture_output=True, text=True).stdout
    cur, al = None, {}
    for ln in aln.splitlines():
        if ln.startswith(">"): cur = ln[1:].strip(); al[cur] = ""
        else: al[cur] += ln.strip()
    out = {}
    for k, a in al.items():
        nt = recode(seqs[k], code); i = 0; s = []
        for ch in a:
            if ch == "-": s.append("---")
            else: s.append(nt[3*i:3*i+3]); i += 1
        out[k] = "".join(s)
    return out

def rnacode_p(seqs, code, query):
    al = codon_msa(seqs, code)
    names = [query] + [k for k in al if k != query]
    with tempfile.TemporaryDirectory() as td:
        f = os.path.join(td, "a.aln")
        w = max(len(n) for n in names) + 2
        L = len(al[query])
        with open(f, "w") as fh:
            fh.write("CLUSTAL W\n\n")
            for i in range(0, L, 60):
                for n in names: fh.write(f"{n[:30]:<32}{al[n][i:i+60]}\n")
                fh.write("\n")
        r = subprocess.run([RNACODE, "-t", "-b", "-n", "200", f], capture_output=True, text=True)
    best = None
    for ln in r.stdout.splitlines():
        c = ln.split("\t")
        if len(c) >= 11 and re.match(r"^\d", c[0]):
            # HSS# strand frame length from to name start end score p
            strand, frame, p = c[1], int(c[2]), float(c[10])
            if best is None or p < best[2]: best = (strand, frame, p, c)
    return best, r.stderr[-300:]
