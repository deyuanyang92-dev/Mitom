"""Per candidate with homologs: current evidence (p_open, p_frame, class) and RNAcode p. True = truth ORFs
(recall set, record left out of the library); false = candidates in the 362 intron-free records."""
import sys, csv, glob, os
sys.path.insert(0, "/media/cody/18T_R1/yzz/Mitom"); sys.path.insert(0, os.path.dirname(__file__))
from multiprocessing import Pool
from Bio import SeqIO
import mito_orf as m
from rnacode_probe import rnacode_p
R = "/media/cody/18T_R1/yzz/_realtest/orf_2026-11-09"
CAP = {}
orig = m.combined_evidence
def cap(o, homologs):
    CAP[id(o)] = [(h["name"], h["s_nt"]) for h in homologs]
    return orig(o, homologs)
m.combined_evidence = cap

def run(args):
    path, targets = args            # targets: list of (start,end,strand) or None (= all candidates)
    r = SeqIO.read(path, "genbank"); code = m.code_of(r, 5)
    cands = m.candidates(r, code)
    if targets is not None:
        sel = []
        for s, e, st in targets:
            cs = [o for o in cands if o.strand == st and max(0, min(o.end, e) - max(o.start, s)) >= 0.5 * (e - s)]
            if cs: sel.append(max(cs, key=lambda o: o.codons))
        cands = sel
    if not cands: return []
    CAP.clear()
    m.search(cands, m.LIBRARY, {r.id, r.id.split(".")[0]}, 1)
    out = []
    for o in cands:
        m.classify(o, code)
        hs = CAP.get(id(o), [])
        if not hs: continue
        seqs = {"query": o.piece_nt}
        for k, (n, nt) in enumerate(hs[:12]): seqs[f"h{k}"] = nt
        try:
            best, err = rnacode_p(seqs, code, "query")
        except Exception as ex:
            best = None
        c = [h for h in o.hits if "p_pos3" in h]
        c = c[0] if c else {}
        out.append([r.id, o.start + 1, o.end, o.strand, o.codons, len(hs), o.cls, f"{c.get('p_open', 1) if 'p_open' in c else (1-m.stop_share(o.piece_nt, code))**c.get('subs',0):.2e}",
                    f"{c.get('p_pos3', 1):.2e}", best[0] if best else "", best[1] if best else "", f"{best[2]:.3g}" if best else "1"])
    return out

if __name__ == "__main__":
    which = sys.argv[1]
    jobs = []
    if which == "true":
        tr = [t for t in csv.DictReader(open(R + "/recall/truth.tsv"), delimiter="\t")
              if t["record"] != "NC_137528.1" and "atpase subunit" not in t["product"].lower()]
        by = {}
        for t in tr: by.setdefault(t["record"], []).append((int(t["start"]) - 1, int(t["end"]), 1 if t["strand"] == "1" else -1))
        jobs = [(f"{R}/recall/in/{a}.gb", v) for a, v in by.items()]
    else:
        jobs = [(f, None) for f in glob.glob(R + "/../intron_algo_2026-11-09/neg/gb/*.gb")]
    w = csv.writer(sys.stdout, delimiter="\t")
    w.writerow(["record", "start", "end", "strand", "codons", "n_hom", "cls", "p_open", "p_frame", "rc_strand", "rc_frame", "rc_p"])
    with Pool(40) as p:
        for rows in p.imap_unordered(run, jobs):
            w.writerows(rows)
