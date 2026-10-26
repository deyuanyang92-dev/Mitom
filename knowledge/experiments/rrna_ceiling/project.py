"""Feasibility test: rRNA ends projected from the closest annotated relatives (leave-species-out / leave-genus-out).

For every NCBI rRNA of a benchmark set: window = MitoM's predicted rRNA (else the truth) +- W nt, in the truth
strand.  blastn (dc-megablast) against the reference library (build_lib.py), drop references of the record's own
species (mode species) or genus (mode genus) and the record itself, align the top K references to the window
(semi-global, Biopython) and map each reference's rRNA start / end into the window.

usage: project.py SET.tsv ROOT LIBDIR OUTDIR {species,genus} [K] [threads]
Writes OUTDIR/proj_<mode>.tsv (one row per truth rRNA and reference rank) and prints a summary.
"""
import csv, os, re, subprocess, sys
from collections import defaultdict
from multiprocessing import Pool

from Bio import SeqIO
from Bio.Align import PairwiseAligner

W = 300
FLANK = 150
COMP = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")


def rc(s):
    return s.translate(COMP)[::-1]


def circ(seq, s, e):
    n = len(seq)
    return (seq * 3)[s + n:e + n]


def parse_loc(txt):
    m = re.match(r"(\d+)\.\.(\d+)([+-])", txt or "")
    return (int(m.group(1)) - 1, int(m.group(2)), 1 if m.group(3) == "+" else -1) if m else None


def binom(org):
    w = re.sub(r"[^A-Za-z0-9. ]", " ", org).split()
    return " ".join(w[:2]).lower()


ALIGNER = None


def aligner():
    global ALIGNER
    if ALIGNER is None:
        a = PairwiseAligner()
        a.mode = "global"
        a.match_score, a.mismatch_score = 2, -3
        a.open_gap_score, a.extend_gap_score = -5, -2
        a.target_end_gap_score = 0.0
        a.query_end_gap_score = 0.0
        ALIGNER = a
    return ALIGNER


def project(job):
    """job = (key, window, ref_id, region, rlen) -> (key, ref_id, pos5, pos3, identity) in window coordinates."""
    key, window, ref_id, region, rlen = job
    a = aligner()
    try:
        aln = a.align(window, region)[0]
    except Exception:
        return key, ref_id, None, None, 0.0
    tblocks, qblocks = aln.aligned
    r5, r3 = FLANK, FLANK + rlen - 1
    p5 = p3 = None
    matches = alen = 0
    for (ts, te), (qs, qe) in zip(tblocks, qblocks):
        if p5 is None and qe > r5:
            p5 = ts + max(0, r5 - qs)
        if qs <= r3:
            p3 = ts + min(qe - 1, r3) - qs
        # identity over the rRNA part of the reference
        lo, hi = max(qs, r5), min(qe, r3 + 1)
        if hi > lo:
            off = ts - qs
            for q in range(lo, hi):
                matches += window[q + off] == region[q]
            alen += hi - lo
    ident = matches / rlen if rlen else 0.0
    return key, ref_id, p5, p3, ident


def main():
    set_tsv, root, libdir, outdir, mode = sys.argv[1:6]
    K = int(sys.argv[6]) if len(sys.argv) > 6 else 5
    threads = int(sys.argv[7]) if len(sys.argv) > 7 else 64
    os.makedirs(outdir, exist_ok=True)
    meta = {r["accession"]: r for r in csv.DictReader(open(set_tsv), delimiter="\t")}
    lib = {r["ref_id"]: r for r in csv.DictReader(open(os.path.join(libdir, "lib.tsv")), delimiter="\t")}
    regions = {rec.id: str(rec.seq) for rec in SeqIO.parse(os.path.join(libdir, "lib.fa"), "fasta")}

    # truth rRNAs and MitoM predictions from the evaluation table
    rows = [r for r in csv.DictReader(open(os.path.join(root, "eval/all/features.tsv")), delimiter="\t")
            if r["tool"] == "mitom_native" and r["kind"] == "rRNA" and r["truth"]]
    queries = {}
    seqs = {}
    for r in rows:
        acc = r["accession"]
        if acc not in meta:
            continue
        if acc not in seqs:
            seqs[acc] = str(SeqIO.read(os.path.join(root, "gb", acc + ".gb"), "genbank").seq).upper()
        seq = seqs[acc]
        n = len(seq)
        t = parse_loc(r["truth"])
        p = parse_loc(r["pred"]) if r["class"] not in ("missed", "strand_error") else None
        s0, e0 = (p[0], p[1]) if p else (t[0], t[1])
        ws, we = s0 - W, e0 + W
        win = circ(seq, ws, we)
        if t[2] == -1:
            win = rc(win)
        key = f"{acc}|{r['name']}|{t[0]}"
        # truth ends in window coordinates (transcript direction)
        if t[2] == 1:
            t5, t3 = t[0] - ws, t[1] - 1 - ws
        else:
            t5, t3 = we - 1 - (t[1] - 1), we - 1 - t[0]
        queries[key] = dict(acc=acc, gene=r["name"], win=win, t5=t5, t3=t3, n=n, cls=r["class"],
                            delta=r["delta"], pred=r["pred"], truth=r["truth"])
    qfa = os.path.join(outdir, "queries.fa")
    with open(qfa, "w") as fh:
        for k, q in queries.items():
            fh.write(f">{k}\n{q['win']}\n")
    db = os.path.join(libdir, "lib.fa")
    if not os.path.exists(db + ".nsq") and not os.path.exists(db + ".00.nsq"):
        subprocess.run(["makeblastdb", "-in", db, "-dbtype", "nucl"], check=True, stdout=subprocess.DEVNULL)
    hits_path = os.path.join(outdir, "hits.tsv")
    if not os.path.exists(hits_path):
        subprocess.run(["blastn", "-task", "dc-megablast", "-query", qfa, "-db", db, "-evalue", "1e-10",
                        "-max_target_seqs", "400", "-num_threads", str(threads),
                        "-outfmt", "6 qseqid sseqid pident length evalue bitscore", "-out", hits_path], check=True)
    best = defaultdict(dict)
    for line in open(hits_path):
        q, s, pid, ln, ev, bits = line.rstrip("\n").split("\t")
        b = float(bits)
        if b > best[q].get(s, 0):
            best[q][s] = b
    jobs = []
    for key, q in queries.items():
        m = meta[q["acc"]]
        sp, genus = binom(m["species"]), m["genus"].lower()
        cands = []
        for s, b in sorted(best.get(key, {}).items(), key=lambda x: -x[1]):
            L = lib[s]
            if L["gene"] != q["gene"] or L["accession"].split(".")[0] == q["acc"].split(".")[0]:
                continue
            if binom(L["organism"]) == sp:
                continue
            if mode == "genus" and L["organism"].split()[0].lower() == genus:
                continue
            cands.append((s, b))
            if len(cands) >= K:
                break
        q["cands"] = cands
        for s, b in cands:
            jobs.append((key, q["win"], s, regions[s], int(lib[s]["length"])))
    with Pool(threads) as pool:
        res = pool.map(project, jobs, chunksize=4)
    proj = defaultdict(dict)
    for key, s, p5, p3, ident in res:
        proj[key][s] = (p5, p3, ident)
    out = os.path.join(outdir, f"proj_{mode}.tsv")
    with open(out, "w") as fh:
        fh.write("key\taccession\tphylum\tgene\tmitom_class\tmitom_delta\trank\tref\tref_organism\tbits\tident"
                 "\trelation\td5\td3\tref_gap5\tref_gap3\n")
        for key, q in queries.items():
            m = meta[q["acc"]]
            for rank, (s, b) in enumerate(q["cands"], 1):
                p5, p3, ident = proj[key].get(s, (None, None, 0))
                L = lib[s]
                tax = L["taxonomy"].split(";")
                rel = ("genus" if L["organism"].split()[0].lower() == m["genus"].lower() else
                       "family" if m.get("family") and m["family"] in tax else
                       "order" if m.get("order") and m["order"] in tax else
                       "class" if m.get("class") and m["class"] in tax else
                       "phylum" if m["phylum"] in tax else "other")
                d5 = "" if p5 is None else p5 - q["t5"]
                d3 = "" if p3 is None else p3 - q["t3"]
                fh.write(f"{key}\t{q['acc']}\t{m['phylum']}\t{q['gene']}\t{q['cls']}\t{q['delta']}\t{rank}\t{s}"
                         f"\t{L['organism']}\t{b:.0f}\t{ident:.3f}\t{rel}\t{d5}\t{d3}\t{L['gap5']}\t{L['gap3']}\n")
            if not q["cands"]:
                fh.write(f"{key}\t{q['acc']}\t{m['phylum']}\t{q['gene']}\t{q['cls']}\t{q['delta']}\t0\t\t\t\t\t\t\t\t\t\n")
    print("queries", len(queries), "jobs", len(jobs), "->", out)


if __name__ == "__main__":
    main()
