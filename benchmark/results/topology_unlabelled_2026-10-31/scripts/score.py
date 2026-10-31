"""Score a run of run_set.sh: topology per sequence vs the truth, by kind; genes found per sequence.
    score.py SET_DIR LABEL [LABEL2 ...] -> SET_DIR/score_<LABEL>.tsv (per sequence) and a summary on stdout"""
import collections
import csv
import sys
from pathlib import Path

from Bio import SeqIO

D = Path(sys.argv[1])
truth = {r["seq_id"]: r for r in csv.DictReader(open(D / "truth.tsv"), delimiter="\t")}
PCG = {"COX1", "COX2", "COX3", "CYTB", "ATP6", "ATP8", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6"}


def genes(rec):
    c = collections.Counter()
    for f in rec.features:
        g = (f.qualifiers.get("gene") or [""])[0].upper()
        if f.type == "CDS" and g in PCG:
            c["pcg"] += 1
        elif f.type == "rRNA":
            c["rrna"] += 1
        elif f.type == "tRNA":
            c["trna"] += 1
    return c


for label in sys.argv[2:]:
    rows = {}
    for dec in D.glob(f"{label}/*/reports/topology_decision.tsv"):
        gb = {r.id: r for r in SeqIO.parse(str(dec.parent.parent / "final_annotation_all.gb"), "genbank")} \
            if (dec.parent.parent / "final_annotation_all.gb").exists() else {}
        for r in csv.DictReader(open(dec), delimiter="\t"):
            g = genes(gb[r["sample_name"]]) if r["sample_name"] in gb else collections.Counter()
            rows[r["sample_name"]] = {**r, **{k: g[k] for k in ("pcg", "rrna", "trna")}}
    out = D / f"score_{label}.tsv"
    with open(out, "w") as fh:
        fh.write("seq_id\tkind\ttruth\tcalled\tok\tpcg\trrna\ttrna\treason\n")
        tab = collections.defaultdict(collections.Counter)
        for sid, t in truth.items():
            r = rows.get(sid)
            called = r["topology"] if r else "missing"
            ok = called == t["truth"]
            amb = " (NCBI 'complete' but linear)" if t["source"].startswith("EU293739") else ""
            key = f"{t['kind']}:{t['truth']}{amb}"
            tab[key]["n"] += 1
            tab[key]["ok"] += ok
            fh.write(f"{sid}\t{t['kind']}\t{t['truth']}\t{called}\t{'yes' if ok else 'NO'}\t"
                     f"{r['pcg'] if r else ''}\t{r['rrna'] if r else ''}\t{r['trna'] if r else ''}\t"
                     f"{r['reason'] if r else ''}\n")
    print(f"== {label}: {len(rows)} sequences annotated")
    print("kind:truth\tn\tcorrect")
    tot = collections.Counter()
    for k in sorted(tab):
        print(f"{k}\t{tab[k]['n']}\t{tab[k]['ok']}")
        if "NCBI 'complete'" not in k:
            tot["n"] += tab[k]["n"]; tot["ok"] += tab[k]["ok"]
    print(f"all (without the ambiguous record)\t{tot['n']}\t{tot['ok']} ({100 * tot['ok'] / max(tot['n'], 1):.1f} %)")
