"""Every de novo intron candidate with its evidence and its truth label.

Input: an intron_only.py run with MITOM_INTRON_WRITE_ALL=1 (every model written that has no internal stop, a start and
no clash) and its eval_introns.py table. Each intron row of <acc>.report.tsv is mapped onto the NCBI record and matched
to the eval table's predictions: exact / le6 / gt6 -> true, fp -> false.

Usage: python3 label_candidates.py <root> <label> > candidates.tsv
"""
import csv
import sys
from pathlib import Path

from Bio import SeqIO

M = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(M / "benchmark"))
from eval_introns import mapper  # noqa: E402

root, label = Path(sys.argv[1]), sys.argv[2]
d = root / "intron_only" / label
ev = list(csv.DictReader(open(root / "eval_io" / f"{label}_introns.tsv"), delimiter="\t"))
by = {}
for r in ev:
    if r["pred"]:
        s, e = (int(x) for x in r["pred"].split(".."))
        by[(r["accession"], s, e)] = r
w = csv.writer(sys.stdout, delimiter="\t")
w.writerow(["accession", "group", "gene", "intron_start", "intron_end", "length", "label", "class", "type", "typed",
            "splice_score", "junction_span", "junction_missing_aa", "junction_extra_aa", "homolog_coverage",
            "splicing_gain", "confidence", "evidence"])
for rep in sorted(d.glob("*.report.tsv")):
    acc = rep.name[: -len(".report.tsv")]
    pred = SeqIO.read(str(d / f"{acc}.gb"), "genbank")
    ref = SeqIO.read(str(root / "gb" / f"{acc}.gb"), "genbank")
    f = mapper(str(ref.seq).upper(), str(pred.seq).upper())
    for r in csv.DictReader(open(rep), delimiter="\t"):
        if not r.get("intron_start") or r.get("written") != "yes":
            continue
        s, e = int(r["intron_start"]), int(r["intron_end"])
        st = 1 if r.get("strand") == "+" else -1
        if s > e or f is None:
            continue
        ms, me, _ = f(s - 1, e, st)
        hit = by.get((acc, ms + 1, me)) or next((v for (a, x, y), v in by.items()
                                                 if a == acc and abs(x - ms - 1) <= 1 and abs(y - me) <= 1), None)
        cls = hit["class"] if hit else "unmatched"
        lab = {"exact": "true", "le6": "true", "gt6": "true", "fp": "false"}.get(cls, "?")
        note = r.get("note", "")
        typ = note.split(";")[0]
        w.writerow([acc, hit["group"] if hit else "", r["gene"], s, e, e - s + 1, lab, cls, typ,
                    "yes" if typ.startswith("group") else "no", r.get("splice_score", ""), r.get("junction_span", ""),
                    r.get("junction_missing_aa", ""), r.get("junction_extra_aa", ""), r.get("homolog_coverage", ""),
                    r.get("splicing_gain", ""), r.get("confidence", ""), r.get("evidence", "")[:120]])
