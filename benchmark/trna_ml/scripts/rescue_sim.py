#!/usr/bin/env python3
"""Replay MitoM-ML with a rescue of the tRNA families it did not call, from the dumps of trna_dump.py and the MiTFi
calls of trna_only.py; writes runs/<variant>/<acc>/final.gbf + status.json for mito_bench_eval.py.

  mitfi:E     a MiTFi call of a family MitoM-ML did not call, E <= E, not overlapping a MitoM-ML call by > OVL nt
              (either strand; the best E per family)
  low:P       the best MitoM-ML candidate (p >= P) of a family it did not call, not overlapping a call by > OVL nt,
              its ends refined by the boundary model as for any call
  both:E:P    low:P first, then mitfi:E for the families still missing

usage: rescue_sim.py --set SET.tsv --root ROOT --variants mitfi:1e-3,low:0.2 [--ovl 10]
"""
import argparse
import csv
import gzip
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(HERE))
import mito_core  # noqa: E402
import trna_ml  # noqa: E402
from Bio import SeqIO  # noqa: E402


def mitfi_calls(gbf: Path, n: int):
    out = []
    for f in SeqIO.read(str(gbf), "genbank").features:
        if f.type != "tRNA":
            continue
        name = f.qualifiers.get("gene", [""])[0]
        notes = " ".join(f.qualifiers.get("note", []))
        m = re.search(r"E ([0-9.eE+-]+)", notes)
        ac = re.search(r"anticodon: ([ACGU]{3})", notes)
        parts = f.location.parts
        strand = f.location.strand
        if len(parts) == 2:                                   # across the origin: one interval with end > n
            tail = next(p for p in parts if int(p.end) == n)
            head = next(p for p in parts if int(p.start) == 0)
            s, e = int(tail.start) + 1, n + int(head.end)
        else:
            s, e = int(f.location.start) + 1, int(f.location.end)
        out.append({"family": name[3:], "start": s, "end": e, "strand": strand, "evalue": float(m.group(1)) if m else 1.0,
                    "anticodon": ac.group(1).lower().replace("u", "t") if ac else "", "src": "MiTFi"})
    return out


def overlap(a, b, n):
    best = 0
    for k in (-n, 0, n):
        best = max(best, min(a["end"], b["end"] + k) - max(a["start"], b["start"] + k) + 1)
    return best


def free(c, calls, n, ovl):
    return all(overlap(c, x, n) <= ovl for x in calls)


def rescue_low(d, genome, p_min, ovl, bclf):
    calls = list(d["calls"])
    have = {c["family"] for c in calls}
    best = {}
    for c in d["cands"]:
        if c["family"] in have or c.get("p", 0) < p_min:
            continue
        if c["family"] not in best or c["p"] > best[c["family"]]["p"]:
            best[c["family"]] = c
    add = []
    for fam, c in sorted(best.items(), key=lambda kv: -kv[1]["p"]):
        c = dict(c, _opts=c["opts"])
        trna_ml.refine_ends([c], genome, bclf)
        if free(c, calls + add, d["n"], ovl):
            add.append(dict(c, src="low"))
    return calls + add


def rescue_mitfi(calls, mitfi, e_max, ovl, n):
    have = {c["family"] for c in calls}
    add = []
    for c in sorted(mitfi, key=lambda c: c["evalue"]):
        if c["family"] in have or c["family"] in {a["family"] for a in add} or c["evalue"] > e_max:
            continue
        if free(c, calls + add, n, ovl):
            add.append(c)
    return calls + add


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--set", required=True)
    p.add_argument("--root", required=True, type=Path)
    p.add_argument("--variants", required=True)
    p.add_argument("--ovl", type=int, default=10)
    a = p.parse_args(argv)
    _, bclf = trna_ml.classifiers()
    rows = list(csv.DictReader(open(a.set), delimiter="\t"))
    for v in a.variants.split(","):
        kind, *par = v.split(":")
        name = "res_" + v.replace(":", "_").replace(".", "p").replace("-", "m")
        done = 0
        for r in rows:
            acc = r["accession"]
            dp = a.root / "cands" / f"{acc}.json.gz"
            mg = a.root / "runs" / "trna_mitfi" / acc / "final.gbf"
            if not dp.exists() or (kind in ("mitfi", "both") and not mg.exists()):
                continue
            d = json.load(gzip.open(dp, "rt"))
            rec = SeqIO.read(str(a.root / "gb" / f"{acc}.gb"), "genbank")
            genome, n, code = str(rec.seq).upper(), len(rec), int(r["genetic_code"] or 5)
            if kind == "mitfi":
                calls = rescue_mitfi(d["calls"], mitfi_calls(mg, n), float(par[0]), a.ovl, n)
            elif kind == "low":
                calls = rescue_low(d, genome, float(par[0]), a.ovl, bclf)
            else:
                calls = rescue_mitfi(rescue_low(d, genome, float(par[1]), a.ovl, bclf), mitfi_calls(mg, n),
                                     float(par[0]), a.ovl, n)
            feats = [mito_core.Feature("tRNA", mito_core.trna_name(c["family"], c.get("anticodon", ""), ""), c["start"] - 1,
                                       c["end"], c["strand"], c.get("src", "MitoM-ML"), anticodon=c.get("anticodon", ""))
                     for c in sorted(calls, key=lambda c: c["start"]) if d["circular"] or c["end"] <= n]
            out = a.root / "runs" / name / acc
            out.mkdir(parents=True, exist_ok=True)
            SeqIO.write(mito_core.genbank_record(genome, feats, code, acc.split(".")[0][:16],
                                                 "circular" if d["circular"] else "linear"), str(out / "final.gbf"), "genbank")
            (out / "status.json").write_text(json.dumps({"tool": name, "accession": acc, "status": "ok", "seconds": 0}))
            done += 1
        print(f"{name}: {done} records", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
