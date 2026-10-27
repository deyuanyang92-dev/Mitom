#!/usr/bin/env python3
"""Every tRNA of the final method (MitoM-ML through drop_wrap_duplicates, the MiTFi rescue) with what the caller knew,
joined to its benchmark class (features.tsv of mito_bench_eval): source (ML / MiTFi rescue), p, whether it is a
weaker second call of a family already called, how many families MitoM-ML called in the record, MiTFi's E.

usage: call_audit.py SET.tsv ROOT FEATURES.tsv TOOL MITFI_RUNS --accs FILE > audit.tsv
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
sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_trna  # noqa: E402
import trna_ml  # noqa: E402
from rescue_sim import mitfi_calls, rescue_mitfi  # noqa: E402
from Bio import SeqIO  # noqa: E402

FOUND = ("exact", "le3", "gt3", "name_unverifiable")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("set")
    p.add_argument("root", type=Path)
    p.add_argument("features")
    p.add_argument("tool")
    p.add_argument("mitfi_runs")
    p.add_argument("--accs", type=Path, required=True)
    a = p.parse_args(argv)
    keep = {x.strip() for x in a.accs.read_text().split() if x.strip()}
    cls = {}
    for r in csv.DictReader(open(a.features), delimiter="\t"):
        if r["tool"] == a.tool and r["kind"] == "tRNA" and r["pred"] and r["accession"] in keep:
            m = re.match(r"(\d+)\.\.(\d+)([+-])", r["pred"])
            cls[(r["accession"], int(m.group(1)), int(m.group(2)), 1 if m.group(3) == "+" else -1)] = r["class"]
    clf, bclf = trna_ml.classifiers()
    w = csv.writer(sys.stdout, delimiter="\t")
    w.writerow(["accession", "family", "start", "end", "strand", "src", "p", "evalue", "second_copy", "n_ml_families",
                "class", "correct"])
    for r in csv.DictReader(open(a.set), delimiter="\t"):
        acc = r["accession"]
        if acc not in keep:
            continue
        d = json.load(gzip.open(a.root / "cands" / f"{acc}.json.gz", "rt"))
        n = d["n"]
        genome = str(SeqIO.read(str(a.root / "gb" / f"{acc}.gb"), "genbank").seq).upper()
        cands = trna_ml.candidates(trna_ml.drop_wrap_duplicates(d["hits"], n),
                                   trna_ml.ext_hits([mito_trna.Pred(**e) for e in d["ext"]]), int(r["genetic_code"] or 5))
        ml = trna_ml.call(cands, clf, trna_ml.THRESHOLD)
        trna_ml.refine_ends(ml, genome, bclf)
        mg = a.root / "runs" / a.mitfi_runs / acc / "final.gbf"
        calls = rescue_mitfi([dict(c, src="ML") for c in ml], mitfi_calls(mg, n), 0.1, 10, n) if mg.exists() else \
            [dict(c, src="ML") for c in ml]
        best_p = {}
        for c in ml:
            best_p[c["family"]] = max(best_p.get(c["family"], 0), c["p"])
        fams = len({c["family"] for c in ml})
        for c in calls:
            klass = cls.get((acc, c["start"], c["end"], c["strand"]), cls.get((acc, c["start"], c["end"] - n, c["strand"]), "?"))
            second = c["src"] == "ML" and c["p"] < best_p.get(c["family"], 0)
            w.writerow([acc, c["family"], c["start"], c["end"], c["strand"], c["src"], f"{c.get('p', 0):.3f}",
                        c.get("evalue", ""), int(second), fams, klass, int(klass in FOUND)])
    return 0


if __name__ == "__main__":
    sys.exit(main())
