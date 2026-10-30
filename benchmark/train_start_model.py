"""Train the start model of mito_curate P11 (mito_start.py) on NCBI records.

1. mito_curate.py --dump_starts --start_model no on the pre-curation annotation of each benchmark record
   (<root>/runs/mitom_lo/<acc>/raw/work/04.final_gb, the leave-out runs: the record's own proteins were not used): every
   start candidate of every protein gene with its features;
2. each candidate labelled 1 when it is the start of the same gene in the NCBI record (the record found again in
   the NCBI sequence, either strand, any rotation);
3. cross-validation by family (GroupKFold: a family's records are never in training and test at once): exact starts
   of P1 alone and of P1 + the model, per phylum;
4. the final model on every record, written as JSON trees (mito_start.StartModel evaluates them; no pickle).

Usage: python3 train_start_model.py --set sets/A.tsv:<root A> [--set sets/B.tsv:<root B>] --out models/start_model
       [--jobs 36]
"""
import argparse
import csv
import json
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import date
from pathlib import Path

import numpy as np
from Bio import SeqIO

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import mito_start  # noqa: E402
import mito_taxa  # noqa: E402
from mito_names import gene_name  # noqa: E402

EXTRA_CLADE = {"Cnidaria": "Cnidaria", "Platyhelminthes": "Platyhelminthes-flatworms",
               "Chaetognatha": "Chaetognatha-arrow-worms"}


def dump(job):
    """Candidates of one record (mito_curate --dump_starts), each with the NCBI start label."""
    root, row, tmp = job
    acc, code = row["accession"], row["genetic_code"]
    pre = root / "runs" / "mitom_lo" / acc / "raw" / "work" / "04.final_gb" / f"{acc}.gbf"
    if not pre.is_file():
        return []
    clade = row["mitoz_clade"] or EXTRA_CLADE.get(row.get("phylum", ""), "Animal")
    out = Path(tmp) / acc
    out.mkdir(parents=True, exist_ok=True)
    (out / "in").mkdir(exist_ok=True)
    src = out / "in" / f"{acc}.gbf"
    if not src.exists():
        src.symlink_to(pre)
    p = subprocess.run([sys.executable, "-W", "ignore", str(HERE.parent / "mito_curate.py"), "-i", str(out / "in"),
                        "-o", str(out / "cur"), "--genetic_code", code, "--no_validate", "--start_model", "no",
                        "--dump_starts", str(out / "starts.tsv"), "--phylum", row.get("phylum", ""),
                        "--aggressive_start_genes", mito_taxa.aggressive_start_genes(clade),
                        *mito_taxa.abbreviated_overlap(clade), *mito_taxa.no_upstream_atg(clade),
                        *mito_taxa.abbreviated_pcg_overlap(clade)], capture_output=True, text=True)
    if p.returncode or not (out / "starts.tsv").is_file():
        return []
    rows = list(csv.DictReader(open(out / "starts.tsv"), delimiter="\t"))
    if not rows:
        return []
    ncbi = SeqIO.read(str(root / "gb" / f"{acc}.gb"), "genbank")
    mine = SeqIO.read(str(pre), "genbank")
    n, seq = len(ncbi.seq), str(ncbi.seq).upper()
    own = str(mine.seq).upper()
    if len(own) != n:
        return []
    k = (seq + seq).find(own)                      # NCBI position of the record's first base (same strand) ...
    rc = k < 0
    if rc:                                         # ... or of its last base on the other strand
        k = (seq + seq).find(str(mine.seq.reverse_complement()).upper())
        if k < 0:
            return []
    starts = defaultdict(set)                      # (gene, strand in the NCBI record) -> 5' ends (1-based)
    for f in ncbi.features:
        if f.type == "CDS":
            g = gene_name(f.qualifiers, "CDS")
            five = int(f.location.parts[0].start) + 1 if f.location.strand == 1 else int(f.location.parts[0].end)
            starts[(g, f.location.strand)].add(five)
    family = next((t for t in ncbi.annotations.get("taxonomy", []) if t.endswith("idae")), row.get("family") or row["phylum"])
    for r in rows:
        pos, strand = int(r["genome_pos"]), 1 if r["strand"] == "+" else -1
        npos = (k + pos - 1) % n + 1 if not rc else (k + n - pos) % n + 1
        nstrand = strand if not rc else -strand
        r["label"] = int(npos in starts.get((gene_name({"gene": [r["gene"]]}, "CDS"), nstrand), ()))
        r["gid"] = f"{root.name}:{acc}:{r['gene']}:{r['strand']}:{int(r['pos']) - int(r['cand'])}"
        r["family"], r["phylum_set"] = family, row["phylum"]
    return rows


def export(model, cols, categories, meta) -> dict:
    """HistGradientBoostingClassifier -> {columns, categories, trees: [{left, right, feature, threshold, value}]};
    leaves have left = -1. (Its trees are scikit-learn internals: read here, at training, only.)"""
    trees = []
    for (pred,) in model._predictors:
        nodes = pred.nodes
        leaf = nodes["is_leaf"].astype(bool)
        trees.append({"left": np.where(leaf, -1, nodes["left"].astype(np.int64)).tolist(),   # (uint32 in the nodes)
                      "right": np.where(leaf, -1, nodes["right"].astype(np.int64)).tolist(),
                      "feature": nodes["feature_idx"].astype(int).tolist(),
                      "threshold": [float(x) for x in nodes["num_threshold"]],
                      "value": [float(x) for x in nodes["value"]]})
    return {**meta, "columns": cols, "categories": categories, "trees": trees}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--set", action="append", required=True, help="benchmark set TSV:root folder (repeatable)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--jobs", type=int, default=36)
    a = ap.parse_args()
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.model_selection import GroupKFold
    import sklearn

    jobs, seen = [], set()
    with tempfile.TemporaryDirectory(prefix="starts.") as tmp:
        for spec in a.set:
            tsv, root = spec.split(":")
            for row in csv.DictReader(open(tsv), delimiter="\t"):
                if row["accession"] not in seen:            # a record in two sets counts once
                    seen.add(row["accession"])
                    jobs.append((Path(root), row, tmp))
        rows = []
        with ProcessPoolExecutor(a.jobs) as ex:
            for rs in ex.map(dump, jobs, chunksize=4):
                rows += rs
    genes = {r["gid"] for r in rows}
    print(f"{len(jobs)} records, {len(genes)} genes, {len(rows)} candidates", file=sys.stderr)
    categories = {c: sorted({r[c] for r in rows}) for c in mito_start.CATEGORICAL}
    cols = mito_start.columns(categories)
    x = mito_start.matrix(rows, categories)
    y = np.array([int(r["label"]) for r in rows])
    gid = np.array([r["gid"] for r in rows])
    fam = np.array([r["family"] for r in rows])
    params = dict(max_iter=300, learning_rate=0.08, max_leaf_nodes=31, random_state=0)

    def exact(scores, mask=None):
        best = {}
        for i in range(len(rows)):
            if mask is not None and not mask[i]:
                continue
            g = gid[i]
            if g not in best or (scores[i], rows[i]["is_p1"]) > (scores[best[g]], rows[best[g]]["is_p1"]):
                best[g] = i
        return sum(y[i] for i in best.values()), len(best)

    cv = np.zeros(len(rows))
    for tr, te in GroupKFold(n_splits=5).split(x, y, groups=fam):
        m = HistGradientBoostingClassifier(**params).fit(x[tr], y[tr])
        cv[te] = m.decision_function(x[te])
    p1 = np.array([float(r["is_p1"]) for r in rows])
    report = {}
    for ph in ["all"] + sorted({r["phylum_set"] for r in rows}):
        mask = None if ph == "all" else np.array([r["phylum_set"] == ph for r in rows])
        (k0, n0), (k1, _) = exact(p1, mask), exact(cv, mask)
        report[ph] = {"genes": int(n0), "p1_exact": int(k0), "model_exact": int(k1)}
        if ph == "all" or n0 >= 100:
            print(f"{ph:16s} genes {n0:6d}  exact starts P1 {k0 / n0:.3f}  P1 + model {k1 / n0:.3f}", file=sys.stderr)
    final = HistGradientBoostingClassifier(**params).fit(x, y)
    meta = {"model": "HistGradientBoostingClassifier " + json.dumps(params), "scikit_learn": sklearn.__version__,
            "trained": date.today().isoformat(), "records": len(jobs), "genes": len(genes), "candidates": len(rows),
            "training_sets": [s.split(":")[0].split("/")[-1] for s in a.set],
            "cross_validation_by_family": report}
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "start_model.json").write_text(json.dumps(export(final, cols, categories, meta)) + "\n")
    # check: the exported trees rank like the trained model
    sm = mito_start.load(str(out / "start_model.json"))
    sample = rows[:2000]
    a1 = np.argsort(sm.score(sample))
    a2 = np.argsort(final.decision_function(mito_start.matrix(sample, categories)))
    print(f"export check: same order {np.array_equal(a1, a2)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
