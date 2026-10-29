"""Evaluate intron annotation against the NCBI records.

Truth (per NCBI record): every `intron` feature (a multi-part one as its span), plus the introns implied by
joined CDS - consecutive parts on one strand separated by >= 142 nt (MFannot's shortest intron; measured across
the origin when the join crosses it, which is 0 nt for a CDS merely cut by the origin) - that no intron feature
already gives (one covering >= 50 % of the gap; an intron nested in it, like a hexacoral cox1 intron inside the nad5
intron, does not count) - many records give an intron only as a CDS join; mito_bench_eval.py counts intron features
only.
Parts on different strands (trans-splicing) and the 1-nt frameshift joins of bird / turtle nad3 are not introns.
Predictions: the same from the predicted records, with the intron type from the note MitoM writes
(group I / group II / type not identified).

A prediction matches a truth intron on the same strand when each covers >= 50 % of the other (so an intron nested in
another one matches neither); pairs are assigned one to one, highest Jaccard index (overlap / union) first, then the
smaller end difference. Classes: exact, le6 (both ends within 6 nt), gt6, missed; unmatched predictions: fp.

Gold-standard records (mitom_work/intron_gold_2026-10-29/gold/gb): an intron feature whose /note has
"gold_class=X" gives the truth row its class (column truth_class; the summary adds one row per class), and a
misc_feature whose /note starts with "eval_introns: ignore" marks an intron whose placement is unknown: truth introns
lying >= 50 % inside it and predictions lying >= 50 % inside it are not scored (rows of class "ignored").

Predicted records: <pred>/<acc>/final.gbf (benchmark runs: MitoM's final record in NCBI coordinates), else
<pred>/<acc>/raw/final_annotation/<acc>.gb or <pred>/<acc>.gb[f] (any GenBank file in its own orientation); a record
is mapped onto the NCBI sequence (rotation, reverse complement) before comparing. When
<pred>/<acc>/raw/reports/intron_report.tsv exists, every missed intron and fp gets the mito_intron.py candidate at
that place (confidence, written, evidence) or, if there is none, the report's notes for that gene.

Usage:
  python3 benchmark/eval_introns.py --gb mitom_work/bench_large/gb --pred mitom_work/bench_large/runs/mitom_native \
      --label bench_large -o eval_introns/ [--acc NC_000933.1,NC_027590.1]
Outputs: <out>/<label>_introns.tsv (one row per truth intron and per fp), <out>/<label>_summary.tsv (per taxon
group and in total); the summary is also printed.
"""
import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Callable, Dict, List, NamedTuple, Optional

from Bio import SeqIO
from Bio.Seq import Seq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from mito_names import gene_name  # noqa: E402

MIN_INTRON = 142
MIN_OVERLAP = 0.5
NEAR = 6
PHYLA = {"Porifera", "Placozoa", "Cnidaria", "Ctenophora", "Annelida", "Mollusca", "Arthropoda", "Chordata",
         "Echinodermata", "Hemichordata", "Nematoda", "Platyhelminthes", "Brachiopoda", "Bryozoa", "Nemertea",
         "Phoronida", "Rotifera", "Acanthocephala", "Chaetognatha", "Tardigrada", "Onychophora", "Priapulida",
         "Kinorhyncha", "Nematomorpha", "Entoprocta", "Cycliophora", "Gastrotricha", "Xenacoelomorpha"}


class Intron(NamedTuple):
    start: int         # 0-based; end exclusive and > n when the intron crosses the origin
    end: int
    strand: int
    gene: str
    source: str        # intron_feature / cds_join
    note: str


def host_gene(rec, s: int, e: int, strand: int) -> str:
    """Normalised name of the gene that contains an intron: joined CDS / rRNA / gene first, then any."""
    for joined in (True, False):
        for kind in ("CDS", "rRNA", "gene"):
            for f in rec.features:
                if f.type == kind and (f.location.strand or 1) == strand and int(f.location.start) <= s \
                        and e <= int(f.location.end) and (len(f.location.parts) > 1) == joined:
                    g = gene_name(f.qualifiers, "rRNA" if kind == "rRNA" else "CDS")
                    if g:
                        return g
    return ""


def introns_of(rec) -> List[Intron]:
    """Intron features plus the introns implied by joined CDS (see module docstring)."""
    n = len(rec)
    feats: List[Intron] = []
    for f in rec.features:
        if f.type != "intron":
            continue
        strand = f.location.strand or 1
        parts = list(f.location.parts)
        s, e = int(f.location.start), int(f.location.end)   # a multi-part intron (NC_022825.1) counts as its span
        if any((int(q.start) < int(p.start)) == (strand == 1) for p, q in zip(parts, parts[1:])):
            s, e = (int(parts[0].start), int(parts[-1].end) + n) if strand == 1 else \
                (int(parts[-1].start), int(parts[0].end) + n)  # written across the origin
        g = gene_name({"gene": f.qualifiers.get("gene", [""])}, "CDS") or host_gene(rec, s, e, strand)
        note = " ".join(f.qualifiers.get("note", []) + f.qualifiers.get("product", []))
        if not any((x.start, x.end, x.strand) == (s, e, strand) for x in feats):
            feats.append(Intron(s, e, strand, g, "intron_feature", note))
    implied: List[Intron] = []
    for f in rec.features:
        if f.type != "CDS" or len(f.location.parts) < 2:
            continue
        g = gene_name(f.qualifiers, "CDS") or (f.qualifiers.get("gene", ["?"])[0])
        for p, q in zip(f.location.parts, f.location.parts[1:]):      # biological order
            if p.strand != q.strand:
                continue                                    # trans-spliced
            plus = (p.strand or 1) == 1
            s, e = (int(p.end), int(q.start)) if plus else (int(q.end), int(p.start))
            if (int(q.start) < int(p.start)) == plus:
                e += n                                      # gap across the origin (0 nt for a CDS cut by it)
            # the same intron given as a feature; an intron nested in this one (corallimorpharian cox1 inside the
            # nad5 intron) is not
            if e - s < MIN_INTRON or any(_overlap((s, e), (x.start, x.end), n) >= MIN_OVERLAP * (e - s)
                                         for x in feats if x.strand == (p.strand or 1)):
                continue
            implied.append(Intron(s, e, p.strand or 1, g, "cds_join", ""))
    return feats + implied


IGNORE = "eval_introns: ignore"


def ignore_regions(rec) -> List[tuple]:
    """(start, end, strand) of the misc_features that mark an intron of unknown placement (gold standard)."""
    out = []
    for f in rec.features:
        if f.type == "misc_feature" and any(x.startswith(IGNORE) for x in f.qualifiers.get("note", [])):
            out.append((int(f.location.start), int(f.location.end), f.location.strand or 1))
    return out


def gold_class(note: str) -> str:
    m = re.search(r"gold_class=([A-Za-z-]+)", note)
    return m.group(1) if m else ""


def _overlap(a, b, n: int) -> int:
    """Overlap of two intervals (0-based, end exclusive, end may exceed n) on a circle of length n."""
    return max(max(0, min(a[1], b[1] + k) - max(a[0], b[0] + k)) for k in (-n, 0, n))


def _jaccard(t: Intron, p: Intron, n: int) -> float:
    ov = _overlap((t.start, t.end), (p.start, p.end), n)
    return ov / ((t.end - t.start) + (p.end - p.start) - ov)


def _deltas(t: Intron, p: Intron, n: int):
    """(start, end) difference prediction - truth, at the circular shift that puts them closest."""
    k = min((-n, 0, n), key=lambda k: abs(p.start + k - t.start) + abs(p.end + k - t.end))
    return p.start + k - t.start, p.end + k - t.end


def mapper(ref: str, seq: str) -> Optional[Callable]:
    """Map 0-based intervals of seq (a rotation of ref or of its reverse complement) onto ref."""
    n = len(ref)
    if len(seq) != n or n < 60:
        return None
    rc = str(Seq(ref).reverse_complement())
    for orient, r in ((1, ref), (-1, rc)):
        i = (r + r).find(seq[:60])
        while i != -1 and i < n:
            if (r + r)[i:i + n] == seq:
                def f(s: int, e: int, strand: int, i=i, orient=orient):
                    ln = e - s
                    if orient == 1:
                        a = (s + i) % n
                        return a, a + ln, strand
                    a = (n - ((s + i) % n) - ln) % n
                    return a, a + ln, -strand
                return f
            i = (r + r).find(seq[:60], i + 1)
    return None


def type_label(note: str) -> str:
    m = re.search(r"group (II|I) intron(?:, subtype (\S+))?(?: \(Infernal/Rfam (\S+?),)?", note)
    if m:
        return f"group {m.group(1)} ({m.group(2) or m.group(3) or '?'})"
    return "type not identified" if "type not identified" in note else ""


def pred_file(pred: Path, acc: str) -> Optional[Path]:
    for p in (pred / acc / "final.gbf", pred / acc / "raw" / "final_annotation" / f"{acc}.gb",
              pred / f"{acc}.gb", pred / f"{acc}.gbf"):
        if p.exists():
            return p
    return None


def report_rows(pred: Path, acc: str, ref: str) -> List[dict]:
    """intron_report.tsv rows of a MitoM run, intron coordinates mapped to the NCBI record (keys s, e, strand)."""
    rep = pred / acc / "raw" / "reports" / "intron_report.tsv"
    own = pred / acc / "raw" / "final_annotation" / f"{acc}.gb"
    if not rep.exists():
        return []
    fmap = mapper(ref, str(SeqIO.read(str(own), "genbank").seq).upper()) if own.exists() else None
    rows = list(csv.DictReader(rep.open(encoding="utf-8"), delimiter="\t"))
    for r in rows:
        r["g"] = gene_name({"gene": [r.get("gene") or ""]}, "CDS")
        if r.get("intron_start") and fmap:
            strand = 1 if r["strand"] == "+" else -1
            r["s"], r["e"], r["strand_n"] = fmap(int(r["intron_start"]) - 1, int(r["intron_end"]), strand)
    return rows


def explain(rows: List[dict], x: Intron, n: int) -> str:
    """The report's candidate at an intron, else its notes for the intron's gene."""
    hits = [r for r in rows if "s" in r and r["strand_n"] == x.strand and _overlap((x.start, x.end), (r["s"], r["e"]), n)]
    if hits:
        r = max(hits, key=lambda r: _overlap((x.start, x.end), (r["s"], r["e"]), n))
        return (f"candidate {r['s'] + 1}..{r['e']} ({r['gene']}, {r['confidence']}, written={r['written']}): "
                f"{r['evidence']} | {r['note']}")
    notes = sorted({f"{r['gene']}: {r['note']}" for r in rows if "s" not in r and r.get("note")
                    and (not x.gene or r["g"] == x.gene or not r["g"])})
    return "; ".join(notes) or "no candidate in intron_report.tsv"


def group_of(rec) -> str:
    lineage = rec.annotations.get("taxonomy", [])
    return next((t for t in lineage if t in PHYLA), lineage[2] if len(lineage) > 2 else "?")


def evaluate(acc: str, gb: Path, pred: Path, label: str) -> List[dict]:
    ncbi_file = gb / f"{acc}.gb"
    pf = pred_file(pred, acc)
    texts = [ncbi_file.read_text(errors="ignore")] + ([pf.read_text(errors="ignore")] if pf else [])
    if not any("intron" in t or "join(" in t for t in texts):
        return []                                           # no intron and no joined feature anywhere
    ncbi = SeqIO.read(str(ncbi_file), "genbank")
    n, ref = len(ncbi), str(ncbi.seq).upper()
    ignore = ignore_regions(ncbi)

    def ignored(x: Intron) -> bool:
        return any(st == x.strand and _overlap((x.start, x.end), (a, b), n) >= MIN_OVERLAP * (x.end - x.start)
                   for a, b, st in ignore)
    truth = [t for t in introns_of(ncbi) if not ignored(t)]
    preds: List[Intron] = []
    status = "no prediction file"
    if pf:
        prec = SeqIO.read(str(pf), "genbank")
        fmap = mapper(ref, str(prec.seq).upper())
        status = "ok" if fmap else "prediction sequence differs from NCBI"
        if fmap:
            for x in introns_of(prec):
                s, e, strand = fmap(x.start, x.end, x.strand)
                preds.append(x._replace(start=s, end=e, strand=strand))
    rep = report_rows(pred, acc, ref) if pf else []
    base = {"set": label, "accession": acc, "group": group_of(ncbi), "organism": ncbi.annotations.get("organism", "")}
    pairs = sorted(((ti, pi) for ti, t in enumerate(truth) for pi, p in enumerate(preds) if p.strand == t.strand
                    and _overlap((t.start, t.end), (p.start, p.end), n) >= MIN_OVERLAP * max(t.end - t.start, p.end - p.start)),
                   key=lambda x: (-_jaccard(truth[x[0]], preds[x[1]], n),
                                  max(map(abs, _deltas(truth[x[0]], preds[x[1]], n)))))
    match: Dict[int, int] = {}
    for ti, pi in pairs:
        if ti not in match and pi not in match.values():
            match[ti] = pi
    out = []
    for ti, t in sorted(enumerate(truth), key=lambda x: x[1].start):
        row = {**base, "class": "missed", "gene": t.gene, "strand": "+" if t.strand == 1 else "-",
               "truth_class": gold_class(t.note),
               "truth": f"{t.start + 1}..{t.end}", "truth_len": t.end - t.start, "truth_source": t.source,
               "truth_note": t.note[:80], "pred": "", "pred_len": "", "d_start": "", "d_end": "", "pred_type": "",
               "pred_gene": "", "report": ""}
        if ti in match:
            p = preds[match[ti]]
            ds, de = _deltas(t, p, n)
            d = max(abs(ds), abs(de))
            row.update({"class": "exact" if d == 0 else "le6" if d <= NEAR else "gt6", "pred": f"{p.start + 1}..{p.end}",
                        "pred_len": p.end - p.start, "d_start": ds, "d_end": de, "pred_type": type_label(p.note),
                        "pred_gene": p.gene})
        else:
            row["report"] = explain(rep, t, n) if status == "ok" else status
        out.append(row)
    for i, p in enumerate(preds):
        if i not in match.values():
            out.append({**base, "class": "ignored" if ignored(p) else "fp", "gene": "",
                        "strand": "+" if p.strand == 1 else "-", "truth_class": "", "truth": "",
                        "truth_len": "", "truth_source": "", "truth_note": "", "pred": f"{p.start + 1}..{p.end}",
                        "pred_len": p.end - p.start, "d_start": "", "d_end": "", "pred_type": type_label(p.note),
                        "pred_gene": p.gene, "report": explain(rep, p, n)})
    return out


def summarise(rows: List[dict], label: str) -> List[dict]:
    by: Dict[str, Counter] = defaultdict(Counter)
    for r in rows:
        for g in (r["group"], "all") + ((f"gold_class {r['truth_class']}",) if r.get("truth_class") else ()):
            by[g][r["class"]] += 1
    out = []
    for g, c in sorted(by.items(), key=lambda kv: (kv[0].startswith("gold_class"), kv[0] == "all", kv[0])):
        found = c["exact"] + c["le6"] + c["gt6"]
        n_truth, n_pred = found + c["missed"], found + c["fp"]
        if g.startswith("gold_class"):               # truth rows only: no fp / precision per class
            n_pred = ""
        out.append({"set": label, "group": g, "n_truth": n_truth, "exact": c["exact"], "le6": c["le6"], "gt6": c["gt6"],
                    "missed": c["missed"], "recall": f"{found / n_truth:.3f}" if n_truth else "", "n_pred": n_pred,
                    "fp": c["fp"] if n_pred != "" else "", "precision": f"{found / n_pred:.3f}" if n_pred else "",
                    **({"ignored": c["ignored"]} if c["ignored"] else {})})
    return out


def write_tsv(path: Path, rows: List[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, list(dict.fromkeys(k for r in rows for k in r)) if rows else ["set"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gb", required=True, help="Folder of NCBI records <acc>.gb (the truth)")
    p.add_argument("--pred", required=True, help="Folder of predictions: <acc>/final.gbf, a MitoM output or <acc>.gb")
    p.add_argument("--label", default="set", help="Set name used in the output file names and rows")
    p.add_argument("--acc", default="", help="Comma-separated accessions (default: every <acc>.gb in --gb)")
    p.add_argument("-o", "--out", required=True, help="Output folder")
    a = p.parse_args(argv)
    gb, pred, out = Path(a.gb), Path(a.pred), Path(a.out)
    accs = a.acc.split(",") if a.acc else sorted(f.name[:-3] for f in gb.glob("*.gb"))
    rows = [r for acc in accs for r in evaluate(acc, gb, pred, a.label)]
    out.mkdir(parents=True, exist_ok=True)
    summ = summarise(rows, a.label)
    write_tsv(out / f"{a.label}_introns.tsv", rows)
    write_tsv(out / f"{a.label}_summary.tsv", summ)
    cols = ["group", "n_truth", "exact", "le6", "gt6", "missed", "recall", "n_pred", "fp", "precision"]
    print(f"[eval_introns] {a.label}: {len(accs)} records")
    print("\t".join(cols))
    for r in summ:
        print("\t".join(str(r[k]) for k in cols))
    return 0


if __name__ == "__main__":
    sys.exit(main())
