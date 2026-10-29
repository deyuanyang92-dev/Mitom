#!/usr/bin/env python3
"""mito_compare.py - compare the original annotation of each record (e.g. as downloaded from NCBI) with the
MitoM annotation of the same sequence, gene by gene.

MitoM rotates / reverse-complements a record so that cox1 comes first; the MitoM genes are mapped back onto the
original sequence first, so both annotations are compared in the original coordinates.

One row per gene (protein genes, rRNAs, tRNAs, introns):
  result  same            same name, strand and both ends
          ends_differ     same gene, at least one end differs (difference_nt = larger end shift)
          name_differs    same place, different gene name
          strand_differs  same place, other strand
          only_original   annotated only in the original record
          only_mitom      annotated only by MitoM
          same_place      tRNA-Leu / tRNA-Ser at the same place, which of the two genes cannot be told (no anticodon)

Usage:  mitom compare -o original_vs_mitom.tsv --original original_annotation/ --mitom final_annotation/ [--genetic_code 5]
        (files are paired by name: original_annotation/<acc>.gb <-> final_annotation/<acc>.gb)
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional

from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, ExactPosition, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mito_names import PCG, gene_name, trna_name, trna_number  # noqa: E402

GB_EXT = {".gb", ".gbk", ".gbf", ".genbank", ".gbff"}
RESULT = {"exact": "same", "le3": "ends_differ", "gt3": "ends_differ", "name_error": "name_differs",
          "strand_error": "strand_differs", "missed": "only_original", "fp": "only_mitom",
          "name_unverifiable": "same_place"}
COLUMNS = ["record", "organism", "type", "gene", "original_location", "mitom_location", "difference_nt", "result"]
KIND = {"PCG": "protein", "rRNA": "rRNA", "tRNA": "tRNA", "intron": "intron"}


# ── the MitoM record on the original sequence (rotation / reverse complement undone) ──────────────────
ANCHOR = 50


class LiftError(Exception):
    pass


def _merge(parts: List[FeatureLocation]) -> List[FeatureLocation]:
    """Join consecutive parts that touch (an origin split undone by rotation)."""
    out: List[FeatureLocation] = []
    for p in parts:
        if out and out[-1].strand == p.strand:
            q = out[-1]
            if p.strand == -1 and int(p.end) == int(q.start):
                out[-1] = FeatureLocation(p.start, q.end, -1)
                continue
            if p.strand != -1 and int(q.end) == int(p.start):
                out[-1] = FeatureLocation(q.start, p.end, p.strand)
                continue
        out.append(p)
    return out


def _shift(part: FeatureLocation, off: int, n: int) -> List[FeatureLocation]:
    """Move one part by off on a circle of length n; split at the origin (biological order kept). Partial markers
    (<start, >end) stay on their position: table2asn accepts a 5'- / 3'-partial CDS without start / stop codon."""
    s = (int(part.start) + off) % n
    e = s + int(part.end) - int(part.start)
    sp = BeforePosition if isinstance(part.start, BeforePosition) else ExactPosition
    ep = AfterPosition if isinstance(part.end, AfterPosition) else ExactPosition
    if e <= n:
        return [FeatureLocation(sp(s), ep(e), part.strand)]
    a, b = FeatureLocation(sp(s), n, part.strand), FeatureLocation(0, ep(e - n), part.strand)
    return [b, a] if part.strand == -1 else [a, b]


def lift_to_ncbi(pred: SeqRecord, ncbi_seq: Seq) -> SeqRecord:
    """Return pred re-expressed on ncbi_seq (same molecule, possibly rotated and/or reverse-complemented)."""
    n = len(ncbi_seq)
    pseq = Seq(str(pred.seq))
    if len(pseq) != n:
        raise LiftError(f"length {len(pseq)} != NCBI {n}")
    ring = str(ncbi_seq).upper() * 2
    for rc in (False, True):
        q = str(pseq.reverse_complement() if rc else pseq).upper()
        hits = [m.start() for m in re.finditer(f"(?={re.escape(q[:ANCHOR])})", ring[:n + ANCHOR - 1])]
        hits = [h for h in hits if ring[h:h + n] == q]
        if len(hits) == 1:
            off = hits[0]
            break
    else:
        raise LiftError("no unique match of the predicted sequence on the NCBI sequence")
    out = SeqRecord(ncbi_seq, id=pred.id, name=pred.name, description=pred.description)
    out.annotations = dict(pred.annotations)
    for f in pred.features:
        parts = list(f.location.parts)
        if rc:                                           # a <start becomes a >end and vice versa
            parts = [FeatureLocation(BeforePosition(n - int(p.end)) if isinstance(p.end, AfterPosition) else n - int(p.end),
                                     AfterPosition(n - int(p.start)) if isinstance(p.start, BeforePosition) else n - int(p.start),
                                     -(p.strand or 1)) for p in parts]
        new = _merge([x for p in parts for x in _shift(p, off, n)])
        loc = new[0] if len(new) == 1 else CompoundLocation(new)
        q = dict(f.qualifiers)
        if "transl_except" in q:                   # positions inside the qualifier move with the feature
            q["transl_except"] = [_lift_except(t, rc, off, n) for t in q["transl_except"]]
        out.features.append(SeqFeature(loc, type=f.type, qualifiers=q))
    return out


def _lift_except(text: str, rc: bool, off: int, n: int) -> str:
    """'(pos:X..Y,aa:TERM)' / '(pos:complement(X),aa:TERM)' re-expressed after reverse-complement and rotation."""
    m = re.match(r"\(pos:(complement\()?(\d+)(?:\.\.(\d+))?\)?,(.*)$", text)
    if not m:
        return text
    def move(x: int) -> int:
        i = n - x if rc else x - 1
        return (i + off) % n + 1
    a, b = move(int(m[2])), move(int(m[3] or m[2]))
    lo, hi = (b, a) if rc else (a, b)              # rc reverses the order; lo > hi only when the span wraps
    if lo <= hi:
        pos = f"{lo}..{hi}" if lo != hi else f"{lo}"
    else:                                          # a 2-nt stop split by the origin
        pos = f"join({lo}{'' if lo == n else f'..{n}'},1{'' if hi == 1 else f'..{hi}'})"
    return f"(pos:complement({pos}),{m[4]}" if bool(m[1]) != rc else f"(pos:{pos},{m[4]}"


# ── gene-by-gene matching (also the benchmark's: benchmark/mito_bench_eval.py) ─────────────────────────
MIN_OVERLAP = 0.5
NEAR = 3
AMBIGUOUS_AA = {"L", "S"}              # two tRNA genes each; told apart only by anticodon
# the two Leu / Ser genes by the codon family they read: L1 = CUN, L2 = UUR, S1 = AGN, S2 = UCN
COMP = str.maketrans("acgtu", "tgcaa")


def ls_class(aa: str, triplet: str, number: str, code: int) -> str:
    """'1' / '2' for the two Leu or Ser tRNA genes, '' when unknown. The triplet from the annotation counts first,
    read as an anticodon when its codon is Leu / Ser in `code`, else as the codon itself (NCBI records write both,
    e.g. trnL1(cta) = codon CUA, trnL(tag) = anticodon); the gene number (trnL1, tRNA-Ser2) only when there is no
    triplet, because submitters number the two genes differently (vertebrate records swap L1 / L2)."""
    if aa not in AMBIGUOUS_AA:
        return ""
    fallback = number if number in ("1", "2") else ""
    t = (triplet or "").lower().replace("u", "t")
    if len(t) != 3 or not set(t) <= set("acgt"):
        return fallback
    table = CodonTable.unambiguous_dna_by_id[code].forward_table
    codon = t.translate(COMP)[::-1].upper()           # anticodon -> codon
    if table.get(codon) != aa:
        codon = t.upper()                              # not an anticodon of this amino acid: it is the codon
        if table.get(codon) != aa:
            return fallback
    if aa == "L":
        return "1" if codon.startswith("CT") else "2"
    return "1" if codon.startswith("AG") else "2"


class Feat(NamedTuple):
    kind: str          # PCG, rRNA, tRNA, intron
    name: str          # cox1 / rrnL / amino acid letter / intron
    strand: int
    start: int         # 0-based; end exclusive and > n when the feature crosses the origin
    end: int
    anticodon: str
    n: int = 0         # sequence length (circular arithmetic)


def span(loc, n: int) -> tuple:
    """(start, end) of a location; parts are in biological order, a wrap past the origin gives end > n."""
    parts = list(loc.parts)
    strand = loc.strand or 1
    wrapped = any((int(b.start) < int(a.start)) if strand != -1 else (int(b.start) > int(a.start))
                  for a, b in zip(parts, parts[1:]))
    if not wrapped:
        return int(loc.start), int(loc.end)
    if strand != -1:
        return int(parts[0].start), int(parts[-1].end) + n
    return int(parts[-1].start), int(parts[0].end) + n


def norm_features(rec: SeqRecord, code: int = 5) -> List[Feat]:
    out = []
    n = len(rec)
    for f in rec.features:
        strand = f.location.strand or 1
        s, e = span(f.location, n)
        if f.type in ("CDS", "rRNA"):
            g = gene_name(f.qualifiers, f.type)
            if g:
                out.append(Feat("PCG" if g in PCG else "rRNA", g, strand, s, e, "", n))
        elif f.type == "tRNA":
            aa, ac = trna_name(f.qualifiers)
            if aa:
                out.append(Feat("tRNA", aa, strand, s, e, ls_class(aa, ac, trna_number(f.qualifiers), code), n))
        elif f.type == "intron":
            out.append(Feat("intron", "intron", strand, s, e, "", n))
    return out


def _shifts(a: Feat) -> tuple:
    return (-a.n, 0, a.n) if a.n else (0,)


def _overlap(a: Feat, b: Feat) -> int:
    """Overlap on the circle (b moved by -n, 0, +n)."""
    return max(max(0, min(a.end, b.end + k) - max(a.start, b.start + k)) for k in _shifts(a))


def _delta(a: Feat, b: Feat) -> int:
    return min(max(abs(b.start + k - a.start), abs(b.end + k - a.end)) for k in _shifts(a))


def _origin_stop(t: Feat, p: Feat) -> bool:
    """The same stop codon written two ways at the record's origin: the truth ends with T / TA at the last base (or, on
    the minus strand, the first base) plus /transl_except, the prediction completes the stop 1-2 nt across the origin
    (a tool that rotated the sequence sees the whole TAA). Same protein, same gene end; only the start is compared."""
    if t.kind != "PCG" or not t.n or (t.end - t.start) % 3 == 0:
        return False
    if t.strand > 0:
        return t.end == t.n and any(p.end + k - t.n in (1, 2) for k in _shifts(t))
    return t.start == 0 and any(p.start + k in (-1, -2) for k in _shifts(t))


def _same_name(t: Feat, p: Feat) -> str:
    """'yes', 'no' or 'unverifiable' (tRNA of an ambiguous amino acid without anticodon on one side)."""
    if t.name != p.name:
        return "no"
    if t.kind != "tRNA" or t.name not in AMBIGUOUS_AA:
        return "yes"
    if not t.anticodon or not p.anticodon:            # here: the Leu / Ser gene class '1' / '2' (ls_class)
        return "unverifiable"
    return "yes" if t.anticodon == p.anticodon else "no"


def match(truth: List[Feat], pred: List[Feat]) -> List[dict]:
    rows, used = [], set()
    for t in sorted(truth, key=lambda x: (x.kind, x.start)):
        cands = [(i, p) for i, p in enumerate(pred) if i not in used and p.kind == t.kind
                 and _overlap(t, p) >= MIN_OVERLAP * (t.end - t.start)]
        row: dict = {"kind": t.kind, "name": t.name, "truth": f"{t.start + 1}..{t.end}{'+' if t.strand > 0 else '-'}",
               "pred": "", "delta": "", "class": "missed"}
        if cands:
            def rank(c):
                same = _same_name(t, c[1]) != "no"
                return (not (same and c[1].strand == t.strand), not same, c[1].strand != t.strand, -_overlap(t, c[1]))
            i, p = min(cands, key=rank)
            used.add(i)
            same = _same_name(t, p)
            delta = _delta(t, p)
            if delta and _origin_stop(t, p):                   # only the 5' end can differ
                k = min(_shifts(t), key=lambda k: abs(p.start + k - t.start) + abs(p.end + k - t.end))
                delta = abs(p.start + k - t.start) if t.strand > 0 else abs(p.end + k - t.end)
            row.update(pred=f"{p.start + 1}..{p.end}{'+' if p.strand > 0 else '-'}", delta=delta)
            if same == "no":
                row["class"] = "name_error"
            elif p.strand != t.strand:
                row["class"] = "strand_error"
            elif same == "unverifiable":
                row["class"] = "name_unverifiable"
            else:
                row["class"] = "exact" if delta == 0 else "le3" if delta <= NEAR else "gt3"
        rows.append(row)
    for i, p in enumerate(pred):
        if i not in used:
            rows.append({"kind": p.kind, "name": p.name, "truth": "",
                         "pred": f"{p.start + 1}..{p.end}{'+' if p.strand > 0 else '-'}", "delta": "", "class": "fp"})
    return rows


def _gene_label(kind: str, name: str) -> str:
    return f"trn{name}" if kind == "tRNA" else name


def compare_record(original, mitom, code: int = 5) -> List[dict]:
    """Rows for one record; `mitom` may be rotated / reverse-complemented relative to `original`."""
    organism = str(original.annotations.get("organism", ""))
    lifted = lift_to_ncbi(mitom, Seq(str(original.seq)))
    rows = []
    for m in match(norm_features(original, code), norm_features(lifted, code)):
        rows.append({"record": original.id, "organism": organism, "type": KIND.get(m["kind"], m["kind"]),
                     "gene": _gene_label(m["kind"], m["name"]), "original_location": m["truth"],
                     "mitom_location": m["pred"], "difference_nt": m["delta"], "result": RESULT[m["class"]]})
    return rows


def gb_by_stem(path: Path) -> Dict[str, Path]:
    files = [path] if path.is_file() else [p for p in path.rglob("*") if p.suffix.lower() in GB_EXT]
    return {p.stem: p for p in sorted(files)}


def compare_dirs(original_dir: Path, mitom_dir: Path, out_tsv: Path, code: int = 5) -> Counter:
    """Write out_tsv for every record present in both directories; returns counts per result
    (+ 'records', 'records_without_original_annotation', 'errors')."""
    orig, mit = gb_by_stem(original_dir), gb_by_stem(mitom_dir)
    rows: List[dict] = []
    n = Counter()
    for stem in sorted(set(orig) & set(mit)):
        try:
            o = SeqIO.read(str(orig[stem]), "genbank")
            m = SeqIO.read(str(mit[stem]), "genbank")
            if not any(f.type in ("CDS", "rRNA", "tRNA") for f in o.features):
                n["records"] += 1                      # nothing to compare: every gene would be 'only_mitom'
                n["records_without_original_annotation"] += 1
                continue
            got = compare_record(o, m, code)
        except (LiftError, ValueError) as exc:
            n["errors"] += 1
            rows.append({"record": stem, "organism": "", "type": "", "gene": "", "original_location": "",
                         "mitom_location": "", "difference_nt": "", "result": f"not_compared: {exc}"})
            continue
        n["records"] += 1
        n.update(r["result"] for r in got)
        rows += got
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return n


def summary_text(n: Counter) -> str:
    genes = sum(n[k] for k in set(RESULT.values()))
    parts = [f"{k} {n[k]}" for k in ("same", "ends_differ", "only_original", "only_mitom", "name_differs",
                                       "strand_differs", "same_place") if n[k]]
    return (f"{n['records']} records, {genes} genes: " + ", ".join(parts)
            + (f"; {n['records_without_original_annotation']} records had no original annotation"
               if n["records_without_original_annotation"] else "")
            + (f"; {n['errors']} records not compared" if n["errors"] else ""))


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--original", required=True, help="Original GenBank file or directory.")
    ap.add_argument("--mitom", required=True, help="MitoM GenBank file or directory.")
    ap.add_argument("-o", "--out", default="original_vs_mitom.tsv")
    ap.add_argument("--genetic_code", type=int, default=5)
    a = ap.parse_args(argv)
    n = compare_dirs(Path(a.original), Path(a.mitom), Path(a.out), a.genetic_code)
    print(summary_text(n))
    print(f"table: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
