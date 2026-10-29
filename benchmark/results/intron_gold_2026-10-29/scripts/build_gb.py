"""Write gold/gb/<acc>.gb: each NCBI record with its introns at the gold placement (criteria: gold/README.md).

Per NCBI intron (classify.py, gold/work/gold.json):
  kept (A, B, C, D-corrected, E)  the host CDS / rRNA join is moved to the gold exon boundaries (/translation redone
                                  when the protein changes) and one intron feature is written at the gold placement,
                                  /note "gold_class=<class>; NCBI <coordinates>; ..." (the NCBI intron feature is
                                  replaced);
  D, excluded (gold "excluded")  intron of unknown placement: the NCBI intron feature is removed and a misc_feature
                                  /note "eval_introns: ignore ..." covers the NCBI placement (benchmark/eval_introns.py
                                  then scores neither truth nor predictions there);
  D, removed (gold "removed")    an NCBI intron feature that is no intron (inside an exon / part of a join gap that is
                                  evaluated on its own): the feature is removed.
Check: benchmark/eval_introns.py's truth of every written record equals the gold introns.
usage: build_gb.py   -> gold/gb/<acc>.gb
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify as V  # noqa: E402
from Bio import SeqIO  # noqa: E402
from Bio.Data import CodonTable  # noqa: E402
from Bio.Seq import Seq  # noqa: E402
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature  # noqa: E402

GOLD = Path(__file__).resolve().parents[1]
E = V.E


def loc(s: int, ln: int, strand: int, n: int):
    """Location of an interval (start mod n, length) that may cross the origin."""
    s %= n
    if s + ln <= n:
        return FeatureLocation(s, s + ln, strand=strand)
    a, b = FeatureLocation(s, n, strand=strand), FeatureLocation(0, s + ln - n, strand=strand)
    return CompoundLocation([a, b] if strand == 1 else [b, a])


def same_interval(f, t: E.Intron, n: int) -> bool:
    parts = list(f.location.parts)
    s, e = int(f.location.start), int(f.location.end)
    strand = f.location.strand or 1
    if any((int(q.start) < int(p.start)) == (strand == 1) for p, q in zip(parts, parts[1:])):
        s, e = (int(parts[0].start), int(parts[-1].end) + n) if strand == 1 else (int(parts[-1].start), int(parts[0].end) + n)
    return (s, e, strand) == (t.start, t.end, t.strand)


def move_join(rec, f, k: int, gold_s: int, gold_len: int, code: int) -> bool:
    """Move the exon boundaries around gap k of feature f so that the gap is the gold placement; True if changed."""
    m = V.model_of(rec, f, code)
    if m.intron(k)[:2] == (gold_s, gold_len):
        return False
    for d5 in range(-62, 63):
        for d3 in range(-62, 63):
            ex = m.moved(k, d5, d3)
            if ex is None or m.intron(k, ex)[:2] != (gold_s, gold_len):
                continue
            parts = list(f.location.parts)
            new = []
            for i, (p, (s, ln, st)) in enumerate(zip(parts, ex)):
                if i in (k, k + 1):
                    q = loc(s, ln, st, m.n)
                    new += list(q.parts)
                else:
                    new.append(p)
            f.location = CompoundLocation(new) if len(new) > 1 else new[0]
            if f.type == "CDS" and "translation" in f.qualifiers:
                table = CodonTable.unambiguous_dna_by_id[m.code]
                nt = str(f.location.extract(rec.seq)).upper()[m.cs - 1:]
                prot = str(Seq(nt[:len(nt) - len(nt) % 3]).translate(table=m.code))
                prot = prot[:-1] if prot.endswith("*") else prot
                if nt[:3] in table.start_codons and m.cs == 1:
                    prot = "M" + prot[1:]
                f.qualifiers["translation"] = [prot.replace("*", "X")]
            return True
    raise ValueError(f"{rec.id}: no placement of gap {k} gives {gold_s}+{gold_len}")


def main() -> int:
    gold = json.load(open(GOLD / "work" / "gold.json"))
    by = defaultdict(list)
    for g in gold:
        by[g["accession"]].append(g)
    out = GOLD / "gb"
    out.mkdir(exist_ok=True)
    bad = 0
    for acc in sorted(p.stem for p in (GOLD.parent / "ncbi_gb").glob("*.gb")):
        rec = SeqIO.read(str(GOLD.parent / "ncbi_gb" / f"{acc}.gb"), "genbank")
        n = len(rec)
        truth = sorted(E.introns_of(rec), key=lambda x: x.start)
        rows = by.get(acc, [])
        assert len(rows) == len(truth), (acc, len(rows), len(truth))
        drop, add, expect, ignore = set(), [], set(), []
        for t, g in zip(truth, rows):
            for i, f in enumerate(rec.features):
                if f.type == "intron" and same_interval(f, t, n):
                    drop.add(i)
            gene_q = next((f.qualifiers["gene"] for f in rec.features if f.type == "intron" and same_interval(f, t, n)
                           and f.qualifiers.get("gene")), None)
            handling = g["gold_handling"]
            if handling == "drop":
                continue
            if handling == "ignore":
                s, ln = (g["gap_s"], g["gap_len"]) if g.get("gap_s") is not None else (t.start % n, t.end - t.start)
                ignore.append(SeqFeature(loc(s, ln, t.strand, n), type="misc_feature", qualifiers={
                    "note": [f"eval_introns: ignore - gold class D, NCBI intron {g['ncbi']}: {g['reason'][:300]}"]}))
                continue
            f, k = V.host_of(rec, t)
            if f is not None and g.get("gold_s") is not None and (g.get("gap_s") is not None):
                move_join(rec, f, k, g["gold_s"], g["gold_len"], code_of(acc))
                gene_q = gene_q or f.qualifiers.get("gene")
            gs, gl = g["gold_s"], g["gold_len"]
            cls = g["class"] + (f" ({g['verdict']})" if g.get("verdict") else "")
            note = (f"gold_class={g['class']}; verdict {cls}; NCBI {g['ncbi']} ({g['ncbi_source']}); "
                    f"{g['type']}; {g['reason'][:400]}")
            add.append(SeqFeature(loc(gs, gl, t.strand, n), type="intron",
                                  qualifiers={**({"gene": gene_q} if gene_q else {}), "note": [note]}))
            expect.add((gs % n, gl, t.strand))
        rec.features = [f for i, f in enumerate(rec.features) if i not in drop] + add + ignore
        rec.features.sort(key=lambda f: (int(f.location.start), f.type != "source", f.type != "gene"))
        SeqIO.write(rec, str(out / f"{acc}.gb"), "genbank")
        back = SeqIO.read(str(out / f"{acc}.gb"), "genbank")
        regions = E.ignore_regions(back)
        got = {(x.start % n, x.end - x.start, x.strand) for x in E.introns_of(back)
               if not any(st == x.strand and E._overlap((x.start, x.end), (a, b), n) >= 0.5 * (x.end - x.start)
                          for a, b, st in regions)}
        if got != expect:
            bad += 1
            print(f"MISMATCH {acc}: truth {sorted(got)} vs gold {sorted(expect)}", file=sys.stderr)
    print(f"{len(by)} records with introns -> {out} ({bad} mismatches)", file=sys.stderr)
    return 1 if bad else 0


_CODE = {}


def code_of(acc: str) -> int:
    if not _CODE:
        import csv
        for r in csv.DictReader(open(GOLD.parent / "records.tsv"), delimiter="\t"):
            _CODE[r["accession"]] = int(r["genetic_code"])
    return _CODE[acc]


if __name__ == "__main__":
    sys.exit(main())
