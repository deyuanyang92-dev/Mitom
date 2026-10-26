"""Build a large order-stratified NC_ benchmark set from a local RefSeq GenBank pool.

The pool is a multi-record GenBank file of RefSeq Metazoa complete mitogenomes (e.g. fetched with
"txid33208[Organism:exp] AND mitochondrion[Filter] AND "complete genome"[Title] AND srcdb_refseq[PROP]").
Ranks come from NCBI Taxonomy (cached in <out>/taxonomy.json). Up to --per_order records per order are
drawn at random (seeded); hybrids, multi-chromosome records and records without cox1 or CDS are skipped.
--all: every usable record of any accession type (e.g. all GenBank records of one phylum); an INSDC record
that an NC_ record of the pool copies is dropped (the NC_ is kept). --accs FILE: exactly the usable records listed
(one accession per line, version optional), e.g. the MitoM-ML tRNA test split.

Writes <out>/gb/<acc>.gb, <out>/bare/<acc>/<acc>.gbf and the record list (mito_sample.py columns).

Usage: python3 make_large_set.py --pool Metazoa.gb --out <root> --tsv sets/mito_bench_large.tsv [--per_order 3 | --all |
       --accs FILE]
"""
import argparse
import csv
import datetime
import json
import random
import re
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

from mito_sample import (COLUMNS, MITOZ_CLADE, annotation_method, evidence_tier, majority_code, split_by_genus,
                         structure_flags, taxonomy, write_bare)


def taxid(rec) -> str:
    for f in rec.features:
        if f.type == "source":
            for x in f.qualifiers.get("db_xref", []):
                if x.startswith("taxon:"):
                    return x[6:]
    return ""


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pool", required=True, type=Path, help="multi-record RefSeq GenBank file")
    p.add_argument("--out", required=True, type=Path, help="data root (gb/, bare/ are created inside)")
    p.add_argument("--tsv", required=True, type=Path)
    p.add_argument("--per_order", type=int, default=3)
    p.add_argument("--all", action="store_true", help="all usable records, any accession type, no per-order cap")
    p.add_argument("--accs", type=Path, help="only these accessions (one per line, version optional)")
    p.add_argument("--seed", type=int, default=20261022)
    a = p.parse_args(argv)
    today = datetime.date.today().isoformat()
    gb_dir, bare_dir = a.out / "gb", a.out / "bare"
    for d in (gb_dir, bare_dir):
        d.mkdir(parents=True, exist_ok=True)

    cand = []                                         # (accession, taxid) of usable records
    copied = set()                                    # INSDC accessions an NC_ record of the pool copies
    wanted = {x.split()[0].split(".")[0] for x in a.accs.read_text().splitlines() if x.strip()} if a.accs else None
    for r in SeqIO.parse(str(a.pool), "genbank"):
        if wanted is not None and r.id.split(".")[0] not in wanted:
            continue
        org = str(r.annotations.get("organism", ""))
        if r.id.startswith("NC_"):
            m = re.search(r"identical to ([A-Z]+\d+)", r.annotations.get("comment", ""))
            copied |= {m.group(1)} if m else set()
        elif not a.all:
            continue
        if " x " in org or not any(f.type == "CDS" for f in r.features):
            continue
        if {"no_cox1", "multi_chromosome"} & set(structure_flags(r)):
            continue
        cand.append((r.id, taxid(r)))
    cache = a.out / "taxonomy.json"
    tax = json.loads(cache.read_text()) if cache.exists() else {}
    need = sorted({t for _, t in cand if t and t not in tax})
    if need:
        tax.update(taxonomy(need))
        cache.write_text(json.dumps(tax))

    by_order = defaultdict(list)
    for acc, t in cand:
        rk = tax.get(t, {})
        if a.all or wanted is not None:
            if rk.get("phylum") and (wanted is not None or acc.split(".")[0] not in copied):
                by_order[(rk["phylum"], rk.get("order", ""))].append(acc)
        elif rk.get("phylum") and rk.get("order"):
            by_order[(rk["phylum"], rk["order"])].append(acc)
    rng = random.Random(a.seed)
    chosen = {}
    for key in sorted(by_order):
        accs = sorted(by_order[key])
        rng.shuffle(accs)
        for acc in accs if a.all or wanted is not None else accs[:a.per_order]:
            chosen[acc] = key
    print(f"[large_set] {len(cand)} usable, {len(copied)} copied by NC_, {len(by_order)} orders, {len(chosen)} chosen", file=sys.stderr)

    rows = []
    for r in SeqIO.parse(str(a.pool), "genbank"):
        if r.id not in chosen:
            continue
        SeqIO.write(r, str(gb_dir / f"{r.id}.gb"), "genbank")
        write_bare(r, bare_dir)
        rk = tax[taxid(r)]
        phylum = rk["phylum"]
        code, assumed = majority_code(r, phylum)
        species = str(r.annotations.get("organism", ""))
        rows.append({"accession": r.id, "source_accession": "", "phylum": phylum, "class": rk.get("class", ""),
                     "order": rk.get("order", ""), "family": rk.get("family", ""),
                     "genus": rk.get("genus", "") or species.split(" ")[0], "species": species,
                     "genetic_code": code, "code_assumed": "yes" if assumed else "",
                     "topology": r.annotations.get("topology", ""), "length": len(r),
                     "evidence_tier": evidence_tier(r), "annotation_method": annotation_method(r, None),
                     "refseq_created": r.annotations.get("date", ""), "structure_flags": ",".join(structure_flags(r)),
                     "mitoz_clade": MITOZ_CLADE.get(phylum) or "", "split": "", "retrieved": today})
    split_by_genus(rows, a.seed)
    a.tsv.parent.mkdir(parents=True, exist_ok=True)
    with a.tsv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["phylum"], r["order"], r["accession"])))
    print(f"[large_set] {len(rows)} records -> {a.tsv}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
