#!/usr/bin/env python3
"""frozen_set.py - leak audit and independent test-set selection for the frozen evaluation (benchmark/frozen/PLAN.md).

"Seen" = anything MitoM was trained, tuned or already evaluated on, collected from the repository itself:
  accessions  every accession written in benchmark/ (sets, tRNA splits, results), models/ (rRNA HMM references,
              intron sites, refseq89m list ref_taxa), knowledge/experiments/ and tests/data/ (version ignored)
  taxa        genus / family of those records: set TSV columns, the lineages of ref_taxa (refseq89m) and
              intron_sites, and the pool records whose accession is seen
  clade CMs   taxon names of the MitoM-ML tRNA clade models (models/trna_ml/clades.tar.gz)
Conservative on purpose: an accession mentioned anywhere counts as seen.

  audit   one row per record of a GenBank pool (or a set TSV): acc_seen, species / genus / family seen, clade CM
          taxon in the lineage (the most specific; clade_cm_rank = its depth from the
          genus, 1 = genus) -> <out>.tsv and a summary on stdout
  select  from a pool: records not seen and of a genus not seen (--level genus, default; family = stricter),
          complete-looking (CDS with cox1, >= 10 kb), at most --per_genus per genus, --per_family per family and
          --per_phylum per phylum, seeded -> accession list for make_large_set.py --accs

Usage:
  python3 frozen_set.py audit  --pool Metazoa_refseq.gb[.gz] -o audit.tsv
  python3 frozen_set.py select --pool Metazoa_refseq.gb[.gz] -o test_accs.txt [--level genus] [--per_phylum 40]
                               [--per_family 2] [--per_genus 1] [--seed 20261108]
"""
import argparse
import csv
import gzip
import io
import random
import re
import sys
import tarfile
from collections import Counter
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set, Tuple

from Bio import SeqIO

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "benchmark"))
from mito_taxa import METAZOAN_PHYLA  # noqa: E402
from mito_sample import source_accession  # noqa: E402

ACC = re.compile(r"\b(?:N[CWTZ]|AC)_\d{6,9}|\b[A-Z]{1,2}\d{5,8}\b")
TEXT_SUFFIXES = (".tsv", ".txt", ".md", ".csv", ".json", ".tsv.gz", ".txt.gz")
SEEN_DIRS = ("benchmark", "models", "knowledge/experiments", "tests/data")
SKIP_DIRS = ("benchmark/frozen",)          # this plan's own outputs are not training data


def _text(p: Path) -> str:
    data = p.read_bytes()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return data.decode("utf-8", errors="replace")


def strip_v(acc: str) -> str:
    return acc.split(".")[0]


def lineage_taxa(lineage: List[str]) -> Tuple[str, str]:
    """(genus, family) from an NCBI lineage list: genus = last element, family = last element ending in -idae."""
    genus = lineage[-1] if lineage else ""
    fam = next((x for x in reversed(lineage) if x.endswith("idae")), "")
    return genus, fam


def seen_accessions(root: Path = ROOT) -> Dict[str, str]:
    """{accession without version: first repository file that names it}."""
    out: Dict[str, str] = {}
    for d in SEEN_DIRS:
        for p in sorted((root / d).rglob("*")):
            rel = p.relative_to(root).as_posix()
            if not p.is_file() or any(rel.startswith(s) for s in SKIP_DIRS) or not p.name.endswith(TEXT_SUFFIXES + (".gb", ".gbf")):
                continue
            for m in ACC.findall(_text(p)):
                out.setdefault(strip_v(m), rel)
    return out


def seen_taxa(root: Path = ROOT) -> Tuple[Set[str], Set[str], Set[str]]:
    """(species, genera, families) of the records in the set TSVs, refseq89m (ref_taxa) and the intron sites."""
    sp, ge, fa = set(), set(), set()
    for t in sorted((root / "benchmark" / "sets").glob("*.tsv")):
        for r in csv.DictReader(io.StringIO(_text(t)), delimiter="\t"):
            sp.add(r.get("species") or "")
            ge.add(r.get("genus") or "")
            fa.add(r.get("family") or "")
    for f, col in ((root / "models" / "ref_taxa.tsv.gz", "lineage"), (root / "models" / "intron_sites.tsv.gz", "lineage")):
        if f.is_file():
            for r in csv.DictReader(io.StringIO(_text(f)), delimiter="\t"):
                g, fam = lineage_taxa((r.get(col) or "").split(";"))
                ge.add(g)
                fa.add(fam)
                if r.get("organism"):
                    sp.add(r["organism"])
    return sp - {""}, ge - {""}, fa - {""}


def clade_cm_taxa(root: Path = ROOT) -> Set[str]:
    """Taxon names of the tRNA clade CMs (clades/<Taxon>__<aa>.cm)."""
    f = root / "models" / "trna_ml" / "clades.tar.gz"
    if not f.is_file():
        return set()
    with tarfile.open(f) as t:
        return {Path(n).name.split("__")[0] for n in t.getnames() if n.endswith(".cm")}


def read_pool(path: Path) -> Iterator[dict]:
    """Records of a (gzipped) GenBank pool: accession, source accession, organism, lineage, phylum, length, cox1."""
    handle = io.StringIO(_text(path))
    for rec in SeqIO.parse(handle, "genbank"):
        lin = list(rec.annotations.get("taxonomy", []))
        genes = {(f.qualifiers.get("gene") or [""])[0].lower() for f in rec.features if f.type == "CDS"}
        genus, fam = lineage_taxa(lin)
        yield {"accession": rec.id, "source_accession": source_accession(rec),
               "organism": rec.annotations.get("organism", ""), "genus": genus, "family": fam,
               "phylum": next((x for x in lin if x in METAZOAN_PHYLA), ""), "lineage": ";".join(lin),
               "length": len(rec.seq), "n_cds": sum(f.type == "CDS" for f in rec.features),
               "cox1": bool(genes & {"cox1", "coi", "co1", "coxi"})}


def audit(rows: List[dict], accs: Dict[str, str], taxa, cms: Set[str]) -> List[dict]:
    sp, ge, fa = taxa
    out = []
    for r in rows:
        hit = accs.get(strip_v(r["accession"])) or (accs.get(strip_v(r["source_accession"])) if r["source_accession"] else None)
        lin = r["lineage"].split(";")
        out.append({**{k: r[k] for k in ("accession", "source_accession", "organism", "phylum", "family", "genus", "length")},
                    "acc_seen": hit or "", "species_seen": "yes" if r["organism"] in sp else "",
                    "genus_seen": "yes" if r["genus"] in ge else "", "family_seen": "yes" if r["family"] in fa else "",
                    "clade_cm": next((x for x in reversed(lin) if x in cms), ""),       # most specific tRNA clade CM
                    "clade_cm_rank": next((str(i) for i, x in enumerate(reversed(lin), 1) if x in cms), "")})
    return out


def select(aud: List[dict], pool: Dict[str, dict], level: str, per_phylum: int, per_family: int, per_genus: int,
           seed: int) -> List[dict]:
    ok = [r for r in aud if not r["acc_seen"] and not r["genus_seen"] and (level == "genus" or not r["family_seen"])
          and r["phylum"] and pool[r["accession"]]["cox1"] and pool[r["accession"]]["n_cds"] >= 10
          and int(r["length"]) >= 10000]
    random.Random(seed).shuffle(ok)
    n_ph, n_fa, n_ge, got = Counter(), Counter(), Counter(), []
    for r in ok:
        if n_ph[r["phylum"]] >= per_phylum or n_ge[r["genus"]] >= per_genus or (r["family"] and n_fa[r["family"]] >= per_family):
            continue
        n_ph[r["phylum"]] += 1
        n_ge[r["genus"]] += 1
        n_fa[r["family"]] += 1
        got.append(r)
    return sorted(got, key=lambda r: (r["phylum"], r["family"], r["accession"]))


def write_tsv(path: Path, rows: List[dict]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]) if rows else ["accession"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    for name in ("audit", "select"):
        s = sub.add_parser(name)
        s.add_argument("--pool", required=True, type=Path, help="GenBank pool (.gb / .gb.gz), e.g. mitom fetch output")
        s.add_argument("-o", "--out", required=True, type=Path)
        if name == "select":
            s.add_argument("--level", choices=["genus", "family"], default="genus")
            s.add_argument("--per_phylum", type=int, default=40)
            s.add_argument("--per_family", type=int, default=2)
            s.add_argument("--per_genus", type=int, default=1)
            s.add_argument("--seed", type=int, default=20261108)
    a = p.parse_args(argv)
    pool = {r["accession"]: r for r in read_pool(a.pool)}
    accs, taxa, cms = seen_accessions(), seen_taxa(), clade_cm_taxa()
    print(f"[frozen_set] seen: {len(accs)} accessions, {len(taxa[1])} genera, {len(taxa[2])} families, "
          f"{len(cms)} clade-CM taxa; pool {len(pool)} records", file=sys.stderr)
    aud = audit(list(pool.values()), accs, taxa, cms)
    c = Counter()
    for r in aud:
        c["acc_seen"] += bool(r["acc_seen"])
        c["genus_seen"] += bool(r["genus_seen"])
        c["family_seen"] += bool(r["family_seen"])
        c["clean_genus"] += not r["acc_seen"] and not r["genus_seen"]
        c["clean_family"] += not r["acc_seen"] and not r["family_seen"]
    print("[frozen_set] " + ", ".join(f"{k} {v}" for k, v in c.items()), file=sys.stderr)
    if a.cmd == "audit":
        write_tsv(a.out, aud)
        return 0
    got = select(aud, pool, a.level, a.per_phylum, a.per_family, a.per_genus, a.seed)
    a.out.write_text("".join(r["accession"] + "\n" for r in got))
    write_tsv(a.out.with_suffix(".audit.tsv"), got)
    per = Counter(r["phylum"] for r in got)
    print(f"[frozen_set] selected {len(got)} records ({a.level} unseen), "
          + ", ".join(f"{k} {v}" for k, v in per.most_common()), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
