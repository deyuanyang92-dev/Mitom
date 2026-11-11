#!/usr/bin/env python3
"""MITOS2 (raw and + MitoM curation) vs MitoM vs MitoM + RNA-seq (mito_rnaseq_annot.py) on the annelid pairs.

Same records for all: MITOS2 = tools_root/runs/mitos2/<acc>/final.gbf, MITOS2+cur = runs/mitos2_cur (mito_bench_run.py
`mitom import` with boundary curation; intron fix, rescue, rRNA HMM off), MITOS2+all = runs/mitos2_full (`mitom
import` defaults: every later MitoM step); NCBI coordinates;
MitoM = mitom_out/annotation/<acc>.gb; MitoM + RNA = rnaseq_annot/<acc>/rnaseq_annotated.gb (both in MitoM's rotated
coordinates; MITOS2 features are mapped onto them by rotation / strand). Pairs whose species check failed are left out.

Per tool, over the 13 protein genes (one row per record in tool_compare.tsv, totals printed):
  found        standard protein genes annotated
  clean        their translation (code 5, joins followed) has no internal stop
  length_ok    protein within 2 aa of the median of annotated relatives of the family (>= 3 in models/pcg_refs)
  end_polya    3' end with >= 3 poly(A) reads from 10 nt before to 30 nt after it (mito_rna sites, not internal
               priming) - an independent check for MITOS2 and MitoM; partly circular for MitoM + RNA, which moves
               ends to poly(A) sites (shown separately for the test pairs)
  introns      the RNA-proven introns (split reads, mito_rnaseq_annot step 1) annotated with exact ends
Usage: compare_tools.py <work> -o tool_compare.tsv
"""
import argparse
import csv
import gzip
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

H = Path(__file__).resolve()
sys.path.insert(0, str(H.parents[3]))
from mito_names import PCG, gene_name  # noqa: E402

TOOLS = ("MITOS2", "MITOS2+cur", "MITOS2+all", "MitoM", "MitoM+RNA")
GENE_REF = {"atp6": "atp6", "atp8": "atp8", "cox1": "cox1", "cox2": "cox2", "cox3": "cox3", "cob": "cob",
            "nad1": "nad1", "nad2": "nad2", "nad3": "nad3", "nad4": "nad4", "nad4l": "nad4l", "nad5": "nad5",
            "nad6": "nad6"}


def mapper(src: str, dst: str):
    """Position map src (0-based) -> (dst position, strand sign) for a rotated and / or reverse-complemented copy."""
    n = len(dst)
    probe = src[:60]
    i = (dst + dst).find(probe)
    if i >= 0:
        return lambda p: ((p + i) % n, 1)
    rc = str(Seq(src).reverse_complement())
    j = (dst + dst).find(rc[:60])                      # rc(src)[0] = src[-1]
    if j >= 0:
        return lambda p: ((j + (len(src) - 1 - p)) % n, -1)
    return None


def genes(rec, mp=None):
    """standard gene -> list of (3' end 1-based in the MitoM frame, strand, protein, parts [(s, e)] in that frame)."""
    out = defaultdict(list)
    for f in rec.features:
        if f.type != "CDS" or f.location is None:
            continue
        g = gene_name(f.qualifiers, "CDS")
        if g not in PCG:
            continue
        code = int(f.qualifiers.get("transl_table", ["5"])[0])
        nt = str(f.location.extract(rec.seq)).upper()
        prot = str(Seq(nt[: len(nt) // 3 * 3]).translate(table=code)).rstrip("*")
        st = f.location.strand or 1
        end0 = int(f.location.end) - 1 if st == 1 else int(f.location.start)
        parts = [(int(p.start), int(p.end)) for p in f.location.parts]
        if mp:
            e, sgn = mp(end0)
            end0, st = e, st * sgn
            parts = sorted(tuple(sorted((mp(s)[0], mp(e2 - 1)[0] + 1))) for s, e2 in parts)
        out[g].append((end0 + 1, st, prot, parts))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("work")
    p.add_argument("-o", "--out", required=True)
    a = p.parse_args()
    w = Path(a.work)
    here = H.parent
    clean = {x.strip() for x in open(here / "clean_pairs.txt") if x.strip()}
    res = {r["mito"]: r for r in csv.DictReader(open(here / "results.tsv"), delimiter="\t")}
    fam = {r["species"]: r["family"] for r in csv.DictReader(
        open(w / "Annelida_genome_inventory/Annelida_species.tsv"), delimiter="\t")}
    refs = defaultdict(list)
    for r in SeqIO.parse(gzip.open(H.parents[3] / "models/pcg_refs.faa.gz", "rt"), "fasta"):
        d = r.description.split("|")
        if len(d) >= 4:
            refs[d[1].lower()].append((d[3], len(r.seq)))
    rows = []
    tot = {g: {t: Counter() for t in TOOLS} for g in ("learning", "test")}
    for mito, r in res.items():
        if r["species_check"] != "ok":
            continue
        acc = mito.split(".")[0]
        grp = "learning" if mito in clean else "test"
        mm = SeqIO.read(w / f"mitom_out/annotation/{mito}.gb", "genbank")
        recs = {"MitoM": (mm, None), "MitoM+RNA": (SeqIO.read(w / f"rnaseq_annot/{acc}/rnaseq_annotated.gb",
                                                               "genbank"), None)}
        for tool, run in (("MITOS2", "mitos2"), ("MITOS2+cur", "mitos2_cur"), ("MITOS2+all", "mitos2_full")):
            m2 = w / f"tools_root/runs/{run}/{mito}/final.gbf"
            if m2.exists():
                mr = SeqIO.read(m2, "genbank")
                recs[tool] = (mr, mapper(str(mr.seq).upper(), str(mm.seq).upper()))
        sites = [(int(s["position"]), s["strand"], int(s["reads"])) for s in csv.DictReader(
            open(w / f"rna_out/{acc}/rna_polya_sites.tsv"), delimiter="\t") if not s["internal_priming"]]
        unstr = "unstranded" in r["library"]
        introns = []  # RNA-proven (split reads): exon 1 ends at a0 (0-based, excl.), exon 2 starts at b0
        for e in csv.DictReader(open(w / f"rnaseq_annot/{acc}/rnaseq_evidence.tsv"), delimiter="\t"):
            if e["step"] == "1" and e["verdict"] == "confirmed" and " intron " in e["feature"]:
                gname, _, span = e["feature"].split(" ")
                i0, i1 = (int(x) for x in span.split(".."))
                introns.append((gname, i0 - 1, i1))
        f_ = fam.get(r["organism"], "")
        for tool in TOOLS:
            if tool not in recs or (tool.startswith("MITOS2") and recs[tool][1] is None):
                continue
            rec, mp = recs[tool]
            gs = genes(rec, mp)
            c = Counter()
            for g, items in gs.items():
                c["found"] += 1
                end1, st, prot, parts = max(items, key=lambda x: len(x[2]))
                c["clean"] += "*" not in prot
                xs = [n for lin, n in refs.get(GENE_REF.get(g, ""), []) if f_ and f_ in lin]
                if len(xs) >= 3:
                    c["length_tested"] += 1
                    c["length_ok"] += abs(len(prot) - statistics.median(xs)) <= 2
                sign = "+" if st == 1 else "-"
                c["end_polya"] += sum(k for pos, s, k in sites if (s == sign or unstr)
                                      and -10 <= (pos - end1) * st <= 30) >= 3
            for gname, a0, b0 in introns:
                c["introns_tested"] += 1
                cands = gs.get(gene_name({"gene": [gname]}, "CDS"), [(0, 1, "", [])])
                parts = max(cands, key=lambda x: len(x[2]))[3]
                c["introns_exact"] += any(s1[1] == a0 and s2[0] == b0 for s1, s2 in zip(parts, parts[1:]))
            tot[grp][tool] += c
            rows.append({"mito": mito, "group": grp, "tool": tool, **{k: c[k] for k in (
                "found", "clean", "length_tested", "length_ok", "end_polya", "introns_tested", "introns_exact")}})
    with open(a.out, "w", newline="") as fh:
        wr = csv.DictWriter(fh, list(rows[0]), delimiter="\t", lineterminator="\n")
        wr.writeheader()
        wr.writerows(rows)
    print("group\ttool\trecords\tfound\tclean\tlength_ok/tested\tend_polya\tintrons exact/tested")
    for grp in ("learning", "test"):
        for t in TOOLS:
            c = tot[grp][t]
            nrec = sum(1 for x in rows if x["group"] == grp and x["tool"] == t)
            print(f"{grp}\t{t}\t{nrec}\t{c['found']}\t{c['clean']}\t{c['length_ok']}/{c['length_tested']}\t"
                  f"{c['end_polya']}\t{c['introns_exact']}/{c['introns_tested']}")


if __name__ == "__main__":
    main()
