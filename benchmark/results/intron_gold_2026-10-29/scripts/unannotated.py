"""Possible introns NCBI did not annotate in the 193 gold-set records (criterion 5) -> gold/unannotated.tsv.

Per record and standard protein gene (NCBI annotation):
  no_cds          a gene / misc_feature of the gene but no CDS (KT989324.1-like: exon 1 given as misc_feature);
  internal_stop   the CDS as NCBI joins it has internal stop codons;
  short / long    the protein is < 85 % or > 115 % of the median of the same gene in the 20 closest other species;
  mito_intron     an intron candidate of mito_intron.py (MitoM's MFannot port, run on the NCBI record) that no gold
                  intron overlaps (written = what MitoM would annotate; low = reported only).
usage: unannotated.py   (after classify.py, and mito_intron.py runs in gold/work/mi_run/out4, out5)
"""
import csv
import json
import pickle
import sys
from pathlib import Path
from statistics import median

sys.path.insert(0, "/media/deyuan/18T_B/yzz/MitoM/benchmark")
from mito_names import PCG, gene_name  # noqa: E402
from Bio import SeqIO  # noqa: E402
from Bio.Data import CodonTable  # noqa: E402

GOLD = Path(__file__).resolve().parents[1]
SET = GOLD.parent


def closest_lengths(hom, gene: str, organism: str, lineage):
    pool = []
    for x in hom:
        v = x["prot"].get(gene)
        if not v or x["organism"] == organism:
            continue
        depth = 0
        for a, b in zip(x["lineage"], lineage):
            if a != b:
                break
            depth += 1
        pool.append((-depth, len(v[0])))
    return [ln for _, ln in sorted(pool)[:20]]


def main() -> int:
    hom = pickle.load(open(GOLD / "work" / "homologs.pkl", "rb"))     # written by homologs.py (this folder)
    info = {r["accession"]: r for r in csv.DictReader(open(SET / "records.tsv"), delimiter="\t")}
    gold = json.load(open(GOLD / "work" / "gold.json"))
    gold_iv = {}
    for g in gold:
        if g.get("gold_s") is not None:
            gold_iv.setdefault(g["accession"], []).append((g["gold_s"], g["gold_s"] + g["gold_len"], g["strand_i"]))
    rows = []
    for acc, inf in sorted(info.items()):
        rec = SeqIO.read(str(SET / "ncbi_gb" / f"{acc}.gb"), "genbank")
        org, lin = rec.annotations.get("organism", ""), rec.annotations.get("taxonomy", [])
        code = int(inf["genetic_code"])
        base = {"accession": acc, "phylum": inf["phylum"], "order": inf["order"], "organism": org}
        cds_genes = set()
        for f in rec.features:
            if f.type != "CDS":
                continue
            g = gene_name(f.qualifiers, "CDS")
            if g not in PCG:
                continue
            cds_genes.add(g)
            t = int(f.qualifiers.get("transl_table", [code])[0])
            cs = int(f.qualifiers.get("codon_start", ["1"])[0])
            try:
                nt = f.location.extract(rec.seq)[cs - 1:]
            except Exception:
                continue
            p = str(nt[:len(nt) - len(nt) % 3].translate(table=t))
            body = p[:-1] if p.endswith("*") else p
            if body.count("*"):
                rows.append({**base, "gene": g, "kind": "internal_stop", "ncbi_feature": f"CDS {f.location}",
                             "detail": f"{body.count('*')} internal stop(s) in {len(body)} codons"})
            lens = closest_lengths(hom, g, org, lin)
            if lens:
                med = median(lens)
                ratio = len(body.rstrip("*")) / med
                if ratio < 0.85 or ratio > 1.15:
                    partial = "<" in str(f.location) or ">" in str(f.location)
                    rows.append({**base, "gene": g, "kind": "short" if ratio < 1 else "long",
                                 "ncbi_feature": f"CDS {f.location}",
                                 "detail": f"{len(body.rstrip('*'))} aa vs median {med:g} of {len(lens)} close species "
                                           f"({ratio:.0%}){'; partial CDS' if partial else ''}"})
        for f in rec.features:
            if f.type in ("gene", "misc_feature"):
                g = gene_name({"gene": f.qualifiers.get("gene", [""])}, "CDS")
                if g in PCG and g not in cds_genes:
                    rows.append({**base, "gene": g, "kind": "no_cds", "ncbi_feature": f"{f.type} {f.location}",
                                 "detail": " ".join(f.qualifiers.get("note", []))[:120]})
                    cds_genes.add(g)
        # mito_intron candidates without a gold intron
        rep = GOLD / "work" / "mi_run" / f"out{code}" / "intron_report.tsv"
        if rep.exists():
            n = len(rec)
            for r in csv.DictReader(open(rep), delimiter="\t"):
                if r["record"] != acc or not r.get("intron_start"):
                    continue
                s, e = int(r["intron_start"]) - 1, int(r["intron_end"])
                st = 1 if r["strand"] == "+" else -1
                ov = any(st == gst and max(0, min(e, b) - max(s, a)) > 0.5 * (e - s)
                         for a, b, gst in gold_iv.get(acc, []))
                if not ov:
                    rows.append({**base, "gene": gene_name({"gene": [r["gene"]]}, "CDS") or r["gene"],
                                 "kind": "mito_intron", "ncbi_feature": "",
                                 "detail": f"{r['intron_start']}..{r['intron_end']} ({r['strand']}, {r['intron_len']} nt), "
                                           f"{r['confidence']}, written={r['written']}: {r['evidence']} | {r['note'][:150]}"})
    with open(GOLD / "unannotated.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, ["accession", "phylum", "order", "organism", "gene", "kind", "ncbi_feature", "detail"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    from collections import Counter
    print(Counter(r["kind"] for r in rows), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
