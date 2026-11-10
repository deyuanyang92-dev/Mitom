"""Protein genes MitoM annotates in complete genome-project mitogenomes (mito_bench_run.py --tools mitom_lo).

Per record: the 13 genes found in MitoM's final annotation, which came from the protein-gene HMM step (CDS /note
'found by MitoM mito_pcg_hmm'), which are missing; where the project annotated the record, its genes for comparison
(same gene within 50 % overlap = agree). Writes per_record.tsv and prints the per-phylum table.
Usage: python3 summarise.py <root> <out dir>
"""
import csv
import glob
import sys
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from mito_names import PCG, gene_name  # noqa: E402

root, out = Path(sys.argv[1]), Path(sys.argv[2])
rows = list(csv.DictReader(open(root / "set.tsv"), delimiter="\t"))
per = []
for r in rows:
    acc = r["accession"]
    gbs = [g for g in glob.glob(str(root / "runs" / "mitom_lo" / acc / "raw" / "annotation" / "*.gb"))
           if "all_records" not in g]
    if not gbs:
        per.append({"accession": acc, "phylum": r["phylum"], "genus": r["genus"], "status": "no output"})
        continue
    rec = SeqIO.read(gbs[0], "genbank")
    got, hmm = {}, set()
    for f in rec.features:
        g = gene_name(f.qualifiers, "CDS") if f.type == "CDS" else ""
        if g in PCG:
            got.setdefault(g, f)
            if any("mito_pcg_hmm" in n for n in f.qualifiers.get("note", [])):
                hmm.add(g)
    proj = SeqIO.read(str(root / "gb" / f"{acc}.gb"), "genbank")
    pg = {gene_name(f.qualifiers, "CDS") for f in proj.features if f.type == "CDS"} - {""}
    per.append({"accession": acc, "phylum": r["phylum"], "genus": r["genus"], "status": "ok", "found": len(got),
                "missing": ",".join(sorted(set(PCG) - set(got))), "by_hmm": ",".join(sorted(hmm)),
                "project_genes": len(pg), "project_only": ",".join(sorted(pg - set(got))),
                "mitom_only": ",".join(sorted(set(got) - pg)) if pg else ""})
cols = ["accession", "phylum", "genus", "status", "found", "missing", "by_hmm", "project_genes", "project_only",
        "mitom_only"]
with open(out / "per_record.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(per)
by = defaultdict(list)
for x in per:
    by[x["phylum"]].append(x)
print("| 门 / phylum | 记录 / records | 13 个全 / all 13 | 缺的基因（次数）/ missing genes (n) | HMM 补回 / added by HMM |")
print("|---|---|---|---|---|")
for ph, xs in sorted(by.items(), key=lambda x: -len(x[1])):
    ok = [x for x in xs if x["status"] == "ok"]
    miss = Counter(g for x in ok for g in x["missing"].split(",") if g)
    hm = Counter(g for x in ok for g in x["by_hmm"].split(",") if g)
    print(f"| {ph} | {len(ok)}/{len(xs)} | {sum(x['found'] == 13 for x in ok)} | "
          f"{', '.join(f'{g} {n}' for g, n in miss.most_common()) or '—'} | "
          f"{', '.join(f'{g} {n}' for g, n in hm.most_common()) or '—'} |")
