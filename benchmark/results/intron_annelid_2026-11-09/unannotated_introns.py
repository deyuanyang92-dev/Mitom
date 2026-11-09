"""Introns MitoM writes in the NCBI annelid mitogenomes that NCBI holds without any annotation (sequence only), and
the candidates it did not write; each gene's protein checked for completeness.

Input: a mito_bench_run.py root with runs/mitom_native/<acc>/raw/ (`mitom run` defaults, genetic code 5).
Per record and protein gene: CDS length (aa), parts, internal stops; introns written (annotation/*.gb intron
features: gene, length, type from the /note); intron_report.tsv rows not written (gene, confidence, evidence).
Usage: python3 unannotated_introns.py <root> <work/set.tsv> > unannotated_introns.tsv
"""
import csv
import glob
import sys
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from mito_names import gene_name  # noqa: E402

root = Path(sys.argv[1])
fam = {}
for f in glob.glob(str(root / "gb" / "*.gb")):
    r = SeqIO.read(f, "genbank")
    tax = r.annotations.get("taxonomy", [])
    fam[r.id] = (next((t for t in tax if t.endswith("idae")), ""), r.annotations.get("organism", ""))
w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
w.writerow(["accession", "family", "organism", "kind", "gene", "start", "end", "length", "detail"])
for d in sorted((root / "runs" / "mitom_native").iterdir()):
    acc = d.name
    gbs = [g for g in glob.glob(str(d / "raw" / "annotation" / "*.gb")) if "all_records" not in g]
    if not gbs:
        w.writerow([acc, *fam.get(acc, ("", "")), "no_output", "", "", "", "", ""])
        continue
    rec = SeqIO.read(gbs[0], "genbank")
    for f in rec.features:
        if f.type == "intron":
            note = " ".join(f.qualifiers.get("note", []))[:120]
            w.writerow([acc, *fam[acc], "intron", f.qualifiers.get("gene", [""])[0], int(f.location.start) + 1,
                        int(f.location.end), len(f), note])
        elif f.type == "CDS":
            g = gene_name(f.qualifiers, "CDS") or f.qualifiers.get("gene", [""])[0]
            code = int(f.qualifiers.get("transl_table", ["5"])[0])
            nt = f.location.extract(rec.seq)
            prot = str(nt[: len(nt) // 3 * 3].translate(table=code)).rstrip("*")
            w.writerow([acc, *fam[acc], "cds", g, int(f.location.start) + 1, int(f.location.end), len(prot),
                        f"parts={len(f.location.parts)}; internal_stops={prot.count('*')}"])
    rep = d / "raw" / "qc" / "intron_report.tsv"
    if rep.exists():
        for r in csv.DictReader(open(rep), delimiter="\t"):
            if r.get("intron_start") and r.get("written") != "yes":
                w.writerow([acc, *fam[acc], "intron_not_written", r.get("gene", ""), r["intron_start"],
                            r["intron_end"], "", f"{r.get('confidence', '')}: {r.get('evidence', '')[:150]}"])
            elif not r.get("intron_start") and "internal stop" in (r.get("note", "") + r.get("evidence", "")):
                w.writerow([acc, *fam[acc], "internal_stop", r.get("gene", ""), "", "", "",
                            (r.get("note", "") or r.get("evidence", ""))[:150]])
