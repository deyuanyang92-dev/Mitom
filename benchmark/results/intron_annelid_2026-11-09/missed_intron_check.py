"""Genes MitoM annotated short (< 70 % of the usual length) or not at all in the unannotated annelid records: is the
rest of the protein elsewhere on the record (an intron MitoM missed), or nowhere (a divergent / missing gene)?

For each such gene: tblastn (genetic code 5, E <= 1e-3) of the annelid proteins of that gene in MITOS2 refseq89m
(models/ref_taxa.tsv.gz lineage Annelida) against the record. Homolog stretches the MitoM CDS covers vs stretches on
the same strand outside it (>= 20 aa of the protein the CDS does not cover, within 6 kb of the CDS): the latter is
the signature of a missed intron (or of a gene annotated at the wrong place).
Usage: python3 missed_intron_check.py <root> <unannotated_introns.tsv> <refseq89m/featureProt> > missed_intron_check.tsv
"""
import csv
import glob
import gzip
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

root, table, fp = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
M = Path(__file__).resolve().parents[3]
ann = {r["accession"] for r in csv.DictReader(gzip.open(M / "models" / "ref_taxa.tsv.gz", "rt"), delimiter="\t")
       if "Annelida" in r["lineage"]}
TYP = {"cox1": 513, "cox2": 228, "cox3": 259, "atp6": 225, "cob": 379, "nad1": 317, "nad2": 335, "nad3": 115,
       "nad4": 449, "nad5": 580, "nad6": 160}
FILE = {"cob": "cob", **{g: g for g in TYP}}
rows = list(csv.DictReader(open(table), delimiter="\t"))
by = defaultdict(list)
for r in rows:
    by[r["accession"]].append(r)
w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
w.writerow(["accession", "family", "organism", "gene", "mitom_aa", "verdict", "outside_hits"])
with tempfile.TemporaryDirectory() as td:
    for acc, rs in sorted(by.items()):
        cds = defaultdict(list)
        for r in rs:
            if r["kind"] == "cds":
                cds[r["gene"]].append((int(r["start"]), int(r["end"]), int(r["length"])))
        todo = [g for g in TYP if g not in cds or sum(x[2] for x in cds[g]) < 0.7 * TYP[g]]
        if not todo or len([g for g in TYP if g not in cds]) >= 6:
            continue
        sub = Path(td) / f"{acc}.fa"
        SeqIO.write(SeqIO.read(str(root / "gb" / f"{acc}.gb"), "genbank"), str(sub), "fasta")
        for g in todo:
            q = Path(td) / f"{g}.faa"
            if not q.exists():
                SeqIO.write([r for r in SeqIO.parse(str(fp / f"{FILE[g]}.fas"), "fasta")
                             if r.id.split(":")[0] in ann], str(q), "fasta")
            out = subprocess.run(["tblastn", "-query", str(q), "-subject", str(sub), "-db_gencode", "5", "-evalue",
                                  "1e-3", "-outfmt", "6 qseqid qstart qend sstart send evalue"],
                                 capture_output=True, text=True).stdout
            have = cds.get(g, [])
            others = [x for k, v in cds.items() if k != g for x in v]
            inside, extends, outside = 0, 0, defaultdict(int)
            for ln in out.splitlines():
                _, qs, qe, ss, se, ev = ln.split("\t")
                s, e = sorted((int(ss), int(se)))
                if any(a - 30 <= s and e <= b + 30 for a, b, _ in have):
                    inside += 1
                elif any(a <= e and s <= b for a, b, _ in have):
                    extends += 1                               # overlaps the CDS and runs past its end
                elif any(a <= e and s <= b for a, b, _ in others) or int(qe) - int(qs) < 20:
                    continue                                   # (inside another gene: cross-homology; tiny)
                elif not have or min(min(abs(s - b), abs(a - e)) for a, b, _ in have) < 6000:
                    outside[(s // 100 * 100, e // 100 * 100)] += 1
            top = sorted(outside.items(), key=lambda x: -x[1])[:3]
            gap = min((min(abs(a - e), abs(s - b)) for (s, e), _ in top for a, b, _ in have), default=None)
            verdict = ("homolog separate from the CDS (>= 142 nt): possible missed intron or misplaced gene"
                       if top and (gap is None or gap >= 142) else
                       "homolog next to the CDS (< 142 nt): CDS too short, not an intron" if top else
                       "homolog runs past the CDS end: CDS too short, not an intron" if extends else
                       "no homolog outside the annotated genes: divergent or absent gene")
            w.writerow([acc, rs[0]["family"], rs[0]["organism"], g, sum(x[2] for x in have), verdict,
                        "; ".join(f"{a}..{b} ({n} refs)" for (a, b), n in top)])
