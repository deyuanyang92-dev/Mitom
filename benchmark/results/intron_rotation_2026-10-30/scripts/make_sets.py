"""Build mito_bench_run roots from the NCBI records in all.gb:
  orig   bare FASTA as in NCBI (no topology label: `mitom run` decides it, as for an assembly)
  exon   circular records rotated so that the origin lies in the middle of the last exon of an intron-bearing gene
  intron circular records rotated so that the origin lies in the middle of that gene's last intron
  near3 / near5  origin 15 nt inside the right / left genomic end of that intron (the CHY records: 15 nt before the
         intron's 3' end), so the gene's part at one end of the record stops right at the junction
The gene: cox1 when it has an intron, else the first gene with one (cis-spliced, parts >= 142 nt apart).
Intron-free records: orig, and exon with the origin in the middle of cox1.
    make_sets.py all.gb groups.tsv records.tsv OUTDIR MITOM_BENCHMARK_DIR"""
import csv
import sys
from pathlib import Path

from Bio import SeqIO

gbfile, groups, records, out, bench = sys.argv[1:6]
sys.path.insert(0, bench)
from mito_sample import COLUMNS, MITOZ_CLADE, DEFAULT_CODE  # noqa: E402
from mito_names import gene_name  # noqa: E402

MIN_INTRON = 142
grp = {r["accession"]: r for r in csv.DictReader(open(groups), delimiter="\t")}
gold = {r["accession"]: r for r in csv.DictReader(open(records), delimiter="\t")}
out = Path(out)


def intron_gene(rec):
    """(gene, last exon (s, e), last intron (s, e)) in 0-based transcript coordinates (may exceed len) or None."""
    n = len(rec)
    found = []
    for f in rec.features:
        if f.type != "CDS" or len(f.location.parts) < 2:
            continue
        parts = list(f.location.parts)
        if len({p.strand for p in parts}) > 1:
            continue                                      # trans-spliced
        g = gene_name(f.qualifiers, "CDS") or ""
        strand = parts[0].strand
        # intron = gap between consecutive parts in transcript order, measured along the strand (mod n)
        gaps = []
        for a, b in zip(parts, parts[1:]):
            gap = (int(b.start) - int(a.end)) % n if strand >= 0 else (int(a.start) - int(b.end)) % n
            gaps.append(gap)
        if max(gaps) < MIN_INTRON:
            continue
        i = max(k for k, x in enumerate(gaps) if x >= MIN_INTRON)   # last intron
        a, b = parts[i], parts[i + 1]
        if strand >= 0:
            intr = (int(a.end), int(a.end) + gaps[i])
        else:
            intr = (int(b.end), int(b.end) + gaps[i])
        found.append((g, (int(b.start), int(b.end)), intr))
    if not found:
        return None
    found.sort(key=lambda x: x[0] != "cox1")
    return found[0]


def cox1(rec):
    for f in rec.features:
        if f.type in ("CDS", "gene") and gene_name(f.qualifiers, "CDS") == "cox1":
            return int(f.location.start), int(f.location.end)
    return None


rows = {k: [] for k in ("orig", "exon", "intron", "near3", "near5")}
log = open(out / "rotations.tsv", "w")
log.write("accession\tgroup\tphylum\tset\tgene\tnew_origin_ncbi_1based\tnote\n")
for rec in SeqIO.parse(gbfile, "genbank"):
    acc = rec.id
    if acc not in grp:
        print("not in groups:", acc, file=sys.stderr)
        continue
    g = grp[acc]
    n = len(rec)
    topo = rec.annotations.get("topology", "circular")
    code = ""
    for f in rec.features:
        if f.type == "CDS" and "transl_table" in f.qualifiers:
            code = f.qualifiers["transl_table"][0]
            break
    code = code or str(DEFAULT_CODE[g["phylum"]])
    plan = {"orig": 0}
    if topo == "circular":
        ig = intron_gene(rec) if g["group"] != "neg" else None
        if ig:
            gene, ex, intr = ig
            plan["exon"] = ((ex[0] + ex[1]) // 2) % n
            plan["intron"] = ((intr[0] + intr[1]) // 2) % n
            plan["near3"] = (intr[1] - 15) % n            # 15 nt before the intron's genomic right end (as CHY)
            plan["near5"] = (intr[0] + 15) % n            # 15 nt after its genomic left end
            what = gene
        elif g["group"] == "neg" and cox1(rec):
            s, e = cox1(rec)
            plan["exon"] = (s + e) // 2
            what = "cox1 (no intron)"
        else:
            what = "no cis-spliced intron gene"
    else:
        what = "linear record: not rotated"
    for k, s in plan.items():
        root = out / k
        (root / "gb").mkdir(parents=True, exist_ok=True)
        (root / "bare" / acc).mkdir(parents=True, exist_ok=True)
        SeqIO.write(rec, str(root / "gb" / f"{acc}.gb"), "genbank")
        seq = str(rec.seq)
        (root / "bare" / acc / f"{acc}.fasta").write_text(f">{acc}\n{seq[s:] + seq[:s]}\n")
        rows[k].append({**{c: "" for c in COLUMNS}, "accession": acc, "phylum": g["phylum"], "genetic_code": code,
                        "topology": topo, "length": n, "mitoz_clade": MITOZ_CLADE.get(g["phylum"]) or "",
                        "split": "test", "retrieved": "efetch"})
        log.write(f"{acc}\t{g['group']}\t{g['phylum']}\t{k}\t{what if k != 'orig' else '-'}\t{s + 1}\t\n")
for k, rs in rows.items():
    with open(out / k / "set.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rs)
    print(k, len(rs), file=sys.stderr)
