"""Which placement gives the better protein where MitoM's intron differs from NCBI's (or NCBI has one MitoM lacks)?

For every truth row of eval_introns.py (classes exact / le6 / gt6 / missed) and every fp: the CDS of that gene as NCBI
splices it and as MitoM splices it (runs/<tool>/<acc>/final.gbf, NCBI coordinates), translated with the record's code;
internal stops; blastp against the refseq89m proteins of the gene (MITOS2 reference set): best bit score, query
coverage of the best hit. The placement with no internal stop and the higher bit score is the better-supported one.

Usage: python3 check_products.py <root> <tool> <refseq89m/featureProt> > products.tsv
"""
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from mito_names import gene_name  # noqa: E402

root, tool, prot = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
BLASTP = "blastp"


def cds_of(rec, gene, s, e):
    """The CDS of `gene` in rec overlapping s..e (1-based), else None."""
    for f in rec.features:
        if f.type == "CDS" and gene_name(f.qualifiers, "CDS") == gene and \
                int(f.location.start) < e and int(f.location.end) >= s:
            return f
    return None


def product(rec, f):
    code = int(f.qualifiers.get("transl_table", ["4"])[0])
    nt = f.location.extract(rec.seq)
    nt = nt[: len(nt) // 3 * 3]
    aa = str(nt.translate(table=code)).rstrip("*")
    return aa, aa.count("*")


def blast(aa, gene):
    db = prot / f"{gene}.fas"
    if not aa or not db.exists():
        return "", ""
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as q:
        q.write(f">q\n{aa.replace('*', 'X')}\n")
    out = subprocess.run([BLASTP, "-query", q.name, "-db", str(db), "-outfmt", "6 bitscore qcovs gaps",
                          "-max_target_seqs", "1", "-seg", "no"], capture_output=True, text=True).stdout.split("\n")[0]
    Path(q.name).unlink()
    return tuple(out.split("\t")[:2]) if out else ("0", "0")


def span(t):
    a, b = t.split("..")
    return int(a), int(b)


rows = list(csv.DictReader(open(Path(__file__).resolve().parent / f"{tool}_introns.tsv"), delimiter="\t"))
w = csv.writer(sys.stdout, delimiter="\t")
w.writerow(["accession", "group", "organism", "class", "gene", "ncbi_intron", "mitom_intron",
            "ncbi_aa", "ncbi_stops", "ncbi_bits", "ncbi_qcov", "mitom_aa", "mitom_stops", "mitom_bits", "mitom_qcov",
            "better"])
for r in rows:
    loc = r["truth"] or r["pred"]
    if not loc:
        continue
    s, e = span(loc)
    ncbi = SeqIO.read(str(root / "gb" / f"{r['accession']}.gb"), "genbank")
    pred = SeqIO.read(str(root / "runs" / tool / r["accession"] / "final.gbf"), "genbank")
    gene = gene_name({"gene": [r["gene"] or r["pred_gene"]]}, "CDS") if (r["gene"] or r["pred_gene"]) else ""
    if not gene:                                     # NCBI intron without /gene: the CDS around it
        for f in ncbi.features + pred.features:
            if f.type == "CDS" and int(f.location.start) < s and int(f.location.end) > e:
                gene = gene_name(f.qualifiers, "CDS")
                break
    res = []
    for rec in (ncbi, pred):
        f = cds_of(rec, gene, s, e) if gene else None
        if f is None:
            res.append(("", "", "", ""))
            continue
        aa, stops = product(rec, f)
        bits, qcov = blast(aa, gene)
        res.append((len(aa), stops, bits, qcov))
    (na, ns, nb, nq), (ma, ms, mb, mq) = res
    better = ""
    if nb and mb:
        sc_n, sc_m = float(nb) - 50 * int(ns), float(mb) - 50 * int(ms)
        better = "same" if abs(sc_n - sc_m) < 1 else ("ncbi" if sc_n > sc_m else "mitom")
    w.writerow([r["accession"], r["group"], r["organism"], r["class"], gene, r["truth"], r["pred"],
                na, ns, nb, nq, ma, ms, mb, mq, better])
