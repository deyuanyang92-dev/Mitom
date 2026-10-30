"""Stress cases for GenBank -> NCBI feature table converters (2026-10-30). in/<case>.gb"""
import copy
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import (AfterPosition, BeforePosition, CompoundLocation, ExactPosition, FeatureLocation,
                            SeqFeature)
from Bio.SeqRecord import SeqRecord

R = Path(__file__).resolve().parent
OUT = R / "in"
OUT.mkdir(exist_ok=True)
DATA = Path("/media/deyuan/18T_B/yzz/MitoM/tests/data")
ann = SeqIO.read(str(DATA / "NC_000931.gb"), "genbank")
intr = SeqIO.read(str(DATA / "NC_000933.gb"), "genbank")
vert = SeqIO.read(str(next(Path("/media/deyuan/18T_B/yzz/mitom_work").glob("*/gb/NC_009732.1.gb"))), "genbank")


def write(name, rec):
    SeqIO.write(rec, str(OUT / f"{name}.gb"), "genbank")


import re


def shift_pos(q, shift, n):
    """/transl_except and /anticodon positions moved by shift (mod n when n): (pos:12..14,...) / complement(...)."""
    out = copy.deepcopy(q)
    for key in ("transl_except", "anticodon"):
        vals = []
        for v in out.get(key, []):
            def mv(m):
                a, b = int(m.group(1)) + shift, int(m.group(2)) + shift if m.group(2) else None
                if n:
                    a, b = (a - 1) % n + 1, None if b is None else (b - 1) % n + 1
                return f"{a}..{b}" if b is not None else f"{a}"
            vals.append(re.sub(r"(\d+)(?:\.\.(\d+))?", mv, v.split(",aa:")[0]) + ",aa:" + v.split(",aa:")[1]
                        if ",aa:" in v else v)
        if key in out:
            out[key] = vals
    return out


def rotate(rec, k):
    """Circular record started at 0-based position k; features moved, those across the new origin joined."""
    n = len(rec)
    new = SeqRecord(rec.seq[k:] + rec.seq[:k], id=rec.id, name=rec.name, description=rec.description,
                    annotations=dict(rec.annotations))
    for f in rec.features:
        if f.type == "source":
            new.features.append(SeqFeature(FeatureLocation(0, n, 1), type="source", qualifiers=f.qualifiers))
            continue
        parts = []
        for p in f.location.parts:
            s, e = (int(p.start) - k) % n, (int(p.end) - k) % n or n
            if s < e:
                parts.append(FeatureLocation(s, e, p.strand))
            else:                                                   # across the new origin
                pp = [FeatureLocation(s, n, p.strand), FeatureLocation(0, e, p.strand)]
                parts += pp if p.strand != -1 else pp[::-1]
        new.features.append(SeqFeature(parts[0] if len(parts) == 1 else CompoundLocation(parts), type=f.type,
                                       qualifiers=shift_pos(f.qualifiers, -k, n)))
    return new


def feat(rec, typ, gene):
    return next(f for f in rec.features if f.type == typ and gene in f.qualifiers.get("gene", [""])[0].upper())


write("c01_annelid_plain", ann)
write("c02_introns", intr)
write("c03_vertebrate_minus_strand", vert)
cox1 = feat(ann, "CDS", "COX1")
write("c04_origin_plus_strand", rotate(ann, int(cox1.location.start) + 700))
nd6 = feat(vert, "CDS", "ND6")
write("c05_origin_minus_strand", rotate(vert, (int(nd6.location.start) + int(nd6.location.end)) // 2))
# c06: a partial sequence, genes cut at both ends (5' partial with codon_start, 3' partial)
a, b = int(cox1.location.start) + 301, int(feat(ann, "CDS", "ND5").location.start) + 500
part = SeqRecord(ann.seq[a:b], id="PARTIAL1", name="PARTIAL1", description="partial", annotations=dict(ann.annotations))
part.annotations["topology"] = "linear"
for f in ann.features:
    s, e = int(f.location.start), int(f.location.end)
    if f.type == "source":
        part.features.append(SeqFeature(FeatureLocation(0, b - a, 1), type="source", qualifiers=f.qualifiers))
    elif len(f.location.parts) == 1 and e > a and s < b:
        ns = BeforePosition(0) if s < a else ExactPosition(s - a)
        ne = AfterPosition(b - a) if e > b else ExactPosition(e - a)
        q = shift_pos(f.qualifiers, -a, 0)
        if f.type == "CDS" and s < a and f.location.strand == 1:
            q["codon_start"] = [str((s - a) % 3 + 1)]
        q.pop("translation", None)
        part.features.append(SeqFeature(FeatureLocation(ns, ne, f.location.strand), type=f.type, qualifiers=q))
write("c06_partial_ends", part)
# c07: qualifiers - pseudo, exception, non-ASCII, quotes, several notes
q = copy.deepcopy(ann)
g = feat(q, "CDS", "ATP8")
g.qualifiers.update(pseudo=[""], note=['frameshift; "check" later', "collected by Müller 测试 µm"])
g.qualifiers.pop("translation", None)
feat(q, "gene", "ATP8").qualifiers["pseudo"] = [""]
feat(q, "CDS", "COX2").qualifiers["exception"] = ["RNA editing"]
next(f for f in q.features if f.type == "source").qualifiers.update(collected_by=["J. Müller"], country=["Germany: Köln"])
write("c07_qualifiers", q)
# c08: exported by Geneious / SnapGene / ApE: their own qualifiers and a non-INSDC feature type
j = copy.deepcopy(ann)
for f in j.features[1:6]:
    f.qualifiers.update(created_by=["User"], modified_by=["User"], label=[f.qualifiers.get("gene", ["x"])[0]],
                        ApEinfo_fwdcolor=["#85dae9"], Name=["thing"])
j.features.append(SeqFeature(FeatureLocation(8000, 8300, 1), type="ORF", qualifiers={"note": ["Geneious type: ORF"]}))
write("c08_editor_junk", j)
# c09: an insertion site between two bases
k = copy.deepcopy(ann)
k.features.append(SeqFeature(FeatureLocation(100, 100, 1), type="misc_feature", qualifiers={"note": ["site 100^101"]}))
write("c09_zero_length_site", k)
# c10: no source feature, no organism
m = copy.deepcopy(ann)
m.features = [f for f in m.features if f.type != "source"]
m.annotations["organism"] = "."
m.annotations.pop("taxonomy", None)
write("c10_no_source_no_organism", m)
# c13: a feature beyond the end of the sequence
o = copy.deepcopy(ann)
o.features.append(SeqFeature(FeatureLocation(len(o) - 50, len(o) + 40, 1), type="misc_feature", qualifiers={"note": ["beyond the end"]}))
try:
    write("c13_beyond_end", o)
except Exception as e:
    print("c13 not writable by Biopython:", e)
# c14: an ID with spaces and slashes, lowercase sequence with N and IUPAC codes
p = copy.deepcopy(ann)
p.id = p.name = "MySample#1/A"
s = str(p.seq).lower()
p.seq = Seq(s[:500] + "nnnnnnnnnn" + s[510:2000] + "rykm" + s[2004:])
write("c14_odd_id_lowercase_iupac", p)
