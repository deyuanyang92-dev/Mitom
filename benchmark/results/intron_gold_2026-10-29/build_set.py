"""All 193 local NCBI records with an intron (2026-10-29), for verifying NCBI's intron annotation (gold standard) and
for testing MitoM's intron annotation. Local files only, no download.

Candidates: intron-ref's site-library scan of the local NCBI data (bench gb/ folders, refpool_metazoa/Metazoa.gb,
annelida/fetch/Annelida.gb, MitoM tests/data); an INSDC record that an NC_ record copies is dropped (the NC_ is kept).
Writes ncbi_gb/<acc>.gb (NCBI record), bare/<acc>/<acc>.gbf (sequence only, MitoM input), set.tsv (mito_bench_run
columns) and records.tsv (phylum, order, organism, genes with introns, dev = used while developing the intron code,
fragment = not a complete mitogenome).
"""
import csv
import io
import json
import sys
from pathlib import Path

from Bio import SeqIO

R = Path(__file__).resolve().parent
sys.path.insert(0, "/media/deyuan/18T_B/yzz/MitoM/benchmark")
from mito_sample import COLUMNS, MITOZ_CLADE, annotation_method, evidence_tier, majority_code, structure_flags, write_bare  # noqa: E402

S = Path("/tmp/claude-1000/-media-deyuan-18T-B-yzz/6bcce157-6435-4fde-a47f-55b33ff6135c/scratchpad")
W = Path("/media/deyuan/18T_B/yzz/mitom_work")
recs = json.loads((S / "intron_recs.json").read_text())
dev = {r["acc"] for r in csv.DictReader(open(S / "intron-algo/sets.tsv"), delimiter="\t") if r["group"].startswith("pos_")}
for d in ("ncbi_gb", "bare"):
    (R / d).mkdir(exist_ok=True)


def keep(r) -> None:
    SeqIO.write(r, str(R / "ncbi_gb" / f"{r.id}.gb"), "genbank")
    write_bare(r, R / "bare")


need = set(recs)
for src in [W / s / "gb" for s in ("bench_large", "bench_metazoa", "bench_annelida", "bench_new", "bench_trna", "bench_vert")]:
    for a in sorted(need):
        if (src / f"{a}.gb").exists():
            keep(SeqIO.read(str(src / f"{a}.gb"), "genbank"))
            need.discard(a)
for p in Path("/media/deyuan/18T_B/yzz/MitoM/tests/data").glob("*.gb"):
    r = SeqIO.read(str(p), "genbank")
    if r.id in need:
        keep(r)
        need.discard(r.id)
for pool in (W / "annelida/fetch/Annelida.gb", W / "refpool_metazoa/Metazoa.gb"):
    if need:
        for r in SeqIO.parse(str(pool), "genbank"):
            if r.id in need:
                keep(r)
                need.discard(r.id)
if need:
    sys.exit(f"not found locally: {sorted(need)}")

rows, info = [], []
for a in sorted(recs):
    r = SeqIO.read(str(R / "ncbi_gb" / f"{a}.gb"), "genbank")
    lin = r.annotations.get("taxonomy", [])
    phylum = next((x for x in lin if x in ("Annelida", "Cnidaria", "Placozoa", "Porifera")), "?")
    code, assumed = majority_code(r, phylum)
    fragment = "complete genome" not in r.description and len(r) < 10000
    rows.append({"accession": a, "source_accession": "", "phylum": phylum, "class": "", "order": recs[a]["order"],
                 "family": "", "genus": lin[-1] if lin else "", "species": r.annotations.get("organism", ""),
                 "genetic_code": code, "code_assumed": "yes" if assumed else "",
                 "topology": r.annotations.get("topology", ""), "length": len(r), "evidence_tier": evidence_tier(r),
                 "annotation_method": annotation_method(r, None), "refseq_created": r.annotations.get("date", ""),
                 "structure_flags": ",".join(structure_flags(r)), "mitoz_clade": MITOZ_CLADE.get(phylum) or "",
                 "split": "test", "retrieved": "local"})
    info.append({"accession": a, "phylum": phylum, "order": recs[a]["order"], "organism": r.annotations.get("organism", ""),
                 "length": len(r), "topology": r.annotations.get("topology", ""), "genetic_code": code,
                 "dev": "yes" if a in dev else "no", "fragment": "yes" if fragment else "no",
                 "intron_genes": ",".join(f"{g}x{n}" for g, n in sorted(recs[a]["genes"].items()))})
for name, data, cols in (("set.tsv", rows, COLUMNS), ("records.tsv", info, list(info[0]))):
    with open(R / name, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t")
        w.writeheader()
        w.writerows(data)
print(f"{len(rows)} records -> {R}", file=sys.stderr)
