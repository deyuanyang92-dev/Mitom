#!/usr/bin/env python3
"""One row per mitogenome / RNA-seq pair: every accession (datasets.tsv), where the reads came from, the species
check and what the RNA-seq showed (mito_rna outputs + --refine).

species_check: the mapped reads' edit rate against the mitogenome. Both come from one genome project, so they should
be one species: <= 2 % ok; > 2 % FLAG (contamination or a mislabelled RNA-seq sample; its expression and 3' ends are
not used). Fewer than MIN_DEPTH protein-gene median depth: too few mitochondrial reads to judge.

Usage: summarise.py <work dir> -o results.tsv
"""
import argparse
import csv
import re
import statistics
from pathlib import Path

DIVERGED = 0.02
MIN_DEPTH = 10
EXTRA = ["reads_source", "reads_used_mapq20", "edit_rate", "species_check", "library", "pcg_median_depth",
         "polya_reads", "polya_sites", "internal_priming_sites", "pcg_3prime_ends_confirmed", "end_corrections",
         "orfs", "short_genes_relative_depth"]
LOG = re.compile(r"\] (\d+) reads \(MAPQ.*?\), (.*?) \(sense share [\d.]+\), protein-gene median depth ([\d.]+), "
                 r"edit rate ([\d.]+), (\d+) poly\(A\) reads, (\d+) sites \(\+ (\d+) internal")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("work")
    p.add_argument("-o", "--out", required=True)
    a = p.parse_args()
    w = Path(a.work)
    pairs = list(csv.DictReader(open(w / "pairs/datasets.tsv"), delimiter="\t"))
    rows = []
    for d in pairs:
        acc = d["mito"].split(".")[0]
        o = w / "rna_out" / acc
        row = dict(d)
        run = d["run"]
        local = any((w / "rna_local" / f"{run}{x}").exists() for x in ("_1.fastq.gz", ".fastq.gz"))
        ncbi = any((w / "rna_sra/fastq" / f"{run}{x}").exists() for x in ("_1.fastq.gz", ".fastq.gz"))
        row["reads_source"] = ("LAN server copy, first 8 M pairs" if local else
                               "NCBI download (organise-sra), whole run" if ncbi else "")
        log = (w / "rna_out" / f"{acc}.log")
        lm = LOG.search(log.read_text()) if log.exists() else None
        if not lm or not (o / "rna_features.tsv").exists():
            row["species_check"] = "not run (reads not available yet)"
            rows.append(row)
            continue
        used, lib, med, er, pa, ns, ip = lm.groups()
        er, med = float(er), float(med)
        row.update(reads_used_mapq20=used, edit_rate=f"{100 * er:.2f} %", library=lib, pcg_median_depth=f"{med:.0f}",
                   polya_reads=pa, polya_sites=ns, internal_priming_sites=ip)
        row["species_check"] = ("too few mitochondrial reads" if med < MIN_DEPTH else
                                f"FLAG: reads differ {100 * er:.1f} % from the mitogenome of the same genome project - "
                                "contamination or mislabelled RNA-seq sample" if er > DIVERGED else "ok")
        feats = list(csv.DictReader(open(o / "rna_features.tsv"), delimiter="\t"))
        pcg = [f for f in feats if f["type"] == "CDS" and not f["feature"].lower().startswith("orf")]
        gb = o / "rna_refined.gb"
        if gb.exists():
            txt = gb.read_text()
            row["pcg_3prime_ends_confirmed"] = f"{txt.count('confirmed by')}/{len(pcg)}"
        corr = list(csv.DictReader(open(o / "rna_end_corrections.tsv"), delimiter="\t")) \
            if (o / "rna_end_corrections.tsv").exists() else []
        row["end_corrections"] = "; ".join(f"{c['gene']} {c['old_aa']}->{c['new_aa']} codons ({c['site_reads']} reads)"
                                           for c in corr)
        row["orfs"] = "; ".join(f"{f['feature']} depth {f['relative_depth']} covered {f['covered_share']}"
                                for f in feats if f["feature"].lower().startswith("orf"))
        short = {f["feature"].upper(): f["relative_depth"] for f in pcg if f["feature"].upper() in ("ATP8", "ND4L", "ND3")}
        row["short_genes_relative_depth"] = ", ".join(f"{k} {v}" for k, v in sorted(short.items()))
        rows.append(row)
    cols = list(pairs[0]) + EXTRA
    with open(a.out, "w", newline="") as fh:
        wr = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        wr.writeheader()
        wr.writerows(rows)
    done = [r for r in rows if r.get("edit_rate")]
    print(f"{len(rows)} pairs, {len(done)} analysed; species check: "
          + ", ".join(f"{k} {sum(r['species_check'].split(':')[0] == k for r in done)}"
                      for k in ("ok", "FLAG", "too few mitochondrial reads")))
    if done:
        print("median edit rate:", f"{statistics.median(float(r['edit_rate'].split()[0]) for r in done):.2f} %")


if __name__ == "__main__":
    main()
