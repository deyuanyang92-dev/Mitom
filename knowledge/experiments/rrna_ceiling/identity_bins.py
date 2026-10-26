"""Agreement of NCBI rRNA ends between a record and its closest annotated relative (other species), by identity.

usage: identity_bins.py proj_species.tsv [...]   (output of project.py; rank-1 reference per truth rRNA)
"""
import csv
import sys

for path in sys.argv[1:]:
    rows = [r for r in csv.DictReader(open(path), delimiter="\t") if r["rank"] == "1" and r["d5"] != "" and r["d3"] != ""]
    for lo, hi in [(0.98, 1.01), (0.97, 0.98), (0.95, 0.97), (0.90, 0.95), (0.85, 0.90), (0.0, 0.85)]:
        sel = [r for r in rows if lo <= float(r["ident"]) < hi]
        if not sel:
            continue
        both = sum(r["d5"] == "0" and r["d3"] == "0" for r in sel)
        e5 = sum(r["d5"] == "0" for r in sel)
        e3 = sum(r["d3"] == "0" for r in sel)
        print(f"{path}\tidentity {lo:.2f}-{min(hi, 1):.2f}\tn={len(sel)}\tboth ends {100 * both / len(sel):.0f}%"
              f"\t5' {100 * e5 / len(sel):.0f}%\t3' {100 * e3 / len(sel):.0f}%")
