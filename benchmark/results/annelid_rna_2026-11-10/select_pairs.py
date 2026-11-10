#!/usr/bin/env python3
"""Pair each annelid genome-project mitogenome with RNA-seq of the same genome project (ncbi-genome-inventory links).

Input: the ncbi-genome-inventory output folder (<taxon>_mitogenomes.tsv, <taxon>_runs.tsv). Only the links that
tool derives are used - no run is paired by species name alone:
  mitogenome -> assemblies: relation_derived in_assembly:<GCA> / same_biosample:<GCA>
  RNA-Seq run -> assemblies: linked_assemblies_derived <GCA>:biosample | bioproject | sibling_project, and
                 specimen_derived ("same specimen: ..." / "different specimen: ToLID a vs b" / "unknown").
Levels (best first): 1 same_specimen, 2 same_project_other_specimen, 3 same_project_specimen_unknown.
Per mitogenome the run chosen: best level, then poly(A)-selected libraries (3' ends), then paired, then the most Gb.

Writes datasets.tsv: every accession of a pair (mitogenome, its BioSample, assembly, genome BioProject; run,
experiment, RNA BioProject, BioSample, library, size) plus the link level and how it was derived; candidates.tsv:
all linked RNA runs per mitogenome. After `mitom rna`, check_species (summarise.py) adds the reads' edit rate:
> 2 % within one genome project = not the same species: contamination or a mislabelled sample, flagged.

Usage: select_pairs.py <inventory dir> <taxon> -o <dir>
"""
import argparse
import csv
from pathlib import Path

LEVELS = {1: "same_specimen", 2: "same_project_other_specimen", 3: "same_project_specimen_unknown"}
COLS = ["mito", "organism", "mito_tax_id", "mito_length", "mito_biosample", "mito_relation", "assembly", "run",
        "experiment", "rna_bioproject", "rna_biosample", "rna_organism", "link", "link_level", "specimen",
        "library_selection", "library_layout", "platform", "bases_Gb", "avg_length", "rna_release", "n_linked_runs"]


def level(link_kind, specimen):
    if specimen.startswith("same specimen") or link_kind == "biosample":
        return 1
    if specimen.startswith("different specimen"):
        return 2
    return 3


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("inventory")
    p.add_argument("taxon")
    p.add_argument("-o", "--out", required=True)
    a = p.parse_args()
    inv = Path(a.inventory)
    mitos = list(csv.DictReader(open(inv / f"{a.taxon}_mitogenomes.tsv"), delimiter="\t"))
    runs = [r for r in csv.DictReader(open(inv / f"{a.taxon}_runs.tsv"), delimiter="\t")
            if r["LibraryStrategy"] == "RNA-Seq"]
    by_gca = {}
    for r in runs:
        for x in filter(None, r["linked_assemblies_derived"].split(";")):
            gca, kind = x.rsplit(":", 1)
            spec = next((s.split(": ", 1)[1] for s in r["specimen_derived"].split(";")
                         if s.startswith(gca + ":")), "unknown")
            by_gca.setdefault(gca, []).append((r, kind, spec.strip()))
    rows, cands = [], []
    for m in mitos:
        gcas = [x.split(":", 1)[1] for x in m["relation_derived"].split(";")
                if x.startswith(("in_assembly:", "same_biosample:"))]
        linked = {}
        for g in gcas:
            for r, kind, spec in by_gca.get(g, []):
                lv = level(kind, spec)
                if r["Run"] not in linked or lv < linked[r["Run"]][3]:
                    linked[r["Run"]] = (r, g, kind, lv, spec)
        if not linked:
            continue
        for r, g, kind, lv, spec in linked.values():
            cands.append({"mito": m["accession"], "run": r["Run"], "assembly": g, "link": kind,
                          "link_level": LEVELS[lv], "specimen": spec, "library_selection": r["LibrarySelection"],
                          "bases_Gb": r["bases_Gb"]})
        r, g, kind, lv, spec = min(linked.values(), key=lambda t: (
            t[3], t[0]["LibrarySelection"] != "PolyA", t[0]["LibraryLayout"] != "PAIRED", -float(t[0]["bases_Gb"] or 0)))
        rows.append({"mito": m["accession"], "organism": m["organism"], "mito_tax_id": m["tax_id"],
                     "mito_length": m["length"], "mito_biosample": m["biosample"],
                     "mito_relation": m["relation_derived"], "assembly": g, "run": r["Run"],
                     "experiment": r["Experiment"], "rna_bioproject": r["BioProject"], "rna_biosample": r["BioSample"],
                     "rna_organism": r["ScientificName"], "link": kind, "link_level": LEVELS[lv], "specimen": spec,
                     "library_selection": r["LibrarySelection"], "library_layout": r["LibraryLayout"],
                     "platform": r["Platform"], "bases_Gb": r["bases_Gb"], "avg_length": r["avgLength"],
                     "rna_release": r["ReleaseDate"], "n_linked_runs": len(linked)})
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rows.sort(key=lambda x: (x["link_level"], x["organism"]))
    with open(out / "datasets.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with open(out / "candidates.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, list(cands[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(cands)
    from collections import Counter
    print(f"{len(rows)} mitogenomes with RNA-seq of their genome project: "
          + ", ".join(f"{k} {v}" for k, v in sorted(Counter(r['link_level'] for r in rows).items())))


if __name__ == "__main__":
    main()
