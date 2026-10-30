"""mitom_topology.py - circular or linear (mitom run / annotate --topology; default auto).

Every record is annotated as its label says: the GenBank LOCUS line, the FASTA header's topology=circular|linear,
else circular for one sequence of >= 10 kb and linear for several sequences in a file or a shorter one
(mitom_inputs.split_or_sanitize_fastas). With --topology auto the core also checks the two ends (mito_core.py
--check_ends: each protein gene at an end aligned to its reference protein across the junction, <id>.ends.tsv), and a
label the sequence contradicts is changed:
  - a gene runs from the end across to the start (its reference aligned >= 30 nt on both sides of the junction):
    circular - a circular genome labelled linear, or an assembly started inside a gene in a multi-sequence FASTA;
  - a circular-labelled record whose gene stops at an end with >= 30 residues of its protein missing, not found at
    the other end: linear - a partial sequence;
  - jellyfish / hydroids and calcareous sponges (linear mitochondrial genomes, mito_taxa.linear_mtdna) labelled
    circular, no gene across the ends: linear.
A changed record is annotated again with its new topology (02.anno_topology/); every record's label, topology and
reason: reports/topology_decision.tsv. --topology label keeps every label; circular / linear set every record
(mitom_inputs)."""
from __future__ import annotations

import argparse
import csv
import dataclasses
import logging
from pathlib import Path
from typing import Callable, Dict, List, Tuple

from Bio import SeqIO

from mito_taxa import linear_mtdna
from mitom_inputs import SampleInfo, write_one_record_fasta
from mitom_log import record_event
from mitom_stage import run_annotation_stage, sample_lineage

TOPOLOGY_PHASE = "02.anno_topology"
DECISION_HEADER = "sample_name\tlength\tlabel\tlabel_from\ttopology\tchanged\treason\n"


def read_ends(workdir: Path, internal_id: str) -> List[Dict[str, str]]:
    """The core's <internal_id>.ends.tsv rows (mito_core.end_evidence); [] when it wrote none."""
    path = workdir / f"{internal_id}.ends.tsv"
    if not path.is_file():
        return []
    with path.open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def decide(label: str, ends: List[Dict[str, str]], lineage: str = "") -> Tuple[str, str]:
    """(topology, reason) of one record from its label, the core's evidence at the ends and its NCBI lineage."""
    joined = next((r for r in ends if r["verdict"] == "joined"), None)
    if joined:
        return "circular", (f"{joined['gene']} runs across the sequence ends (its reference protein aligns "
                            f"{joined['before_nt']} nt before and {joined['after_nt']} nt after the junction)")
    cut = next((r for r in ends if r["verdict"] == "cut"), None)
    if label == "circular" and cut:
        return "linear", (f"{cut['gene']} stops at the junction of the sequence ends with residues {cut['residues']} "
                          f"of its protein, the rest is not on the other side: the ends do not join")
    taxon = linear_mtdna(lineage)
    if label == "circular" and taxon:
        return "linear", f"{taxon}: linear mitochondrial genome, no gene across the ends"
    return label, "as labelled"


def auto_mode(args: argparse.Namespace) -> bool:
    """--topology auto with MitoM's core, whose evidence it needs (MitoZ / MITOS2 / imported annotations keep their
    labels)."""
    return getattr(args, "topology", "auto") == "auto" and getattr(args, "engine", "native") == "native"


def relabel(si: SampleInfo, topology: str) -> SampleInfo:
    """si with `topology` and an internal FASTA that says so (<internal_id>.<topology>.fasta beside the first)."""
    new_fa = si.internal_fasta.with_name(f"{si.internal_id}.{topology}.fasta")
    write_one_record_fasta(SeqIO.read(str(si.internal_fasta), "fasta"), new_fa, si.internal_id, topology)
    return dataclasses.replace(si, topology=topology, internal_fasta=new_fa.resolve())


def linear_lineages_first(samples: List[SampleInfo], args: argparse.Namespace) -> List[SampleInfo]:
    """--topology auto before round 1: a circular-labelled record of a lineage with a linear mitochondrial genome
    (mito_taxa.linear_mtdna, from its GenBank taxonomy) is annotated as linear from the start; apply_topology
    annotates it again only when a gene runs across its ends."""
    if not auto_mode(args):
        return samples
    return [relabel(si, "linear") if si.topology == "circular" and linear_mtdna(sample_lineage(si)) else si
            for si in samples]


def apply_topology(samples: List[SampleInfo], results: Dict[str, Dict], args: argparse.Namespace, out_root: Path,
                   reports: Path, info: Callable[[str], None] = logging.info) \
        -> Tuple[List[SampleInfo], Dict[str, Dict]]:
    """--topology auto after round 1: decide each annotated record's topology (decide, from its label), annotate the
    ones whose round-1 topology differs again, and write reports/topology_decision.tsv. Returns the samples (those
    annotated again with their new topology and FASTA) and the round-1 results (those: the new annotation)."""
    auto = auto_mode(args)
    decided: Dict[str, Tuple[str, str]] = {}
    changed: List[SampleInfo] = []
    for si in samples:
        r = results.get(si.internal_id, {})
        if not r.get("ok"):
            continue
        label = si.label or si.topology
        if auto:
            decided[si.internal_id] = decide(label, read_ends(Path(r["workdir"]), si.internal_id), sample_lineage(si))
        else:
            decided[si.internal_id] = (si.topology, f"set by --topology {si.topology}" if si.topology != label
                                       else "as labelled")
        topology, why = decided[si.internal_id]
        if topology != si.topology:
            changed.append(relabel(si, topology))
            logging.info(f"[TOPOLOGY] {si.sample_name}: annotated again as {topology} - {why}")
    new: Dict[str, SampleInfo] = {}
    if changed:
        redo = run_annotation_stage(samples=changed, phase_dirname=TOPOLOGY_PHASE, out_root=out_root, args=args,
                                    fasta_selector=lambda s: s.internal_fasta)
        new = {si.internal_id: si for si in changed if redo.get(si.internal_id, {}).get("ok")}
        results = {iid: redo[iid] if iid in new else r for iid, r in results.items()}
    rows = []
    for si in samples:
        if si.internal_id not in decided:
            continue
        topology, why = decided[si.internal_id]
        if topology != si.topology and si.internal_id not in new:
            record_event(si.sample_name, f"annotation as {topology} failed: the record stays {si.topology}",
                         logging.WARNING)
            topology, why = si.topology, f"{why}; annotation as {topology} failed, kept {si.topology}"
        label = si.label or si.topology
        length = sum(len(x) for x in SeqIO.parse(str(si.internal_fasta), "fasta"))
        rows.append(f"{si.sample_name}\t{length}\t{label}\t{si.topology_source}\t{topology}\t"
                    f"{'yes' if topology != label else 'no'}\t{why}\n")
    reports.mkdir(parents=True, exist_ok=True)
    (reports / "topology_decision.tsv").write_text(DECISION_HEADER + "".join(rows), encoding="utf-8")
    if auto:
        n_lin = sum(1 for si in samples if si.internal_id in decided and decided[si.internal_id][0] != (si.label or si.topology))
        info(f"topology: {n_lin} record(s) not as labelled ({len(new)} annotated again), the others as labelled: "
             "reports/topology_decision.tsv")
    return [new.get(si.internal_id, si) for si in samples], results
