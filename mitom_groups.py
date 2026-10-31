"""mitom_groups.py - --clade auto per record: the phylum and genetic code of every input record (GenBank: its taxonomy
and CDS /transl_table; no taxonomy, FASTA: mito_classify.py, blastx to the reference proteins), and a run over
several phyla / codes as one run per group, merged into one output folder.

  decide(args)            -> Plan: one row per record (reports/clade_decision.tsv) and the groups
  run_groups(args, plan)  -> exit code: input split per group (out/_groups/<group>/input), `mitom run` per group
                             (in parallel, --jobs shared), the outputs merged into out/ (final_annotation/,
                             summary.tsv, ncbi_submission/, plots/, reports/ ..., run.log)
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import logging
import os
import re
import shutil
import subprocess
import sys
import time
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO

from mito_common import safe_name
from mito_taxa import METAZOAN_PHYLA, PHYLUM_CLADE, PHYLUM_CODE, PHYLUM_SYNONYMS
from mitom_inputs import _is_fasta_file, _is_gb_file, _sniff_input_type

HERE = Path(__file__).resolve().parent
FALLBACK = ("Annelida-segmented-worms", 5, "")     # nothing could be decided (the default before --clade auto)
DECISION_HEADER = ("record\tfile\tsource\tphylum\tclade\tgenetic_code\tcode_from\tsupport\tcodon_check\tclosest\t"
                   "group\tnote")


@dataclass
class Row:
    record: str
    file: str
    source: str = ""             # GenBank taxonomy / blastx (mito_classify) / undecided
    phylum: str = ""
    clade: str = ""
    code: int = 0
    code_from: str = ""
    support: str = ""
    codon_check: str = ""
    closest: str = ""
    note: str = ""
    group: str = ""

    def tsv(self) -> str:
        return "\t".join([self.record, self.file, self.source, self.phylum or "-", self.clade, str(self.code),
                          self.code_from, self.support, self.codon_check, self.closest, self.group, self.note])


@dataclass
class Plan:
    rows: List[Row] = field(default_factory=list)
    error: str = ""
    warnings: List[str] = field(default_factory=list)

    def groups(self) -> Dict[Tuple[str, int, str], List[Row]]:
        out: Dict[Tuple[str, int, str], List[Row]] = {}
        for r in self.rows:
            out.setdefault((r.clade, r.code, r.phylum), []).append(r)
        return dict(sorted(out.items(), key=lambda kv: -len(kv[1])))


LAST_PLAN: Optional[Plan] = None      # the plan of this run (workflow_run writes reports/clade_decision.tsv)


def group_label(key: Tuple[str, int, str]) -> str:
    clade, code, phylum = key
    return safe_name(f"{phylum or clade.split('-')[0]}_gc{code}")


def _phylum(tax: List[str]) -> str:
    return next((t for t in tax if t in METAZOAN_PHYLA), "") or next((PHYLUM_SYNONYMS[t] for t in tax
                                                                      if t in PHYLUM_SYNONYMS), "")


def input_files(in_path: Path) -> Tuple[str, List[Path]]:
    """('genbank' | 'fasta', the files the run reads) - as iter_inputs: a folder with both kinds uses the GenBank."""
    kind = _sniff_input_type(in_path)
    if in_path.is_file():
        return ("fasta" if kind == "fasta" else "genbank"), [in_path.resolve()]
    pick = _is_fasta_file if kind == "fasta" else _is_gb_file
    return ("fasta" if kind == "fasta" else "genbank"), sorted(p.resolve() for p in in_path.rglob("*")
                                                                if p.is_file() and pick(p))


def decide(args: argparse.Namespace) -> Plan:
    """The phylum, clade and genetic code of every input record (--clade auto). A taxonomy that is not animal is an
    error; an explicit --genetic_code is used for every record."""
    import mito_classify
    plan = Plan()
    if not Path(args.input).exists():
        plan.error = f"input not found: {args.input}"
        return plan
    kind, files = input_files(Path(args.input))
    if not files:
        plan.error = f"no GenBank or FASTA file in {args.input}"
        return plan
    todo: List[Tuple[Row, str]] = []                      # rows to classify by blastx, with the sequence
    for f in files:
        try:
            recs = list(SeqIO.parse(str(f), "genbank" if kind == "genbank" else "fasta"))
        except Exception as exc:                          # noqa: BLE001  (named in the error)
            plan.error = f"--clade auto: cannot read {f}: {exc}; fix the file or give --clade and --genetic_code"
            return plan
        for rec in recs:
            row = Row(rec.id, f.name)
            tax = [str(t) for t in rec.annotations.get("taxonomy", []) if t] if kind == "genbank" else []
            if tax:
                ph = _phylum(tax)
                if not ph and "Metazoa" not in tax:
                    plan.error = (f"--clade auto: {rec.id}: the taxonomy names no animal phylum (lineage ends "
                                  f"'{'; '.join(tax[-3:])}'); MitoM annotates animal mitochondrial genomes - give "
                                  f"--clade and --genetic_code to annotate it anyway")
                    return plan
                if ph:
                    row.source, row.phylum = "GenBank taxonomy", ph
                    codes = [str(x.qualifiers["transl_table"][0]) for x in rec.features
                             if x.type == "CDS" and "transl_table" in x.qualifiers]
                    if codes and Counter(codes).most_common(1)[0][0].isdigit():
                        row.code, row.code_from = int(Counter(codes).most_common(1)[0][0]), "CDS transl_table"
                    else:
                        row.code = mito_classify.code_of_lineage(";".join(tax)) or PHYLUM_CODE.get(ph, 5)
                        row.code_from = "lineage (reference genomes of the same group)"
                else:
                    row.note = f"taxonomy without a known phylum ('{'; '.join(tax[-2:])}')"
            if not row.phylum:
                todo.append((row, str(rec.seq)))
            plan.rows.append(row)
    if todo:
        msg = (f"[CLADE] {len(todo)} record(s) without a taxonomy: phylum and genetic code by blastx to the reference "
               f"proteins (mito_classify.py, a few seconds per record) ...")
        logging.info(msg)
        print(msg, file=sys.stderr, flush=True)
        try:
            got = mito_classify.classify_records([(r.record, s) for r, s in todo], getattr(args, "mitos_refdir", None),
                                                 threads=max(1, int(args.threads) * int(args.max_tasks)))
        except (FileNotFoundError, RuntimeError) as exc:
            got = []
            plan.warnings.append(f"no blastx classification ({exc})")
        for (row, _), c in zip(todo, got):
            row.support, row.codon_check, row.closest = str(c.support), c.codon_check, c.closest
            if c.genetic_code:                            # phylum '' = uncertain: clade Animal, code from codons
                row.source, row.phylum, row.code = "blastx (mito_classify)", c.phylum, c.genetic_code
                row.code_from, row.note = c.code_from, c.note
                if not c.phylum:
                    plan.warnings.append(f"{row.record}: {c.note}: --clade Animal --genetic_code {c.genetic_code} "
                                         f"({c.code_from}); give --clade if you know the group")
                elif c.codon_check.startswith("conflict"):
                    plan.warnings.append(f"{row.record}: {c.note} ({c.codon_check}) - check, or give --genetic_code")
            else:
                row.note = "; ".join(x for x in (row.note, c.note) if x)
    undecided = [r for r in plan.rows if not r.code]
    decided = [r for r in plan.rows if r.code]
    for r in decided:
        r.clade = PHYLUM_CLADE.get(r.phylum, "Animal")
    if undecided:
        if decided:                                       # short contigs of a sample: with the most common group
            ph, clade, code = Counter((r.phylum, r.clade, r.code) for r in decided).most_common(1)[0][0]
            why = "with the most common group of the input"
        else:
            clade, code, ph = FALLBACK
            why = "the previous default"
        for r in undecided:
            r.source, r.phylum, r.clade, r.code, r.code_from = "undecided", ph, clade, code, why
        plan.warnings.append(f"{len(undecided)} record(s) without a phylum ({undecided[0].note or 'no hits'}"
                             f"{', ...' if len(undecided) > 1 else ''}): --clade {clade} --genetic_code {code} ({why})")
    if getattr(args, "genetic_code", None):
        for r in plan.rows:
            r.code, r.code_from = int(args.genetic_code), "--genetic_code"
    for key, rows in plan.groups().items():
        for r in rows:
            r.group = group_label(key)
    return plan


def write_decision(plan: Plan, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join([DECISION_HEADER] + [r.tsv() for r in plan.rows]) + "\n", encoding="utf-8")


# ── several groups: one run per group, merged ────────────────────────────────────────────────────────────────────────
def _split_input(plan: Plan, args: argparse.Namespace, groups_dir: Path) -> Dict[str, Path]:
    """Each group's records into <groups_dir>/<group>/input/<file name> (GenBank records copied as text)."""
    kind, files = input_files(Path(args.input))
    group_of = [r.group for r in plan.rows]           # decide() read the records in this order
    dirs = {g: groups_dir / g / "input" for g in dict.fromkeys(group_of)}
    for d in dirs.values():
        if d.parent.is_dir():
            shutil.rmtree(d.parent)
        d.mkdir(parents=True)
    k = 0
    for f in files:
        if kind == "genbank":
            text = f.read_text(encoding="utf-8", errors="replace")
            chunks = [c for c in re.split(r"(?m)^//[ \t]*$\n?", text) if c.strip()]
            ids = [r.id for r in SeqIO.parse(str(f), "genbank")]
            if len(chunks) != len(ids):
                raise RuntimeError(f"{f}: {len(ids)} records but {len(chunks)} '//' blocks - cannot split it by group")
            pieces = [c.lstrip("\n") + "//\n" for c in chunks]
        else:
            recs = list(SeqIO.parse(str(f), "fasta"))
            ids = [r.id for r in recs]
            pieces = [f">{r.description}\n{str(r.seq)}\n" for r in recs]
        for rid, piece in zip(ids, pieces):
            g = group_of[k]
            k += 1
            with (dirs[g] / f.name).open("a", encoding="utf-8") as fh:
                fh.write(piece)
    return {g: (d / files[0].name if len(files) == 1 else d) for g, d in dirs.items()}


def _concat_tsv(srcs: List[Path], dst: Path) -> None:
    """Tables with one header line: the header once, then every table's rows."""
    head, rows = None, []
    for s in srcs:
        lines = s.read_text(encoding="utf-8", errors="replace").splitlines()
        if not lines:
            continue
        head = head or lines[0]
        rows += lines[1:] if lines[0] == head else lines
    if head is not None:
        dst.write_text("\n".join([head] + rows) + "\n", encoding="utf-8")


def merge_outputs(out_root: Path, group_outs: Dict[str, Path]) -> None:
    """The outputs of the group runs, merged into out_root in the layout of one run."""
    from mitom_output import OUTPUT_LAYOUT
    outs = [p for p in group_outs.values() if p.is_dir()]
    for name, _ in OUTPUT_LAYOUT:                          # a previous run's results in this folder
        old = out_root / name.rstrip("/")
        if name != "run.log" and old.is_dir() and not old.is_symlink():
            shutil.rmtree(old)
        elif name != "run.log" and old.is_file():
            old.unlink()
    for sub in ("final_annotation", "original_annotation", "ncbi_submission", "plots", "work"):
        for o in outs:
            src = o / sub
            if not src.is_dir():
                continue
            if sub == "work":                              # --keep_intermediate yes: work/<group>/
                shutil.copytree(src, out_root / "work" / next(g for g, x in group_outs.items() if x == o),
                                dirs_exist_ok=True)
                continue
            for p in sorted(src.rglob("*")):
                if p.is_file() and not (sub == "ncbi_submission" and p.suffix == ".tsv" and p.parent == src):
                    t = out_root / sub / p.relative_to(src)
                    t.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(p, t)
    for name in ("submission_summary.tsv", "validation_messages.tsv"):
        srcs = [o / "ncbi_submission" / name for o in outs if (o / "ncbi_submission" / name).is_file()]
        if srcs:
            _concat_tsv(srcs, out_root / "ncbi_submission" / name)
    for name in ("summary.tsv", "original_vs_mitom.tsv"):
        srcs = [o / name for o in outs if (o / name).is_file()]
        if srcs:
            _concat_tsv(srcs, out_root / name)
    with (out_root / "final_annotation_all.gb").open("w", encoding="utf-8") as fh:
        for o in outs:
            if (o / "final_annotation_all.gb").is_file():
                fh.write((o / "final_annotation_all.gb").read_text(encoding="utf-8", errors="replace"))
    rep = out_root / "reports"
    names = sorted({p.relative_to(o / "reports") for o in outs if (o / "reports").is_dir()
                    for p in (o / "reports").rglob("*") if p.is_file()})
    for rel in names:
        srcs = [o / "reports" / rel for o in outs if (o / "reports" / rel).is_file()]
        t = rep / rel
        t.parent.mkdir(parents=True, exist_ok=True)
        if rel.suffix == ".tsv":
            _concat_tsv(srcs, t)
        elif rel.suffix == ".json" or rel.parts[0] == "failed_logs":
            for g, o in group_outs.items():
                if (o / "reports" / rel).is_file():
                    dst = rep / "groups" / f"{g}.{rel.name}" if rel.suffix == ".json" else t
                    dst.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(o / "reports" / rel, dst)
        else:
            t.write_text("".join(s.read_text(encoding="utf-8", errors="replace") for s in srcs), encoding="utf-8")
    first = next((o / "README.txt" for o in outs if (o / "README.txt").is_file()), None)
    if first:
        txt = first.read_text(encoding="utf-8")
        (out_root / "README.txt").write_text(
            txt + "\nThis run had several phyla / genetic codes (--clade auto): each group was annotated with its own "
                  "clade and code (reports/clade_decision.tsv) and the results merged here; the log of each group "
                  "is reports/groups/<group>.run.log (mitom debug reads them all).\n", encoding="utf-8")


def run_groups(args: argparse.Namespace, plan: Plan, argv: List[str], say=print) -> int:
    """`mitom run` once per group (in parallel; --jobs shared by the number of records), then merge_outputs."""
    out_root = Path(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    groups_dir = out_root / "_groups"
    groups = plan.groups()
    inputs = _split_input(plan, args, groups_dir)
    total, jobs = len(plan.rows), int(args.max_tasks)
    say("")
    say(f"MitoM - {total} records of {len(groups)} groups (phylum / genetic code; --clade auto, "
        f"reports/clade_decision.tsv):")
    tasks = []
    for key, rows in groups.items():
        clade, code, phylum = key
        g = group_label(key)
        n_jobs = max(1, min(len(rows), round(jobs * len(rows) / total)))
        cmd = [sys.executable, str(HERE / "mitom_pipeline.py"), *argv, "-i", str(inputs[g]),
               "-o", str(groups_dir / g / "out"), "--clade", clade, "--genetic_code", str(code), "--jobs", str(n_jobs)]
        if phylum:
            cmd += ["--phylum", phylum]
        tasks.append((g, rows, cmd))
        say(f"    {g:<28} {len(rows):>4} record(s)   --clade {clade} --genetic_code {code}   ({n_jobs} at a time)")
    logging.info("[GROUP] " + "; ".join(f"{g}: {len(rows)} records: {' '.join(cmd)}" for g, rows, cmd in tasks))
    say("")

    def one(t):
        g, rows, cmd = t
        t0 = time.time()
        with (groups_dir / g / "console.txt").open("w", encoding="utf-8") as fh:
            rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, env={**os.environ, "MITOM_GROUP_RUN": g}).returncode
        say(f"    [{g}] finished: exit {rc}, {len(rows)} record(s), {time.time() - t0:.0f} s"
            + (f"  - see {groups_dir / g / 'console.txt'}" if rc else ""))
        return g, rc

    with cf.ThreadPoolExecutor(max_workers=len(tasks)) as ex:
        rcs = dict(ex.map(one, tasks))
    outs = {g: groups_dir / g / "out" for g, _, _ in tasks}
    merge_outputs(out_root, outs)
    write_decision(plan, out_root / "reports" / "clade_decision.tsv")
    for g, _, _ in tasks:                    # each group's log: reports/groups/<group>.run.log (`mitom debug` reads them)
        lg = outs[g] / "run.log"
        if lg.is_file():
            (out_root / "reports" / "groups").mkdir(parents=True, exist_ok=True)
            shutil.copy2(lg, out_root / "reports" / "groups" / f"{g}.run.log")
        logging.info(f"[GROUP] {g}: exit {rcs[g]}; log reports/groups/{g}.run.log")
    if getattr(args, "keep_intermediate", "no") != "yes" and all(rc == 0 for rc in rcs.values()):
        shutil.rmtree(groups_dir, ignore_errors=True)
        shutil.rmtree(out_root / "00.input", ignore_errors=True)        # unpacked / converted input (mitom_inputs)
    rc = max(rcs.values()) if rcs else 1
    say(f"\nMitoM finished ({len(groups)} groups merged, exit {rc}): {out_root}")
    if rc:
        failed = [g for g, r in rcs.items() if r]
        say(f"    ! groups with errors: {', '.join(failed)} (their files: {groups_dir}/<group>/console.txt, out/run.log)")
    return rc
