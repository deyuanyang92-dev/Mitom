#!/usr/bin/env python3
"""mito_engines.py - MitoZ and MITOS2 as annotation engines of MitoM, and their existing results as MitoM input.

`mitom run --engine mitoz|mitos2` runs that program on each sample instead of MitoM's own core (mito_core.py);
`mitom import` reads results made earlier (MitoZ result folders or .gbf files, MITOS2 result folders, GenBank files
of any annotator). Either way the annotation becomes the record the later MitoM steps read - reorientation, round 2,
rescue, rRNA HMM, introns, curation with NCBI names, NCBI submission files, genome maps - written like the core's
<workdir>/<prefix>.gbf (MITOS2 names translated to the core's: cox1 -> COX1, rrnS -> s-rRNA; tRNAs trnL1 ...).

Neither program is installed with MitoM; mito_env finds them on PATH or in a conda environment, or give the folder:
  MitoZ 3.6      mamba create -n mitoz -c bioconda -c conda-forge mitoz=3.6        (--mitoz_bin <env>/bin)
  MITOS2 2.1.10  mamba create -n mitos -c bioconda -c conda-forge mitos=2.1.10     (--mitos_bin <env>/bin;
                 reference data refseq89m: mitom setup)

Usage (mitom_pipeline.py runs this per sample; the output is <workdir>/<prefix>.gbf):
  mito_engines.py mitoz  -i sample.fasta -o workdir --prefix S00001 --clade Mollusca --genetic_code 5 [--threads 4]
  mito_engines.py mitos2 -i sample.fasta -o workdir --prefix S00001 --genetic_code 5
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_env  # noqa: E402
import mito_names  # noqa: E402
from mito_common import GB_SUFFIXES, contig_record_id, run_tool, topology_from_header  # noqa: E402
from mito_common import read_fasta as _read_fasta  # noqa: E402

FA_SUFFIXES = (".fa", ".fasta", ".fna", ".fas")
MITOZ_INSTALL = "mamba create -n mitoz -c bioconda -c conda-forge mitoz=3.6"
MITOS2_INSTALL = "mamba create -n mitos -c bioconda -c conda-forge mitos=2.1.10"


def read_fasta(fa: Path) -> List[Tuple[str, str, str]]:
    """(id, sequence, topology) of every sequence: topology=linear in the header, else circular (as mito_core.py)."""
    return [(r.id, str(r.seq).upper(), topology_from_header(r.description)) for r in _read_fasta(fa)]


record_id = contig_record_id            # record name as mito_core.py gives it (mito_common)


def engine_env(bin_dir: str, isolated: bool = False) -> Dict[str, str]:
    """Environment for an engine: its own folder first on PATH (it calls its helper programs by name); isolated
    keeps only the system folders behind it (MITOS2 must not meet another BLAST / Infernal)."""
    rest = "/usr/bin:/bin" if isolated else os.environ.get("PATH", "")
    return dict(os.environ, PATH=f"{bin_dir}{os.pathsep}{rest}")


# ── MitoZ ────────────────────────────────────────────────────────────────────
def run_mitoz(fa: Path, workdir: Path, prefix: str, clade: str, code: int, threads: int = 1,
              bin_dir: Optional[str] = None, species: Optional[str] = None) -> List[SeqRecord]:
    """MitoZ `annotate` on one sample (in workdir/mitoz); its records with the topology of the input FASTA (MitoZ
    writes `linear` for every record) and no definition line (MitoZ writes `topology=linear.`)."""
    exe = mito_env.find("mitoz", bin_dir)
    if not exe:
        raise RuntimeError(f"MitoZ not found: install it ({MITOZ_INSTALL}) or give its folder with --mitoz_bin")
    wd = workdir / "mitoz"
    wd.mkdir(parents=True, exist_ok=True)
    cmd = [exe, "annotate", "--workdir", str(wd), "--outprefix", prefix, "--thread_number", str(threads),
           "--clade", clade, "--genetic_code", str(code), "--fastafiles", str(fa)]
    if species:
        cmd += ["--species_name", species]
    print("[engine] " + " ".join(cmd), file=sys.stderr, flush=True)
    rc = run_tool(cmd, cwd=wd, env=engine_env(str(Path(exe).parent)), capture=False).returncode
    if rc:
        raise RuntimeError(f"MitoZ failed (exit {rc})")
    gbfs = sorted(wd.glob("*.result/*.gbf"))
    if not gbfs:
        raise RuntimeError(f"MitoZ wrote no GenBank file ({wd}/*.result/*.gbf)")
    topology = {i: t for i, _, t in read_fasta(fa)}
    recs = [r for g in gbfs for r in SeqIO.parse(str(g), "genbank")]
    for r in recs:
        r.annotations["topology"] = topology.get(r.id, topology.get(r.name, "circular"))
        r.description = ""
    return recs


# ── MITOS2 ───────────────────────────────────────────────────────────────────
RUNMITOS_FAST = Path(__file__).resolve().parent / "scripts" / "runmitos_fast.py"


def mitos2_command(fa: Path, out: Path, code: int, bin_dir: str, refdir: str, ref: str = mito_env.MITOS_REF,
                   mode: str = "full", skip: Sequence[str] = (), linear: bool = False, noplots: bool = True) -> List[str]:
    """The MITOS2 command of every MitoM module (run it with engine_env(bin_dir, isolated=True)): mode full = runmitos
    as installed; sensitive / fast = scripts/runmitos_fast.py (the Infernal searches in parallel, same results; fast
    also skips the sensitive re-search of missing tRNAs); skip = MITOS2 parts not annotated (prot, rrna, trna,
    intron, oril, orih)."""
    exe = [str(Path(bin_dir) / "runmitos")] if mode == "full" else \
        [str(Path(bin_dir) / "python"), str(RUNMITOS_FAST), *(["--no_sensitive"] if mode == "fast" else [])]
    return [*exe, "-i", str(fa), "-c", str(code), "-o", str(out), "-r", ref, "--refdir", str(refdir),
            *[x for part in skip for x in (f"--{part}", "0")], *(["--noplots"] if noplots else []),
            *(["--linear"] if linear else [])]


class BedRow(NamedTuple):
    seq: str          # sequence name
    start: int        # 0-based; start > end for a feature across the origin
    end: int
    name: str         # MITOS2 name: nad5, nad5_1 (2nd copy), cox2-a (fragment), trnL1(tag), OH ...
    score: str
    strand: int


def parse_result_bed(text: str) -> List[BedRow]:
    """The features of a MITOS2 result.bed text (lines with fewer than 6 columns skipped)."""
    rows = []
    for line in text.splitlines():
        c = line.split("\t")
        if len(c) >= 6:
            rows.append(BedRow(c[0], int(c[1]), int(c[2]), c[3], c[4], -1 if c[5].strip() == "-" else 1))
    return rows


def read_result_bed(result_dir: Path) -> List[BedRow]:
    """The features of result_dir/result.bed (MITOS2's result folder)."""
    return parse_result_bed((Path(result_dir) / "result.bed").read_text())


MITOS_TRNA = re.compile(r"^trn([A-Z])(\d?)(?:_(\d+))?(?:-([a-z]))?\(([a-z-]{3})\)", re.I)
MITOS_GENE = re.compile(r"^(.+?)(?:_(\d+))?(?:-([a-z]))?$")       # nad5, nad5_1 (2nd copy), cox2-a (fragment)
MITOS_SKIP = re.compile(r"^(OH|OL)", re.I)
MAX_FRAG_GAP = 3000      # CDS fragments (-a/-b) further apart than this are not joined


def _gap(a: FeatureLocation, b: FeatureLocation) -> int:
    """Nucleotides between two intervals (< 0 when they overlap)."""
    return max(int(b.start) - int(a.end), int(a.start) - int(b.end))


def _chains(frags: List[tuple]) -> List[List[FeatureLocation]]:
    """Fragments (-a, -b ...) in letter order -> chains of parts to join: only parts on the same strand
    with a gap of 0..MAX_FRAG_GAP nt to the previous part are joined; others start a new feature."""
    chains: List[List[FeatureLocation]] = []
    for _, parts in sorted(frags, key=lambda x: x[0]):
        last = chains[-1][-1] if chains else None
        if last is not None and last.strand == parts[0].strand and 0 <= _gap(last, parts[0]) <= MAX_FRAG_GAP:
            chains[-1].extend(parts)
        else:
            chains.append(list(parts))
    return chains


def mitos_to_record(result_dir: Path, seq: Seq, acc: str) -> SeqRecord:
    """GenBank-like record from MITOS2 result.bed (0-based). Copies (_0, _1) stay separate features of the
    same gene; CDS fragments (-a, -b) are joined in letter order when on the same strand and at most
    MAX_FRAG_GAP nt apart (otherwise separate features); rRNA fragments are never joined (rescue then takes
    the longest piece); start > end wraps the origin; OH/OL skipped."""
    n = len(seq)
    groups: Dict[tuple, List[tuple]] = {}
    for i, row in enumerate(read_result_bed(result_dir)):
        if MITOS_SKIP.match(row.name):
            continue
        strand = row.strand
        m = MITOS_TRNA.match(row.name)
        if m:
            key = ("tRNA", f"trn{m.group(1).upper()}{m.group(2)}", m.group(5).lower().strip("-"), m.group(3) or "")
            frag = m.group(4) or ""
        else:
            g = MITOS_GENE.match(row.name)
            base = g.group(1) if g else row.name
            kind = "rRNA" if base.lower().startswith("rrn") else "CDS"
            # rRNA pieces are separate features: the row number makes each key unique
            key = (kind, base, "", f"line{i}" if kind == "rRNA" else (g.group(2) or "") if g else "")
            frag = (g.group(3) or "") if g else ""
        s0, e0 = row.start, row.end
        if s0 <= e0:
            parts = [FeatureLocation(s0, e0, strand)]
        else:                                    # circular feature across the origin
            a, b = FeatureLocation(s0, n, strand), FeatureLocation(0, e0, strand)
            parts = [b, a] if strand == -1 else [a, b]
        groups.setdefault(key, []).append((frag, parts))
    r = SeqRecord(seq, id=acc)
    r.annotations["molecule_type"] = "DNA"
    for (kind, name, ac, _), frags in groups.items():
        for parts in _chains(frags):
            loc = parts[0] if len(parts) == 1 else CompoundLocation(parts)
            q = {"gene": [name]}
            if kind == "tRNA":
                q["anticodon"] = [ac]
            r.features.append(SeqFeature(loc, type="gene", qualifiers={"gene": [name]}))
            r.features.append(SeqFeature(loc, type=kind, qualifiers=q))
    return r


def run_mitos2(fa: Path, workdir: Path, prefix: str, code: int, refdir: Optional[str] = None,
               bin_dir: Optional[str] = None) -> List[SeqRecord]:
    """MITOS2 runmitos (unchanged, refseq89m) on each sequence of one sample (in workdir/mitos2/<record>)."""
    b, d = mito_env.mitos_bin(bin_dir), mito_env.mitos_refdir(refdir)
    if not b:
        raise RuntimeError(f"MITOS2 (runmitos) not found: install it ({MITOS2_INSTALL}) or give --mitos_bin")
    if not d:
        raise RuntimeError("MITOS2 reference data refseq89m not found: run `mitom setup`, or give --mitos_refdir")
    seqs = read_fasta(fa)
    recs = []
    for k, (sid, seq, topology) in enumerate(seqs):
        rid = record_id(prefix, sid, k, len(seqs))
        out = workdir / "mitos2" / rid
        (out / "result").mkdir(parents=True, exist_ok=True)       # runmitos refuses a missing -o folder
        (out / "in.fa").write_text(f">{rid}\n{seq}\n")
        cmd = mitos2_command(out / "in.fa", out / "result", code, b, d, linear=topology == "linear")
        print("[engine] " + " ".join(cmd), file=sys.stderr, flush=True)
        rc = run_tool(cmd, cwd=out, env=engine_env(b, isolated=True), capture=False).returncode
        if rc or not (out / "result" / "result.bed").is_file():
            raise RuntimeError(f"MITOS2 failed for {rid} (exit {rc}, no result.bed)")
        recs.append(mitos2_record(out / "result", seq, rid, code, topology))
    return recs


def mitos2_record(result_dir: Path, seq: str, rid: str, code: Optional[int], topology: str) -> SeqRecord:
    """MITOS2 result.bed as the record mito_core.py writes: the core's names and qualifiers (CDS /codon_start
    /transl_table, tRNA /product and `anticodon: XXX` note), a gene feature before each; fragments of a gene joined
    and features across the origin written as joins (mito_rescue.mitos_to_record); OH / OL and names without a core
    name (introns) left out."""
    mitos = mitos_to_record(Path(result_dir), Seq(seq), rid)
    rec = SeqRecord(Seq(seq), id=rid, name=rid[:16], description="")
    rec.annotations.update(molecule_type="DNA", topology=topology, data_file_division="INV")
    rec.features.append(SeqFeature(FeatureLocation(0, len(seq), 1), type="source",
                                   qualifiers={"organism": [rid], "organelle": ["mitochondrion"],
                                               "mol_type": ["genomic DNA"]}))
    feats = []
    for f in mitos.features:
        name = f.qualifiers["gene"][0]
        if f.type == "tRNA":
            ac = (f.qualifiers.get("anticodon") or [""])[0].upper().replace("T", "U")
            q = {"gene": [name], "product": [f"tRNA-{mito_names.AA1_TO_3.get(name[3:4].upper(), 'Xaa')}"],
                 "note": ([f"anticodon: {ac}"] if ac else []) + ["MITOS2"]}
        elif f.type in ("CDS", "rRNA") and name in mito_names.CORE_NAMES:
            gene, product = mito_names.CORE_NAMES[name]
            q = {"gene": [gene], "product": [product]}
            if f.type == "CDS":
                q.update(codon_start=["1"], **({"transl_table": [str(code)]} if code else {}))
            q["note"] = ["MITOS2"]
        else:
            continue
        feats.append(SeqFeature(f.location, type=f.type, qualifiers=q))
    for f in sorted(feats, key=lambda x: int(x.location.start)):
        rec.features.append(SeqFeature(f.location, type="gene", qualifiers={"gene": list(f.qualifiers["gene"])}))
        rec.features.append(f)
    return rec


def mitos2_topology(result_dir: Path) -> Optional[str]:
    """circular / linear from MITOS2's result.gff (region line, Is_circular=True / False); None when not stated."""
    gff = Path(result_dir) / "result.gff"
    m = re.search(r"Is_circular=(True|False)", gff.read_text(errors="ignore")) if gff.is_file() else None
    return None if not m else "circular" if m.group(1) == "True" else "linear"


def mitos2_seq_ids(result_dir: Path) -> List[str]:
    """Sequence names in MITOS2's result.bed (first column)."""
    return sorted({row.seq for row in read_result_bed(result_dir)})


# ── existing results (mitom import) ───────────────────────────────────────────
def find_results(src: Path) -> Tuple[List[Path], List[Path]]:
    """(GenBank files, MITOS2 result folders) under src. A folder with MitoZ result folders (*.result/*.gbf) gives
    only those (MitoZ also leaves a draft .gbf in its tmp_* folder); MITOS2 result folders are the folders with a
    result.bed."""
    src = Path(src)
    if src.is_file():
        return ([src], []) if src.suffix.lower() in GB_SUFFIXES else ([], [src.parent] if src.name == "result.bed" else [])
    mitos = sorted({p.parent for p in src.rglob("result.bed")})
    mitoz = sorted(src.rglob("*.result/*.gbf"))
    if mitoz:
        return mitoz, mitos
    gbs = sorted(p for p in src.rglob("*") if p.is_file() and p.suffix.lower() in GB_SUFFIXES
                 and not any(m in p.parents for m in mitos))
    return gbs, mitos


def written_by_mitoz(rec: SeqRecord) -> bool:
    """MitoZ writes `linear` and the definition `topology=linear.` for every record, whatever the genome."""
    return rec.description.strip().rstrip(".") == "topology=linear"


def import_results(src: Path, fastas: Sequence[Path] = (), code: Optional[int] = None,
                   topology: str = "auto") -> Tuple[List[SeqRecord], List[str]]:
    """Annotated records from existing results under src and notes for the log. topology auto: MITOS2 from its
    result.gff, GenBank from the LOCUS line, MitoZ output circular (it cannot tell; give --topology linear for
    contigs). A MITOS2 result needs its genome sequence: --fasta, or a FASTA file in the result folder or the
    folder above it, with the sequence named as in result.bed."""
    gbs, mitos = find_results(src)
    notes: List[str] = []
    recs: List[SeqRecord] = []
    n_mitoz = 0
    for g in gbs:
        for r in SeqIO.parse(str(g), "genbank"):
            if topology != "auto":
                r.annotations["topology"] = topology
            elif written_by_mitoz(r):
                r.annotations["topology"] = "circular"
                n_mitoz += 1
            if written_by_mitoz(r):
                r.description = ""
            recs.append(r)
    if n_mitoz:
        notes.append(f"{n_mitoz} MitoZ record(s) taken as circular (MitoZ writes `linear` for all; "
                     "use --topology linear for contigs)")
    if mitos:
        if not code:
            raise ValueError("MITOS2 results: give --genetic_code (the code MITOS2 ran with)")
        seqs: Dict[str, str] = {}
        for fa in [*fastas, *(p for m in mitos for d in (m, m.parent) for p in sorted(d.glob("*"))
                               if p.is_file() and p.suffix.lower() in FA_SUFFIXES)]:
            for r in SeqIO.parse(str(fa), "fasta"):
                seqs.setdefault(r.id, str(r.seq).upper())
        for m in mitos:
            for sid in mitos2_seq_ids(m):
                if sid not in seqs:
                    raise ValueError(f"MITOS2 result {m}: genome sequence '{sid}' not found; give it with --fasta")
                topo = topology if topology != "auto" else (mitos2_topology(m) or "circular")
                recs.append(mitos2_record(m, seqs[sid], sid, code, topo))
        notes.append(f"{len(mitos)} MITOS2 result folder(s) converted (genetic code {code})")
    seen: Dict[str, int] = {}
    for r in recs:                                          # one record per name
        seen[r.id] = seen.get(r.id, 0) + 1
        if seen[r.id] > 1:
            r.id = r.name = f"{r.id}_{seen[r.id]}"
    return recs, notes


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="engine", required=True)
    for name in ("mitoz", "mitos2"):
        s = sub.add_parser(name)
        s.add_argument("-i", "--input", required=True, help="FASTA of one sample (several sequences = its contigs)")
        s.add_argument("-o", "--out", required=True, help="work folder; the annotation is <out>/<prefix>.gbf")
        s.add_argument("--prefix", required=True)
        s.add_argument("--genetic_code", type=int, required=True)
        s.add_argument("--threads", type=int, default=1)
        if name == "mitoz":
            s.add_argument("--clade", required=True)
            s.add_argument("--mitoz_bin", default=None, help="folder with the mitoz program")
            s.add_argument("--species_name", default=None)
        else:
            s.add_argument("--refdir", default=None, help="folder that holds refseq89m")
            s.add_argument("--mitos_bin", default=None, help="folder with runmitos")
    a = p.parse_args(argv)
    out = Path(a.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    fa = Path(a.input).resolve()
    try:
        if a.engine == "mitoz":
            recs = run_mitoz(fa, out, a.prefix, a.clade, a.genetic_code, a.threads, a.mitoz_bin, a.species_name)
        else:
            recs = run_mitos2(fa, out, a.prefix, a.genetic_code, a.refdir, a.mitos_bin)
    except RuntimeError as e:
        print(f"[engine] {e}", file=sys.stderr)
        return 1
    SeqIO.write(recs, str(out / f"{a.prefix}.gbf"), "genbank")
    n = {t: sum(f.type == t for r in recs for f in r.features) for t in ("CDS", "rRNA", "tRNA")}
    print(f"[engine] {a.engine}: {len(recs)} record(s), {n['CDS']} protein genes, {n['rRNA']} rRNAs, "
          f"{n['tRNA']} tRNAs -> {out / (a.prefix + '.gbf')}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
