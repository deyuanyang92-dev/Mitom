#!/usr/bin/env python3
"""mito_trna.py - multi-tool tRNA annotation of mitochondrial genomes, GeSeq style.

Several tRNA predictors are run on each record and ALL their predictions are written into one
GenBank file for review (e.g. in Geneious), as GeSeq does with its third-party tRNA tools
(Tillich et al. 2017, Nucleic Acids Res 45:W6, doi:10.1093/nar/gkx391). A second GenBank keeps one
tRNA per locus for submission. The number of tRNAs is never imposed: it varies widely in animal
mtDNA (e.g. 1-2 in most anthozoans, 23-24 with duplications in some annelids), so every
prediction is judged on its own evidence.

Sources (any subset, --tools):
  mitos2    MITOS2 runmitos, tRNA module only (metazoan reference data, e.g. refseq89m)
  arwen     ARWEN (metazoan mitochondrial tRNA; recommended by GeSeq for metazoan mtDNA)
  aragorn   ARAGORN (-mt)
  trnascan  tRNAscan-SE 2 (-O organellar mode)
  input     tRNA features already present in the input GenBank
  mitom     MitoM-ML (trna_ml.py): MITOS2 tRNA CMs + clade CMs on the record's NCBI lineage, a learned caller
            using ARWEN / ARAGORN as extra evidence, learned ends. Unseen-family benchmark (1,069 RefSeq invertebrate
            mitogenomes, 22,099 tRNAs): recall 98.5 % at 0.17 extra tRNAs per genome vs 94.9 % / 0.19 for the
            MiTFi CMs alone; 0 calls on 150 shuffled genomes (projects/invert_trna/REPORT.md). Lineage: GenBank
            taxonomy of the record, else --lineage (without one only the pan-metazoan CMs are used)

Support is counted in independent method families (MITOS2; ARWEN+ARAGORN as one, same author and
algorithm; tRNAscan-SE; input annotation), not in tools.

Locus groups: predictions on the same strand overlapping >= 50 % of the shorter one, with the
same amino acid AND anticodon (so trnL1/trnL2 and trnS1/trnS2 stay separate). Overlapping
predictions with different anticodons form separate groups flagged as a conflict.

Selection --select mitom (default when MitoM-ML is installed: models/trna_ml, cmsearch, cmalign, refseq89m): the
selected tRNAs are the MitoM-ML calls (its boundaries and anticodons); the other sources stay in the review GenBank.
Selection --select mitos2_or_support (default otherwise): a group is selected when MITOS2 found it, or when
>= --min_support (2) independent method families agree. Overlapping groups with different anticodons
(conflict): the one with stronger evidence (MITOS2, then more families) is kept, ties drop both.
--select mitos2: the input tRNAs are replaced by the MITOS2 tRNAs (tools mitos2,input;
the input tRNAs stay visible in the review GenBank). Boundaries otherwise come from the first source
in --priority (MITOS2 first: most exact boundaries vs NCBI).

MITOS2 runs in fast mode by default (scripts/runmitos_fast.py): runmitos re-searches every tRNA model it did not
find with a sensitive cmsearch (--max, E 0.1, single CPU): 640 s vs 8 s on one anthozoan (Metridium), plus spurious
tRNAs; fast mode skips that pass, --mitos_mode sensitive runs it in parallel (fast mode missed 28/28 of the weak tRNAs
found only by that pass, 2026-10-20 dev benchmark), --mitos_mode full is the unmodified runmitos. `mitom run` does
not run this script: its tRNAs come from the annotation core (mito_core.py, MitoM-ML with MiTFi).

The annotation core's tRNA step is here too (core_trnas, called by mito_core.py for each record): MitoM-ML calls
kept when sure enough or confirmed by MiTFi, MiTFi's ends for vertebrates, MiTFi's calls of the families MitoM-ML does
not call; MiTFi alone (cmsearch with the MiTFi models of refseq89m, a sensitive pass outside the protein gene / rRNA
interiors, mitfi_import.jar) when MitoM-ML is not installed.

Outputs per input file: review/<name>.gbf (all predictions; /note = source, score, anticodon,
group, support), selected/<name>.gbf (selected tRNAs only; input for mito_curate.py -i out/selected),
and trna_report.tsv (one row per group).

Usage:
  python3 mito_trna.py -i gb_or_fasta -o out_dir --genetic_code 5 \
      --mitos_refdir DIR --mitos_ref refseq89m --mitos_bin ENV/bin --bin ENV/bin
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Dict, List, Optional, Tuple

from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mito_env  # noqa: E402  (finds external programs: config, PATH, conda environments)
from mito_names import AA3, trna_symbol as trna_name  # noqa: E402  (trnL1/L2, trnS1/S2 by anticodon)
from mito_names import AA1_TO_3 as ONE2THREE  # noqa: E402
from mito_common import GB_SUFFIXES, old_spelling, run_tool, write_tsv  # noqa: E402
from mito_common import cli_guard, count_seq_records, InputError, prepare_input, require_records  # noqa: E402
from mito_common import CoreFeature, positions, run_or_raise  # noqa: E402  (the annotation core's tRNA step)

GB_EXT = set(GB_SUFFIXES)
TIMEOUT = 3600


@dataclass
class Pred:
    source: str
    start: int            # 1-based inclusive
    end: int
    strand: int
    aa: str               # one-letter
    anticodon: str        # 5'->3', lowercase dna (e.g. 'tag'); '' if unknown
    score: str = ""
    group: int = -1
    extra: Dict[str, str] = field(default_factory=dict)

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def run(cmd: List[str], cwd: Optional[Path] = None, env: Optional[Dict[str, str]] = None) -> str:
    try:
        return run_tool(cmd, timeout=TIMEOUT, cwd=cwd, env=env).stdout
    except (subprocess.TimeoutExpired, OSError) as exc:
        print(f"[mito_trna] {cmd[0]} failed: {exc}", file=sys.stderr)
        return ""


# ── predictors ───────────────────────────────────────────────────────────────
LASLETT = re.compile(r"m?tRNA-([A-Z][a-z]{2}|\?\?\?)\S*\s+(c?)\[(\d+),(\d+)\]\s+(\S+)?\s*\(([acgtnACGTN]{3})\)")


def parse_laslett(text: str, source: str) -> List[Pred]:
    """ARWEN / ARAGORN text output ('mtRNA-Leu  c[100,165]  33  (tag)'). A tRNA across the origin of a
    circular sequence is reported as [start,end] with start > end; it is kept with end = end + L
    (extended coordinates, written as a join in GenBank)."""
    out = []
    for m in LASLETT.finditer(text):
        aa = AA3.get(m.group(1), "")
        if aa:
            s, e = int(m.group(3)), int(m.group(4))
            out.append(Pred(source, s, e, -1 if m.group(2) else 1, aa, m.group(6).lower(), m.group(5) or "",
                            extra={"wrap": "1"} if s > e else {}))
    return out


def unwrap(preds: List[Pred], L: int) -> List[Pred]:
    """Origin-spanning predictions (start > end) -> extended coordinates end + L."""
    for p in preds:
        if p.start > p.end:
            p.end += L
    return preds


def arwen(fa: Path, code: int, b: Dict[str, str], tmp: Path, a) -> List[Pred]:
    return parse_laslett(run([b["arwen"], f"-gc{code}", "-mtx", "-w", str(fa)]), "ARWEN")


def aragorn(fa: Path, code: int, b: Dict[str, str], tmp: Path, a) -> List[Pred]:
    return parse_laslett(run([b["aragorn"], "-mt", f"-gc{code}", "-w", str(fa)]), "ARAGORN")


def trnascan(fa: Path, code: int, b: Dict[str, str], tmp: Path, a) -> List[Pred]:
    res = tmp / "trnascan.out"
    res.unlink(missing_ok=True)
    run([b["tRNAscan-SE"], "-O", "-q", "--brief", "-o", str(res), str(fa)])
    out = []
    for line in res.read_text().splitlines() if res.exists() else []:
        c = [x.strip() for x in line.split("\t")]
        if len(c) < 9 or c[4][:3] not in AA3:
            continue
        s, e = int(c[2]), int(c[3])
        out.append(Pred("tRNAscan-SE", min(s, e), max(s, e), 1 if s < e else -1, AA3[c[4][:3]],
                        c[5].lower().replace("u", "t"), c[8]))
    return out


# trnL1(tag), trnS2(tga), extra copies of one model as trnM_0(cat) / trnM_1(cat), and trnQ(---) when the
# anticodon could not be read (kept, anticodon '')
MITOS_NAME = re.compile(r"^trn([A-Z])\d?(?:_\d+)?\(([acgtnu]{3}|---)\)", re.I)


def mitos2(fa: Path, code: int, b: Dict[str, str], tmp: Path, a) -> List[Pred]:
    """MITOS2, tRNA prediction only (proteins, rRNA, introns and OL skipped for speed)."""
    if a.mitos_mode not in ("fast", "sensitive", "full"):          # any other value used to act as 'sensitive'
        raise ValueError(f"mitos_mode must be fast, sensitive or full, not {a.mitos_mode!r}")
    import mito_engines
    od = tmp / "mitos"
    shutil.rmtree(od, ignore_errors=True)
    od.mkdir()
    env_bin = str(Path(b["runmitos"]).parent)
    run(mito_engines.mitos2_command(fa, od, code, env_bin, a.mitos_refdir, ref=a.mitos_ref, mode=a.mitos_mode,
                                    skip=("prot", "rrna", "intron", "oril"), noplots=False),
        env=mito_engines.engine_env(env_bin, isolated=True))
    bed = od / "result.bed"
    if not bed.exists():                              # failed run, not "no tRNA found"
        raise RuntimeError("MITOS2 produced no result.bed")
    return parse_mitos_bed(bed.read_text())


def parse_mitos_bed(text: str) -> List[Pred]:
    """tRNA lines of MITOS2 result.bed (mito_engines.parse_result_bed; 0-based start)."""
    import mito_engines
    out = []
    for row in mito_engines.parse_result_bed(text):
        m = MITOS_NAME.match(row.name)
        if m:
            out.append(Pred("MITOS2", row.start + 1, row.end, row.strand, m.group(1).upper(),
                            m.group(2).lower().replace("u", "t").strip("-"), row.score, extra={"mitos_name": row.name}))
    return out


TOOLS = {"mitos2": mitos2, "arwen": arwen, "aragorn": aragorn, "trnascan": trnascan}


def mitom_ready(b: Dict[str, str], refdir: Optional[str]) -> bool:
    import importlib.util
    import trna_ml
    return bool(refdir and b.get("cmsearch") and b.get("cmalign") and (trna_ml.MODELS / "caller_train.tsv.gz").exists()
                and ((trna_ml.MODELS / "clades.tar.gz").exists() or (trna_ml.MODELS / "clades").is_dir())
                and importlib.util.find_spec("sklearn") is not None)     # else the first call fails, not MiTFi taken


def mitom(seq: str, lineage: str, code: int, b: Dict[str, str], a, preds: List[Pred], tmp: Path) -> List[Pred]:
    """MitoM-ML calls; ARWEN / ARAGORN predictions of this record are used as evidence (run here when not
    among --tools)."""
    import trna_ml
    ext = [p for p in preds if p.source in ("ARWEN", "ARAGORN")]
    if not ext:
        fa = tmp / "g.fa"
        for t in ("arwen", "aragorn"):
            if t in b:
                ext += unwrap(TOOLS[t](fa, code, b, tmp, a), len(seq))
    else:
        ext = [Pred(p.source, p.start, p.end + (len(seq) if p.start > p.end else 0), p.strand, p.aa, p.anticodon)
               for p in ext]
    out = []
    for c in trna_ml.predict(seq, lineage, code, Path(a.mitos_refdir) / a.mitos_ref / "ncRNA", b["cmsearch"],
                             b["cmalign"], ext, cpu=a.threads):
        out.append(Pred("MitoM-ML", c["start"], c["end"], c["strand"], c["aa"], c["anticodon"], f"{c['p']:.3f}",
                        extra={"family": c["family"], "bound_p": str(c["bound_p"])}))
    return out


def from_input(rec) -> List[Pred]:
    """tRNA features of the input annotation (amino acid from /product, anticodon from /note or name). One of three
    readers of a tRNA's amino acid / anticodon, kept apart on purpose (see mito_names.trna_name; mito_curate N1)."""
    out = []
    for f in rec.features:
        if f.type != "tRNA":
            continue
        txt = " ".join(f.qualifiers.get("product", []) + f.qualifiers.get("gene", []) + f.qualifiers.get("note", []))
        m = re.search(r"tRNA-([A-Z][a-z]{2})", txt)
        aa = AA3.get(m.group(1), "") if m else ""
        if not aa:
            g = re.match(r"trn([A-Z])", f.qualifiers.get("gene", [""])[0])
            aa = g.group(1) if g else ""
        ac = re.search(r"\(([acgtuACGTU]{3})\)", txt)
        if aa:
            out.append(Pred("input", int(f.location.start) + 1, int(f.location.end), f.location.strand or 1, aa,
                            ac.group(1).lower().replace("u", "t") if ac else ""))
    return out


# ── grouping and selection ──────────────────────────────────────────────────
CIRC_L = 0                                            # length of the current (circular) record


def overlap(p: "Pred", q: "Pred") -> int:
    """Overlap on a circular sequence: origin-spanning predictions use extended coordinates (end > L),
    so compare also with one of them shifted by -L."""
    best = min(p.end, q.end) - max(p.start, q.start) + 1
    if CIRC_L:
        for dp, dq in ((-CIRC_L, 0), (0, -CIRC_L)):
            best = max(best, min(p.end + dp, q.end + dq) - max(p.start + dp, q.start + dq) + 1)
    return best


def group(preds: List[Pred]) -> List[List[Pred]]:
    groups: List[List[Pred]] = []
    for p in sorted(preds, key=lambda x: (x.start, x.end)):
        for g in groups:
            q = g[0]
            ov = overlap(p, q)
            same = p.aa == q.aa and (not p.anticodon or not q.anticodon or p.anticodon == q.anticodon)
            if p.strand == q.strand and same and ov >= 0.5 * min(p.length, q.length):
                g.append(p)
                break
        else:
            groups.append([p])
    for i, g in enumerate(groups):
        for p in g:
            p.group = i + 1
    return groups


def conflicts(groups: List[List[Pred]]) -> Dict[int, List[int]]:
    """Groups overlapping another group on the same strand with a different aa/anticodon."""
    out: Dict[int, List[int]] = {}
    for i, g in enumerate(groups):
        for j, h in enumerate(groups):
            if i == j or g[0].strand != h[0].strand:
                continue
            if overlap(g[0], h[0]) >= 0.5 * min(g[0].length, h[0].length):
                out.setdefault(i + 1, []).append(j + 1)
    return out


# Independent method families: ARWEN and ARAGORN share author and algorithm (Laslett & Canback),
# so their agreement is one line of evidence, not two (Metridium: both report trnF/trnS2 where NCBI
# and the covariance-model methods have no tRNA).
FAMILY = {"MitoM-ML": "MitoM-ML (MITOS2 + clade covariance models, learned caller)",
          "MITOS2": "MITOS2 (covariance models, metazoan mt tRNAs)", "ARWEN": "ARWEN/ARAGORN (heuristic)",
          "ARAGORN": "ARWEN/ARAGORN (heuristic)", "tRNAscan-SE": "tRNAscan-SE (covariance models)",
          "input": "input annotation"}


def families(g: List[Pred]) -> List[str]:
    return sorted({FAMILY.get(p.source, p.source) for p in g})


def strength(g: List[Pred], select: str) -> tuple:
    """Evidence rank used to settle conflicts: MITOS2 found it (mitos2_or_support mode), then families."""
    return (select == "mitos2_or_support" and any(p.source == "MITOS2" for p in g), len(families(g)))


def anticodon_of(g: List[Pred]) -> str:
    c = Counter(p.anticodon for p in g if p.anticodon)
    return c.most_common(1)[0][0] if c else ""


def choose(g: List[Pred], priority: List[str]) -> Pred:
    order = {s: i for i, s in enumerate(priority)}
    return min(g, key=lambda p: order.get(p.source, 99))


def feature(p: Pred, ac: str, note: str, L: int) -> SeqFeature:
    name = trna_name(p.aa, ac)
    loc: Any
    if p.end > L:                                     # across the origin of a circular sequence
        parts = [FeatureLocation(p.start - 1, L, strand=p.strand), FeatureLocation(0, p.end - L, strand=p.strand)]
        loc = CompoundLocation(parts if p.strand == 1 else parts[::-1])
    else:
        loc = FeatureLocation(p.start - 1, p.end, strand=p.strand)
    return SeqFeature(loc, type="tRNA",
                      qualifiers={"gene": [name + (f"({ac.replace('t', 'u')})" if ac else "")],   # the core's style trnN(guu)
                                  "product": [f"tRNA-{ONE2THREE[p.aa]}"], "note": [note]})


# ── the annotation core's tRNA step (mito_core.py): MitoM-ML checked by MiTFi; MiTFi alone without MitoM-ML
CM_EVALUE = 0.01           # cmsearch fast pass and MiTFi threshold of that pass (MITOS2 runmitos default --ncev)
NO_SENSITIVE_TRNA = {"cnidaria"}   # phyla whose missing tRNA models are not searched again (sensitive pass):
                           # cnidarian mitogenomes encode 1-6 tRNAs, the ~20-model --max pass found junk (18 records,
                           # leave-out: 13 of 17 extra tRNAs, 2 real ones; 344 -> 50 s median; a false trnE inside nad5
                           # was an NCBI error in the stress test)
CM_EVALUE_SENSITIVE = 0.1  # cmsearch pass for models without a hit (MITOS2 sensitive pass)
CM_SENSITIVE = ["--nohmm", "--F6", "3e-4", "--notrunc"]   # ... run without the HMM filters, as MITOS2's --max, but
                           # with the CM's own CYK filter (P <= 3e-4), outside the protein gene / rRNA interiors
                           # (gene_free): MITOS2's --max on the whole record took minutes (NC_014885.1, 17 models, 2
                           # threads: 297 s, now 18 s; 542 benchmark records: 11.9 -> 1.0 h CPU; final tRNAs the same
                           # except 14 weak MiTFi calls inside genes / rRNAs, none in NCBI, and 2 other wrong ones). HMM
                           # filters (--mid) lost degenerate tRNAs: recall 95.6 -> 94.5 % on the 187 unseen-family
                           # records that have this pass
TRNA_SENSITIVE_FLANK = 100  # nt at each end of a protein gene / rRNA still searched (a tRNA plus 30 nt)
TRNA_NAMES = ["trn" + f for f in "A C D E F G H I K L1 L2 M N P Q R S1 S2 T V W Y".split()]
TRNA_RESCUE_OVL = 10       # a MiTFi tRNA of a family MitoM-ML did not call is added when it overlaps no call by more
                           # than this many nt (rescue_trnas)
TRNA_SURE_P = 0.9          # MitoM-ML calls below this p need MiTFi's call of the same tRNA, and a weaker second call
                           # of a family needs at least this p (trna_second_opinion; tuning halves: 37 % of the calls
                           # below 0.9 without MiTFi were wrong, 2 % above, 30 % of the weaker second copies below)
COPY_EVALUE = 1e-3         # MiTFi copies (second hits of a model) only this good (MITOS2 runmitos.py)
MITFI_TRNA_OVL = 50        # tRNAs may overlap this much in MiTFi: runmitos --maxtrnaovl default (its function default 10 is
                           # not what the command runs); with 10, NC_060484.1 lost trnL2 overlapping trnP by 12 nt (MITOS2
                           # keeps both); MiTFi -overlap 50 -method l reproduces MITOS2's own tRNAout.nc line for line


def mitfi_files(mitos_bin: Optional[str]) -> Tuple[Path, Path]:
    """mitfi_import.jar and mitfiGeneticCodes of the mitos package (searched next to the MITOS2 python)."""
    bins = [Path(mitos_bin)] if mitos_bin else []
    mb = mito_env.mitos_bin(None)
    if mb:
        bins.append(Path(mb))
    for b in bins:
        for p in sorted((b.parent / "lib").glob("python3*/site-packages/mitos/mitfi/parser/mitfi_import.jar")):
            return p, p.parent / "mitfiGeneticCodes"
    for p in Path(__file__).resolve().parent.glob("models/mitfi/mitfi_import.jar"):
        return p, p.parent / "mitfiGeneticCodes"
    raise SystemExit("[core] MiTFi (mitfi_import.jar of the mitos package) not found")


def cmsearch(cm: str, model: Path, genome_fa: Path, out: Path, sensitive: bool, log, rrna: bool = False,
             z: Optional[float] = None) -> Path:
    """cmsearch as MITOS2 runs it (local mode, --cpu 1 --notextw --smxsize 80000); fast pass --noF4b -E 0.01,
    sensitive pass -E 0.1 with CM_SENSITIVE (rRNA: --nohmm). Local mode finds the truncated tRNAs of nematodes. z =
    search space (Mb) for the E-values: a circular record is searched with ORIGIN_EXT nt appended, MITOS2 searches it at
    its length."""
    cmd = [cm, "--cpu", "1", "--notextw", "--smxsize", "80000", "--tblout", str(out.with_suffix(".tbl"))]
    if z:
        cmd += ["-Z", f"{z:.6f}"]
    if rrna:                                   # same rrnL hit, 536 -> 463 MB, 3.8 -> 2.7 s (alignment matrix cap)
        cmd += ["--mxsize", "128"]
    if sensitive:
        cmd += (["--nohmm"] if rrna else CM_SENSITIVE) + ["-E", str(CM_EVALUE_SENSITIVE)]
    else:
        cmd += ["--noF4b", "-E", str(CM_EVALUE)]
    p = run_or_raise(cmd + [str(model), str(genome_fa)])
    out.write_text(p.stdout)
    return out


def gene_free(genes: List[CoreFeature], n: int, length: int, flank: int = TRNA_SENSITIVE_FLANK) -> List[Tuple[int, int]]:
    """Stretches (0-based, half-open) of the searched sequence of this length outside the interiors of the protein genes
    / rRNAs (a gene minus `flank` nt at each end: the core's gene ends may cover a tRNA, also where two genes abut);
    n = record length (a circular record carries its first nt again after n)."""
    cov = bytearray(length)
    for f in genes:
        for k in (-n, 0, n):
            a, b = max(0, f.start + k + flank), min(length, f.end + k - flank)
            if a < b:
                cov[a:b] = b"\x01" * (b - a)
    return [(m.start(), m.end()) for m in re.finditer(b"\x00+", bytes(cov))]


CM_ROW = re.compile(r"^(\s+\(\d+\)\s+[!?](?:\s+\S+){3}\s+)s(\d+)(\s+)(\d+)(\s+)(\d+)(\s.*)$")        # hit table
CM_HIT = re.compile(r"^(\s+\(\d+\)\s+[!?](?:\s+\S+){3}\s+cm(?:\s+\S+){3}\s+)(\d+)(\s+)(\d+)(\s.*)$")  # alignment header
CM_ALN = re.compile(r"^(\s+)s(\d+)(\s+)(\d+)(\s.*\s)(\d+)(\s*)$")                                     # alignment row


def shift_cmout(text: str, name: str) -> str:
    """cmsearch output for stretches named s<offset> (gene_free) rewritten as a search of the whole sequence `name`
    (what MiTFi reads): target name and coordinates of each hit in the hit table, in its alignment header and in its
    alignment row."""
    out, off = [], 0
    for line in text.splitlines():
        row, hit, aln = CM_ROW.match(line), CM_HIT.match(line), CM_ALN.match(line)
        if line.startswith(">> s"):
            off = int(line.split()[1][1:])
            line = f">> {name}" + line[len(f">> s{off}"):]
        elif row:
            o = int(row.group(2))
            line = (f"{row.group(1)}{name}{row.group(3)}{int(row.group(4)) + o}{row.group(5)}{int(row.group(6)) + o}"
                    f"{row.group(7)}")
        elif hit:
            line = f"{hit.group(1)}{int(hit.group(2)) + off}{hit.group(3)}{int(hit.group(4)) + off}{hit.group(5)}"
        elif aln and int(aln.group(2)) == off:
            line = (f"{aln.group(1)}{name}{aln.group(3)}{int(aln.group(4)) + off}{aln.group(5)}{int(aln.group(6)) + off}"
                    f"{aln.group(7)}")
        out.append(line)
    return "\n".join(out) + "\n"


def annotate_trnas(genome: str, genome_fa: Path, refdir: Path, code: int, threads: int, tools: dict, tmp: Path, log,
                   n: Optional[int] = None, phylum: Optional[str] = None,
                   genes: Optional[List[CoreFeature]] = None) -> List[CoreFeature]:
    """MITOS2's tRNA step: cmsearch per MiTFi model (fast pass, then a sensitive pass for the tRNA names without a hit,
    CM_SENSITIVE outside the interiors of the protein genes / rRNAs `genes`), MiTFi for selection, structure and
    anticodon. genome / genome_fa = the searched sequence, n = record length (circular records carry ORIGIN_EXT appended
    nt; positions are mapped back modulo n)."""
    n = n or len(genome)
    models = sorted((refdir / "ncRNA").glob("trn*.cm"))
    seqname = genome_fa.name

    def search(model: Path, sensitive: bool, fa: Path = genome_fa) -> Path:
        # E-values for the record's length on both strands, as MITOS2 (which searches rotations of the record, not the
        # record + ORIGIN_EXT): the 3,000 appended nt made every E ~1.2x larger on a 15-kb record, so hits MITOS2 keeps
        # at E 0.08-0.1 were lost (NC_060482.1 / NC_060708.1 acanthocephalan trnD, E 0.096 -> 0.12)
        return cmsearch(tools["cmsearch"], model, fa, tmp / f"{seqname}_{model.stem}.cmout", sensitive, log,
                        z=2 * n / 1e6)

    jar, codes = tools["mitfi"]

    def mitfi(outs: List[Path], evalue: float) -> List[CoreFeature]:
        cmd = [tools["java"], "-jar", str(jar), "-code", str(code), "-codefile", str(codes), "-structure", "-genome",
               str(genome_fa), "-glength", str(n), "-evalue", str(evalue), "-method", "l", "-overlap",
               str(MITFI_TRNA_OVL)]
        for o in outs:
            cmd += ["-in", str(o)]
        return parse_mitfi(run_or_raise(cmd).stdout, log, n)

    # MITOS2's procedure (runmitos.py): fast pass judged at E <= 0.001 - a model whose only hits are weaker counts
    # as empty and is searched again (MITOS2: --max on the whole record); the final MiTFi run accepts E <= 0.1
    with ThreadPoolExecutor(max(1, threads)) as ex:
        outs = list(ex.map(lambda m: search(m, False), models))
    feats = mitfi(outs, CM_EVALUE)
    accepted = {f.name for f in feats}          # MITOS2 re-runs the models whose NAME is missing (MiTFi names by anticodon:
    missing = [m for m in models if m.stem not in accepted]   # an L1 gene is often found by trnL2.cm)
    if (phylum or "").lower() in NO_SENSITIVE_TRNA:
        log.info(f"[core] tRNA sensitive pass skipped for {phylum}")
        missing = []
    free = gene_free(genes or [], n, len(genome))
    if missing and free:
        with open(genome_fa) as fh:
            name = fh.readline()[1:].split()[0]
        free_fa = tmp / f"{seqname}.free.fa"
        free_fa.write_text("".join(f">s{s}\n{genome[s:e]}\n" for s, e in free))
        log.info(f"[core] tRNA sensitive pass for {len(missing)} model(s) on {sum(e - s for s, e in free)} of "
                 f"{len(genome)} nt (outside protein gene / rRNA interiors): " + " ".join(m.stem for m in missing))

        def sensitive(model: Path) -> None:
            out = search(model, True, free_fa)
            out.write_text(shift_cmout(out.read_text(), name))
        with ThreadPoolExecutor(max(1, threads)) as ex:
            list(ex.map(sensitive, missing))
        feats = mitfi(outs, CM_EVALUE_SENSITIVE)
    log.info(f"[core] tRNAs: {len(feats)}")
    return feats


def annotate_trnas_ml(genome: str, n: int, code: int, threads: int, tools: dict, lineage: str, circular: bool,
                      log) -> List[CoreFeature]:
    """MitoM-ML tRNAs of one record (genome = the record without the appended ORIGIN_EXT nt)."""
    import trna_ml
    ext: list = []
    with tempfile.TemporaryDirectory(prefix="mito_core_trna.") as td:
        fa = Path(td) / "g.fa"
        fa.write_text(f">g\n{genome}\n")
        for t in ("arwen", "aragorn"):
            if tools.get(t):
                ext += unwrap(TOOLS[t](fa, code, {t: tools[t]}, Path(td), None), n)
    calls = trna_ml.predict(genome, lineage, code, Path(tools["refdir"]) / "ncRNA", tools["cmsearch"], tools["cmalign"], ext,
                            cpu=max(1, threads))
    feats = []
    for c in calls:
        if c["end"] > n and not circular:
            continue
        feats.append(CoreFeature("tRNA", mitfi_trna_name(c["family"], c["anticodon"], ""), c["start"] - 1, c["end"], c["strand"],
                             f"MitoM-ML p {c['p']:.3f} ends p {c['bound_p']:.3f}", anticodon=c["anticodon"],
                             score=c["p"], model="MitoM-ML"))
    log.info(f"[core] tRNAs (MitoM-ML, lineage {lineage or 'none'}, {'ARWEN/ARAGORN evidence' if ext else 'no ARWEN/ARAGORN'}):"
             f" {len(feats)}")
    return feats


def rescue_trnas(ml: List[CoreFeature], mitfi: List[CoreFeature], n: int) -> List[CoreFeature]:
    """MitoM-ML calls plus MiTFi's calls of the families MitoM-ML did not call (lowest E first, one per family), each
    when it overlaps no call by more than TRNA_RESCUE_OVL nt (either strand, across the origin). The learned caller
    rejects weak hits in genomes with few strong tRNAs, where MiTFi still finds some (unseen-family benchmark,
    benchmark/trna_ml/scripts/rescue_sim.py: F1 98.41 -> 98.66 % on the tuning half, 98.77 -> 98.83 % on the other
    half; Acanthocephala recall 17 -> 43 %; no cnidarian record changed)."""
    def evalue(f: CoreFeature) -> float:
        m = re.search(r" E ([0-9.eE+-]+)", f.note)
        return float(m.group(1)) if m else 1.0

    def overlap(a: CoreFeature, b: CoreFeature) -> int:
        return max(min(a.end, b.end + k) - max(a.start, b.start + k) for k in (-n, 0, n))
    out, have = list(ml), {f.name for f in ml}
    for f in sorted(mitfi, key=evalue):
        if f.name not in have and all(overlap(f, g) <= TRNA_RESCUE_OVL for g in out):
            out.append(f)
            have.add(f.name)
    return out


def same_trna(a: CoreFeature, b: CoreFeature, n: int) -> bool:
    """Same tRNA name and strand, overlapping >= 50 % of the shorter (across the origin too)."""
    o = max(min(a.end, b.end + k) - max(a.start, b.start + k) for k in (-n, 0, n))
    return a.name == b.name and a.strand == b.strand and o >= 0.5 * min(a.end - a.start, b.end - b.start)


def trna_second_opinion(ml: List[CoreFeature], mitfi: Optional[List[CoreFeature]], n: int) -> List[CoreFeature]:
    """MitoM-ML calls kept when sure enough: a weaker second call of a family needs p >= TRNA_SURE_P, and any call with
    p < TRNA_SURE_P needs MiTFi's call of the same tRNA (not applied without MiTFi). Chosen on the tuning halves of
    the unseen-family set and the GenBank-only annelids (precision >= MITOS2 + 0.1 points, most recall), then
    checked on the held-out halves (benchmark/trna_ml/README.md)."""
    best: Dict[str, float] = {}
    for f in ml:
        best[f.name] = max(best.get(f.name, 0.0), f.score)
    out = [f for f in ml if f.score >= best[f.name] or f.score >= TRNA_SURE_P]
    if mitfi is not None:
        out = [f for f in out if f.score >= TRNA_SURE_P or any(same_trna(f, g, n) for g in mitfi)]
    return out


def mitfi_ends(ml: List[CoreFeature], mitfi: List[CoreFeature], n: int) -> List[CoreFeature]:
    """MitoM-ML calls with MiTFi's ends wherever MiTFi calls the same tRNA (name and strand, >= 50 % of the shorter).
    For vertebrates: the learned ends never saw one and are 1 nt off NCBI more often than MiTFi's (434 large-set
    chordates: both ends exact 90.2 -> 91.4 %, MITOS2 91.1 %; 188 held-out vertebrates 91.3 -> 92.7 %)."""
    out = []
    for f in ml:
        pf = positions(f, n)                     # the MiTFi call that overlaps most (not the first in start order)
        m = max((g for g in mitfi if same_trna(f, g, n)), key=lambda g: (len(pf & positions(g, n)), g.score),
                default=None)
        out.append(f if m is None else CoreFeature(f.kind, f.name, m.start, m.end, f.strand, f.note + ", ends by MiTFi",
                                               anticodon=f.anticodon, score=f.score, model=f.model))
    return out


def vertebrate(lineage: str) -> bool:
    """Lineage of a vertebrate; a bare 'Metazoa;Chordata' (FASTA input, --phylum only) counts as one."""
    tokens = {t.strip() for t in lineage.split(";")}
    return bool(tokens & {"Vertebrata", "Craniata"}) or lineage.strip() == "Metazoa;Chordata"


def mitfi_trna_name(aa: str, anticodon: str, model: str) -> str:
    """MITOS2 naming: the amino acid MiTFi derives from the anticodon; the two Leu / Ser genes by codon family
    (L1 = CUN: anticodon UAG; L2 = UUR: UAA; S1 = AGN: GCU / UCU; S2 = UCN: UGA), else the model's number."""
    aa = (aa or "").strip().upper()
    if re.fullmatch(r"[LS][12]", aa):                  # MiTFi names the two Leu / Ser genes itself (codon family)
        return "trn" + aa
    aa = aa[:1]
    ac = anticodon.replace("u", "t")
    if aa == "L":
        return "trnL1" if ac == "tag" else "trnL2" if ac == "taa" else ("trnL" + (re.search(r"trnL([12])", model) or [None, ""])[1])
    if aa == "S":
        return "trnS1" if ac in ("gct", "tct") else "trnS2" if ac == "tga" else ("trnS" + (re.search(r"trnS([12])", model) or [None, ""])[1])
    if aa and aa in ONE2THREE:
        return "trn" + aa
    m = re.search(r"(trn[A-Z](?:[12])?)", model)
    return m.group(1) if m else ""


def parse_mitfi(text: str, log, n: Optional[int] = None) -> List[CoreFeature]:
    """MiTFi output (what mitos.mitfi.parse reads): three lines per feature, the first
    '>acc|code|start|stop|strand|qstart|qstop|bits|evalue|pvalue|gc|acpos|anticodon|aa|model|tophit|copy|g/l',
    then the sequence and the structure."""
    feats: List[CoreFeature] = []
    for line in text.splitlines():
        line = line.strip()
        if not line.startswith((">", "|")):
            continue
        c = line[1:].split("|")
        try:
            start, stop, strand = int(c[2]) - 1, int(c[3]), int(c[4])
            bits, evalue = float(c[7]), float(c[8])
            acpos = int(c[11]) - 1
            anticodon = c[12].strip().lower().replace("u", "t") if acpos >= 0 else ""
            aa, model = c[13], c[14]
            m_aa = re.match(r"trn(\w\d*)", aa)       # no anticodon found: MiTFi writes the model name here ("trnL1";
            aa = m_aa.group(1) if m_aa else aa        # MITOS2 mitfi parser: "workaround"); [:1] made it Thr (NC_058285.1)
            tophit, copy = int(c[15]), int(c[16])
        except (IndexError, ValueError) as exc:
            log.warning(f"[core] MiTFi line not understood ({exc}): {line[:80]}")
            continue
        if n and start >= n:                     # found again in the appended first ORIGIN_EXT nt of a circular record
            start, stop = start - n, stop - n
        if not (tophit or (copy and evalue <= COPY_EVALUE)):  # MITOS2: a copy (second hit of a model) only when E <= 0.001
            continue
        name = mitfi_trna_name(aa, anticodon, model)
        if not name or start >= stop:
            continue
        mm = re.search(r"(trn[A-Z](?:[12])?)\.cm", model)
        feats.append(CoreFeature("tRNA", name, start, stop, strand, f"MiTFi {model} bits {bits:.1f} E {evalue:.1e}"
                             + (" copy" if copy and not tophit else ""), anticodon=anticodon, acpos=acpos, score=bits,
                             model=mm.group(1) if mm else model))
    # tRNAs overlapping by more than MiTFi allows (the same tRNA found at i and at i + n in a circular record: complete
    # overlap): the better score is kept
    feats.sort(key=lambda f: -f.score)
    kept: List[CoreFeature] = []
    L = n or max((f.end for f in feats), default=1)
    for f in feats:
        pf = positions(f, L)
        if all(len(pf & positions(k, L)) <= MITFI_TRNA_OVL
               # a piece of a better hit of the same tRNA (its part past the origin of a circular record, <= 50 nt,
               # found again at the sequence start: 2026-10-31 review, trnE across the origin cut to 32 nt)
               and not (k.name == f.name and k.strand == f.strand and len(pf & positions(k, L)) >= 0.9 * len(pf))
               for k in kept):
            kept.append(f)
    return sorted(kept, key=lambda f: f.start)


def core_trnas(genome: str, searched: str, genome_fa: Path, refdir: Path, code: int, threads: int, tools: dict,
               tmp: Path, log, n: int, phylum: Optional[str], lineage: str, circular: bool, use_ml: bool,
               genes: List[CoreFeature]) -> List[CoreFeature]:
    """The tRNAs of one record for the annotation core (mito_core.py): MitoM-ML calls (annotate_trnas_ml) kept when
    sure enough or confirmed by MiTFi (trna_second_opinion), with MiTFi's ends for vertebrates (mitfi_ends) and MiTFi's
    calls of the families MitoM-ML did not call (rescue_trnas); without MitoM-ML (use_ml False) MiTFi's tRNAs
    (annotate_trnas). genome = the record, searched / genome_fa = the record as searched (ORIGIN_EXT nt appended when
    circular), n = its length, genes = the protein genes / rRNAs found (the sensitive MiTFi pass searches outside them)."""
    if not use_ml:
        return annotate_trnas(searched, genome_fa, refdir, code, threads, tools, tmp, log, n, phylum, genes=genes)
    try:
        trnas = annotate_trnas_ml(genome, n, code, threads, tools, lineage, circular, log)
    except RuntimeError as exc:                       # cmsearch / cmalign failed: MiTFi's tRNAs, said on the console
        log.warning(f"[core] MitoM-ML failed ({str(exc)[:300]}): tRNAs by MiTFi only")
        return annotate_trnas(searched, genome_fa, refdir, code, threads, tools, tmp, log, n, phylum, genes=genes)
    mitfi = annotate_trnas(searched, genome_fa, refdir, code, threads, tools, tmp, log, n, phylum,
                           genes=genes) if tools.get("mitfi") else None
    kept = trna_second_opinion(trnas, mitfi, n)
    if len(kept) < len(trnas):
        log.info(f"[core] MitoM-ML calls not sure enough without MiTFi: "
                 f"{' '.join(f.name for f in trnas if f not in kept)}")
    trnas = kept
    if mitfi is not None:
        if vertebrate(lineage):
            trnas = mitfi_ends(trnas, mitfi, n)
            log.info(f"[core] vertebrate: tRNA ends of MiTFi for {sum('ends by MiTFi' in f.note for f in trnas)}")
        missing = [x for x in TRNA_NAMES if x not in {f.name for f in trnas}]
        if missing:
            found = rescue_trnas(trnas, mitfi, n)
            log.info(f"[core] MitoM-ML did not call {' '.join(missing)}; added from MiTFi: "
                     f"{' '.join(f.name for f in found[len(trnas):]) or 'none'}")
            trnas = found
    return trnas



# ── driver ───────────────────────────────────────────────────────────────────
def find(name: str, bin_dir: Optional[str], mitos_bin: Optional[str]) -> Optional[str]:
    return mito_env.find(name, mitos_bin if name == "runmitos" else bin_dir)


def main(argv=None) -> int:
    from Bio import SeqIO
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True, help="GenBank or FASTA file / directory")
    p.add_argument("-o", "--out", required=True)
    p.add_argument("--genetic_code", type=int, default=5)
    p.add_argument("--tools", default=None,
                   help="default: mitom,arwen,aragorn,input with --select mitom; mitos2,arwen,aragorn,trnascan,input "
                        "with mitos2_or_support / support; mitos2,input with --select mitos2")
    p.add_argument("--priority", default="MitoM-ML,MITOS2,ARWEN,ARAGORN,tRNAscan-SE,input",
                   help="source whose boundaries are used for the selected tRNA of a group")
    p.add_argument("--select", choices=["auto", "mitom", "mitos2_or_support", "support", "mitos2"], default="auto",
                   help="selected tRNAs: auto (default) = mitom when MitoM-ML is installed, else mitos2_or_support; "
                        "mitom = the MitoM-ML calls (a record whose run fails keeps its input tRNAs); mitos2 = replace the input tRNAs with the MITOS2 tRNAs "
                        "(MITOS2 boundaries and anticodons; a record whose MITOS2 run fails is kept unchanged); "
                        "mitos2_or_support (default) = found by MITOS2, or supported by >= "
                        "--min_support independent method families; support = families only. Benchmark "
                        "(39 annelid + 71 anthozoan/sponge/placozoan NCBI mitogenomes): MITOS2 alone had the "
                        "fewest false positives; ARWEN/ARAGORN alone gave hundreds in tRNA-poor genomes")
    p.add_argument("--min_support", type=int, default=2,
                   help="independent method families needed for a group to be selected (ARWEN+ARAGORN count once)")
    p.add_argument("--bin", default=None, help="directory with arwen, aragorn, tRNAscan-SE")
    p.add_argument("--mitos_bin", default=None, help="directory with runmitos (MITOS2 environment bin)")
    p.add_argument("--refdir", dest="mitos_refdir", metavar="DIR", default=None, help="folder that holds the MITOS2 reference data "
                   "refseq89m")
    old_spelling(p, "--mitos_refdir", "mitos_refdir")
    p.add_argument("--mitos_ref", default="refseq89m")
    p.add_argument("--lineage", default="",
                   help="NCBI lineage (';'-separated) for MitoM-ML clade models when the input has no GenBank taxonomy")
    p.add_argument("--threads", type=int, default=2, help="CPUs for MitoM-ML cmsearch / cmalign")
    p.add_argument("--mitos_mode", choices=["fast", "sensitive", "full"], default="fast",
                   help="fast (default): MITOS2 without its sensitive-mode fallback (scripts/runmitos_fast.py); "
                        "full: unmodified runmitos (much slower when tRNAs are missing, e.g. anthozoans)")
    a = p.parse_args(argv)
    a.mitos_refdir = mito_env.mitos_refdir(a.mitos_refdir, a.mitos_ref)

    found = {n: find(n, a.bin, a.mitos_bin) for n in ("arwen", "aragorn", "tRNAscan-SE", "runmitos")}
    found.update({n: mito_env.find(n, a.bin) or mito_env.find(n) for n in ("cmsearch", "cmalign")})
    ml_ok = mitom_ready({k: v for k, v in found.items() if v}, a.mitos_refdir)
    if a.select == "auto":
        a.select = "mitom" if ml_ok and (a.tools is None or "mitom" in a.tools.split(",")) else "mitos2_or_support"
    if a.tools is None:
        a.tools = {"mitos2": "mitos2,input", "mitom": "mitom,arwen,aragorn,input"}.get(
            a.select, "mitos2,arwen,aragorn,trnascan,input")
    want = [t.strip() for t in a.tools.split(",") if t.strip()]
    if "mitom" in want and not ml_ok:
        print("[mito_trna] mitom: MitoM-ML needs models/trna_ml, cmsearch, cmalign and MITOS2 refseq89m - skipped",
              file=sys.stderr)
        want.remove("mitom")
    need = {"mitos2": "runmitos", "arwen": "arwen", "aragorn": "aragorn", "trnascan": "tRNAscan-SE"}
    for t in list(want):
        if t in need and not found[need[t]]:
            print(f"[mito_trna] {t}: {need[t]} not found - skipped", file=sys.stderr)
            want.remove(t)
        if t == "mitos2" and t in want and not a.mitos_refdir:
            print("[mito_trna] mitos2: MITOS2 reference data (refseq89m) not found - skipped "
                  "(give --mitos_refdir, or run install.sh)", file=sys.stderr)
            want.remove(t)
    if a.select == "mitos2" and "mitos2" not in want:
        print("[mito_trna] --select mitos2 needs MITOS2 (runmitos via --mitos_bin, and --mitos_refdir)", file=sys.stderr)
        return 1
    if a.select == "mitom" and "mitom" not in want:
        print("[mito_trna] --select mitom needs MitoM-ML (models/trna_ml, cmsearch, cmalign, refseq89m)", file=sys.stderr)
        return 1
    b: Dict[str, str] = {k: v for k, v in found.items() if v}
    src, out = Path(prepare_input(a.input)), Path(a.out)
    require_records(count_seq_records(src), "sequence records", a.input)
    files = sorted(x for x in src.rglob("*") if x.is_file() and x.suffix.lower() in GB_EXT | {".fa", ".fasta", ".fna"}) \
        if src.is_dir() else [src]
    if a.select == "mitos2":
        # MITOS2 runs on the sequence of each input GenBank record, so its coordinates always match that
        # record; the input must therefore be the final, reoriented annotation (mitom_pipeline 04.final_gb),
        # never a FASTA or the round-1 output from before reorientation
        bad = [x for x in files if x.suffix.lower() not in GB_EXT or any(
            part.startswith(("00.", "01.anno_round1", "02.anno_topology", "01.mitoz_anno1")) for part in x.parts)]   # (01.mitoz_anno1: older runs)
        if bad:
            print(f"[mito_trna] --select mitos2 needs the final reoriented GenBank (e.g. 04.final_gb/), "
                  f"not: {', '.join(map(str, bad[:3]))}", file=sys.stderr)
            return 2
    out.mkdir(parents=True, exist_ok=True)
    priority = [x.strip() for x in a.priority.split(",")]
    report: List[Dict[str, str]] = []
    for f in files:
        fmt = "genbank" if f.suffix.lower() in GB_EXT else "fasta"
        recs = list(SeqIO.parse(str(f), fmt))
        review, selected = [], []
        for rec in recs:
            rec.annotations.setdefault("molecule_type", "DNA")
            with tempfile.TemporaryDirectory() as td:
                tmp = Path(td)
                fa = tmp / "g.fa"
                fa.write_text(f">{rec.id}\n{rec.seq}\n")
                preds: List[Pred] = []
                failed = []
                for t in want:
                    if t in ("input", "mitom"):
                        preds += from_input(rec) if t == "input" else []
                        continue
                    try:
                        preds += TOOLS[t](fa, a.genetic_code, b, tmp, a)
                    except RuntimeError as exc:
                        failed.append(t)
                        print(f"[mito_trna] {rec.id}: {t} failed ({exc})", file=sys.stderr)
                if "mitom" in want:                       # after ARWEN / ARAGORN: their predictions are its evidence
                    lineage = a.lineage or ";".join(rec.annotations.get("taxonomy", []))
                    try:
                        preds += mitom(str(rec.seq), lineage, a.genetic_code, b, a, preds, tmp)
                    except (RuntimeError, OSError, ValueError) as exc:
                        failed.append("mitom")
                        print(f"[mito_trna] {rec.id}: mitom failed ({exc})", file=sys.stderr)
            sel_tool = {"mitos2": "mitos2", "mitom": "mitom"}.get(a.select)
            if sel_tool and sel_tool in failed:
                # never drop the existing tRNAs because the selecting method crashed: record written unchanged
                print(f"[mito_trna] {rec.id}: {sel_tool} failed - input tRNAs kept unchanged", file=sys.stderr)
                review.append(rec)
                selected.append(rec)
                report.append({"file": f.name, "record": rec.id, "group": "", "tRNA": "", "anticodon": "", "strand": "",
                               "start": "", "end": "", "boundary_source": "", "n_families": "", "sources": "",
                               "coords_by_source": "", "conflict": "", "selected": f"{sel_tool} failed; input kept"})
                continue
            L = len(rec.seq)
            global CIRC_L
            CIRC_L = L if rec.annotations.get("topology", "circular") == "circular" else 0
            preds = unwrap(preds, L)
            groups = group(preds)
            conf = conflicts(groups)
            sources_run = sorted({x for x in ("MITOS2" if "mitos2" in want else "", "ARWEN" if "arwen" in want else "",
                                              "MitoM-ML" if "mitom" in want else "",
                                              "ARAGORN" if "aragorn" in want else "",
                                              "tRNAscan-SE" if "trnascan" in want else "",
                                              "input" if "input" in want else "") if x})
            # old tRNAs and their gene features are replaced (the core writes a gene feature per tRNA)
            trna_genes = {g for x in rec.features if x.type == "tRNA" for g in x.qualifiers.get("gene", [])}
            others = [x for x in rec.features if x.type != "tRNA" and not (
                x.type == "gene" and (set(x.qualifiers.get("gene", [])) & trna_genes
                                      or x.qualifiers.get("gene", [""])[0].lower().startswith("trn")))]
            rev_feats, sel_feats = list(others), list(others)
            for gi, g in enumerate(groups, 1):
                ac = anticodon_of(g)
                srcs = sorted({x.source for x in g})
                fams = families(g)
                n_fam_run = len({FAMILY.get(x, x) for x in sources_run})
                support = (f"{len(fams)}/{n_fam_run} independent method families; sources: {', '.join(srcs)}")
                flag = f"; CONFLICT with group(s) {conf[gi]}" if gi in conf else ""
                for x in g:
                    rev_feats.append(feature(x, ac, f"source: {x.source}; score {x.score or 'NA'}; anticodon "
                                                    f"{x.anticodon or 'NA'}; group {gi}; support {support}{flag}", L))
                best = choose(g, priority)
                # a conflict removes a group only when a conflicting group has equal or stronger
                # evidence (ARWEN-only hits overlapping a MITOS2 tRNA must not remove it)
                beaten = any(strength(groups[j - 1], a.select) >= strength(g, a.select) for j in conf.get(gi, []))
                if a.select in ("mitos2", "mitom"):
                    sel_src = "MITOS2" if a.select == "mitos2" else "MitoM-ML"
                    keep = sel_src in srcs
                    best = next(x for x in g if x.source == sel_src) if keep else best
                elif a.select == "mitos2_or_support":
                    keep = ("MITOS2" in srcs or len(fams) >= a.min_support) and not beaten
                else:
                    keep = len(fams) >= a.min_support and not beaten
                if keep:
                    ac_sel = (best.anticodon or ac) if a.select in ("mitos2", "mitom") else ac
                    tf = feature(best, ac_sel, f"tRNA by mito_trna.py: boundaries from {best.source}; "
                                               f"support {support}", L)
                    sel_feats += [SeqFeature(tf.location, type="gene", qualifiers={"gene": tf.qualifiers["gene"]}), tf]
                report.append({"file": f.name, "record": rec.id, "group": str(gi), "tRNA": trna_name(best.aa, ac),
                               "anticodon": ac, "strand": "+" if best.strand == 1 else "-",
                               "start": str(best.start), "end": str(best.end), "boundary_source": best.source,
                               "n_families": str(len(fams)), "sources": ",".join(srcs),
                               "coords_by_source": "; ".join(f"{x.source}:{x.start}..{x.end}" for x in g),
                               "conflict": ",".join(map(str, conf.get(gi, []))), "selected": "yes" if keep else "no"})
            key = lambda x: (int(x.location.start), x.type != "source", x.type != "gene")  # noqa: E731
            r1, r2 = rec[:], rec[:]
            r1.annotations, r2.annotations = dict(rec.annotations), dict(rec.annotations)
            r1.features, r2.features = sorted(rev_feats, key=key), sorted(sel_feats, key=key)
            review.append(r1)
            selected.append(r2)
        for sub, recs_out in (("review", review), ("selected", selected)):
            (out / sub).mkdir(exist_ok=True)
            SeqIO.write(recs_out, str(out / sub / f"{f.stem}.gbf"), "genbank")
    cols = ["file", "record", "group", "tRNA", "anticodon", "strand", "start", "end", "boundary_source",
            "n_families", "sources", "coords_by_source", "conflict", "selected"]
    write_tsv(out / "trna_report.tsv", cols, report)
    sel = sum(r["selected"] == "yes" for r in report)
    rule = "MITOS2 tRNAs replace input tRNAs" if a.select == "mitos2" else \
        "MitoM-ML tRNAs replace input tRNAs" if a.select == "mitom" else \
        ("found by MITOS2 or " if a.select == "mitos2_or_support" else "") + \
        f">= {a.min_support} independent method families"
    print(f"[mito_trna] {len(files)} file(s); tools {', '.join(want)}; {len(report)} tRNA loci, {sel} selected "
          f"({rule}{'' if a.select in ('mitos2', 'mitom') else '; conflicts -> stronger evidence'}) -> {out}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(cli_guard(main))
