#!/usr/bin/env python3
"""mito_nuc.py - nuclear marker genes in contigs: rDNA (18S, ITS1, 5.8S, ITS2, 28S) and histones (H1, H2A, H2B, H3, H4).

Reference-based, like MitoFinder for mitogenomes: the contigs are searched with complete NCBI annotations and
each gene boundary is projected from the reference alignments onto the contig.

rDNA: blastn (megablast, word 16) of the candidate regions against the reference rDNA units; per reference the best
collinear chain of HSPs is kept (several chains = several copies on one contig, e.g. a tandem array); every
annotated reference boundary (18S/5.8S/28S start and end) is mapped through the alignment column by column
(or projected <= EXTRAP nt past the end of an HSP) and the TOPN best references vote, weighted by chain bit
score; the position with the most weight wins. A boundary that no reference reaches is projected without
alignment from the chain end; when that lands outside the contig the gene is partial (< / >). ITS1 / ITS2 are
the spacers between the rRNAs.
Histones: tblastn of the reference proteins, HSPs grouped into loci, exonerate protein2genome and protein2dna
(intronless) with the TOPK best proteins of the locus, best score wins (spliced models x INTRON_PEN per
intron); the CDS is extended in frame to the next stop codon and, when the model does not start with ATG at
protein position 1, to the farthest upstream in-frame ATG before a stop. Nuclear code (transl_table 1). Each copy
is classed replication-dependent (tandem cluster) / replacement (H3.3, spliced) / single-copy from H3 residues
31/87/89/90, introns and its coverage relative to the single-copy coverage (histone_class); the per-gene consensus
takes a cluster copy first (pick_histone).
Massive assemblies: both searches first run with one reference per genus (rDNA; low-complexity stretches do not
seed) / per gene and phylum (histones) against a BLAST database of all contigs; only the regions hit (+-PAD nt) go
to the steps above, each gene set in its own regions, and at most MAX_RDNA_REGIONS rDNA regions and
MAX_HIST_REGIONS histone regions per gene of them (best first; the time is bounded whatever the assembly).
Time and memory: docs/METHODS.md, nuclear markers.

Per-gene consensus (multi-copy rDNA split over many contigs of one assembly): see consensus(); output
consensus/<gene>.fasta and consensus.tsv (which contigs were used, which were dropped and why).

Accuracy (benchmark/nuc_bench.py, nuc_gold_bench.py; knowledge/topics/nuclear_markers.md): docs/METHODS.md,
nuclear markers.

Usage:
  mitom nuc fetch -o nuc_refs                           # NCBI GenBank records: rdna.gb, histone.gb
  mitom nuc build --gb nuc_refs/rdna.gb nuc_refs/histone.gb [my_geneious.gb] -o models/nuc
  mitom nuc annotate -i contigs.fasta -o out [--genes rdna,histone] [--threads 16] [--lib models/nuc]
Outputs of annotate in -o: nuc_features.tsv (one row per feature), nuc_annotation.gb (contigs with features),
genes/<gene>.fasta (every copy; + histone .faa), consensus/<gene>.fasta + consensus.tsv (one sequence per
gene), nuc.log.
"""
from __future__ import annotations

import argparse
import csv
import logging
import os
import re
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq, UndefinedSequenceError
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqFeature import BeforePosition, AfterPosition
from Bio.SeqRecord import SeqRecord

import mito_env  # noqa: E402
import mito_common  # noqa: E402
from mito_nuc_lib import (  # noqa: E402, F401  (re-exported: _label, _partial, histone_name, rrna_name)
    _contig,
    _cov,
    _label,
    _partial,
    _ref_gene,
    build,
    DB_PIECE,
    END_TOL,
    fetch,
    histone_name,
    HISTONES,
    LIB,
    Library,
    MAX_CONS_SEED,
    MAX_HIST_REGIONS,
    MAX_RDNA_REGIONS,
    PAD,
    PREFILTER_NT,
    PRODUCT,
    program_output,
    RRNA,
    rrna_name,
    SEED_COV,
    STOPS,
    tool,
    UNIT,
)
from mito_nuc_rdna import (  # noqa: E402, F401  (re-exported: _cm_rescue_58s, _cmsearch, _tidy, cm_ends, cm_refine, rdna_loci, RFAM)
    _cm_rescue_58s,
    _cmsearch,
    _tidy,
    annotate_rdna,
    cm_ends,
    cm_refine,
    rdna_loci,
    RFAM,
)
from mito_nuc_histone import (  # noqa: E402, F401  (re-exported: _exonerate, _exonerate_all, _orf, _tblastn, H3_SITES)
    _exonerate,
    _exonerate_all,
    _orf,
    _tblastn,
    annotate_histones,
    CLASS_RANK,
    H3_SITES,
    histone_class,
    single_copy_cov,
)
from mito_common import cli_guard, count_seq_records, InputError, prepare_input, require_records, run_tool, write_gb  # noqa: E402


# ── prefilter: candidate regions of the input ──────────────────────────────
def blast_db(fa: Path, db: Path, threads: int) -> str:
    """nucleotide BLAST database of a large FASTA: makeblastdb on `threads` pieces at once, joined by an alias
    database (same total length and sequence count, so the same e-values as one database); one makeblastdb of a
    460 Mb SPAdes assembly took 25-35 s. BLAST threads share the work by volume: one volume per thread"""
    if any(c.isspace() for c in str(db)):            # makeblastdb / blastdb_aliastool cut paths at spaces (D18)
        raise InputError(f"BLAST cannot write a database into a path with spaces: {db.parent} (use another -o)")
    size = fa.stat().st_size
    n = max(1, min(threads, size // DB_PIECE))
    if n == 1:
        program_output([tool("makeblastdb"), "-in", str(fa), "-dbtype", "nucl", "-out", str(db)])
        return str(db)
    parts = [db.parent / f"{db.name}_{i}" for i in range(n)]
    with open(fa, "rb") as f:
        cuts = [0]                      # pieces start at the first record after each n-th of the file
        for k in range(1, n):
            f.seek(max(cuts[-1], k * size // n))
            pos, buf = f.tell(), b""
            while b"\n>" not in buf:
                more = f.read(1 << 20)
                if not more:
                    break
                buf += more
            if b"\n>" in buf:
                cuts.append(pos + buf.index(b"\n>") + 1)
        cuts.append(size)
        for p, a, b in zip(parts, cuts, cuts[1:]):
            f.seek(a)
            with open(f"{p}.fa", "wb") as o:
                while a < b:
                    chunk = f.read(min(b - a, 1 << 24))
                    if not chunk:                # the file shrank after stat() (review 11-03 #4, 2026-11-05)
                        raise InputError(f"{fa} changed while it was split for BLAST: run again")
                    o.write(chunk)
                    a += len(chunk)
    parts = parts[:len(cuts) - 1]

    def mk(p: Path) -> None:
        program_output([tool("makeblastdb"), "-in", f"{p}.fa", "-dbtype", "nucl", "-out", str(p)])
        Path(f"{p}.fa").unlink()
    with ThreadPoolExecutor(len(parts)) as pool:
        list(pool.map(mk, parts))
    program_output([tool("blastdb_aliastool"), "-dblist", " ".join(map(str, parts)), "-dbtype", "nucl", "-out", str(db),
         "-title", db.name])
    return str(db)


def _windows(hits: List[List[str]]) -> Dict[str, List[Tuple[int, int]]]:
    """contig -> merged (start, end) windows around the hits (+- PAD, clipped); hits: [sseqid, sstart, send, slen, ...]"""
    win = defaultdict(list)
    for s, a, b, n, *_ in hits:
        a, b, n = sorted((int(a), int(b))) + [int(n)]
        win[s].append((max(1, a - PAD), min(n, b + PAD)))
    merged = {}
    for s, ws in win.items():
        ws.sort()
        m = [list(ws[0])]
        for a, b in ws[1:]:
            if a <= m[-1][1] + 1:
                m[-1][1] = max(m[-1][1], b)
            else:
                m.append([a, b])
        merged[s] = [tuple(x) for x in m]
    return merged


def _best_windows(windows: Dict[str, List[Tuple[int, int]]], hits: List[List[str]], score, cap: int,
                  what: str) -> Dict[str, List[Tuple[int, int]]]:
    """the `cap` windows with the highest score(hits inside the window) (ties: contig coverage), all when cap is 0;
    a fixed amount of work whatever the number of candidates (2,800 rDNA contigs of 469 species: the consensus alone
    ran > 10 min)"""
    inside = defaultdict(list)
    for h in hits:
        a = min(int(h[1]), int(h[2]))
        ws = windows[h[0]]
        inside[(h[0], next(w for w in reversed(ws) if w[0] <= a))].append(h)
    if not cap or len(inside) <= cap:
        return windows
    keep = sorted(inside, key=lambda k: (score(inside[k]), float(_cov(k[0]) or 0), k), reverse=True)[:cap]
    logging.info(f"[NUC] {what}: {len(inside)} candidate regions, the {cap} best kept ({len(inside) - cap} skipped)")
    out: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for s, w in sorted(keep):
        out[s].append(w)
    return dict(out)


def prefilter(contigs: Path, lib: Library, genes: List[str], work: Path, threads: int,
              n_seq: int) -> Tuple[Dict[str, Dict[str, List[Tuple[int, int]]]], set]:
    """gene set -> contig -> windows hit by a representative reference, and the gene sets whose work is bounded
    (assembly: candidates < 10 % of the n_seq input sequences); each set is searched in its own windows.
    rDNA: low-complexity stretches do not seed (dust): the flanks of the reference units (IGS / ETS) hold
    microsatellites ((CA)n, (CATA)n, (CT)n) that hit ~10k contigs of a 460 Mb assembly, none of them with an rRNA gene;
    searching them took 20 of 22 min and 16 GB (SRR27947581). At most MAX_RDNA_REGIONS rDNA regions (most aligned
    bits of one reference first) and MAX_HIST_REGIONS histone regions per gene (closest to a reference protein
    first) of an assembly go on, so the time is bounded"""
    db = blast_db(contigs, work / "contigs", threads)
    out, bounded = {}, set()
    if "rdna" in genes:
        hits = [x.split("\t") for x in program_output([tool("blastn"), "-task", "megablast", "-word_size", "20", "-query",
                                             str(lib.rdna_rep), "-db", db, "-evalue", "1e-20", "-dust", "yes",
                                             "-max_target_seqs", "100000", "-num_threads", str(threads), "-outfmt",
                                             "6 sseqid sstart send slen qseqid bitscore"]).splitlines()]

        def bits(hs: List[List[str]]) -> float:         # aligned bits of the best reference in the region
            per: Dict[str, float] = defaultdict(float)
            for h in hs:
                per[h[4]] += float(h[5])
            return max(per.values())
        win = _windows(hits)
        if 10 * sum(map(len, win.values())) < n_seq:
            bounded.add("rdna")
        out["rdna"] = _best_windows(win, hits, bits, MAX_RDNA_REGIONS if "rdna" in bounded else 0, "rDNA")
    if "histone" in genes:
        hits = [x.split("\t") for x in program_output([tool("tblastn"), "-query", str(lib.hist_rep), "-db", db, "-evalue", "1e-10",
                                             "-seg", "no", "-max_target_seqs", "100000", "-num_threads", str(threads),
                                             "-outfmt", "6 sseqid sstart send slen qseqid bitscore pident"]).splitlines()]
        win, keep = _windows(hits), set()
        if 10 * sum(map(len, win.values())) < n_seq:
            bounded.add("histone")
        for g in HISTONES:
            gh = [h for h in hits if lib.hist[h[4]]["gene"] == g]
            for s, ws in _best_windows(_windows(gh), gh, lambda hs: max((float(h[6]), float(h[5])) for h in hs),
                                       MAX_HIST_REGIONS if "histone" in bounded else 0, f"histone {g}").items():
                # a window of one gene lies inside a window of all genes, which is kept
                keep.update((s, w) for a, b in ws for w in win[s] if w[0] <= a and b <= w[1])
        out["histone"] = {s: [w for w in ws if (s, w) in keep] for s, ws in win.items()
                          if any((s, w) in keep for w in ws)}
    return out, bounded


def _offset(region_id: str) -> int:
    m = re.search(r":(\d+)-\d+$", region_id)
    return int(m.group(1)) - 1 if m else 0


class FastaIndex:
    """records of a FASTA file by id (first word of the title), read on demand like SeqIO.index (which took 5 s on
    the 900k contigs of a 460 Mb assembly); the file is scanned in blocks, with Python work per record title only.
    plain: every title is a BLAST-safe id alone, the ids are unique, every record has sequence and there is no space,
    tab or CR (read_input would write the same records); nt: number of nucleotides of a plain file"""
    BLOCK = 1 << 26
    TITLE = re.compile(b"\n>([^\n]*)")

    def __init__(self, path: Path):
        self.path, self.nt, self.plain = path, 0, True
        titles: List[bytes] = []
        starts: List[int] = []
        with open(path, "rb") as f:
            base, rest = 0, b"\n"           # a record starts after a newline: the file gets one in front
            while True:
                b = f.read(self.BLOCK)
                data = rest + b
                cut = data.rfind(b"\n>") if b else len(data)
                if cut <= 0:                        # a record longer than the block
                    rest = data
                    continue
                t = self.TITLE.findall(data, 0, cut)
                starts += [base + m.start() for m in re.finditer(b"\n>", data[:cut])]
                if any(data.find(c, 0, cut + 1) >= 0 for c in (b" ", b"\t", b"\r", b"\n\n")):
                    self.plain = False
                    t = [(x.split(None, 1) or [b""])[0] for x in t]
                self.nt += cut - sum(map(len, t)) - len(t) - data.count(b"\n", 0, cut)
                titles += t
                base += cut
                rest = data[cut:]
                if not b:
                    break
        ends = starts[1:] + [base - 1]
        joined = b"\n".join(titles)
        # every title a BLAST-safe id alone, unique, followed by sequence (a title line is len + 2 bytes)
        if self.plain and (not starts or starts[0] != 0 or b"" in titles or
                           re.search(b"[^A-Za-z0-9_.\n-]", joined) or
                           any(e - a <= len(x) + 2 for a, e, x in zip(starts, ends, titles))):
            self.plain = False
        self.at: Dict[str, Tuple[int, int]] = dict(zip(joined.decode(errors="replace").split("\n"),
                                                       zip(starts, ends)))
        self.plain = self.plain and len(titles) == len(self.at)
        self.fd = os.open(path, os.O_RDONLY)

    def __getitem__(self, k: str) -> SeqRecord:
        a, b = self.at[k]
        title, _, seq = os.pread(self.fd, b - a, a).partition(b"\n")
        return SeqRecord(Seq(seq.translate(None, b" \t\r\n").decode()), id=k, name=k,
                         description=title[1:].decode(errors="replace").rstrip())

    def keys(self):
        return self.at.keys()

    def close(self) -> None:
        """close the file descriptor; calling it again does nothing"""
        fd, self.fd = self.fd, None
        if fd is not None:
            os.close(fd)

    def __enter__(self) -> "FastaIndex":
        return self

    def __exit__(self, *exc) -> None:
        self.close()


def open_input(path: Path, work: Path) -> Tuple[Path, Optional[FastaIndex], Dict[str, str], int, int, set]:
    """read_input, and its index when it has one already: a plain FASTA file (SPAdes scaffolds) is used as it is
    (rewriting a 460 Mb assembly took 15 s)"""
    with open(path, "rb") as f:
        first = f.readline()
    if path.suffix.lower() not in (".gb", ".gbk", ".gbf", ".genbank", ".gbff") and not re.search(b"[ \t\r]", first):
        idx = FastaIndex(path)
        if idx.plain:
            return path, idx, {}, len(idx.at), idx.nt, set()
        idx.close()
    fa, org, n, nt, circ = read_input(path, work)
    return fa, None, org, n, nt, circ


def read_input(path: Path, work: Path) -> Tuple[Path, Dict[str, str], int, int, set]:
    """FASTA or GenBank -> FASTA with unique BLAST-safe ids; returns it, organism by id (GenBank input),
    number of sequences and of nucleotides, and the ids of circular sequences (GetOrganelle paths named
    '...(circular)', GenBank topology circular). Records are read back on demand (FastaIndex), so an assembly of
    millions of contigs is never held in memory."""
    fmt = "genbank" if path.suffix.lower() in (".gb", ".gbk", ".gbf", ".genbank", ".gbff") else "fasta"
    ids, org, nt, circ = set(), {}, 0, set()
    fa = work / "input.fa"
    with open(fa, "w") as o:
        for r in SeqIO.parse(str(path), fmt):
            rid = re.sub(r"[^A-Za-z0-9_.-]", "_", r.id)
            try:
                bytes(r.seq)
            except UndefinedSequenceError:     # CON / WGS master records list pieces, no sequence
                logging.warning(f"[NUC] {r.id}: no sequence in the record, skipped")
                continue
            if rid in ids:
                raise SystemExit(f"duplicate sequence id {rid} in {path}")
            ids.add(rid)
            nt += len(r.seq)
            if "organism" in r.annotations:
                org[rid] = r.annotations["organism"]
            if r.annotations.get("topology") == "circular" or r.id.endswith("(circular)"):
                circ.add(rid)
            desc = r.description[len(r.id):].strip() if r.description.startswith(r.id) else r.description
            o.write(f">{rid} {desc}\n{r.seq}\n" if desc else f">{rid}\n{r.seq}\n")
    return fa, org, len(ids), nt, circ


def _cut(n: int, spans: List[Tuple[int, int]]) -> Optional[int]:
    """0-based new start of a circular sequence of n nt: the middle of the longest stretch no feature covers.
    spans: 1-based feature ends, positions taken modulo n (features of the sequence repeated); None without
    features or without a free stretch."""
    cov = bytearray(n)
    for a, b in spans:
        if b - a + 1 >= n:
            return None
        for k in range(a - 1, b):
            cov[k % n] = 1
    if not spans or all(cov):
        return None
    s0, run, best = cov.index(1), 0, (0, 0)     # best: (length, start) of the longest free run
    for i in range(1, n + 1):
        k = (s0 + i) % n
        run = 0 if cov[k] else run + 1
        if run > best[0]:
            best = (run, (k - run + 1) % n)
    return (best[1] + best[0] // 2) % n


def rotate_circular(fa: Path, circ: set, lib: "Library", genes: List[str], work: Path, threads: int,
                    leave_out: Dict[str, str]) -> None:
    """A circular sequence (a GetOrganelle rDNA unit) starts anywhere, also inside a gene, which then came out cut
    in two partial pieces. Each is searched three times over (no gene cut by its ends) and fa rewritten with it
    turned to start in the middle of its longest stretch without genes (the IGS of an rDNA unit)."""
    d = work / "circular"
    d.mkdir()
    seqs = {r.id: r.seq for r in SeqIO.parse(str(fa), "fasta") if r.id in circ}
    logging.info(f"[NUC] {len(seqs)} circular sequence(s): genes searched first to start them between genes")
    with open(d / "x3.fa", "w") as o:
        for k, s in seqs.items():
            o.write(f">{k}\n{s * 3}\n")
    spans: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    if "rdna" in genes:
        loci = annotate_rdna(d / "x3.fa", lib, d, threads, leave_out, {k: 3 * len(s) for k, s in seqs.items()}, False)
        for q, ls in loci.items():
            spans[q] += [(f["lo"], f["hi"]) for L in ls for f in L["feats"].values()]
    if "histone" in genes:
        hist = annotate_histones(d / "x3.fa", lib, d, threads, leave_out, {k: s * 3 for k, s in seqs.items()})
        for q, ms in hist.items():
            spans[q] += [(m["exons"][0][0], m["exons"][-1][1]) for m in ms]
    cut = {}
    for k, s in seqs.items():
        c = _cut(len(s), spans[k])
        if c is None and spans[k]:
            logging.warning(f"[NUC] {k}: circular, genes all around it: not turned")
        elif c:
            cut[k] = c
            logging.info(f"[NUC] {k}: circular, turned to start at nt {c + 1}, between genes")
    if not cut:
        return
    tmp = work / "input_turned.fa"
    with open(tmp, "w") as o:
        for r in SeqIO.parse(str(fa), "fasta"):
            if r.id in cut:
                c = cut[r.id]
                r = SeqRecord(r.seq[c:] + r.seq[:c], id=r.id,
                              description=f"{r.description} (circular, starts at nt {c + 1} of the input)")
            SeqIO.write(r, o, "fasta")
    tmp.replace(fa)


def write_regions(windows: Dict[str, List[Tuple[int, int]]], recs, fa: Path) -> Dict[str, Seq]:
    """whole contig when one window covers it, else contig:start-end pieces"""
    seqs = {}
    with open(fa, "w") as o:
        for c, ws in windows.items():
            s = recs[c].seq
            for a, b in ws:
                rid = c if (a == 1 and b == len(s)) else f"{c}:{a}-{b}"
                seqs[rid] = s[a - 1:b]
                o.write(f">{rid}\n{seqs[rid]}\n")
    return seqs


def _loc(a: int, b: int, st: int, p_lo: bool, p_hi: bool) -> FeatureLocation:
    return FeatureLocation(BeforePosition(a - 1) if p_lo else a - 1, AfterPosition(b) if p_hi else b, strand=st)


def collect(recs, rdna, hist, base_cov: Optional[float] = None) -> Tuple[List[dict], Dict[str, List[SeqFeature]]]:
    rows, feats = [], defaultdict(list)

    def ratio(c: str) -> Optional[float]:
        cov = _cov(c)
        return float(cov) / base_cov if cov and base_cov else None
    for rid, loci in (rdna or {}).items():
        c, off = _contig(rid), _offset(rid)
        for L in loci:
            for g in UNIT:
                if g not in L["feats"]:
                    continue
                f = L["feats"][g]
                a, b = f["lo"] + off, f["hi"] + off
                n = len(recs[c].seq)
                p_lo, p_hi = f["p_lo"], f["p_hi"]
                p5, p3 = (p_lo, p_hi) if L["st"] == 1 else (p_hi, p_lo)
                ftype = "rRNA" if g in RRNA else "misc_RNA"
                note = f"boundaries {f['how']}; best reference {L['ref']} {L['ref_org']} ({L['identity']} % identity)"
                feats[c].append(SeqFeature(_loc(a, b, L["st"], p_lo, p_hi), type=ftype,
                                           qualifiers={"product": [PRODUCT[g]], "note": [note]}))
                rows.append({"contig": c, "contig_len": n, "contig_cov": _cov(c), "type": ftype, "gene": g,
                             "start": a, "end": b, "strand": "+" if L["st"] == 1 else "-",
                             "partial_5": int(p5), "partial_3": int(p3), "length": b - a + 1, "exons": f"{a}-{b}",
                             "best_ref": L["ref"], "best_ref_organism": L["ref_org"], "identity": L["identity"],
                             "evidence": f["how"], "aligned_nt": f.get("aligned", ""),
                             "copy_ratio": f"{ratio(c):.1f}" if ratio(c) is not None else ""})
    for rid, models in (hist or {}).items():
        c, off = _contig(rid), _offset(rid)
        for m in models:
            ex = [(x + off, y + off) for x, y in m["exons"]]
            lo_p, hi_p = (m["p5"], m["p3"]) if m["st"] == 1 else (m["p3"], m["p5"])
            parts = [_loc(x, y, m["st"], lo_p and i == 0, hi_p and i == len(ex) - 1) for i, (x, y) in enumerate(ex)]
            if m["st"] == -1:
                parts = parts[::-1]
            loc = parts[0] if len(parts) == 1 else CompoundLocation(parts)
            seq = loc.extract(recs[c].seq)
            cs = m.get("codon_start", 1)
            seq = seq[cs - 1:]
            prot = str(seq[: len(seq) // 3 * 3].translate(table=1)).rstrip("*")
            rt = ratio(c)
            cls, why = histone_class(m["gene"], prot, len(ex), rt)
            q = {"gene": [m["gene"]], "product": [PRODUCT[m["gene"]]], "codon_start": [str(cs)], "transl_table": ["1"],
                 "note": [f"exonerate model from {m['ref']} {m['ref_org']}"
                          + (f"; {cls} histone ({why})" if cls else "")]}
            if "*" not in prot:
                q["translation"] = [prot]
            feats[c].append(SeqFeature(loc, type="CDS", qualifiers=q))
            rows.append({"contig": c, "contig_len": len(recs[c].seq), "contig_cov": _cov(c), "type": "CDS",
                         "gene": m["gene"], "start": ex[0][0], "end": ex[-1][1], "strand": "+" if m["st"] == 1 else "-",
                         "partial_5": int(m["p5"]), "partial_3": int(m["p3"]), "length": sum(y - x + 1 for x, y in ex),
                         "exons": ";".join(f"{x}-{y}" for x, y in ex), "best_ref": m["ref"],
                         "best_ref_organism": m["ref_org"], "identity": "",
                         "evidence": "internal stop codon" if "*" in prot else "exonerate", "aligned_nt": "",
                         "copy_ratio": f"{rt:.1f}" if rt is not None else "", "class": cls, "class_evidence": why, "codon_start": cs})
    return rows, feats


def write_outputs(out: Path, recs, org: Dict[str, str], rows: List[dict], feats: Dict[str, List[SeqFeature]]) -> None:
    rows.sort(key=lambda r: (r["contig"], r["start"]))
    cols = ["contig", "contig_len", "contig_cov", "type", "gene", "start", "end", "strand", "partial_5", "partial_3",
            "length", "exons", "best_ref", "best_ref_organism", "identity", "evidence", "aligned_nt", "copy_ratio",
            "class", "class_evidence", "codon_start"]
    with open(out / "nuc_features.tsv", "w", newline="") as o:
        w = csv.DictWriter(o, cols, delimiter="\t", restval="")
        w.writeheader()
        w.writerows(rows)
    gbrecs = []
    for c, fs in sorted(feats.items()):          # fixed order: identical files whatever the thread count
        r = recs[c]
        g = SeqRecord(r.seq, id=c[:16] if len(c) > 16 else c, name=c[:16], description=r.description,
                      annotations={"molecule_type": "DNA", "organism": org.get(c, "unknown"),
                                   "topology": "linear"})
        g.features = sorted(fs, key=lambda f: int(f.location.start))
        gbrecs.append(g)
    write_gb(gbrecs, out / "nuc_annotation.gb")
    gd = out / "genes"
    gd.mkdir(exist_ok=True)
    by_gene = defaultdict(list)
    for c, fs in sorted(feats.items()):
        for f in fs:
            g = f.qualifiers.get("gene", [None])[0] or \
                next(k for k, v in PRODUCT.items() if v == f.qualifiers["product"][0])
            s = f.location.extract(recs[c].seq)
            tag = f"{c}:{int(f.location.start) + 1}-{int(f.location.end)}({'+' if f.location.strand == 1 else '-'})"
            by_gene[g].append(SeqRecord(s, id=tag, description=""))
            if "translation" in f.qualifiers:
                by_gene[g + ".faa"].append(SeqRecord(Seq(f.qualifiers["translation"][0]), id=tag, description=""))
    for g, rs in by_gene.items():
        SeqIO.write(rs, str(gd / (g if g.endswith(".faa") else f"{g}.fasta")), "fasta")


# ── NCBI submission files ──────────────────────────────────────────────────
# source qualifiers carried from GenBank input to the table2asn FASTA definition line
SRC_MODS = ("isolate", "specimen_voucher", "strain", "clone", "culture_collection", "bio_material", "country",
            "geo_loc_name", "lat_lon", "collection_date", "collected_by", "identified_by", "host", "tissue_type",
            "dev_stage", "sex", "altitude")
DEF_NAME = {"18S": "18S ribosomal RNA gene", "5.8S": "5.8S ribosomal RNA gene", "28S": "28S ribosomal RNA gene",
            "ITS1": "internal transcribed spacer 1", "ITS2": "internal transcribed spacer 2"}


def _join(items: List[str]) -> str:
    return items[0] if len(items) == 1 else f"{items[0]} and {items[1]}" if len(items) == 2 else \
        ", ".join(items[:-1]) + ", and " + items[-1]


def definition(organism: str, src: Dict[str, str], feats: List[SeqFeature]) -> str:
    """NCBI rRNA / ITS definition line: parts in sequence order, consecutive parts of equal completeness grouped
    ('18S ribosomal RNA gene, internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed
    spacer 2, complete sequence; and 28S ribosomal RNA gene, partial sequence'); histone CDS as '... gene(s),
    complete / partial cds'"""
    groups: List[list] = []          # [kind, status, names]
    for f in sorted(feats, key=lambda f: int(f.location.start)):
        part = isinstance(f.location.start, BeforePosition) or isinstance(f.location.end, AfterPosition) or \
            any(isinstance(p.start, BeforePosition) or isinstance(p.end, AfterPosition) for p in f.location.parts)
        if f.type == "CDS":
            kind, name = "cds", f.qualifiers["product"][0]
        elif f.type == "gene":                       # /pseudo
            kind, name = "rna", f"{PRODUCT[f.qualifiers['gene'][0]]} pseudogene"
        else:
            g = next(k for k, v in PRODUCT.items() if v == f.qualifiers["product"][0])
            kind, name = "rna", DEF_NAME[g]
        status = "partial" if part else "complete"
        if groups and groups[-1][:2] == [kind, status]:
            groups[-1][2].append(name)
        else:
            groups.append([kind, status, [name]])
    parts = []
    for kind, status, names in groups:
        if kind == "cds":
            parts.append(f"{_join(names)} gene{'s' if len(names) > 1 else ''}, {status} cds")
        else:
            parts.append(f"{_join(names)}, {status} sequence")
    text = parts[0] if len(parts) == 1 else "; ".join(parts[:-1]) + "; and " + parts[-1]
    tag = next((f"{lab} {src[k]}" for k, lab in (("specimen_voucher", "voucher"), ("isolate", "isolate"),
                                                 ("strain", "strain"), ("clone", "clone")) if src.get(k)), "")
    return " ".join(x for x in (organism, tag, text) if x)


def _source_quals(inp: Path, ids: set) -> Dict[str, Dict[str, str]]:
    if inp.suffix.lower() not in (".gb", ".gbk", ".gbf", ".genbank", ".gbff"):
        return {}
    out = {}
    for r in SeqIO.parse(str(inp), "genbank"):
        rid = re.sub(r"[^A-Za-z0-9_.-]", "_", r.id)
        if rid in ids:
            src = next((f for f in r.features if f.type == "source"), None)
            out[rid] = {k: src.qualifiers[k][0] for k in SRC_MODS if src is not None and k in src.qualifiers}
    return out


def write_submission(out: Path, inp: Path, recs, org: Dict[str, str], feats: Dict[str, List[SeqFeature]],
                     organism: Optional[str] = None, src_table: Optional[Path] = None) -> Optional[Path]:
    """submission/ : nuc.fsa (table2asn source modifiers + NCBI definition line), nuc.tbl (5-column feature
    table without the tool's notes), sequences with every rRNA on the minus strand reverse-complemented; with
    table2asn installed nuc.sqn is made and validated (-V vb -Z: nuc.val, nuc.dr, validation.txt); runs of >= 10 N
    become gaps there (-gaps-min 10, 100 N = unknown length), as NCBI asks"""
    import ncbi_validate as nv
    d = out / "submission"
    d.mkdir(exist_ok=True)
    src = _source_quals(inp, set(feats))
    if src_table:
        for r in csv.DictReader(open(src_table), delimiter="\t"):
            sid = r.pop("Sequence_ID", None) or r.pop("sequence_id", None)
            if sid:
                src.setdefault(sid, {}).update({k: v for k, v in r.items() if v})
    fsa, tbl, missing = [], [], 0
    for c in sorted(feats):
        fs = sorted(feats[c], key=lambda f: int(f.location.start))
        rec = SeqRecord(recs[c].seq, id=c, description="")
        rec.features = []
        for f in fs:
            q = {k: v for k, v in f.qualifiers.items() if k not in ("note", "translation")}
            if f.type == "CDS":
                # linked to its gene by overlap: a /gene name on the CDS makes table2asn take the first gene of
                # that name in a cluster (CDSgeneRange, GeneXrefStrandProblem)
                q.pop("gene", None)
                cds = str(f.location.extract(recs[c].seq)[int(q["codon_start"][0]) - 1:]).upper()
                last = f.location.parts[-1]                  # 3'-most part (Biopython keeps transcript order)
                p3 = isinstance(last.end, AfterPosition) if f.location.strand != -1 else \
                    isinstance(last.start, BeforePosition)
                if not p3 and (len(cds) % 3 or cds[-3:] not in STOPS):
                    f.qualifiers.pop("translation", None)     # complete 3' end without a stop codon: pseudogene
            if f.type != "CDS":
                rec.features.append(SeqFeature(f.location, type=f.type, qualifiers=q))
                continue
            # one gene span per CDS (table2asn MultiIntervalGene otherwise); a model with an internal stop codon
            # is a pseudogene: gene /pseudo without CDS
            loc = f.location
            span = FeatureLocation(loc.parts[0].start if loc.strand != -1 else loc.parts[-1].start,
                                   loc.parts[-1].end if loc.strand != -1 else loc.parts[0].end, strand=loc.strand)
            gq = {"gene": f.qualifiers["gene"]}
            if "translation" not in f.qualifiers:
                gq.update(pseudo=[""], note=["nonfunctional due to internal stop codon"])
            rec.features.append(SeqFeature(span, type="gene", qualifiers=gq))
            if "translation" in f.qualifiers:
                rec.features.append(SeqFeature(loc, type="CDS", qualifiers=q))
        rna = [f for f in rec.features if f.type in ("rRNA", "misc_RNA")]
        if rna and all(f.location.strand == -1 for f in rna):
            rec = rec.reverse_complement(id=c, description="", features=True)
        sq = dict(src.get(c, {}))
        name = sq.pop("organism", None) or org.get(c) or organism
        if not name:
            missing += 1
            name = "unknown"
        mods = [f"[organism={name}]", "[moltype=genomic DNA]"] + \
            [f"[{k}={v}]" for k, v in sq.items() if k != "title"]
        if any(f.type == "CDS" for f in rec.features):
            mods.append("[gcode=1]")
        title = sq.get("title") or definition(name, sq, [f for f in rec.features if f.type != "gene"] +
                                              [f for f in rec.features if "pseudo" in f.qualifiers])
        fsa.append(f">{c} {' '.join(mods)} {title}\n{rec.seq}\n")
        tbl.append("\n".join(nv.feature_lines(rec)) + "\n")
    (d / "nuc.fsa").write_text("".join(fsa))
    (d / "nuc.tbl").write_text("".join(tbl))
    if missing:
        logging.warning(f"[NUC] submission: {missing} sequence(s) without organism (give --organism or --src)")
    t2a = mito_env.find("table2asn")
    if not t2a or not fsa:
        return d
    val = d / "nuc.val"
    val.unlink(missing_ok=True)                        # never report the last run's messages
    try:
        p = run_tool([t2a, "-i", str(d / "nuc.fsa"), "-f", str(d / "nuc.tbl"), "-V", "vb", "-Z", "-a", "s",
                      "-gaps-min", "10", "-gaps-unknown", "100",
                      "-o", str(d / "nuc.sqn")], timeout=3600)
    except subprocess.TimeoutExpired:                  # (table2asn 1.28 can hang on some SeqIDs)
        logging.warning("[NUC] submission: table2asn did not finish in 1 h: not validated")
        return d
    if not val.exists():                               # a crash is 'not validated', never '0 error(s)' (2026-11-05)
        logging.warning(f"[NUC] submission: table2asn exit {p.returncode}, no nuc.val: not validated - "
                        f"{' '.join(p.stderr.split())[-200:]}")
        (d / "validation.txt").write_text(f"NOT VALIDATED: table2asn exit {p.returncode}\n")
        return d
    msgs = Counter()
    for line in val.read_text(errors="ignore").splitlines() if val.exists() else []:
        m = nv.MSG_RE.match(line.strip())
        if m and not nv.IGNORE_CODES.search(m.group(2)):
            msgs[(m.group(1).upper(), m.group(2))] += 1
    lines = [f"table2asn -V vb -Z on {len(fsa)} sequence(s); submitter / publication messages not counted "
             "(add your template with table2asn -t template.sbt)"]
    lines += [f"{sev}\t{code}\t{n}" for (sev, code), n in msgs.most_common()] or ["no sequence or feature messages"]
    (d / "validation.txt").write_text("\n".join(lines) + "\n")
    err = sum(n for (sev, _), n in msgs.items() if sev in ("ERROR", "REJECT", "CRITICAL", "FATAL"))
    logging.info(f"[NUC] submission: {d / 'nuc.fsa'} + nuc.tbl; table2asn {err} error(s), "
                 f"{sum(msgs.values()) - err} warning(s)/info  ({d / 'validation.txt'}, discrepancy report nuc.dr)")
    return d


def _gene_seq(r: dict, recs) -> Seq:
    s = recs[r["contig"]].seq
    seq = sum((s[int(a) - 1:int(b)] for a, b in (x.split("-") for x in r["exons"].split(";"))), Seq(""))
    return seq if r["strand"] == "+" else seq.reverse_complement()


def _pileup(frags: List[str], wt: List[float], backbone: str, work: Path) -> Tuple[str, Dict[int, Tuple[float, int]]]:
    """consensus of the fragments placed on the backbone by their best blastn HSP (plus strand): per backbone
    position the weighted majority base or deletion; an insertion is kept when most of the spanning weight has
    it. Uppercase backbone letters are sample sequence, lowercase ones only a reference placeholder: positions no
    fragment spans keep the backbone letter in lowercase. Returns the consensus and, per fragment, (identity,
    aligned columns) over the uppercase backbone positions."""
    (work / "bb.fa").write_text(f">bb\n{backbone}\n")
    with open(work / "fr.fa", "w") as o:
        for i, f in enumerate(frags):
            o.write(f">{i}\n{f}\n")
    out = program_output([tool("blastn"), "-task", "blastn", "-query", str(work / "fr.fa"), "-subject", str(work / "bb.fa"),
               "-strand", "plus", "-dust", "no", "-evalue", "1e-5",
               "-outfmt", "6 qseqid qstart qend sstart send bitscore qseq sseq"])
    n = len(backbone)
    best: Dict[int, tuple] = {}
    for line in out.splitlines():
        q, qs, qe, ss, se, bs, qa, sa = line.split("\t")
        q, qs, qe, ss, se = int(q), int(qs), int(qe), int(ss), int(se)
        if q in best and float(bs) <= best[q][2]:
            continue
        # local alignments stop before terminal mismatches: extend ungapped to the end of the overlap (<= 20 nt)
        f = frags[q]
        e5, e3 = min(qs - 1, ss - 1, 20), min(len(f) - qe, n - se, 20)
        qa = f[qs - 1 - e5:qs - 1] + qa + f[qe:qe + e3]
        sa = backbone[ss - 1 - e5:ss - 1] + sa + backbone[se:se + e3]
        best[q] = (ss - e5, se + e3, float(bs), qa.upper(), sa)
    base = [Counter() for _ in range(n + 1)]
    ins = [Counter() for _ in range(n + 1)]       # insertion after backbone position k
    span = [0.0] * (n + 2)
    ident: Dict[int, Tuple[float, int]] = {}
    for q, (ss, se, _, qa, sa) in best.items():
        w, k, buf, m, c = wt[q], ss - 1, "", 0, 0
        for x, y in zip(qa, sa):
            if y == "-":
                buf += x
                continue
            if buf:
                ins[k][buf] += w
            buf = ""
            k += 1
            base[k][x] += w
            if y.isupper():
                c += 1
                m += x == y
        for k in range(ss, se + 1):
            span[k] += w
        ident[q] = (m / c if c else 0.0, c)
    cons = []
    for k in range(1, n + 1):
        if span[k] == 0:
            cons.append(backbone[k - 1].lower())
            continue
        x = base[k].most_common(1)[0][0]
        cons.append("" if x == "-" else x)
        if k < n and ins[k]:
            s_, w_ = ins[k].most_common(1)[0]
            if w_ > min(span[k], span[k + 1]) / 2:
                cons.append(s_)
    return "".join(cons), ident


def _assemble(frags: List[str], wt: List[float], ref: Optional[str], work: Path) -> str:
    """backbone from fragments of one genome: contained fragments are removed and dovetailing ones merged
    (blastn, >= 98 % identity, >= 40 nt overlap, a copy variant within END_TOL nt of the ends allowed) until nothing
    changes; a piece end is joined to the best-covered piece that continues it, and only when that piece is best
    continued by it (where copy variants fork, the better-covered branch is joined); several pieces left are
    ordered on the reference gene; neighbours sharing an exact end overlap of >= 15 nt are joined on it, others
    with lowercase n for the reference distance between them (a piece placed over a heavier one is dropped).
    Without a reference only the heaviest piece is kept."""
    pieces = sorted(zip(frags, wt), key=lambda x: -x[1])
    changed = True
    while changed and len(pieces) > 1:
        changed = False
        with open(work / "pc.fa", "w") as o:
            for i, (f, _) in enumerate(pieces):
                o.write(f">{i}\n{f}\n")
        out = program_output([tool("blastn"), "-task", "blastn", "-query", str(work / "pc.fa"), "-subject", str(work / "pc.fa"),
                   "-strand", "plus", "-dust", "no", "-perc_identity", "98", "-evalue", "1e-10",
                   "-outfmt", "6 qseqid sseqid qstart qend sstart send length"])
        hsps = [h for h in (tuple(map(int, line.split("\t"))) for line in out.splitlines())
                if h[0] != h[1] and h[6] >= 40]

        def inside(h) -> bool:          # q inside t (a variant at q's ends allowed where t has sequence opposite it)
            q, t, qs, qe, ss, se = h[:6]
            lq, lt = len(pieces[q][0]), len(pieces[t][0])
            return lq <= lt and qs <= END_TOL + 1 and qe >= lq - END_TOL and ss >= qs and lt - se >= lq - qe

        def dovetail(h) -> bool:        # q's end overlaps t's start
            q, t, _, qe, ss, se = h[:6]
            return not inside(h) and qe >= len(pieces[q][0]) - END_TOL and ss <= END_TOL + 1 and \
                se < len(pieces[t][0]) - 2
        nxt: Dict[int, int] = {}        # best-covered (lowest index) piece continuing each end
        prv: Dict[int, int] = {}
        for h in filter(dovetail, hsps):
            nxt[h[0]], prv[h[1]] = min(nxt.get(h[0], h[1]), h[1]), min(prv.get(h[1], h[0]), h[0])
        # every containment / join between pieces not yet touched in this round is applied (one pair per
        # round needed a blastn per fragment: hours for a batch of hundreds of GenBank records)
        used, gone = set(), set()
        for h in hsps:
            q, t, qs, qe, ss, se, ln = h
            if q in used or t in used:
                continue
            if inside(h):
                pieces[t] = (pieces[t][0], pieces[t][1] + pieces[q][1])
                gone.add(q)
            elif dovetail(h) and nxt[q] == t and prv[t] == q:
                pieces[q] = (pieces[q][0][:qe] + pieces[t][0][se:], pieces[q][1] + pieces[t][1])
                gone.add(t)
            else:
                continue
            used.update((q, t))
            changed = True
        pieces = sorted((x for i, x in enumerate(pieces) if i not in gone), key=lambda x: -x[1])
    if len(pieces) == 1 or not ref:
        return pieces[0][0]
    with open(work / "pc.fa", "w") as o:
        for i, (f, _) in enumerate(pieces):
            o.write(f">{i}\n{f}\n")
    (work / "rf.fa").write_text(f">r\n{ref}\n")
    out = program_output([tool("blastn"), "-task", "blastn", "-query", str(work / "pc.fa"), "-subject", str(work / "rf.fa"),
               "-strand", "plus", "-dust", "no", "-evalue", "1e-5",
               "-outfmt", "6 qseqid qstart qend sstart send bitscore"])
    place: Dict[int, Tuple[int, int, float]] = {}
    for line in out.splitlines():
        q, qs, qe, ss, se, bs = line.split("\t")
        q, qs, qe, ss, se = int(q), int(qs), int(qe), int(ss), int(se)
        if q not in place or float(bs) > place[q][2]:
            place[q] = (ss - (qs - 1), se + (len(pieces[q][0]) - qe), float(bs))   # piece ends in reference positions
    kept: List[int] = []
    for q in sorted(place, key=lambda q: -pieces[q][1]):
        a, b, _ = place[q]
        if all(min(b, place[k][1]) - max(a, place[k][0]) < 50 for k in kept):
            kept.append(q)
    if not kept:                # pieces too short for blastn to place
        return pieces[0][0]
    kept.sort(key=lambda q: place[q][0])
    seq = pieces[kept[0]][0]
    for prev, q in zip(kept, kept[1:]):
        nxt = pieces[q][0]
        # assembler contigs often share an exact end overlap shorter than the blastn one (k-mer overlap)
        k = next((k for k in range(min(len(seq), len(nxt), 300), 14, -1) if seq[-k:] == nxt[:k]), 0)
        seq += nxt[k:] if k else "n" * max(1, place[q][0] - place[prev][1] - 1) + nxt
    return seq


def consensus(rows: List[dict], recs, lib: "Library", work: Path, min_ident: float,
              bounded: bool = False) -> List[dict]:
    """One sequence per gene from all contigs. rDNA genes / spacers: fragments are cut at scaffold gaps (runs of
    >= 10 N); the pieces of the fragments with >= SEED_COV of the best contig coverage (contig cov_ of SPAdes names,
    else 1 = all fragments) are joined into a backbone by their overlaps and ordered on the reference gene of the
    best fragment (_assemble); a spacer has no reference, so its complete copies among them make the backbone when
    there are any, ordered on the best copy. Every fragment is placed on it and kept when it has >= min_ident over
    >= 50 nt of backbone sequence; the others (other organisms, another rDNA type, contamination) are listed as
    dropped. Each position takes the coverage-weighted majority of the kept fragments (_pileup); stretches no kept
    fragment covers are N (ends stripped). Histones: the CDS of the best-covered contig, complete and >= 80 % of the
    median reference length first (copies are paralogs; a consensus would be a chimera)."""
    out = []
    for g in UNIT:
        frags = [r for r in rows if r["gene"] == g]
        if not frags:
            continue
        seqs = [str(_gene_seq(r, recs)).upper() for r in frags]
        wt = [float(r["contig_cov"] or 1) for r in frags]
        top = max(range(len(frags)), key=lambda i: (wt[i], int(frags[i]["length"])))
        ref = _ref_gene(lib, frags[top]["best_ref"], g) if g in RRNA else None
        if len(frags) == 1:          # one copy (a GenBank record, a PCR product): the copy itself, no alignment
            ident, keep, cons = {0: (1.0, len(seqs[0]))}, [0], seqs[0]
        else:
            # scaffold gaps are not sequence: SPAdes puts N between contigs that overlap by k as well as over
            # stretches a short contig of its own holds; the pieces between the N join like contigs
            parts = [(i, p) for i, s in enumerate(seqs) for p in re.split("N{10,}", s) if p]
            # seed: the well-covered fragments (same genome as the best one in an assembly of mixed samples)
            seed = {i for i in range(len(frags)) if wt[i] >= SEED_COV * wt[top]}
            # assembly: the best-covered ones (the all-against-all blastn of _assemble grows with their number squared)
            if bounded:
                seed = set(sorted(seed, key=lambda i: (wt[i], len(seqs[i]), -i), reverse=True)[:MAX_CONS_SEED])
            if not ref:
                seed = {i for i in seed if frags[i]["partial_5"] == frags[i]["partial_3"] == "0"} or seed
                # the best copy orders the pieces it was cut into at scaffold gaps
                ref = seqs[max(seed, key=lambda i: (wt[i], len(seqs[i])))]
            backbone = _assemble([p for i, p in parts if i in seed], [wt[i] for i, _ in parts if i in seed], ref, work)
            _, idt = _pileup([p for _, p in parts], [wt[i] for i, _ in parts], backbone, work)
            m, c = defaultdict(float), defaultdict(int)       # identity of each fragment over all its pieces
            for k, (x, n) in idt.items():
                m[parts[k][0]] += x * n
                c[parts[k][0]] += n
            ident = {i: (m[i] / c[i] if c[i] else 0.0, c[i]) for i in range(len(frags))}
            # the best-covered fragment only when nothing else qualifies: it may be the other rDNA type of a
            # mixed sample, and its votes made a chimera of the two
            keep = [i for i in range(len(frags)) if ident[i][0] >= min_ident and ident[i][1] >= 50] or [top]
            kept = [(i, p) for i, p in parts if i in keep]
            cons, _ = _pileup([p for _, p in kept], [wt[i] for i, _ in kept], backbone, work)
        seq = re.sub("[a-z]", "N", cons).strip("N")
        out.append({"gene": g, "length": len(seq), "n_fragments": len(keep), "n_dropped": len(frags) - len(keep),
                    "internal_N": seq.count("N"),
                    "partial_5": int(all(frags[i]["partial_5"] == "1" for i in keep)),
                    "partial_3": int(all(frags[i]["partial_3"] == "1" for i in keep)),
                    "contigs": ";".join(frags[i]["contig"] for i in keep),
                    "dropped": ";".join(f"{frags[i]['contig']}({ident[i][0] * 100:.1f}% over {ident[i][1]} nt)"
                                        for i in range(len(frags))
                                        if i not in keep),
                    "seq": seq})
    for h in HISTONES:
        cds = [r for r in rows if r["gene"] == h]
        if not cds:
            continue
        lens = sorted(int(r["length_aa"]) for r in lib.hist.values() if r["gene"] == h)
        r = pick_histone(cds, 3 * lens[len(lens) // 2] if lens else 0)
        n_rep = sum(c.get("class") == "replacement" for c in cds)
        out.append({"gene": h, "length": int(r["length"]), "n_fragments": 1, "n_dropped": 0, "internal_N": 0,
                    "partial_5": int(r["partial_5"]), "partial_3": int(r["partial_3"]), "contigs": r["contig"],
                    "class": r.get("class", ""),
                    "dropped": f"{len(cds) - 1} other copies" + (f" ({n_rep} replacement)" if n_rep else ""),
                    "codon_start": int(r.get("codon_start") or 1), "seq": str(_gene_seq(r, recs))})
    return out


def pick_histone(cds: List[dict], typical: int) -> dict:
    """The copy that stands for a histone gene: replication-dependent (cluster) copies first, then undecided, then
    single-copy, replacement variants (H3.3, spliced genes) last - they are other paralogs; within a class complete
    and >= 80 % of the typical length, then coverage, then length. A cluster copy that SPAdes left as a fragment
    is reported partial rather than replaced by a complete variant."""
    return max(cds, key=lambda r: (CLASS_RANK.get(r.get("class", ""), 2),
                                   r["partial_5"] == r["partial_3"] == "0" and int(r["length"]) >= 0.8 * typical,
                                   float(r["contig_cov"] or 0), int(r["length"])))


def write_consensus(out: Path, cons: List[dict]) -> None:
    d = out / "consensus"
    d.mkdir(exist_ok=True)
    multi = any(c.get("sample") for c in cons)
    cols = (["sample"] if multi else []) + ["gene", "length", "n_fragments", "n_dropped", "internal_N", "partial_5",
                                            "partial_3", "class", "contigs", "dropped"]
    with open(out / "consensus.tsv", "w", newline="") as o:
        w = csv.DictWriter(o, cols, delimiter="\t", extrasaction="ignore", restval="")
        w.writeheader()
        w.writerows(cons)
    fa, faa = defaultdict(str), defaultdict(str)     # one file per gene, one sequence per sample
    for c in cons:
        part = "".join(x for x, y in (("5'-partial ", c["partial_5"]), ("3'-partial ", c["partial_3"])) if y)
        name = re.sub(r"[^A-Za-z0-9_.-]", "_", f"{c['sample']}_{c['gene']}") if multi else c["gene"]
        fa[c["gene"]] += f">{name} {part}fragments={c['n_fragments']} length={c['length']}\n{c['seq']}\n"
        if c["gene"] in HISTONES:
            nt = c["seq"][int(c.get("codon_start") or 1) - 1:]          # a partial 5' end: the frame is codon_start
            faa[c["gene"]] += f">{name}\n{str(Seq(nt[:len(nt) // 3 * 3]).translate(table=1)).rstrip('*')}\n"
    for g, txt in fa.items():
        (d / f"{g}.fasta").write_text(txt)
    for g, txt in faa.items():
        (d / f"{g}.faa").write_text(txt)


def annotate(inp: Path, out: Path, lib_dir: Path, genes: List[str], threads: int,
             leave_out: Optional[Path] = None, consensus_identity: float = 0.99,
             organism: Optional[str] = None, src_table: Optional[Path] = None, refine: bool = True) -> List[dict]:
    out.mkdir(parents=True, exist_ok=True)
    work = out / "tmp"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir()
    lo = {}
    if leave_out:
        for line in open(leave_out):
            if line.strip():
                k, v = line.rstrip("\n").split("\t")[:2]
                lo[k] = v
    lib = Library(lib_dir, work)
    t0 = time.time()
    fa, recs, org, n_seq, n_nt, circ = open_input(inp, work)
    try:
        if circ:
            rotate_circular(fa, circ, lib, genes, work, threads, lo)
        recs = recs or FastaIndex(fa)
        logging.info(f"[NUC] {n_seq} sequences, {n_nt:,} nt ({time.time() - t0:.0f} s)")
        # the prefilter (one reference per genus) only pays off on large assemblies; smaller inputs are searched whole
        t0 = time.time()
        if n_nt > PREFILTER_NT:
            windows, bounded = prefilter(fa, lib, genes, work, threads, n_seq)
        else:
            whole = {r.id: [(1, len(r))] for r in SeqIO.parse(str(fa), "fasta")}
            windows, bounded = {g: whole for g in genes}, set()
        regions, seqs = {}, {}
        for g in genes:
            regions[g] = work / f"regions_{g}.fa"
            seqs[g] = write_regions(windows[g], recs, regions[g])
            logging.info(f"[NUC] prefilter {g}: {len(windows[g])} contig(s), {len(seqs[g])} region(s), "
                         f"{sum(len(s) for s in seqs[g].values()):,} nt ({time.time() - t0:.0f} s)")
        rdna = hist = None
        if seqs.get("rdna"):
            rdna = annotate_rdna(regions["rdna"], lib, work, threads, lo, {k: len(v) for k, v in seqs["rdna"].items()},
                                 refine, "rdna" in bounded)
        if seqs.get("histone"):
            t0 = time.time()
            hist = annotate_histones(regions["histone"], lib, work, threads, lo, seqs["histone"])
            logging.info(f"[NUC] histones (tblastn + exonerate): {time.time() - t0:.0f} s")
        base = single_copy_cov(recs.keys()) if hist else None
        if base:
            logging.info(f"[NUC] single-copy coverage {base:.1f}x (length-weighted median of SPAdes contigs >= 1 kb)")
        rows, feats = collect(recs, rdna, hist, base)
        write_outputs(out, recs, org, rows, feats)
        if feats:
            write_submission(out, inp, recs, org, feats, organism, src_table)
        if rows and consensus_identity > 0:
            rows_s = [{k: str(v) for k, v in r.items()} for r in rows]
            # GenBank records of several organisms are several samples: a consensus over them would join conserved
            # 18S of different species (and an all-against-all blastn of hundreds of fragments took hours)
            groups: Dict[str, List[dict]] = defaultdict(list)
            for r in rows_s:
                groups[org.get(r["contig"], "")].append(r)
            if len(groups) == 1:
                groups = {"": rows_s}
            items = sorted(groups.items())
            if len(items) == 1:
                res = [consensus(items[0][1], recs, lib, work, consensus_identity, "rdna" in bounded)]
            else:
                seqs_by_id = {r.id: r for r in (recs[c] for c in {r["contig"] for r in rows_s})}
                lib._seqs = SeqIO.to_dict(SeqIO.parse(str(lib.rdna_fa), "fasta"))     # SeqIO.index is not thread-safe

                def one(i: int) -> List[dict]:
                    wd = work / f"cons{i}"      # each sample its own blastn scratch files
                    wd.mkdir(exist_ok=True)
                    return consensus(items[i][1], seqs_by_id, lib, wd, consensus_identity)
                with ThreadPoolExecutor(max(1, threads)) as pool:
                    res = list(pool.map(one, range(len(items))))
            cons = [dict(c, sample=sample) for (sample, _), cs in zip(items, res) for c in cs]
            write_consensus(out, cons)
            if len(groups) > 1:
                logging.info(f"[NUC] consensus: {len(cons)} sequences of {len(groups)} organisms  ({out / 'consensus.tsv'})")
            else:
                logging.info("[NUC] consensus: " + ", ".join(f"{c['gene']} {c['length']} nt from {c['n_fragments']}"
                                                             for c in cons) + f"  ({out / 'consensus.tsv'})")
    finally:
        if recs is not None:
            recs.close()
    shutil.rmtree(work)
    counts = Counter(r["gene"] for r in rows)
    found = ", ".join(f"{g} {counts[g]}" for g in (*UNIT, *HISTONES) if counts[g])
    logging.info(f"[NUC] features: {found or 'none'}")
    logging.info(f"[NUC] {out / 'nuc_features.tsv'}  {out / 'nuc_annotation.gb'}  {out / 'genes'}/")
    return rows


def build_parser() -> argparse.ArgumentParser:
    """The command line (main); options() gives the same options without one."""
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fetch", help="download the NCBI reference records (rdna.gb, histone.gb)")
    f.add_argument("-o", "--out", type=Path, default=Path("nuc_refs"))
    f.add_argument("--email")
    f.add_argument("--api_key")
    b = sub.add_parser("build", help="build the reference library from annotated GenBank files")
    b.add_argument("--gb", nargs="+", type=Path, required=True)
    b.add_argument("-o", "--out", type=Path, default=LIB)
    b.add_argument("--threads", type=int, default=8)
    b.add_argument("--curated", nargs="*", type=Path, default=[],
                   help="GenBank files of trusted references (curated.gb of mitom nuc fetch, own checked records)")
    a = sub.add_parser("annotate", help="annotate rDNA and histone genes in contigs (FASTA or GenBank)")
    a.add_argument("-i", "--input", type=Path, required=True)
    a.add_argument("-o", "--out", type=Path, required=True)
    a.add_argument("--lib", type=Path, default=LIB, help="reference library (mitom nuc build)")
    a.add_argument("--genes", default="rdna,histone", help="rdna, histone or both (default)")
    a.add_argument("--threads", type=int, default=8)
    a.add_argument("--consensus_identity", type=float, default=0.99,
                   help="fragments below this identity to a gene's consensus are left out of it (default 0.99; "
                        "0 = no consensus)")
    a.add_argument("--no_cm", action="store_true",
                   help="rRNA ends from the reference projection only (default: placed with the Rfam models)")
    a.add_argument("--organism", help="organism of FASTA input, for the NCBI submission files (GenBank input: /organism)")
    a.add_argument("--src", type=Path, help="NCBI source table (TSV: Sequence_ID, organism, isolate, specimen_voucher, "
                                            "country, ...; title = own definition line) for the submission files")
    a.add_argument("--leave_out", type=Path, help=argparse.SUPPRESS)   # benchmark: TSV sequence id -> genus to skip
    return ap


def options(cmd: str, **given) -> argparse.Namespace:
    """Options for run(): cmd "fetch", "build" (gb) or "annotate" (input, out); paths as Path."""
    return mito_common.options(build_parser(), cmd, **given)


def run(args: argparse.Namespace) -> int:
    """args.cmd fetch / build / annotate (the functions fetch(), build(), annotate() do the work). Returns 0."""
    handlers = [logging.StreamHandler(sys.stdout)]
    if args.cmd == "annotate":
        args.out.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(args.out / "nuc.log", mode="w"))
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", handlers=handlers)
    if args.cmd == "fetch":
        fetch(args.out, args.email, args.api_key)
    elif args.cmd == "build":
        n_r, n_h = build(args.gb, args.out, args.threads, args.curated)
        logging.info(f"[NUC] library {args.out}: {n_r} rDNA records, {n_h} histone proteins")
    else:
        genes = [g.strip() for g in args.genes.split(",") if g.strip()]
        bad = set(genes) - {"rdna", "histone"}
        if bad:
            raise InputError(f"--genes: unknown {', '.join(sorted(bad))} (rdna, histone)")
        src = Path(prepare_input(args.input))
        require_records(count_seq_records(src), "sequence records", args.input)
        annotate(src, args.out, args.lib, genes, args.threads, args.leave_out, args.consensus_identity,
                 args.organism, args.src, not args.no_cm)
    return 0


def main(argv=None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    sys.exit(cli_guard(main))
