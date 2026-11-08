"""mito_nuc_histone.py - the histone part of mito_nuc.py: H1-H4 genes by tblastn + exonerate on the reference proteins,
replication-dependent / variant class, single-copy coverage (annotate_histones)."""
from __future__ import annotations
import logging
import re
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from Bio.Seq import Seq
from mito_nuc_lib import (  # noqa: E402
    _contig,
    INTRON_PEN,
    Library,
    MIN_H34_PARTIAL_ID,
    MIN_PROT_COVER,
    program_output,
    STOPS,
    tool,
    TOPK,
)
from mito_common import run_tool


# ── histones ────────────────────────────────────────────────────────────────
# Replication-dependent histones (the tandem cluster; the H3 / H4 markers of phylogenetics) vs replacement
# variants. Canonical H3 (H3.1 / H3.2) and H3.3 differ at four residues in animals (numbered without the initiator
# Met): A31 S87 V89 M90 vs S31 A87 I89 G90; replication-dependent genes have no introns; H2A.X ends in SQ(E/D)Phi.
# Found on MJ20260820041 (2026-10-28): SPAdes splits the histone cluster into 84-134 nt contigs at ~4000x, so the
# "complete, best covered" H3 of the consensus was a single-copy H3.3 (14x) in 3 of 6 annelid assemblies and
# canonical H3 in 2 - different paralogs under one gene name.
H3_SITES = (31, 87, 89, 90)


H3_SIGNATURE = {"canonical H3": "ASVM", "H3.3": "SAIG"}


H2AX_END = re.compile(r"SQ[ED][ACFGILMVWY]$")


CLUSTER_RATIO = 5.0     # contig coverage / single-copy coverage at or above which a copy counts as a cluster copy


SINGLE_RATIO = 2.0      # below: a single-copy gene (MJ20260820041: H2A.X, extra H2B / H3 at 0.8-1.5x; cluster 100-270x)


CLASS_RANK = {"replication-dependent": 3, "": 2, "single-copy": 1, "replacement": 0}


def histone_class(gene: str, prot: str, n_exons: int, ratio: Optional[float]) -> Tuple[str, str]:
    """('replication-dependent' | 'single-copy' | 'replacement' | '', evidence) for one histone CDS (protein with or
    without Met). Sequence evidence (H3 residues, introns) outranks coverage."""
    why, cls = [], ""
    if gene == "H3":
        i = prot.find("ARTKQTAR")
        if i >= 0 and len(prot) - i >= max(H3_SITES):
            r = "".join(prot[i + p - 1] for p in H3_SITES)
            score = {k: sum(a == b for a, b in zip(r, v)) for k, v in H3_SIGNATURE.items()}
            (best, hi), (_, lo) = sorted(score.items(), key=lambda x: -x[1])
            if hi >= 3 and hi > lo:
                why.append(f"H3 residues 31/87/89/90 {r} ({best})")
                cls = "replication-dependent" if best == "canonical H3" else "replacement"
    if gene == "H2A" and H2AX_END.search(prot):
        why.append(f"H2A.X C-terminus {prot[-4:]}")
    if n_exons > 1:
        why.append(f"{n_exons - 1} intron(s)")
        cls = "replacement"
    if ratio is not None:
        why.append(f"copy ratio {ratio:.1f}")
        if not cls and ratio >= CLUSTER_RATIO:
            cls = "replication-dependent"
        elif not cls and ratio < SINGLE_RATIO:
            cls = "single-copy"
    return cls, "; ".join(why)


def single_copy_cov(ids) -> Optional[float]:
    """Coverage of the single-copy genome from SPAdes contig names (NODE_1_length_5000_cov_12.3): the
    length-weighted median coverage of contigs >= 1 kb; None for other names."""
    pairs = []
    for c in ids:
        m = re.search(r"length_(\d+)_cov_([\d.]+)", c)
        if m and int(m.group(1)) >= 1000:
            pairs.append((float(m.group(2)), int(m.group(1))))
    if not pairs:
        return None
    pairs.sort()
    half, acc = sum(n for _, n in pairs) / 2, 0
    for cov, n in pairs:
        acc += n
        if acc >= half:
            return cov
    return None


def _exonerate_all(target: str, prots: List[str], model: str, td: Path) -> List[Optional[tuple]]:
    """the best alignment of each protein (--bestn 1 counts per query) in one exonerate run: (score, exons, protein
    start, end, length, % identity) or None. One run per protein and model was 22,700 runs of ~20 ms for the 3,700 H3 loci of
    SRR2006110; 2 per locus now"""
    (td / "p.fa").write_text("".join(f">p{i}\n{p}\n" for i, p in enumerate(prots)))
    p = run_tool([tool("exonerate"), "--model", model, "--showalignment", "no", "--showvulgar", "no",
                  "--showtargetgff", "yes", "--bestn", "1", "--geneticcode", "1", "-q", str(td / "p.fa"),
                  "-t", target, "--ryo", "SCORE %qi %s %qab %qae %ql %pi\n"])
    if p.returncode:                             # a crash is not 'no alignment' (GLM 4.2, 2026-11-05)
        logging.warning(f"exonerate failed (exit {p.returncode}) on {target}: {' '.join(p.stderr.split())[-200:]}")
    out = p.stdout
    res: List[list] = [[None, [], 0, 0, 0, 0.0] for _ in prots]
    ex: List[tuple] = []
    for line in out.splitlines():
        if line.startswith("SCORE"):             # follows the GFF of its alignment
            _, q, sc, qab, qae, ql, pi = line.split()
            r = res[int(q[1:])]
            r[0], r[2], r[3], r[4], r[5] = int(sc), int(qab), int(qae), int(ql), float(pi)
            r[1] += ex
            ex = []
            continue
        x = line.split("\t")
        # protein2genome lists exons; protein2dna (no introns) only one similarity line
        if len(x) > 8 and (x[2] == "exon" or x[2] == "similarity" and model == "protein2dna"):
            ex.append((int(x[3]), int(x[4]), x[6]))
    return [(sc, e, qab, qae, ql, pi) if sc and e else None for sc, e, qab, qae, ql, pi in res]


def _exonerate(target: str, prot: str, model: str, td: Path) -> Optional[tuple]:
    return _exonerate_all(target, [prot], model, td)[0]


def _canonical_introns(tseq: str, ex: List[tuple]) -> bool:
    """every intron of an exonerate model (forward coordinates of tseq) has GT/GC ... AG splice sites"""
    ex = sorted(ex)
    for (a1, b1, st), (a2, b2, _) in zip(ex, ex[1:]):
        left, right = tseq[b1:b1 + 2].upper(), tseq[a2 - 3:a2 - 1].upper()
        if (st == "+" and (left not in ("GT", "GC") or right != "AG")) or \
                (st == "-" and (left != "CT" or right not in ("AC", "GC"))):
            return False
    return True


def _orf(seq: Seq, ex: List[Tuple[int, int]], st: int,
         at_start: bool) -> Tuple[List[Tuple[int, int]], bool, bool, int]:
    """extend the exon model (forward coordinates) to a stop codon and an ATG; returns exons, 5' and 3' partial
    (of the gene, not of the forward coordinates) and codon_start. An end left open runs to the sequence end
    (NCBI: a partial CDS end must be at the sequence end; the first codon then starts at codon_start 1-3)"""
    n = len(seq)
    cs = seq if st == 1 else seq.reverse_complement()
    tr = (lambda a, b: (a, b)) if st == 1 else (lambda a, b: (n - b + 1, n - a + 1))
    E = sorted(tr(a, b) for a, b in ex)
    p, p3 = E[-1][1], True
    for _ in range(300):
        c = str(cs[p:p + 3]).upper()
        if len(c) < 3:
            break
        p += 3
        if c in STOPS:
            p3 = False
            break
    E[-1] = (E[-1][0], n if p3 and n - p < 3 else min(p, n))
    a0, p5, codon_start = E[0][0], False, 1
    if not (at_start and str(cs[a0 - 1:a0 + 2]).upper() == "ATG"):
        q, cand = a0, []
        while q - 3 >= 1:
            c = str(cs[q - 4:q - 1]).upper()
            if c in STOPS:
                break
            q -= 3
            if c == "ATG":
                cand.append(q)
        if cand:
            E[0] = (cand[-1], E[0][1])
        else:
            p5 = str(cs[a0 - 1:a0 + 2]).upper() != "ATG"
            if p5 and q - 3 < 1:
                E[0], codon_start = (1, E[0][1]), q   # open to the contig end
    return sorted(tr(a, b) for a, b in E), p5, p3, codon_start


def _tblastn(lib: Library, db: str, work: Path, threads: int) -> List[str]:
    """tblastn of the reference proteins against db, lines as one tblastn run writes them. 677 proteins are 269
    different sequences (one H3 sequence 87 times): each sequence is searched once (its hits are the same for every
    copy) in single-threaded runs on `threads` pieces at once; one tblastn of all on 16 threads took 40-50 s on the
    3,700 histone regions of SRR2006110 (460 Mb assembly)"""
    first: Dict[str, str] = {}
    for ref, p in lib.prot.items():
        first.setdefault(p, ref)
    uniq = list(first.items())
    n = max(1, min(threads, len(uniq)))
    for i in range(n):
        with open(work / f"hist_q{i}.faa", "w") as o:
            o.writelines(f">{ref}\n{p}\n" for p, ref in uniq[i::n])

    def search(i: int) -> str:
        return program_output([tool("tblastn"), "-query", str(work / f"hist_q{i}.faa"), "-db", db, "-evalue", "1e-5", "-outfmt",
                    "6 qseqid sseqid sstart send bitscore", "-max_target_seqs", "100000", "-num_threads", "1",
                    "-seg", "no"])
    by_query: Dict[str, List[str]] = defaultdict(list)
    with ThreadPoolExecutor(n) as pool:
        for out in pool.map(search, range(n)):
            for line in out.splitlines():
                by_query[line.split("\t", 1)[0]].append(line.split("\t", 1)[1])
    return [f"{ref}\t{x}" for ref, p in lib.prot.items() for x in by_query[first[p]]]


def annotate_histones(regions: Path, lib: Library, work: Path, threads: int, leave_out: Dict[str, str],
                      seqs: Dict[str, Seq]) -> Dict[str, List[dict]]:
    db = str(work / "regions")
    program_output([tool("makeblastdb"), "-in", str(regions), "-dbtype", "nucl", "-out", db])
    hits = defaultdict(list)
    for line in _tblastn(lib, db, work, threads):
        ref, s, a, b, bs = line.split("\t")
        r = lib.hist[ref]
        if leave_out and leave_out.get(_contig(s)) == r["genus"] or _contig(s) == r["accession"]:
            continue
        a, b = int(a), int(b)
        hits[s].append({"ref": ref, "lo": min(a, b), "hi": max(a, b), "st": 1 if b >= a else -1, "bs": float(bs)})

    def job(s: str) -> List[dict]:
        seq, res = seqs[s], []
        loci: List[dict] = []
        for h in sorted(hits[s], key=lambda h: (h["st"], h["lo"])):
            if loci and loci[-1]["st"] == h["st"] and h["lo"] <= loci[-1]["hi"]:
                loci[-1]["hi"] = max(loci[-1]["hi"], h["hi"])
                loci[-1]["h"].append(h)
            else:
                loci.append({"st": h["st"], "lo": h["lo"], "hi": h["hi"], "h": [h]})
        with tempfile.TemporaryDirectory(dir=work) as td:
            td = Path(td)
            for L in loci:
                a, b = max(0, L["lo"] - 1 - 3000), min(len(seq), L["hi"] + 3000)
                (td / "t.fa").write_text(f">t\n{seq[a:b]}\n")
                refs = []
                for h in sorted(L["h"], key=lambda h: -h["bs"]):
                    if h["ref"] not in refs:
                        refs.append(h["ref"])
                # each different protein once (the best references of a locus often share theirs), all in one run
                prots = list(dict.fromkeys(lib.prot[ref] for ref in refs[:TOPK]))
                done = {(p, model): m for model in ("protein2genome", "protein2dna")
                        for p, m in zip(prots, _exonerate_all(str(td / "t.fa"), prots, model, td))}
                models = []
                for ref in refs[:TOPK]:
                    for model in ("protein2genome", "protein2dna"):
                        m = done[(lib.prot[ref], model)]
                        # a spliced model needs real splice sites (H1 genes have no introns: a gap in the
                        # alignment there is a frameshift or a divergent tail, OU230952.1)
                        if m and not _canonical_introns(str(seq[a:b]), m[1]):
                            m = None
                        if m:
                            sc = m[0] * (INTRON_PEN ** (len(m[1]) - 1) if model == "protein2genome" else 1)
                            models.append((sc, ref, m))
                if not models:
                    continue
                sc, ref, (_, ex, qab, qae, ql, pi) = max(models, key=lambda x: x[0])
                st = 1 if ex[0][2] == "+" else -1
                ex = [(x + a, y + a) for x, y, _ in ex]
                exons, p5, p3, cstart = _orf(seq, ex, st, qab == 0)
                at_edge = exons[0][0] <= 3 or exons[-1][1] >= len(seq) - 2
                if (qae - qab) / ql < MIN_PROT_COVER and not at_edge:
                    continue
                # H3 and H4 are >= ~85 % identical across animals: a partial model far below that is a piece of a
                # degenerate copy (SRR2006110: 3,688 'H3' of 61 % median identity, retrocopies with frameshifts and
                # stop codons, ~87 % nt identity to each other); complete divergent genes exist (MT043459.1 H3 50 %)
                if lib.hist[ref]["gene"] in ("H3", "H4") and pi < MIN_H34_PARTIAL_ID and (p5 or p3):
                    continue
                res.append({"gene": lib.hist[ref]["gene"], "st": st, "exons": exons, "p5": p5, "p3": p3,
                            "codon_start": cstart, "ref": lib.hist[ref]["accession"], "ref_org": lib.hist[ref]["organism"], "score": sc})
        # same CDS found from two overlapping loci: keep the better one
        res.sort(key=lambda r: -r["score"])
        kept = []
        for r in res:
            if not any(k["st"] == r["st"] and k["exons"][0][0] <= r["exons"][-1][1] and
                       r["exons"][0][0] <= k["exons"][-1][1] for k in kept):
                kept.append(r)
        return kept

    with ThreadPoolExecutor(max(1, threads)) as pool:
        return dict(zip(list(hits), pool.map(job, list(hits))))
