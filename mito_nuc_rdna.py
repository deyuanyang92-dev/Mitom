"""mito_nuc_rdna.py - the rDNA part of mito_nuc.py: 18S-ITS1-5.8S-ITS2-28S loci from blastn chains of the reference units,
boundaries by reference votes and the Rfam covariance models, a missing 5.8S searched again (annotate_rdna)."""
from __future__ import annotations
import logging
import re
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from mito_nuc_lib import (  # noqa: E402
    _contig,
    _ref_gene,
    _strand,
    BLAST_TARGETS,
    chains,
    CHUNK_NT,
    CURATED_MIN_ID,
    EXTRAP,
    LIB,
    Library,
    map_pos,
    MAX_CM_JOBS,
    MAX_ITS,
    MAX_RESCUE,
    MIN_CHAIN_BITS,
    MIN_GENE_OVERLAP,
    PIECE_MIN_ID,
    program_output,
    RRNA,
    tool,
    tool_ok,
    TOPN,
)
from mito_common import run_tool


def project(h, b: int, maxd: Optional[int]) -> Optional[int]:
    """reference position b outside the HSP -> query position, ungapped from the nearer HSP end"""
    qs, qe, ss, se = h[:4]
    st = _strand(h)
    d1, d2 = (ss - b) * st, (b - se) * st
    if d1 > 0 and (maxd is None or d1 <= maxd):
        return qs - d1
    if d2 > 0 and (maxd is None or d2 <= maxd):
        return qe + d2
    return None


def rdna_loci(hits: Dict[str, list], lib: Library, qlen: int, min_overlap: int = MIN_GENE_OVERLAP,
              min_bits: float = MIN_CHAIN_BITS) -> List[dict]:
    """one locus per rDNA copy: chains of all references grouped by query overlap; boundaries by weighted vote"""
    allc = []
    for ref, hs in hits.items():
        for c in chains(hs, min_bits):
            allc.append((sum(h[4] for h in c), ref, c))
    allc.sort(key=lambda x: -x[0])
    loci: List[dict] = []
    for sc, ref, c in allc:
        lo, hi, st = c[0][0], c[-1][1], _strand(c[0])
        for L in loci:
            ov = min(hi, L["hi"]) - max(lo, L["lo"])
            if L["st"] == st and ov > 0.5 * min(hi - lo, L["hi"] - L["lo"]):
                if ref not in {r for _, r, _ in L["chains"]}:
                    L["chains"].append((sc, ref, c))
                break
        else:
            loci.append({"lo": lo, "hi": hi, "st": st, "chains": [(sc, ref, c)]})
    out = []
    for L in loci:
        top = L["chains"][:TOPN]
        st = L["st"]
        feats = {}
        for g in RRNA:
            votes = {"start": Counter(), "end": Counter()}
            proj = {"start": [], "end": []}
            covered = 0
            for sc, ref, c in top:
                if g in lib.rdna[ref]:
                    a, b = lib.rdna[ref][g][:2]
                    cov = set()
                    for h in c:
                        cov.update(range(max(a, min(h[2], h[3])), min(b, max(h[2], h[3])) + 1))
                    covered = max(covered, len(cov))
            trusted, trusted_id = {}, {}
            for side, k in (("start", 0), ("end", 1)):
                usable = [x for x in L["chains"] if g in lib.rdna[x[1]] and lib.rdna[x[1]][g][2 + k]]
                # a curated reference aligned over the boundary (or <= EXTRAP nt past its alignment) decides it
                for sc, ref, c in (x for x in usable if x[1] in lib.curated and _identity(x[2]) >= CURATED_MIN_ID):
                    pos = lib.rdna[ref][g][k]
                    m = next((x for x in (map_pos(h, pos) for h in c) if x is not None), None)
                    if m is None:
                        m = next((x for x in (project(h, pos, EXTRAP) for h in c) if x is not None), None)
                    if m is not None:
                        trusted[side], trusted_id[side] = int(m), _identity(c)
                        break
                # a side no top reference may vote on (partial or failed QC) takes the next references that may
                for sc, ref, c in [x for x in usable if x in top] or usable[:TOPN]:
                    pos = lib.rdna[ref][g][k]
                    m = next((x for x in (map_pos(h, pos) for h in c) if x is not None), None)
                    if m is None:
                        m = next((x for x in (project(h, pos, EXTRAP) for h in c) if x is not None), None)
                    if m is not None:
                        votes[side][m] += sc
                    else:
                        near = min(c, key=lambda h: min(abs(h[2] - pos), abs(h[3] - pos)))
                        proj[side].append(project(near, pos, None))
            ends = {}
            for side in ("start", "end"):
                if side in trusted:
                    ends[side] = (trusted[side], "curated")
                elif votes[side]:
                    ends[side] = (int(votes[side].most_common(1)[0][0]), "aligned")
                elif [p for p in proj[side] if p is not None]:
                    p = sorted(x for x in proj[side] if x is not None)
                    ends[side] = (p[len(p) // 2], "projected")
                else:
                    ends[side] = None
            if ends["start"] is None or ends["end"] is None:
                continue
            g5, g3 = ends["start"], ends["end"]           # gene 5' / 3' in query coordinates
            lo_, hi_ = (g5, g3) if st == 1 else (g3, g5)
            p_lo, p_hi = lo_[0] < 1, hi_[0] > qlen
            lo_p, hi_p = max(1, lo_[0]), min(qlen, hi_[0])
            if hi_p < lo_p or (hi_p - lo_p < 20 and not (p_lo or p_hi)):      # short pieces only at a contig end
                continue
            # aligned query nt inside the feature itself (an alignment of a spacer or of another gene does not
            # carry a gene projected kilobases away) and the bits of that part (picks between overlapping loci)
            qcov, sup, aln = 0, 0.0, set()
            for sc, ref, c in top:
                if g not in lib.rdna[ref]:
                    continue
                ga, gb = lib.rdna[ref][g][:2]
                pos, bits = set(), 0.0
                for h in c:
                    # the part of the HSP inside the reference gene, on the query (ungapped offset)
                    r1, r2 = max(ga, min(h[2], h[3])), min(gb, max(h[2], h[3]))
                    if r2 < r1:
                        continue
                    q1, q2 = sorted(((h[0] + (r - h[2])) if _strand(h) == 1 else (h[0] + (h[2] - r))) for r in (r1, r2))
                    a_, b_ = max(q1, lo_p, h[0]), min(q2, hi_p, h[1])
                    if b_ >= a_:
                        pos.update(range(a_, b_ + 1))
                        bits += h[4] * (b_ - a_ + 1) / (h[1] - h[0] + 1)
                aln |= pos
                qcov, sup = max(qcov, len(pos)), max(sup, bits)
            feats[g] = {"lo": lo_p, "hi": hi_p, "p_lo": p_lo, "p_hi": p_hi, "sup": sup, "aligned": qcov, "aln": aln,
                        "cur_id": {5: trusted_id.get("start", 0), 3: trusted_id.get("end", 0)},
                        "weak": covered < min_overlap or qcov < min_overlap,
                        "how": f"{lo_[1]}/{hi_[1]}" if st == 1 else f"{hi_[1]}/{lo_[1]}"}
        # a gene needs MIN_GENE_OVERLAP aligned nt, except a piece cut by the contig end next to a supported gene
        if all(f["weak"] for f in feats.values()):
            continue
        feats = {g: f for g, f in feats.items() if not f["weak"] or f["p_lo"] or f["p_hi"]}
        best = top[0][1]
        ident = _identity(top[0][2])
        out.append({"st": st, "feats": feats, "ref": best, "ref_org": lib.org[best][0], "identity": ident})
    _dedup(out)
    for L in out:
        _spacers(L["feats"], L["st"], qlen)
    return [L for L in out if L["feats"]]


def _dedup(loci: List[dict]) -> None:
    """one rRNA per place: a chain spanning parts of two units (28S - spacer - 18S) forms its own locus beside the
    unit's locus, and both report the shared gene; keep the copy with more aligned bits inside it"""
    fs = sorted(((f["sup"], i, g) for i, L in enumerate(loci) for g, f in L["feats"].items() if g in RRNA),
                reverse=True)
    kept: List[tuple] = []
    for _, i, g in fs:
        f, st = loci[i]["feats"][g], loci[i]["st"]
        if any(g == g2 and st == st2 and min(f["hi"], b) - max(f["lo"], a) > 0.5 * min(f["hi"] - f["lo"], b - a)
               for g2, st2, a, b in kept):
            del loci[i]["feats"][g]
        else:
            kept.append((g, st, f["lo"], f["hi"]))


def _spacers(feats: dict, st: int, qlen: int) -> None:
    """ITS1 / ITS2 between the left and right rRNA in forward coordinates; a missing rRNA is allowed when the
    contig ends within MAX_ITS nt (partial spacer), not when 18S and 28S are there and only 5.8S was not found"""
    for its in ("ITS1", "ITS2"):
        feats.pop(its, None)
    if "5.8S" not in feats and "18S" in feats and "28S" in feats:
        return
    for its, (g1, g2) in (("ITS1", ("18S", "5.8S")), ("ITS2", ("5.8S", "28S"))):
        L, R = (feats.get(g1), feats.get(g2)) if st == 1 else (feats.get(g2), feats.get(g1))
        if (L and L["p_hi"]) or (R and R["p_lo"]) or not (L or R):
            continue
        a = L["hi"] + 1 if L else (1 if R["lo"] - 1 <= MAX_ITS else None)
        b = R["lo"] - 1 if R else (qlen if qlen - L["hi"] <= MAX_ITS else None)
        if a is not None and b is not None and b - a >= 9:
            feats[its] = {"lo": a, "hi": b, "p_lo": L is None, "p_hi": R is None, "how": "between rRNAs"}


def _rescue_58s(loci: Dict[str, List[dict]], regions: Path, lib: Library, work: Path, leave_out: Dict[str, str],
                threads: int, qlens: Dict[str, int]) -> None:
    """5.8S (~155 nt) of divergent taxa is often below the megablast seed. Searched again with blastn word 11 and
    voted on like a whole unit: the spacer where 18S and 28S were found but no 5.8S, and the whole sequence
    (<= MAX_RESCUE nt) where no 5.8S was found at all (5.8S RNA records, ITS amplicons, short contigs)"""
    spacers = {}
    for q, ls in loci.items():
        if any("5.8S" in L["feats"] for L in ls):
            continue
        n = qlens[q]
        for i, L in enumerate(ls):
            f = L["feats"]
            if "18S" in f and "28S" in f:
                a, b = (f["18S"]["hi"] + 1, f["28S"]["lo"] - 1) if L["st"] == 1 else \
                    (f["28S"]["hi"] + 1, f["18S"]["lo"] - 1)
                if b - a >= 100:
                    spacers[f"{q}|{i}|{a}"] = (q, a, b)
                continue
            # one rRNA: the stretch on the spacer side of it up to the contig end (ITS amplicons, contigs ending in ITS)
            for g, side in (("18S", 1), ("28S", -1)):
                if g in f:
                    a, b = (f[g]["hi"] + 1, min(n, f[g]["hi"] + MAX_RESCUE)) if side * L["st"] == 1 else \
                        (max(1, f[g]["lo"] - MAX_RESCUE), f[g]["lo"] - 1)
                    if b - a >= 100:
                        spacers[f"{q}|{i}|{a}"] = (q, a, b)
        if not ls and n <= MAX_RESCUE:      # rDNA hits but no locus: 5.8S RNA records, short ITS contigs
            spacers[f"{q}|-1|1"] = (q, 1, n)
    if not spacers:
        return
    seqs = SeqIO.index(str(regions), "fasta")
    # one single-threaded blastn per chunk of the stretches, `threads` at a time: one blastn over thousands of
    # stretches kept few CPUs busy (41 of 57 min on a 475 Mb SPAdes assembly)
    keys = list(spacers)
    n = max(1, min(threads, len(keys)))
    for i in range(n):
        with open(work / f"rdna_58s{i}.fa", "w") as o:
            for k in keys[i::n]:
                q, a, b = spacers[k]
                o.write(f">{k}\n{seqs[q].seq[a - 1:b]}\n")

    def blast(i: int) -> str:
        return program_output([tool("blastn"), "-task", "blastn", "-word_size", "11", "-query", str(work / f"rdna_58s{i}.fa"),
                    "-db", lib.rdna_db, "-evalue", "1e-5", "-dust", "no", "-num_threads", "1",
                    "-max_target_seqs", "500",
                    "-outfmt", "6 qseqid sseqid qstart qend sstart send bitscore qseq sseq"])
    hits: Dict[str, Dict[str, list]] = defaultdict(lambda: defaultdict(list))
    with ThreadPoolExecutor(n) as pool:
        for line in "".join(pool.map(blast, range(n))).splitlines():
            k, s, qs, qe, ss, se, bs, qa, sa = line.split("\t")
            q = spacers[k][0]
            if leave_out and leave_out.get(_contig(q)) == lib.org[s][1] or _contig(q) == s:
                continue
            hits[k][s].append((int(qs), int(qe), int(ss), int(se), float(bs), qa, sa))
    for k in (k for k in keys if k in hits):
        h = hits[k]
        q, a, b = spacers[k]
        i = int(k.split("|")[1])
        for x in rdna_loci(h, lib, b - a + 1, min_overlap=40, min_bits=50):
            f = x["feats"].get("5.8S")
            if not f or (i >= 0 and x["st"] != loci[q][i]["st"]):
                continue
            # a piece cut by the searched stretch is only kept at the contig end
            if (f["p_lo"] and a > 1) or (f["p_hi"] and b < qlens[q]):
                continue
            f = dict(f, lo=f["lo"] + a - 1, hi=f["hi"] + a - 1, p_lo=f["p_lo"] and a == 1, p_hi=f["p_hi"] and b == qlens[q],
                     aln={x + a - 1 for x in f["aln"]})
            if i < 0:
                L = dict(x, feats={})
                loci[q].append(L)
            else:
                L = loci[q][i]
            L["feats"]["5.8S"] = f
            _spacers(L["feats"], L["st"], qlens[q])
            break


def _identity(chain) -> float:
    m = n = 0
    for h in chain:
        m += sum(1 for x, y in zip(h[5], h[6]) if x == y)
        n += len(h[5])
    return round(100 * m / n, 1) if n else 0.0


def annotate_rdna(regions: Path, lib: Library, work: Path, threads: int, leave_out: Dict[str, str],
                  qlens: Dict[str, int], refine: bool = True, bounded: bool = False) -> Dict[str, List[dict]]:
    fmt = "6 qseqid sseqid qstart qend sstart send bitscore qseq sseq"
    # BLAST threads scale poorly over many short queries: single-threaded blastn on chunks of the regions, `threads`
    # at a time; each chunk's hits (with the alignment strings) are turned into loci and dropped right away, so memory
    # depends on the chunk size, not on the input size
    recs = sorted(SeqIO.parse(str(regions), "fasta"), key=lambda r: -len(r))
    n_chunks = max(1, threads, -(-sum(len(r) for r in recs) // CHUNK_NT))
    chunks: List[list] = [[] for _ in range(n_chunks)]
    size = [0] * n_chunks
    for r in recs:
        i = size.index(min(size))
        chunks[i].append(r)
        size[i] += len(r)
    del recs

    def blast(i: int) -> Dict[str, List[dict]]:
        fa, tsv = work / f"rdna_q{i}.fa", work / f"rdna_blast{i}.tsv"
        SeqIO.write(chunks[i], str(fa), "fasta")
        with open(tsv, "w") as o:
            run_tool([tool("blastn"), "-task", "megablast", "-word_size", "16", "-query", str(fa),
                      "-db", lib.rdna_db,
                      "-evalue", "1e-10", "-outfmt", fmt, "-num_threads", "1",
                      "-max_target_seqs", str(BLAST_TARGETS),
                      "-dust", "no"], stdout=o, check=True)
            # curated references again with blastn word 11: their conserved ends must align in distant taxa too
            if lib.curated_db:
                run_tool([tool("blastn"), "-task", "blastn", "-word_size", "11", "-query", str(fa),
                          "-db", lib.curated_db, "-evalue", "1e-5", "-outfmt", fmt, "-num_threads", "1",
                          "-dust", "no"], stdout=o, check=True)
        hits: Dict[str, Dict[str, list]] = defaultdict(lambda: defaultdict(list))
        for line in open(tsv):
            q, s, qs, qe, ss, se, bs, qa, sa = line.rstrip("\n").split("\t")
            if leave_out and leave_out.get(_contig(q)) == lib.org[s][1] or _contig(q) == s:
                continue
            h = (int(qs), int(qe), int(ss), int(se), float(bs), qa, sa)
            if h not in hits[q][s]:      # a curated reference found by both searches
                hits[q][s].append(h)
        fa.unlink()
        tsv.unlink()
        return {q: rdna_loci(h, lib, qlens[q]) for q, h in hits.items()}

    t0 = time.time()
    loci: Dict[str, List[dict]] = {}
    with ThreadPoolExecutor(max(1, threads)) as pool:
        for part in pool.map(blast, [i for i in range(n_chunks) if chunks[i]]):
            loci.update(part)
    logging.info(f"[NUC] rDNA blastn + loci ({n_chunks} chunks): {time.time() - t0:.0f} s")
    t0 = time.time()
    _rescue_58s(loci, regions, lib, work, leave_out, threads, qlens)
    logging.info(f"[NUC] 5.8S rescue: {time.time() - t0:.0f} s")
    seqs = SeqIO.index(str(regions), "fasta")
    for q, ls in loci.items():
        _tidy(ls, str(seqs[q].seq).upper(), qlens[q], final=not refine, ref_gene=lambda r, g: _ref_gene(lib, r, g))
    if refine:
        rs = {q: seqs[q].seq for q in loci}
        t0 = time.time()
        n = cm_refine(loci, rs, lib.dir, work, threads, bounded)
        logging.info(f"[NUC] rRNA ends placed by the Rfam models: {n} ({time.time() - t0:.0f} s)")
        for q, ls in loci.items():
            _tidy(ls, str(rs[q]).upper(), qlens[q], final=True, ref_gene=lambda r, g: _ref_gene(lib, r, g))
    return {q: [L for L in ls if L["feats"]] for q, ls in loci.items()}


# Rfam covariance models: the mature ends of the curated RefSeq rRNAs (human, rat, macaque, rabbit, cow, hamster,
# Drosophila, C. elegans) sit at column 1 / the last column, with CM_OFFSET nt outside (5', 3'), measured with
# cmalign on the 28 curated records 2026-10-25 (C. elegans 5.8S 3' +2, rat 28S 3' +0: the mammal value is used)
RFAM = {"18S": ("RF01960", 1831), "5.8S": ("RF00002", 154), "28S": ("RF02543", 3401)}


CM_OFFSET = {"18S": (0, 0), "5.8S": (1, 3), "28S": (2, 1)}


CM_WIN = 600      # nt on each side of the predicted end searched with the model


CM_CONF = 5       # end columns that must hold residues of posterior probability >= 0.95


HMM_MAX_EXTEND = 3     # HMM hits that stop at most this many columns short of the model end are taken ...


EDGE_MARGIN = 20       # nt an end must leave to the outer edge of its search window


END_MARGIN = 20        # an unconfirmed end this close to the sequence end is partial


HMM_MAX_SHIFT = 30     # ... when they move the end at most this far from the reference projection


SAME_SPECIES_ID = 99.5  # a curated reference this close (same species) decides its boundary even against the model


def _read_sto(path: Path) -> Tuple[str, Dict[str, str], Dict[str, str]]:
    seq: Dict[str, str] = defaultdict(str)
    pp: Dict[str, str] = defaultdict(str)
    rf = ""
    for line in open(path):
        if line.startswith("#=GC RF"):
            rf += line.split()[-1]
        elif line.startswith("#=GR") and line.split()[2] == "PP":
            _, k, _, v = line.split()
            pp[k] += v
        elif line.strip() and not line.startswith("#") and not line.startswith("//"):
            k, v = line.split()
            seq[k] += v
    return rf, seq, pp


def _cmsearch(cm: Path, seqs: List[Tuple[str, str]], work: Path, tag: str, threads: int, opts: List[str],
              align: bool = False, toponly: bool = True) -> List[tuple]:
    """cmsearch of (name, sequence) pairs in `threads` pieces at once: one cmsearch reads a file of short windows as
    one block for one worker (38 windows of the 28S model: 14 s at --cpu 16, as at --cpu 1). -Z is the size of all
    pieces together (nt, x2 for both strands), so the E-values and filter thresholds are those of one search of all
    (Infernal truncates -Z to whole nt). Returns per piece (tblout rows, -A alignment: RF line, sequences, PP)."""
    n = max(1, min(threads, len(seqs)))
    z = sum(len(s) for _, s in seqs) * (1 if toponly else 2)

    def one(i: int) -> tuple:
        fa, tbl, sto = (work / f"cm_{tag}_{i}.{x}" for x in ("fa", "tbl", "sto"))
        with open(fa, "w") as o:
            o.writelines(f">{k}\n{s}\n" for k, s in seqs[i::n])
        program_output([tool("cmsearch"), "--cpu", "1", "--notextw", "-Z", f"{(z + 0.5) / 1e6:.9f}", *opts, "--tblout", str(tbl),
             *(["-A", str(sto)] if align else []), str(cm), str(fa)])
        rows = [line.split() for line in open(tbl) if not line.startswith("#")]
        return (rows, *(_read_sto(sto) if align and sto.exists() and sto.stat().st_size else ("", {}, {})))
    with ThreadPoolExecutor(n) as pool:
        return list(pool.map(one, range(n)))


def cm_ends(jobs: List[tuple], seqs, gene: str, cm: Path, work: Path, threads: int, win: int = CM_WIN,
            hmm_only: bool = False, max_extend: int = 0, edge: Optional[set] = None,
            trunc: Optional[set] = None) -> Dict[int, int]:
    """jobs: (contig, strand, side 5|3, predicted end); returns job index -> end position homologous to the curated
    mature end (cmsearch on a window of win nt around the end; the hit must reach the model's first (5') / last
    (3') column with confident residues and the window must hold sequence beyond the end).
    hmm_only: the model's profile HMM only (no structure, ~25x faster); its local hits stop a few columns short of
    the model ends, which are added ungapped when at most max_extend columns are missing (no posterior check).
    edge: filled with the jobs whose hit runs into the outer window edge without reaching the model end (the end
    is further away than win); trunc: the same where that edge is the sequence end (the gene is cut there)"""
    tag = f"{gene}_{win}{'_hmm' if hmm_only else ''}"
    meta, windows = {}, []
    for i, (c, st, side, p) in enumerate(jobs):
        n = len(seqs[c])
        lo, hi = max(1, p - win), min(n, p + win)
        s = seqs[c][lo - 1:hi]
        # the window's outer edge (beyond the end) is the sequence end itself
        meta[f"j{i}"] = (i, st, side, lo, hi, (lo == 1) if (side == 5) == (st == 1) else (hi == n))
        windows.append((f"j{i}", str(s if st == 1 else s.reverse_complement())))
    o5, o3 = CM_OFFSET[gene]
    clen = RFAM[gene][1]
    out: Dict[int, int] = {}
    best: Dict[str, tuple] = {}
    aln: Dict[str, tuple] = {}             # window -> the alignment of its piece (columns, sequences, PP)
    for rows, rf, al, pp in _cmsearch(cm, windows, work, tag, threads, ["-E", "1e-5", "--toponly"] +
                                      (["--hmmonly"] if hmm_only else []), align=not hmm_only):
        cons = [i for i, ch in enumerate(rf) if ch not in ".-~"]
        for x in rows:
            if x[9] == "+" and (x[0] not in best or float(x[14]) > best[x[0]][4]):
                best[x[0]] = (int(x[5]), int(x[6]), int(x[7]), int(x[8]), float(x[14]))
                aln[x[0]] = (cons, al, pp)
    for k, (mf, mt, sf, se, _) in best.items():
        cons, al, pp = aln[k]
        i, st, side, lo, hi, at_seq_end = meta[k]
        wl = hi - lo + 1

        def room(t: int) -> bool:
            """an end placed close to the outer window edge is an artefact of the window (OY756919.1: CM end 4 nt
            from a +-150 window edge, not confident in +-600), unless that edge is the sequence end"""
            ok = at_seq_end or (t - 1 if side == 5 else wl - t) >= EDGE_MARGIN
            if not ok and edge is not None:
                edge.add(i)
            return ok
        runs_out = (side == 5 and mf > 1 and sf <= 3) or (side == 3 and mt < clen and se >= wl - 2)
        if runs_out and at_seq_end and trunc is not None:
            trunc.add(i)
        if hmm_only:
            if side == 5 and mf - 1 <= max_extend and sf - (mf - 1) - o5 >= 1:
                t = sf - (mf - 1) - o5
            elif side == 3 and clen - mt <= max_extend and se + (clen - mt) + o3 <= wl:
                t = se + (clen - mt) + o3
            else:
                continue
            if not room(t):
                continue
            pos = lo + t - 1 if st == 1 else hi - t + 1
            # without posterior probabilities a move far from the reference projection is left to the full model
            # (OY756919.1 Polypodium 28S 3': HMM +146 nt, the covariance model not confident)
            if abs(pos - jobs[i][3]) <= HMM_MAX_SHIFT:
                out[i] = pos
            continue
        if edge is not None and not at_seq_end and runs_out:
            edge.add(i)
        s, p = al.get(f"{k}/{sf}-{se}"), pp.get(f"{k}/{sf}-{se}")
        if s is None or p is None:        # no posterior probabilities for the hit: not confident
            continue
        cols = cons[:CM_CONF] if side == 5 else cons[::-1][:CM_CONF]
        if not all(s[c] not in ".-~" and p[c] == "*" for c in cols):
            continue
        # the end may be the sequence end itself (mature rRNA records, e.g. PDB chains) when the offset nt are there;
        # lineages have up to 1 nt fewer outside the model (C. elegans 5.8S 3' +2 vs +3), so a sequence end 1 nt
        # short of the offset is that end (9BH5_A8, C. elegans 5.8S chain)
        if side == 5 and mf == 1 and sf - o5 >= 1:
            t = sf - o5
        elif side == 5 and mf == 1 and at_seq_end and o5 and sf - o5 == 0:
            t = 1
        elif side == 3 and mt == clen and se + o3 <= wl:
            t = se + o3
        elif side == 3 and mt == clen and at_seq_end and o3 and se + o3 == wl + 1:
            t = wl
        else:
            continue
        if room(t):
            out[i] = lo + t - 1 if st == 1 else hi - t + 1
    return out


def _cm_rescue_58s(loci: Dict[str, List[dict]], seqs, d: Path, work: Path, threads: int) -> int:
    """a 5.8S cut by the sequence start / end next to its ITS2 (ITS2 amplicons: the primer sits in the 5.8S) is often
    too short for blastn; the stretch between the sequence end and the 28S is searched with the 5.8S model, and a
    hit ending at the model's last column becomes a 5' partial 5.8S"""
    jobs, where = [], []
    for q, ls in loci.items():
        n = len(seqs[q])
        for L in ls:
            f, st = L["feats"], L["st"]
            if "5.8S" in f or "28S" not in f or "18S" in f:
                continue
            a, b = (1, f["28S"]["lo"] - 1) if st == 1 else (f["28S"]["hi"] + 1, n)
            if b - a < 20 or (st == 1 and a != 1) or (st == -1 and b != n):
                continue
            jobs.append((q, st, a, b))
            where.append(L)
    if not jobs:
        return 0
    stretches = [(f"r{i}", str(seqs[q][a - 1:b] if st == 1 else seqs[q][a - 1:b].reverse_complement()))
                 for i, (q, st, a, b) in enumerate(jobs)]
    o3, done = CM_OFFSET["5.8S"][1], 0
    for x in (x for rows, *_ in _cmsearch(d / f"{RFAM['5.8S'][0]}.cm", stretches, work, "58s_rescue", threads,
                                          ["-E", "1e-3"], toponly=False) for x in rows):
        i, mf, mt, sf, se = int(x[0][1:]), int(x[5]), int(x[6]), int(x[7]), int(x[8])
        q, st, a, b = jobs[i]
        L = where[i]
        if x[9] != "+" or mt != RFAM["5.8S"][1] or sf > 3 or mf == 1 or "5.8S" in L["feats"]:
            continue        # must be the 3' part of the 5.8S running from the sequence start
        t = min(se + o3, b - a + 1)
        lo, hi = (a, a + t - 1) if st == 1 else (b - t + 1, b)
        L["feats"]["5.8S"] = {"lo": lo, "hi": hi, "p_lo": st == 1, "p_hi": st == -1, "sup": 0.0, "aligned": t,
                              "aln": set(range(lo, hi + 1)), "how": "projected/cm"}
        _spacers(L["feats"], st, len(seqs[q]))
        done += 1
    return done


def cm_refine(loci: Dict[str, List[dict]], seqs, lib_dir: Path, work: Path, threads: int,
              bounded: bool = False) -> int:
    """move every complete rRNA end to the position homologous to the curated mature end (Rfam CM); ends the model
    does not place confidently keep the reference projection. Returns the number of ends placed by the model.
    bounded (assembly): the covariance model searches only the MAX_CM_JOBS best supported ends per gene, in +-150 nt
    windows"""
    d = lib_dir / "rfam"
    if not d.is_dir():
        d = LIB / "rfam"
    if not tool_ok("cmsearch") or not d.is_dir():
        logging.warning("[NUC] cmsearch or the Rfam models not found: rRNA ends from the reference projection only")
        return 0
    placed = _cm_rescue_58s(loci, seqs, d, work, threads)
    for g, (acc, _) in RFAM.items():
        jobs, where = [], []
        for q, ls in loci.items():
            for L in ls:
                f = L["feats"].get(g)
                if not f:
                    continue
                st = L["st"]
                n = len(seqs[q])
                for side in (5, 3):
                    lo_side = (side == 5) == (st == 1)
                    pos = f["lo"] if lo_side else f["hi"]
                    # a partial end at the sequence end may be the mature end itself (rRNA records, PDB chains)
                    if (f["p_lo"] if lo_side else f["p_hi"]) and min(pos - 1, n - pos) > CM_WIN:
                        continue
                    # an end decided by a curated reference of (nearly) the same species keeps it: the nt outside
                    # the model columns differ between lineages (C. elegans 5.8S 3' +2, mammals +3)
                    if f["how"].split("/")[0 if side == 5 else 1] == "curated" and \
                            f.get("cur_id", {}).get(side, 0) >= SAME_SPECIES_ID:
                        continue
                    jobs.append((q, st, side, pos))
                    where.append((f, lo_side, st))
        if not jobs:
            continue
        # fast first: the model's profile HMM on +-150 nt agrees with the full model whenever its hit stops <= 3
        # columns short of the model end (4366 ends, benchmark/nuc_cm_speed.py, 2026-10-25); the other ends go to
        # the covariance model on +-150 nt, then +-600 nt
        placed_at: Dict[int, Tuple[int, str]] = {}
        todo, log, cut = list(range(len(jobs))), [], set()
        for how, kw in (("hmm", dict(win=150, hmm_only=True, max_extend=HMM_MAX_EXTEND)),
                        ("cm", dict(win=150)), ("cm", dict(win=CM_WIN))):
            if kw["win"] == CM_WIN:          # the wide window only where the hit ran out of the narrow one
                todo = [i for i in todo if i in cut]
            # assembly: a covariance model search takes 0.5-4 s per +-150 nt window and 15-45 s per +-600 nt window
            # (which placed 1 end in 30 assemblies): the narrow windows of the best supported genes only
            cap = MAX_CM_JOBS if kw["win"] < CM_WIN else 0
            if bounded and how == "cm" and len(todo) > cap:
                log.append(f"{len(todo) - cap} of {len(todo)} not searched")
                todo = sorted(sorted(todo, key=lambda i: -where[i][0].get("sup", 0))[:cap])
            if not todo:
                break
            edge: set = set()
            trunc: set = set()
            t0 = time.time()
            res = cm_ends([jobs[i] for i in todo], seqs, g, d / f"{acc}.cm", work, threads, edge=edge, trunc=trunc,
                          **kw)
            log.append(f"{how} {kw['win']}: {len(res)}/{len(todo)} {time.time() - t0:.0f} s")
            placed_at.update({todo[k]: (pos, how) for k, pos in res.items()})
            cut = {todo[k] for k in edge}
            # a gene the profile already shows cut by the sequence end cannot get a complete end from the full model
            todo = [i for k, i in enumerate(todo) if k not in res and k not in trunc]
        logging.info(f"[NUC] {g} ends: " + ", ".join(log))
        for i, (pos, how_) in placed_at.items():
            f, lo_side, st = where[i]
            f["lo" if lo_side else "hi"] = pos
            f["p_lo" if lo_side else "p_hi"] = False
            how = f["how"].split("/")                   # 5' / 3' evidence
            how[0 if jobs[i][2] == 5 else 1] = how_
            f["how"] = "/".join(how)
            placed += 1
        # an end the models could not place that lies at the sequence end is not shown to be complete: partial,
        # run to the sequence end (LC765217.1 18S "complete" at nt 2, OK209792.1 at nt 1: both start inside the 18S)
        for i, (q, st, side, pos) in enumerate(jobs):
            if i in placed_at:
                continue
            f, lo_side, _ = where[i]
            n = len(seqs[q])
            if lo_side and pos <= END_MARGIN:
                f.update(lo=1, p_lo=True)
            elif not lo_side and n - pos < END_MARGIN:
                f.update(hi=n, p_hi=True)
    return placed


def _ends(f: dict, st: int) -> Tuple[str, str]:
    """evidence of the low / high end in forward coordinates ("how" is 5'/3')"""
    a, b = f["how"].split("/")
    return (a, b) if st == 1 else (b, a)


def _piece_matches_ref(f: dict, st: int, seq: str, ref_seq: Optional[str], min_id: float = PIECE_MIN_ID) -> bool:
    """an rRNA piece cut by the contig end is the reference gene's end: >= min_id ungapped identity with the
    reference's 3' end (5'-partial piece) or 5' end (3'-partial piece)"""
    if not ref_seq:
        return False
    piece = seq[f["lo"] - 1:f["hi"]].upper()
    if st == -1:
        piece = str(Seq(piece).reverse_complement())
    p5, p3 = (f["p_lo"], f["p_hi"]) if st == 1 else (f["p_hi"], f["p_lo"])
    if p5 == p3 or len(piece) > len(ref_seq):
        return False
    other = ref_seq[-len(piece):] if p5 else ref_seq[:len(piece)]
    return sum(a == b for a, b in zip(piece, other.upper())) >= min_id * len(piece)


def _tidy(loci: List[dict], seq: str, qlen: int, final: bool = False, ref_gene=None) -> None:
    """(0) final pass (after the Rfam models had their chance): an rRNA with no aligned nt inside it and both ends
    only projected stays only when its sequence is the reference gene's end (_piece_matches_ref; ref_gene(ref, gene)
    gives it) - MJ20260820041 CHY007: 20 nt of ITS1 at a contig end became a "5'-partial 18S" projected from the
    5.8S across the reference's shorter ITS1 (no likeness to the conserved 18S 3' end), cutting the sample's
    519-nt ITS1 to a "complete" 499 nt; a real end too short to align (7 nt of 5.8S) matches and stays;
    (1) no end inside a run of >= 10 N; a projected end does not cross a run of >= 10 N (scaffold gap of unknown length) when no part of the
    gene beyond the gap aligns: the gene stops at the gap, partial; (2) two different rRNAs never overlap, on either strand (a reverse copy projected over a unit
    is false): the one with fewer aligned bits is dropped when both its ends are projected or it is on the other
    strand, else cut back to the other one (partial on that side)"""
    gaps = [(m.start() + 1, m.end()) for m in re.finditer("N{10,}", seq)]
    if final:
        for L in loci:
            for g in [g for g, f in L["feats"].items() if g in RRNA and not f.get("aligned", 1)
                      and _ends(f, L["st"]) == ("projected", "projected")]:
                ref = ref_gene(L.get("ref"), g) if ref_gene and L.get("ref") else None
                if not _piece_matches_ref(L["feats"][g], L["st"], seq, ref):
                    del L["feats"][g]
    for L in loci:
        for g, f in L["feats"].items():
            if g not in RRNA:
                continue
            lo_how, hi_how = _ends(f, L["st"])
            aln = f.get("aln", set())
            for a, b in gaps:
                if a <= f["lo"] <= b < f["hi"]:          # an end inside the gap: the gene starts after it
                    f.update(lo=b + 1, p_lo=True)
                if f["lo"] < a <= f["hi"] <= b:
                    f.update(hi=a - 1, p_hi=True)
                if not (f["lo"] < a and b < f["hi"]):
                    continue
                # only an end with no alignment of the gene beyond the gap (a gap inside an aligned gene is kept,
                # as most submitters do)
                if hi_how == "projected" and not any(x > b for x in aln) and any(x < a for x in aln):
                    f.update(hi=a - 1, p_hi=True)
                elif lo_how == "projected" and not any(x < a for x in aln) and any(x > b for x in aln):
                    f.update(lo=b + 1, p_lo=True)
    fs = sorted(((f["sup"], i, g) for i, L in enumerate(loci) for g, f in L["feats"].items() if g in RRNA),
                reverse=True)
    kept: List[tuple] = []
    for _, i, g in fs:
        L = loci[i]
        f = L["feats"][g]
        for a, b, st in kept:
            if min(f["hi"], b) - max(f["lo"], a) < 20:
                continue
            if _ends(f, L["st"]) == ("projected", "projected") or st != L["st"]:
                del L["feats"][g]
                break
            if f["lo"] < a:
                f.update(hi=min(f["hi"], a - 1), p_hi=True)
            else:
                f.update(lo=max(f["lo"], b + 1), p_lo=True)
            if f["hi"] - f["lo"] < 20:
                del L["feats"][g]
                break
        else:
            kept.append((f["lo"], f["hi"], L["st"]))
    for L in loci:
        _spacers(L["feats"], L["st"], qlen)
    # (3) a spacer neither overlaps an rRNA (of any copy) nor starts / ends in a gap (table2asn
    # BadRRNAcomponentOverlapRRNA, FeatureInsideGap)
    rr = [(f["lo"], f["hi"]) for L in loci for g, f in L["feats"].items() if g in RRNA]
    for L in loci:
        for its in ("ITS1", "ITS2"):
            f = L["feats"].get(its)
            if not f:
                continue
            for a, b in rr + gaps:
                if a <= f["lo"] <= b:
                    f.update(lo=b + 1, p_lo=True)
                if a <= f["hi"] <= b:
                    f.update(hi=a - 1, p_hi=True)
            if f["hi"] - f["lo"] < 9 or any(f["lo"] < a and b < f["hi"] for a, b in rr):
                del L["feats"][its]
