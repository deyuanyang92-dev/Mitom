#!/usr/bin/env python3
"""trna_ml.py - learned mitochondrial tRNA caller for invertebrates (MitoM-ML), used by MitoM's annotation core
(mito_core.py) and by mito_trna.py --tools mitom.

Built and benchmarked in projects/invert_trna (REQUIREMENTS #53). For one sequence:
  1. scan   cmsearch with the 22 MITOS2 refseq89m tRNA CMs (MiTFi models, Juhling 2012 NAR 40:2833) plus the
            clade CMs on the record's NCBI lineage (models/trna_ml/clades: 220 clades x 22 families, built from
            RefSeq invertebrate mitogenomes, 200 seqs per model, cmcalibrate -L 0.02); circular sequence + 300 nt;
            a --mid second pass for families without a hit at E <= 1e-3; anticodon read from the alignment of
            every hit to the pan CM of its family; L1/L2 and S1/S2 named by that anticodon.
  2. call   loci (same-strand hits overlapping >= 50 %) x families -> candidates with 28 features (pan / clade
            scores, model coverage, competing families, genome context, ARWEN / ARAGORN support, opposite-strand
            "shadow" of a stronger tRNA); a gradient-boosting classifier keeps the best candidate per locus with
            p >= threshold; same-strand calls overlapping > 50 % keep the more probable.
  3. ends   candidate (start, end) pairs from all source ends +-1 scored by a second model (source agreement,
            CM end reached, acceptor-stem pairing).
The two models are trained on first use from the training tables shipped in models/trna_ml (genomes of
held-out families, i.e. scanned with clade models that never saw their family) and cached as
models/trna_ml/cache/*.pkl, so no pickled model is downloaded or shared between scikit-learn versions.

Benchmark (1,069 RefSeq invertebrate mitogenomes of families unseen by the clade models and absent from MITOS2
refseq89m; 22,099 RefSeq tRNAs; projects/invert_trna/REPORT.md): see that report for the numbers.

Usage (normally through mito_trna.py):
  python3 trna_ml.py -i genome.fa --lineage "Metazoa;Arthropoda;..." --genetic_code 5
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import pickle   # models trained locally from the shipped tables and read back by this module only
import statistics
import sys
import tarfile
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

HERE = Path(__file__).resolve().parent
MODELS = HERE / "models" / "trna_ml"
FAMS = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L1", "L2", "M", "N", "P", "Q", "R", "S1", "S2", "T", "V", "W", "Y"]
WRAP = 300
RESCUE_E = 1e-3
RESCUE_MODE = "--mid"            # --max (MITOS2 sensitive pass) was > 10 min per genome, --mid 13 s
THRESHOLD = 0.5
COMP = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")
# anticodon (5'->3', DNA) -> Leu / Ser family; invertebrate convention L1 = CUN, L2 = UUR, S1 = AGN, S2 = UCN
LS_FAM = {"tag": "L1", "taa": "L2", "tct": "S1", "gct": "S1", "tga": "S2", "cct": "S1", "tcc": "S1", "act": "S1"}
EXT_LS = {"tag": "L1", "taa": "L2", "tct": "S1", "tga": "S2", "gct": "S1", "tcc": "S1"}
FEATS = ["pan_bits", "pan_logE", "clade_bits", "clade_logE", "best_logE", "depth", "n_models", "n_clade_models",
         "len", "cov", "trunc", "delta_bits_next_fam", "fam_rank", "fam_n_cands", "fam_gap_to_best", "best_fam_at_locus",
         "rel_to_fam_best", "is_LS", "overlap_stronger", "n_strong_genome", "gc", "ac_known", "ac_match", "arwen", "aragorn",
         "only_ext", "overlap_stronger_opp", "opp_bits_ratio"]
BFEATS = ["s_n", "s_bits", "s_pan", "s_clade", "s_arwen", "s_aragorn", "s_full", "s_shift", "s_off",
          "e_n", "e_bits", "e_pan", "e_clade", "e_arwen", "e_aragorn", "e_full", "e_shift", "e_off",
          "len", "acc0", "acc1", "acc2", "wc1", "first_nt_T", "last_nt_A"]
PAIR = {("G", "C"), ("C", "G"), ("A", "T"), ("T", "A"), ("G", "T"), ("T", "G")}
WC = {("G", "C"), ("C", "G"), ("A", "T"), ("T", "A")}


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def sh(cmd: List[str]) -> None:
    from mito_common import run_tool
    run_tool(cmd)


# ── models ───────────────────────────────────────────────────────────────────
def clade_dir(models: Path = MODELS) -> Path:
    """Clade CMs, unpacked from clades.tar.gz on first use."""
    d = models / "clades"
    if not (d / "models.json").exists():
        tgz = models / "clades.tar.gz"
        if not tgz.exists():
            raise FileNotFoundError(f"{tgz} missing (MitoM-ML clade models)")
        with tarfile.open(tgz) as t:
            t.extractall(models, filter="data")
    return d


def _train(table: Path, feats: List[str], out: Path, **kw):
    from sklearn.ensemble import HistGradientBoostingClassifier
    with gzip.open(table, "rt", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    clf = HistGradientBoostingClassifier(random_state=0, **kw)
    clf.fit([[float(r[k]) for k in feats] for r in rows], [int(r["label"]) for r in rows],
            sample_weight=[float(r.get("weight") or 1) for r in rows])   # shuffled-genome negatives: weight 0.1
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "wb") as fh:
        pickle.dump(clf, fh)
    return clf


def classifiers(models: Path = MODELS):
    """(caller, boundary model): cached pickles, trained here from the shipped tables when missing or older."""
    import sklearn
    cache = models / "cache" / f"sklearn{sklearn.__version__}"
    out = []
    for name, feats, kw in (("caller", FEATS, dict(max_iter=400, learning_rate=0.05, max_leaf_nodes=31,
                                                   l2_regularization=1.0)),
                            ("bounds", BFEATS, dict(max_iter=500, learning_rate=0.05, max_leaf_nodes=31))):
        table, pkl = models / f"{name}_train.tsv.gz", cache / f"{name}.pkl"
        if pkl.exists() and pkl.stat().st_mtime >= table.stat().st_mtime:
            with open(pkl, "rb") as fh:
                out.append(pickle.load(fh))
        else:
            out.append(_train(table, feats, pkl, **kw))
    return out[0], out[1]


def model_paths(pan: Path, cdir: Optional[Path], lineage: str) -> List[Tuple[Path, str, int]]:
    """[(cm, clade, depth)]: the 22 pan CMs (depth 0) and every clade CM whose clade is on the lineage."""
    out = [(pan / f"trn{f}.cm", "pan", 0) for f in FAMS]
    if cdir is None or not lineage:
        return out
    lin = [x.strip() for x in lineage.replace("; ", ";").split(";")]
    meta = json.load(open(cdir / "models.json"))["models"]
    for name, m in sorted(meta.items()):
        if m["clade"] in lin and m.get("status") in ("ok", "cached"):
            out.append((cdir / f"{name}.cm", m["clade"], lin.index(m["clade"])))
    return out


# ── 1. scan ──────────────────────────────────────────────────────────────────
def _hairpin_loops(ss: str) -> List[Tuple[int, int]]:
    out, i, n = [], 0, len(ss)
    while i < n:
        if ss[i] in "<([{":
            j = i
            while j < n and ss[j] in "<([{":
                j += 1
            k = j
            while k < n and ss[k] not in "<([{>)]}":
                k += 1
            if k < n and ss[k] in ">)]}":
                out.append((j, k))
            i = k
        else:
            i += 1
    return out


def _parse_sto(text: str) -> Tuple[Dict[str, str], str]:
    seqs: Dict[str, str] = {}
    ss = ""
    for line in text.splitlines():
        if line.startswith("#=GC SS_cons"):
            ss += line.split()[-1]
        elif line and not line.startswith("#") and not line.startswith("//"):
            k, v = line.split()[:2]
            seqs[k] = seqs.get(k, "") + v
    return seqs, ss


def anticodon_from_alignment(aln: str, ss: str) -> str:
    """Middle 3 nt of the 7-nt anticodon loop (2nd of the 3 cloverleaf hairpins)."""
    loops = _hairpin_loops(ss)
    if len(loops) < 2:
        return ""
    lo = loops[1] if len(loops) >= 3 else loops[0]
    cols = list(range(lo[0], lo[1]))
    nts = [aln[c] for c in cols if aln[c] not in ".-~"]
    if len(nts) != 7 or len([c for c in cols if ss[c] != "~"]) < 5:
        return ""
    return "".join(nts[2:5]).lower().replace("u", "t")


def _clen(cm: Path) -> int:
    for line in open(cm):
        if line.startswith("CLEN"):
            return int(line.split()[1])
    return 0


def _parse_tbl(tbl: Path, info, L: int, mode: str) -> List[Dict[str, Any]]:
    rows = []
    for line in tbl.read_text().splitlines() if tbl.exists() else []:
        if line.startswith("#"):
            continue
        c = line.split()
        clade, depth, clen, fam = info.get(c[2], ("?", -1, 0, "?"))
        s, e = int(c[7]), int(c[8])
        lo, hi = min(s, e), max(s, e)
        if lo > L:                                   # duplicate in the appended part
            continue
        rows.append({"model": c[2], "clade": clade, "depth": depth, "family": fam, "start": lo, "end": hi,
                     "strand": 1 if c[9] == "+" else -1, "mfrom": int(c[5]), "mto": int(c[6]), "clen": clen,
                     "trunc": c[10], "bits": float(c[14]), "evalue": float(c[15]), "gc": float(c[12]), "mode": mode})
    return rows


def drop_wrap_duplicates(rows: List[Dict[str, Any]], L: int) -> List[Dict[str, Any]]:
    """A gene at the sequence start is found twice in seq + seq[:WRAP]: at 1.. and, with its ends +-1 nt, again from
    position L in the appended copy (a hit starting at L is kept as one across the origin). Of a hit across the origin
    and a hit of the same model and strand inside the sequence that cover >= 50 % of the shorter once shifted by L,
    the higher-scoring one stays (the one inside on a tie): NC_067910.1 trnM 1..64 was called twice, as 1..63 and as
    14404..14468 (21 duplicate calls in 1,057 test records). A real gene across the origin outscores the truncated
    hit of its part after the origin."""
    drop = set()
    for i, h in enumerate(rows):
        if h["end"] <= L:
            continue
        for j, x in enumerate(rows):
            if x["end"] > L or x["model"] != h["model"] or x["strand"] != h["strand"] or j in drop:
                continue
            o = min(x["end"], h["end"] - L) - max(x["start"], h["start"] - L) + 1
            if o >= 0.5 * min(x["end"] - x["start"] + 1, h["end"] - h["start"] + 1):
                drop.add(i if h["bits"] <= x["bits"] else j)
                break
    return [h for i, h in enumerate(rows) if i not in drop]


def scan(seq: str, lineage: str, pan: Path, cdir: Optional[Path], cmsearch: str, cmalign: str, tmp: Path,
         cpu: int = 1, rescue: bool = True) -> List[Dict[str, Any]]:
    seq = seq.upper()
    L = len(seq)
    q = tmp / "q.fa"
    ext = seq + seq[:WRAP]
    q.write_text(f">q\n{ext}\n")
    paths = model_paths(pan, cdir, lineage)
    info = {}
    with open(tmp / "all.cm", "w") as fh:
        for p, clade, depth in paths:
            txt = p.read_text()
            fh.write(txt)
            nm = next(x.split()[1] for x in txt.splitlines() if x.startswith("NAME"))
            info[nm] = (clade, depth, _clen(p), p.stem.split("__")[-1].replace("trn", ""))
    tbl = tmp / "hits.tbl"
    base = [cmsearch, "--cpu", str(cpu), "--tblout", str(tbl), "-E", "10", "--notextw"]
    sh(base + [str(tmp / "all.cm"), str(q)])
    rows = _parse_tbl(tbl, info, L, "")
    if rescue:
        good = {r["family"] for r in rows if r["evalue"] <= RESCUE_E}
        miss = [f for f in FAMS if f not in good]
        if miss:
            with open(tmp / "miss.cm", "w") as fh:
                for p, _, _ in paths:
                    if p.stem.split("__")[-1].replace("trn", "") in miss:
                        fh.write(p.read_text())
            sh(base + [RESCUE_MODE, str(tmp / "miss.cm"), str(q)])
            seen = {(r["model"], r["start"], r["end"]) for r in rows}
            rows += [r for r in _parse_tbl(tbl, info, L, "max") if (r["model"], r["start"], r["end"]) not in seen]
    rows = drop_wrap_duplicates(rows, L)
    by: Dict[str, List[int]] = {}
    for i, r in enumerate(rows):
        r["anticodon"] = ""
        if r["evalue"] <= 1.0:
            by.setdefault(r["family"], []).append(i)
    for fam, idx in by.items():
        fa, sto = tmp / f"ac_{fam}.fa", tmp / f"ac_{fam}.sto"
        with open(fa, "w") as fh:
            for i in idx:
                r = rows[i]
                sub = ext[r["start"] - 1:r["end"]]
                fh.write(f">h{i}\n{sub if r['strand'] == 1 else rc(sub)}\n")
        sh([cmalign, "--cpu", str(cpu), "--outformat", "Stockholm", "-o", str(sto), str(pan / f"trn{fam}.cm"), str(fa)])
        if sto.exists():
            al, ss = _parse_sto(sto.read_text())
            for k, v in al.items():
                rows[int(k[1:])]["anticodon"] = anticodon_from_alignment(v, ss)
    return rows


# ── 2. call ──────────────────────────────────────────────────────────────────
def ac_aa(ac: str, code: int) -> str:
    from Bio.Data.CodonTable import unambiguous_dna_by_id
    if len(ac) != 3 or set(ac) - set("acgt"):
        return ""
    codon = ac[::-1].translate(str.maketrans("acgt", "TGCA"))
    t = unambiguous_dna_by_id.get(int(code or 5), unambiguous_dna_by_id[5])
    return t.forward_table.get(codon, "*")


def lg(e: float) -> float:
    return math.log10(max(e, 1e-300))


def _ov(a, b) -> int:
    return min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1


def _loci(hits):
    out: List[List[Dict[str, Any]]] = []
    for h in sorted(hits, key=lambda h: h["evalue"]):
        for L in out:
            r = L[0]
            if r["strand"] == h["strand"] and _ov(r, h) >= 0.5 * min(r["end"] - r["start"] + 1, h["end"] - h["start"] + 1):
                L.append(h)
                break
        else:
            out.append([h])
    return out


def ext_hits(preds) -> List[Dict[str, Any]]:
    """ARWEN / ARAGORN predictions (mito_trna.Pred: source, start, end, strand, aa, anticodon) as pseudo-hits."""
    out = []
    for p in preds:
        fam = EXT_LS.get(p.anticodon, "") if p.aa in "LS" else p.aa
        if fam and p.source in ("ARWEN", "ARAGORN"):
            out.append({"model": p.source, "clade": "ext", "depth": -1, "family": fam, "start": p.start, "end": p.end,
                        "strand": p.strand, "evalue": 1.0, "bits": 0.0, "gc": 0.0, "mfrom": 1, "mto": 1, "clen": 1,
                        "trunc": "no", "anticodon": "", "ac_known": 1, "ac_match": 1})
    return out


def candidates(hits: List[Dict[str, Any]], ext: List[Dict[str, Any]], code: int) -> List[Dict[str, Any]]:
    """Candidates (locus x family) with the FEATS features; same computation as projects/invert_trna learn.py."""
    for h in hits:
        a = ac_aa(h.get("anticodon", ""), code)
        h["ac_known"] = int(bool(a))
        if a and h["family"][0] == a and a in "LS" and LS_FAM.get(h["anticodon"], h["family"]) != h["family"]:
            h["family"] = LS_FAM[h["anticodon"]]
        h["ac_match"] = int(a == h["family"][0]) if a else 0
    hs = hits + ext
    cands = []
    for li, L in enumerate(_loci(hs)):
        byfam = defaultdict(list)
        for h in L:
            byfam[h["family"]].append(h)
        fam_best = {fa: max(v, key=lambda h: h["bits"]) for fa, v in byfam.items()}
        order = sorted(fam_best, key=lambda fa: -fam_best[fa]["bits"])
        for fa, v0 in byfam.items():
            ex = {h["model"] for h in v0 if h["clade"] == "ext"}
            v = [h for h in v0 if h["clade"] != "ext"] or v0
            pan = [h for h in v if h["clade"] == "pan"]
            cl = [h for h in v if h["clade"] not in ("pan", "ext")]
            deep = max(cl, key=lambda h: (h["depth"], h["bits"])) if cl else None
            best = min(v, key=lambda h: (h["evalue"], h["model"] != "ARWEN"))
            others = [fam_best[o]["bits"] for o in order if o != fa]
            acs = Counter(h["anticodon"] for h in v if h.get("anticodon"))
            cands.append({
                "locus": li, "family": fa, "strand": best["strand"], "start": best["start"], "end": best["end"],
                "anticodon": acs.most_common(1)[0][0] if acs else "",
                "pan_bits": max((h["bits"] for h in pan), default=-20),
                "pan_logE": min((lg(h["evalue"]) for h in pan), default=2),
                "clade_bits": deep["bits"] if deep else -20, "clade_logE": lg(deep["evalue"]) if deep else 2,
                "best_logE": lg(best["evalue"]), "depth": deep["depth"] if deep else 0, "n_models": len(v),
                "n_clade_models": len(cl), "len": best["end"] - best["start"] + 1,
                "cov": (best["mto"] - best["mfrom"] + 1) / max(best["clen"], 1), "trunc": int(best["trunc"] != "no"),
                "delta_bits_next_fam": fam_best[fa]["bits"] - (max(others) if others else 0),
                "best_fam_at_locus": int(order[0] == fa), "is_LS": int(fa[0] in "LS"), "gc": best["gc"],
                "ac_known": best["ac_known"], "ac_match": best["ac_match"],
                "arwen": int("ARWEN" in ex), "aragorn": int("ARAGORN" in ex), "only_ext": int(best["clade"] == "ext"),
                "_bits": fam_best[fa]["bits"], "_opts": v0})
    fam_c = defaultdict(list)
    for c in cands:
        fam_c[c["family"]].append(c)
    for v in fam_c.values():
        v.sort(key=lambda c: -c["_bits"])
        for r, c in enumerate(v):
            c["fam_rank"], c["fam_n_cands"] = r, len(v)
            c["fam_gap_to_best"] = v[0]["_bits"] - c["_bits"]
            c["rel_to_fam_best"] = c["_bits"] / max(v[0]["_bits"], 1)
    strong = [c for c in cands if c["best_logE"] < -5 and c["best_fam_at_locus"]]
    n_strong = len({c["family"] for c in strong})
    for c in cands:
        c["n_strong_genome"] = n_strong
        c["overlap_stronger"] = max([_ov(c, s) for s in strong if s is not c and s["strand"] == c["strand"]
                                     and s["_bits"] > c["_bits"] and s["locus"] != c["locus"]] or [0])
        opp = [s for s in strong if s["strand"] != c["strand"] and _ov(c, s) > 0.5 * c["len"]]
        c["overlap_stronger_opp"] = max([_ov(c, s) for s in opp] or [0])
        c["opp_bits_ratio"] = max([s["_bits"] for s in opp] or [0]) / max(c["_bits"], 1)
    return cands


def call(cands, clf, thr: float = THRESHOLD):
    """Best candidate per locus with p >= thr; same-strand calls overlapping > 50 % keep the more probable."""
    if not cands:
        return []
    for c, q in zip(cands, clf.predict_proba([[float(c[k]) for k in FEATS] for c in cands])[:, 1]):
        c["p"] = float(q)
    best: Dict[int, Dict[str, Any]] = {}
    for c in cands:
        if c["p"] >= thr and (c["locus"] not in best or c["p"] > best[c["locus"]]["p"]):
            best[c["locus"]] = c
    final: List[Dict[str, Any]] = []
    for c in sorted(best.values(), key=lambda c: -c["p"]):
        if any(x["strand"] == c["strand"] and _ov(x, c) > 0.5 * min(x["end"] - x["start"], c["end"] - c["start"])
               for x in final):
            continue
        final.append(c)
    return final


# ── 3. ends ──────────────────────────────────────────────────────────────────
def _stem(a: str, b: str) -> Tuple[int, int]:
    p = w = 0
    for x, y in zip(a, reversed(b)):
        p += (x, y) in PAIR
        w += (x, y) in WC
    return p, w


def end_pairs(c, gen: str) -> List[Dict[str, Any]]:
    """Candidate (start, end) pairs of one call with the BFEATS features (as projects/invert_trna bounds.py)."""
    strand, L = c["strand"], len(gen)
    g2 = gen + gen
    opts = [{"m": h["model"], "clade": h["clade"], "s": h["start"], "e": h["end"], "bits": h["bits"],
             "full5": h["mfrom"] == 1, "full3": h["mto"] == h["clen"]} for h in c["_opts"]]
    S, E = defaultdict(list), defaultdict(list)
    for o in opts:
        S[o["s"]].append(o)
        E[o["e"]].append(o)
    med_s, med_e = statistics.median(o["s"] for o in opts), statistics.median(o["e"] for o in opts)
    cs = {s + d: (s, d) for s in S for d in (-1, 0, 1)}
    ce = {e + d: (e, d) for e in E for d in (-1, 0, 1)}
    for s in list(cs):
        if s in S:
            cs[s] = (s, 0)
    for e in list(ce):
        if e in E:
            ce[e] = (e, 0)
    out = []
    for s, (s0, ds) in cs.items():
        for e, (e0, de) in ce.items():
            n = e - s + 1
            if not 40 <= n <= 110:
                continue
            sub = g2[s - 1:e] if e <= 2 * L and s >= 1 else ""
            seq = sub if strand == 1 else rc(sub)
            if len(seq) != n:
                continue
            accs = [_stem(seq[:7], seq[-7 - d:len(seq) - d])[0] for d in (0, 1, 2)]
            so, eo = S[s0], E[e0]
            f5, f3 = ("full5", "full3") if strand == 1 else ("full3", "full5")
            out.append({"s": s, "e": e,
                        "s_n": len(so), "s_bits": max(o["bits"] for o in so),
                        "s_pan": int(any(o["clade"] == "pan" for o in so)),
                        "s_clade": int(any(o["clade"] not in ("pan", "ext") for o in so)),
                        "s_arwen": int(any(o["m"] == "ARWEN" for o in so)),
                        "s_aragorn": int(any(o["m"] == "ARAGORN" for o in so)),
                        "s_full": int(any(o[f5] for o in so if o["clade"] != "ext")), "s_shift": ds, "s_off": s - med_s,
                        "e_n": len(eo), "e_bits": max(o["bits"] for o in eo),
                        "e_pan": int(any(o["clade"] == "pan" for o in eo)),
                        "e_clade": int(any(o["clade"] not in ("pan", "ext") for o in eo)),
                        "e_arwen": int(any(o["m"] == "ARWEN" for o in eo)),
                        "e_aragorn": int(any(o["m"] == "ARAGORN" for o in eo)),
                        "e_full": int(any(o[f3] for o in eo if o["clade"] != "ext")), "e_shift": de, "e_off": e - med_e,
                        "len": n, "acc0": accs[0], "acc1": accs[1], "acc2": accs[2],
                        "wc1": _stem(seq[:7], seq[-8:-1])[1], "first_nt_T": int(seq[0] == "T"),
                        "last_nt_A": int(seq[-1] == "A")})
    return out


def refine_ends(calls, gen: str, bclf) -> None:
    for c in calls:
        pairs = end_pairs(c, gen)
        if pairs:
            p = bclf.predict_proba([[float(r[k]) for k in BFEATS] for r in pairs])[:, 1]
            i = max(range(len(pairs)), key=lambda k: p[k])
            c["start"], c["end"], c["bound_p"] = pairs[i]["s"], pairs[i]["e"], float(p[i])


# ── driver ───────────────────────────────────────────────────────────────────
def predict(seq: str, lineage: str, code: int, pan: Path, cmsearch: str, cmalign: str, ext_preds=(),
            models: Path = MODELS, cpu: int = 1, thr: float = THRESHOLD, use_clades: bool = True) -> List[Dict[str, Any]]:
    """tRNA calls of one (circular) sequence: dicts with family (A..Y, L1/L2, S1/S2), aa, anticodon, start, end
    (1-based; end > len(seq) for a tRNA across the origin), strand, p (classifier probability), bound_p."""
    from threadpoolctl import threadpool_limits      # installed with scikit-learn
    with threadpool_limits(limits=max(1, cpu)):     # its OpenMP pool otherwise takes every core: 16 parallel records
        clf, bclf = classifiers(models)              # ran 99 threads each (load 160 on 72 cores)
        cdir = clade_dir(models) if use_clades else None
        with tempfile.TemporaryDirectory() as td:
            hits = scan(seq, lineage, pan, cdir, cmsearch, cmalign, Path(td), cpu)
        calls = call(candidates(hits, ext_hits(ext_preds), code), clf, thr)
        refine_ends(calls, seq.upper(), bclf)
    return [{"family": c["family"], "aa": c["family"][0], "anticodon": c["anticodon"], "start": c["start"],
             "end": c["end"], "strand": c["strand"], "p": round(c["p"], 4), "bound_p": round(c.get("bound_p", 0), 4)}
            for c in sorted(calls, key=lambda c: c["start"])]


def main(argv=None) -> int:
    sys.path.insert(0, str(HERE))
    import mito_env
    from Bio import SeqIO
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", required=True, help="FASTA or GenBank (lineage then read from the record)")
    ap.add_argument("--lineage", default="", help="NCBI lineage, ';'-separated (GenBank input: from the record)")
    ap.add_argument("--genetic_code", type=int, default=5)
    ap.add_argument("--mitos_refdir", default=None)
    ap.add_argument("--cpu", type=int, default=1)
    a = ap.parse_args(argv)
    ref = mito_env.mitos_refdir(a.mitos_refdir)
    cms, cma = mito_env.find("cmsearch"), mito_env.find("cmalign")
    if not ref or not cms or not cma:
        print("[trna_ml] needs MITOS2 refseq89m (--mitos_refdir), cmsearch and cmalign", file=sys.stderr)
        return 2
    fmt = "fasta" if Path(a.input).suffix.lower() in (".fa", ".fasta", ".fna") else "genbank"
    w = csv.writer(sys.stdout, delimiter="\t")
    w.writerow(["record", "family", "anticodon", "start", "end", "strand", "p", "bound_p"])
    for rec in SeqIO.parse(a.input, fmt):
        lin = a.lineage or ";".join(rec.annotations.get("taxonomy", []))
        for c in predict(str(rec.seq), lin, a.genetic_code, Path(ref) / "ncRNA", cms, cma, cpu=a.cpu):
            w.writerow([rec.id, c["family"], c["anticodon"], c["start"], c["end"], c["strand"], c["p"], c["bound_p"]])
    return 0


if __name__ == "__main__":
    sys.exit(main())
