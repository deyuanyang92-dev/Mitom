"""Classify every NCBI intron from the evidence of verify.py and choose the gold placement (criteria: gold/README.md).

usage: classify.py [--evidence work/evidence.pkl]  -> gold/introns.tsv, gold/summary.tsv, gold/work/gold.json
"""
import argparse
import csv
import json
import pickle
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, "/media/deyuan/18T_B/yzz/MitoM/benchmark")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from mito_names import gene_name  # noqa: E402
from Bio import SeqIO  # noqa: E402

GOLD = Path(__file__).resolve().parents[1]
CM_SURE = 1e-5          # a covariance-model hit at least this good fixes the intron type
MAJ = 2 / 3             # "better than NCBI" / "better than the runner-up": on >= 2/3 of the homologs
SUBST = 0.8             # a correction that only changes residues (NCBI junction clean): better on >= 80 % of homologs
STRONG_WINS, STRONG_GAIN = 11, 20.0
# class "consensus": NCBI's product and junction are clean and the change rests on a canonical / clade-consensus
# placement that aligns better; the gold moves only on strong evidence: better for >= STRONG_WINS of 12 homologs and
# either a gain >= STRONG_GAIN or a fully canonical placement of a confidently typed intron (CM E <= CM_SURE) where
# NCBI's placement is not canonical. NCBI's junction region counts as wrong (D-corrected) when one side aligns at
# < 0.3 identity and the correction brings it to >= 0.5.


def ncbi_region_wrong(junc, alt_junc) -> bool:
    return min(junc["id_before"], junc["id_after"]) < 0.3 <= 0.5 <= min(alt_junc["id_before"], alt_junc["id_after"])


def strong(a, m: str, typed: bool, ncbi_canonical: bool) -> bool:
    gain = a["score"] - a["base_score"]
    return a["wins"] >= STRONG_WINS and (gain >= STRONG_GAIN or (m == "full" and typed and not ncbi_canonical))


def coords(s: int, ln: int, n: int) -> str:
    e = s + ln
    return f"{s + 1}..{e if e <= n else e - n}" + (" (across origin)" if e > n else "")


def rules(r) -> tuple:
    """Allowed canonical rules (set), type text, type note."""
    kind, ev, model = r.get("cm") or ("", None, "")
    nt = r["note_type"]
    if kind and ev is not None and ev <= CM_SURE:
        flag = f"; NCBI note says group {nt}" if nt and nt != kind else ""
        return {kind}, f"group {kind} ({model}, E={ev:.1g}){flag}", ""
    if kind:
        allowed = {kind, nt} if nt else {"I", "II"}
        return allowed, f"group {kind}? (weak: {model}, E={ev:.1g})" + (f"; NCBI note group {nt}" if nt else ""), "weak type"
    if nt:
        return {nt}, f"group {nt} (NCBI note; no CM hit)", "no CM hit"
    return {"I", "II"}, "untyped (no CM hit, no NCBI note)", "untyped"


def meets(kinds, allowed) -> str:
    """'full' when a placement meets the rule of an allowed type, '5p' for the group II 5' motif alone, else ''."""
    if any(k in allowed for k in kinds if k in ("I", "II")):
        return "full"
    if "II5" in kinds and "II" in allowed:
        return "5p"
    return ""


def junction_text(j) -> str:
    return f"{j['indel']}/{j['n']} homologs with indel at +-5 aa (net {j['net']:+g}); identity {j['id_before']}/{j['id_after']}"


def clean(j) -> bool:
    return j["indel"] < 0.5 * j["n"]


def decide(r) -> dict:
    """Verdict A / B / C / D / D-corrected with the gold placement (genome start, length) and the reason."""
    n, out = r["n"], {}
    allowed, type_text, type_flag = rules(r)
    out["type"], out["type_flag"] = type_text, type_flag
    s, ln = r["gap_s"], r["gap_len"]
    ncbi = (s, ln)
    flags = []
    if (r["ncbi_s"] % n, r["ncbi_e"] - r["ncbi_s"]) != (s, ln):
        flags.append(f"NCBI intron feature {coords(r['ncbi_s'] % n, r['ncbi_e'] - r['ncbi_s'], n)} differs from the "
                     f"CDS join gap {coords(s, ln, n)}; the join is used")
    prod, junc = r["product"], r["junction"]
    if prod["start"].startswith("bad") or prod["stop"].startswith("none"):
        flags.append(f"CDS end problem not caused by the intron: start {prod['start']}, stop {prod['stop']}")
    if prod["translation"] not in ("same", "no /translation"):
        flags.append(f"/translation {prod['translation']}")
    ncbi_meets = meets(r["sig_kinds"], allowed)
    sign = r["strand"]
    win_can = [(g, meets(k, allowed)) for g, *_, k in r["window_sig"] if meets(k, allowed)]
    alts = r["alts"]
    nh = junc["n"]
    syn = [(a, p, meets(p["kinds"], allowed)) for a in alts if a["same_protein"] for p in a["placements"]
           if meets(p["kinds"], allowed)]
    better = []
    need = (MAJ if (prod["stops"] or not clean(junc)) else SUBST) * nh
    for a in alts:
        if a["same_protein"] or not clean(a["junction"]) or a["wins"] < need or a["score"] <= r["score"]:
            continue
        pl = [(p, meets(p["kinds"], allowed)) for p in a["placements"] if meets(p["kinds"], allowed)]
        if pl:
            pl.sort(key=lambda x: (x[1] != "full", abs(x[0]["d5"]) + abs(x[0]["d3"])))
            better.append((a, pl[0][0], pl[0][1]))
    better.sort(key=lambda x: (x[2] != "full", -x[0]["score"]))
    unique = None
    if better:
        top = better[0][0]
        rest = [b for b in better[1:] if b[2] == better[0][2]]
        if not rest or sum(x > y for x, y in zip(top["scores"], rest[0][0]["scores"])) >= MAJ * nh:
            unique = better[0]
    out["ncbi_signals"] = f"{r['sig_last']}|{r['sig_5']}...{r['sig_3']} ({','.join(r['sig_kinds']) or 'none'})"
    out["window"] = f"{r['window'][0]:+d}..{r['window'][1]:+d}" + \
        (f"; canonical at {','.join(f'{g:+d}({m})' for g, m in win_can)}" if win_can else "; no canonical placement")

    def alt_text(a, p, m):
        return (f"5' {p['d5']:+d} / 3' {p['d3']:+d} nt ({','.join(p['kinds'])}{', 5p only' if m == '5p' else ''}): "
                f"{a['aa_diff']}, score {a['score'] - r['score']:+.1f}, better for {a['wins']}/{nh}, "
                f"junction {junction_text(a['junction'])}")

    def moved(p):                       # gold placement of an alternative (genome start, length)
        return p["start"], p["len"]

    def substitution(u, what: str) -> dict:
        """A correction of a clean NCBI product: D-corrected when NCBI's junction region does not align, else class
        consensus (gold moved only on strong evidence)."""
        a, p, m = u
        a = {**a, "base_score": r["score"]}
        if ncbi_region_wrong(junc, a["junction"]):
            return {"verdict": "D-corrected", "gold": moved(p),
                    "reason": f"NCBI junction region does not align ({junction_text(junc)}); {what}: {alt_text(a, p, m)}"}
        if strong(a, m, type_flag == "", ncbi_meets == "full"):
            return {"verdict": "consensus", "gold": moved(p),
                    "reason": f"NCBI product and junction clean but NCBI placement not canonical; {what} (strong: gold "
                              f"moved): {alt_text(a, p, m)}"}
        return {"verdict": "consensus", "gold": ncbi,
                "reason": f"NCBI product and junction clean; {what} but the evidence is not strong, so the gold keeps "
                          f"NCBI's placement; alternative: {alt_text(a, p, m)}"}

    stops = prod["stops"]
    if stops or not clean(junc):
        why = f"{stops} internal stop(s) in the spliced CDS" if stops else f"junction indel vs homologs: {junction_text(junc)}"
        if unique:
            a, p, m = unique
            out.update(verdict="D-corrected", gold=moved(p), reason=f"{why}; corrected: {alt_text(a, p, m)}")
        elif better:
            out.update(verdict="D", gold=None, reason=f"{why}; {len(better)} canonical placements remove it, none "
                                                       f"clearly best: {alt_text(*better[0])}")
        elif stops:
            out.update(verdict="D", gold=None, reason=f"{why}; no canonical placement within 60 nt gives a clean CDS")
        elif min(junc["id_before"], junc["id_after"]) < 0.6 <= 0.9 <= max(junc["id_before"], junc["id_after"]):
            out.update(verdict="D", gold=None, reason=f"{why}; one side of the junction stops aligning and no canonical "
                                                       "placement within 60 nt removes the indel")
        else:
            flags.append(f"{why}; no canonical placement within 60 nt removes it (lineage-specific indel?)")
            stops = None
    if "verdict" not in out:
        if ncbi_meets == "full":
            out.update(verdict="A", gold=ncbi, reason="canonical NCBI placement, product and junction clean")
        elif any(m == "full" for _, m in win_can):
            g = min((g for g, m in win_can if m == "full"), key=abs)
            out.update(verdict="B", gold=((s + g) % n, ln),
                       reason=f"same spliced sequence; canonical placement {g * sign:+d} nt (transcript) from NCBI's")
        elif any(m == "full" for *_, m in syn):
            a, p, m = min(((a, p, m) for a, p, m in syn if m == "full"), key=lambda x: abs(x[1]["d5"]) + abs(x[1]["d3"]))
            out.update(verdict="B", gold=moved(p),
                       reason=f"same protein (synonymous difference), canonical placement 5' {p['d5']:+d} / 3' {p['d3']:+d} nt")
        elif unique and unique[2] == "full":
            out.update(substitution(unique, "the canonical placement aligns better"))
        elif ncbi_meets == "5p":
            out.update(verdict="A", gold=ncbi, reason="NCBI placement has the group II 5' motif (3' end not AY)")
            flags.append("group II 3' end not AY")
        elif any(m == "5p" for _, m in win_can):
            g = min((g for g, m in win_can if m == "5p"), key=abs)
            out.update(verdict="B", gold=((s + g) % n, ln),
                       reason=f"same spliced sequence; group II 5' motif {g * sign:+d} nt (transcript) from NCBI's (3' not AY)")
        elif any(m == "5p" for *_, m in syn):
            a, p, m = min(((a, p, m) for a, p, m in syn if m == "5p"), key=lambda x: abs(x[1]["d5"]) + abs(x[1]["d3"]))
            out.update(verdict="B", gold=moved(p),
                       reason=f"same protein (synonymous), group II 5' motif at 5' {p['d5']:+d} / 3' {p['d3']:+d} nt (3' not AY)")
        elif unique:
            out.update(substitution(unique, "the group II 5' motif placement aligns better"))
        else:
            out.update(verdict="C", gold=ncbi, reason="product and junction clean; no canonical placement in the window")
            if better:
                flags.append(f"better-scoring canonical placements exist but none clearly best: {alt_text(*better[0])}")
    if out["verdict"] in ("A", "B", "C") and unique and out.get("gold") != moved(unique[1]):
        flags.append(f"note: a canonical placement aligns better: {alt_text(*unique)}")
    out["best_alt"] = alt_text(*better[0]) if better else ""
    out["flags"] = flags
    return out


def placement_error(r, cds, n: int):
    """An intron feature that matches no join gap of its gene: ('keep', why) for the half of a trans-splicing intron
    next to a piece of a two-strand CDS, ('drop', why) when it lies inside an exon or covers only part of a join gap
    (the join gap is the intron and is evaluated on its own row)."""
    s, e = r["ncbi_s"], r["ncbi_e"]
    for f in cds:
        parts = list(f.location.parts)
        if any(E_overlap((int(p.start), int(p.end)), (s, e), n) for p in parts):
            return "drop", (f"NCBI error: the intron feature lies inside exon {f.location} of the CDS; the CDS join gap "
                            "is the intron (own row)")
        for p_, q in zip(parts, parts[1:]):
            if p_.strand == q.strand:
                gs, ge = (int(p_.end), int(q.start)) if (p_.strand or 1) == 1 else (int(q.end), int(p_.start))
                if E_overlap((gs, ge), (s, e), n):
                    return "drop", (f"NCBI error: the intron feature covers only {e - s} nt of the CDS join gap "
                                    f"{gs + 1}..{ge}, which is the intron (own row)")
        if len({p.strand for p in parts}) > 1:
            ends = [int(p.end) for p in parts if p.strand == r["strand"]] + [int(p.start) for p in parts if p.strand == r["strand"]]
            if any(min(abs(x - s), abs(x - e)) <= 3 for x in ends):
                return "keep", (f"half of a trans-splicing intron next to a piece of the two-strand CDS {f.location}: "
                                f"CM type {r['cm'][0] or '-'} ({r['cm'][2]}, E={r['cm'][1]}); not checkable: product, "
                                "junction, the partner half")
    return "keep", f"intron feature without a join gap (contained in: {r.get('containing', '')}); not checkable"


def E_overlap(a, b, n):
    return max(max(0, min(a[1], b[1] + k) - max(a[0], b[0] + k)) for k in (-n, 0, n))


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--evidence", default=str(GOLD / "work" / "evidence.pkl"))
    p.add_argument("--out", default=str(GOLD))
    a = p.parse_args()
    rows = pickle.load(open(a.evidence, "rb"))         # written by verify.py (this folder)
    out_rows, gold = [], []
    recs = {}
    for r in rows:
        n = r["n"]
        if r["accession"] not in recs:
            recs[r["accession"]] = SeqIO.read(str(GOLD.parent / "ncbi_gb" / f"{r['accession']}.gb"), "genbank")
        rec = recs[r["accession"]]
        base = {"accession": r["accession"], "phylum": r["phylum"], "order": r["order"], "organism": r["organism"],
                "gene": r["gene"], "strand": "+" if r["strand"] == 1 else "-",
                "ncbi": coords(r["ncbi_s"] % n, r["ncbi_e"] - r["ncbi_s"], n), "ncbi_source": r["source"]}
        kind, ev, model = r.get("cm") or ("", None, "")
        cds = [f for f in rec.features if f.type == "CDS" and gene_name(f.qualifiers, "CDS") == r["gene"]]
        split = len(cds) > 1 or any(len({x.strand for x in f.location.parts}) > 1 for f in cds)
        if r["host_kind"] == "none":
            handling, why = placement_error(r, cds, n)
            row = {**base, "gold": "removed" if handling == "drop" else base["ncbi"],
                   "class": "D" if handling == "drop" else "E", "verdict": "",
                   "type": f"{kind or '-'} {model} E={ev:.1g}" if ev is not None else "no CM hit",
                   "product": "", "junction": "", "ncbi_signals": "", "window": "", "best_alt": "", "clade": "",
                   "homologs": "", "reason": why, "gold_handling": handling}
            out_rows.append(row)
            gold.append({**row, "gold_s": r["ncbi_s"] % n if handling != "drop" else None,
                         "gold_len": r["ncbi_e"] - r["ncbi_s"], "strand_i": r["strand"], "n": n})
            continue
        if r["host_kind"] == "rRNA":
            rj = r["rrna_junction"]
            sig = f"{r['sig_last']}|{r['sig_5']}...{r['sig_3']} ({','.join(r['sig_kinds']) or 'none'})"
            row = {**base, "gold": base["ncbi"], "class": "E", "verdict": "",
                   "type": f"group {kind} ({model}, E={ev:.1g})" if kind else "no CM hit",
                   "product": "", "junction": f"rRNA spliced vs {rj['n']} homolog rRNAs: {rj['indel']} with an indel within +-15 nt",
                   "ncbi_signals": sig, "window": f"{r['window'][0]:+d}..{r['window'][1]:+d}", "best_alt": "", "clade": "",
                   "homologs": f"{len(r['homologs'])} rRNAs",
                   "reason": "intron of an rRNA gene: covariance-model type and splice signals checked, rRNA homology "
                             "across the junction checked; not checkable: protein product / junction", "gold_handling": "keep"}
            out_rows.append(row)
            gold.append({**row, "gold_s": r["gap_s"], "gold_len": r["gap_len"], "strand_i": r["strand"], "n": n})
            continue
        e_reason = ("fragment (LC670763.1, 1.5-kb cox1 piece): start / stop of the gene not checkable" if r["fragment"] else
                    "intron of a trans-spliced (split) gene: checked as a cis intron; the trans-splicing itself and "
                    "the gene ends of the pieces are not checkable" if split else "")
        d = decide(r)
        pr, j = r["product"], r["junction"]
        cls = "E" if e_reason else d["verdict"]
        g = d["gold"]
        handling = "keep" if g else "ignore"
        row = {**base, "gold": coords(g[0], g[1], n) if g else "excluded", "class": cls,
               "verdict": d["verdict"] if e_reason else "", "type": d["type"],
               "product": f"{pr['aa']} aa, stops {pr['stops']}, start {pr['start']}, stop {pr['stop']}, /translation {pr['translation']}",
               "junction": junction_text(j), "ncbi_signals": d["ncbi_signals"], "window": d["window"],
               "best_alt": d["best_alt"], "clade": "",
               "homologs": f"{len(r['homologs'])} ({sum(not h for _, h in r['homologs'])} without intron): "
                           + ",".join(a_ for a_, _ in r["homologs"][:4]) + ("..." if len(r["homologs"]) > 4 else ""),
               "reason": "; ".join(x for x in [e_reason, d["reason"]] + d["flags"] if x), "gold_handling": handling}
        out_rows.append(row)
        gold.append({**row, "gold_s": g[0] if g else None, "gold_len": g[1] if g else None, "exons": r["exons"],
                     "k": r["k"], "strand_i": r["strand"], "n": n, "gap_s": r["gap_s"], "gap_len": r["gap_len"],
                     "ref_pos": r.get("ref_pos"), "phase": r["phase"],
                     "gold_ref_pos": next((a.get("ref_pos") for a in r["alts"] for p_ in a["placements"]
                                           if g and (p_["start"], p_["len"]) == tuple(g)), r.get("ref_pos"))})
    clade(out_rows, gold, rows)
    for row, g in zip(out_rows, gold):          # the clade pass changes rows: keep gold.json in step
        g.update({k: row[k] for k in ("class", "verdict", "reason", "gold", "gold_handling", "type", "clade")})
        if row["gold"] not in ("excluded", "removed"):
            assert row["gold"] == coords(g["gold_s"], g["gold_len"], g["n"]), (row["accession"], row["gold"])
    out = Path(a.out)
    with open(out / "introns.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, list(out_rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    summ = Counter((r["phylum"], r["class"]) for r in out_rows)
    classes = ["A", "B", "C", "consensus", "D-corrected", "D", "E"]
    with open(out / "summary.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["phylum"] + classes + ["total", "in_gold", "excluded_ignore", "removed_not_intron",
                                          "gold_moved_from_ncbi"])
        for ph in sorted({p_ for p_, _ in summ}) + ["all"]:
            sel = [r for r in out_rows if ph == "all" or r["phylum"] == ph]
            c = [sum(1 for r in sel if r["class"] == cl) for cl in classes]
            excl = sum(1 for r in sel if r["gold"] == "excluded")
            rem = sum(1 for r in sel if r["gold"] == "removed")
            moved = sum(1 for r in sel if r["gold"] not in ("excluded", "removed") and r["gold"] != r["ncbi"])
            w.writerow([ph] + c + [len(sel), len(sel) - excl - rem, excl, rem, moved])
    json.dump(gold, open(out / "work" / "gold.json", "w"), indent=0, default=str)
    print(open(out / "summary.tsv").read())
    return 0


def clade(rows, gold, evid) -> None:
    """Clade consistency (criterion 4): the junction of each intron as (reference residue, phase) in the Metridium
    senile protein; introns of one gene and phylum whose gold junctions (after the corrections above) lie within 3 nt
    of each other (reference CDS) form a site (single linkage); the site's consensus is the most
    common (residue, phase) of its members with a clean junction (>= 3 members, >= half of them). A C intron whose
    NCBI placement differs from the consensus moves to the same-product placement that matches it (class B, clade);
    an intron with an unresolved junction indel moves to the unique placement that matches it, gives a clean
    junction and aligns better (D-corrected, clade); other differences are flagged."""
    import verify as V
    V.HOM = pickle.load(open(GOLD / "work" / "homologs.pkl", "rb"))      # written by homologs.py (this folder)
    ref = {g: p for x in V.HOM if x["acc"] == V.REF_ACC for g, (p, _) in x["prot"].items()}
    info = []
    for row, g, r in zip(rows, gold, evid):
        if r.get("host_kind") != "CDS" or r["gene"] not in ref:
            info.append(None)
            continue
        rec = SeqIO.read(str(GOLD.parent / "ncbi_gb" / f"{r['accession']}.gb"), "genbank")
        t = V.E.Intron(r["ncbi_s"], r["ncbi_e"], r["strand"], r["gene"], r["source"], r["ncbi_note"])
        f, k = V.host_of(rec, t)
        m = V.model_of(rec, f, r["code"])
        prot = m.protein().rstrip("*")
        aln = V.AL.align(prot, ref[r["gene"]])[0]
        pmap = {}
        for (p0, p1), (h0, _) in zip(*aln.aligned):
            for i in range(p1 - p0):
                pmap[p0 + i] = h0 + i

        def site(jnt, pm=pmap, n_=len(prot)):
            j = min(max(int(jnt // 3), 0), n_ - 1)
            near = min(pm, key=lambda x: abs(x - j)) if pm else j
            return pm.get(j, pm.get(near, -1)), jnt % 3

        def placement_site(gs, gl):
            """(ref residue, phase) of a placement (genome start, length) that keeps the product."""
            for d in range(-61, 62):
                ex = m.moved(k, d, d)
                if ex is not None and m.intron(k, ex)[:2] == (gs, gl):
                    return site(m.junction_nt(k, ex))
            return None
        cur = (g.get("gold_s"), g.get("gold_len"))
        here = placement_site(*cur) if cur[0] is not None else None
        if here is None and cur[0] is not None:                     # a correction changed the product
            for d5 in range(-61, 62):
                for d3 in range(-61, 62):
                    ex = m.moved(k, d5, d3)
                    if ex is not None and m.intron(k, ex)[:2] == cur:
                        p2 = m.protein(ex).rstrip("*")
                        al2 = V.AL.align(p2, ref[r["gene"]])[0]
                        pm2 = {p0 + i: h0 + i for (p0, p1), (h0, _) in zip(*al2.aligned) for i in range(p1 - p0)}
                        here = site(m.junction_nt(k, ex), pm2, len(p2))
                        break
                if here is not None:
                    break
        pat = None
        if cur[0] is not None:
            last, _, i3 = V.signals(m, cur[0], cur[1], r["strand"])
            pat = (last, i3[-1:])                                     # exon-1 last base, intron last base
        info.append({"m": m, "k": k, "prot": prot, "site": site, "here": here, "ncbi": site(m.junction_nt(k)),
                     "clean": clean(r["junction"]) and not r["product"]["stops"], "r": r, "pat": pat})
    by = defaultdict(list)
    for i, x in enumerate(info):
        if x and x["ncbi"][0] >= 0:
            x["key"] = x["here"] if x["here"] and x["here"][0] >= 0 else x["ncbi"]   # gold site (after corrections)
            by[(rows[i]["gene"], rows[i]["phylum"])].append(i)
    nt = lambda st: 3 * st[0] + st[1]                         # junction as a nucleotide offset in the reference CDS
    for (gene, phylum), idx in by.items():
        idx.sort(key=lambda i: nt(info[i]["key"]))
        clusters, cur = [], [idx[0]]
        for i in idx[1:]:
            if nt(info[i]["key"]) - nt(info[cur[-1]]["key"]) <= 3:
                cur.append(i)
            else:
                clusters.append(cur)
                cur = [i]
        clusters.append(cur)
        sites = []
        for cl in clusters:
            voters = [i for i in cl if info[i]["clean"] and info[i]["here"] and rows[i]["class"] != "D"
                      and rows[i]["gold"] not in ("excluded", "removed")]
            votes = Counter(info[i]["here"] for i in voters)
            cons, nv = votes.most_common(1)[0] if votes else (None, 0)
            ok = cons is not None and nv >= 3 and nv >= 0.5 * sum(votes.values())
            # family rule: the terminal bases of the consensus members when most of them break the generic rule
            at_cons = [info[i]["pat"] for i in voters if info[i]["here"] == cons]
            generic = sum(1 for la, w in at_cons if (la == "T" and w == "G")) if at_cons else 0
            fam = Counter(at_cons).most_common(1)[0][0] if at_cons and generic < 0.5 * len(at_cons) else None
            sites.append({"cl": cl, "cons": cons, "nv": nv, "nvotes": sum(votes.values()), "ok": ok, "fam": fam,
                          "generic": generic, "n_cons": len(at_cons), "near": ""})
        for st in sites:                    # a lone outlier (misplaced by NCBI) joins the nearest site within 20 aa
            if st["ok"]:
                continue
            best = min((o for o in sites if o["ok"]), default=None,
                       key=lambda o: min(abs(nt(info[i]["key"]) - nt(o["cons"])) for i in st["cl"]))
            if best and min(abs(nt(info[i]["key"]) - nt(best["cons"])) for i in st["cl"]) <= 60:
                st.update({k_: best[k_] for k_ in ("cons", "nv", "nvotes", "ok", "fam", "generic", "n_cons")})
                st["near"] = f"; assigned to the nearest site (consensus aa {best['cons'][0] + 1})"
        for st in sites:
            cl, cons, nv, ok, fam = st["cl"], st["cons"], st["nv"], st["ok"], st["fam"]
            orders = Counter(rows[i]["order"] for i in cl)
            at_cons, generic = [None] * st["n_cons"], st["generic"]
            for i in cl:
                x, row = info[i], rows[i]
                desc = (f"{gene} site at reference aa {x['ncbi'][0] + 1}, phase {x['ncbi'][1]}; site of {len(cl)} "
                        f"{phylum} introns (" + ", ".join(f"{o} {c}" for o, c in orders.most_common(3)) + ")")
                if not ok:
                    row["clade"] = desc + "; no consensus (too few clean members)"
                    continue
                desc += f"{st['near']}; consensus aa {cons[0] + 1} phase {cons[1]} ({nv}/{st['nvotes']})"
                if fam:
                    desc += f"; family rule exon ...{fam[0]} | intron ...{fam[1]} ({len(at_cons) - generic}/{len(at_cons)} " \
                            "consensus members lack exon-U + omega-G)"
                verdict = row["verdict"] or row["class"]
                unresolved = "no canonical placement within 60 nt removes it" in row["reason"] or row["gold"] == "excluded"
                moved = ""
                if unresolved:
                    moved = clade_move(x, cons, row, gold[i], V, "fix")
                elif x["here"] != cons and verdict == "C":
                    moved = clade_move(x, cons, row, gold[i], V, "same") or clade_move(x, cons, row, gold[i], V, "fix")
                elif x["here"] != cons and fam:            # the family lacks omega-G: its consensus is the rule
                    moved = clade_move(x, cons, row, gold[i], V, "same", override=True) or \
                        clade_move(x, cons, row, gold[i], V, "fix", override=True)
                if moved:
                    row["clade"] = desc + f"; {moved}"
                elif x["here"] == cons:
                    row["clade"] = desc + "; gold placement agrees"
                    if fam and verdict == "C" and x["pat"] and x["pat"][0] == fam[0]:
                        set_verdict(row, "A")
                        row["reason"] = row["reason"].replace(
                            "product and junction clean; no canonical placement in the window",
                            f"product and junction clean; NCBI placement is the clade consensus: exon ...{x['pat'][0]} | "
                            f"intron ...{x['pat'][1]} (family rule: exon-U at the consensus junction; the family lacks omega-G)")
                else:
                    row["clade"] = desc + f"; gold placement differs (aa {x['here'][0] + 1 if x['here'] else '?'} " \
                                          f"phase {x['here'][1] if x['here'] else '?'})"


def set_verdict(row, v: str) -> None:
    if row["class"] == "E":
        row["verdict"] = v
    else:
        row["class"] = v


def clade_move(x, cons, row, g, V, mode: str, override: bool = False) -> str:
    """Move a C or D intron to the clade-consistent placement; returns the description or ''.
    mode 'same': a placement with the same product (identical spliced sequence or protein) -> B;
    mode 'fix': a placement with another product that gives a clean junction and aligns better -> D-corrected."""
    m, k, r = x["m"], x["k"], x["r"]
    n = r["n"]
    unresolved = "no canonical placement within 60 nt removes it" in row["reason"] or row["gold"] == "excluded"
    if mode == "same":
        for d in sorted(range(-60, 61), key=abs):            # same product: identical spliced sequence or protein
            ex = m.moved(k, d, d)
            if ex is None or d == 0:
                continue
            if m.protein(ex).rstrip("*") != x["prot"] or x["site"](m.junction_nt(k, ex)) != cons:
                continue
            gs, gl, _ = m.intron(k, ex)
            same_nt = m.spliced(ex) == m.spliced()
            g.update(gold_s=gs, gold_len=gl)
            row["gold"] = coords(gs, gl, n)
            why = (f"same {'spliced sequence' if same_nt else 'protein (synonymous)'}; placement {d:+d} nt (transcript) "
                   + ("matches the clade consensus (the family lacks omega-G; the generic canonical placement is overridden)"
                      if override else "matches the clade consensus (no canonical placement)"))
            set_verdict(row, "B")
            if override:
                row["reason"] = why + "; was: " + row["reason"]
            else:
                row["reason"] = row["reason"].replace("product and junction clean; no canonical placement in the window", why)
            return f"moved to the consensus: {why}"
        return ""
    homs = V.homologs(r["gene"], r["organism"], V_lineage(V, r))       # the same homologs as verify.py
    base = V.scores(x["prot"], homs)
    need = (MAJ if unresolved else SUBST) * len(homs)
    found = {}
    for d5 in range(-60, 61):
        for d3 in range(-60, 61):
            if (d5 - d3) % 3:
                continue
            ex = m.moved(k, d5, d3)
            if ex is None or m.intron(k, ex)[1] < 50:
                continue
            p2 = m.protein(ex)
            if "*" in p2.rstrip("*"):
                continue
            q = m.junction_nt(k, ex)
            if x["site"](q)[1] != cons[1]:
                continue
            p2 = p2.rstrip("*")
            al2 = V.AL.align(p2, V_ref(V, r["gene"]))[0]
            pm2 = {p0 + i: h0 + i for (p0, p1), (h0, _) in zip(*al2.aligned) for i in range(p1 - p0)}
            j = min(int(q // 3), len(p2) - 1)
            if pm2.get(j, -1) != cons[0]:
                continue
            if p2 not in found:
                found[p2] = (d5, d3, ex, q)
    cands = []
    for p2, (d5, d3, ex, q) in found.items():
        jc = V.junction_check(p2, q / 3, homs)
        sc = V.scores(p2, homs)
        if clean(jc) and sum(a > b for a, b in zip(sc, base)) >= need:
            cands.append((sum(sc) / len(sc), d5, d3, ex, jc, sc))
    cands.sort(key=lambda c: -c[0])
    if not cands or (len(cands) > 1 and sum(a > b for a, b in zip(cands[0][5], cands[1][5])) < MAJ * len(homs)):
        return ""
    sc, d5, d3, ex, jc, scs = cands[0]
    gain = sc - sum(base) / len(base)
    wins = sum(a > b for a, b in zip(scs, base))
    ncbi_wrong = bool(r["product"]["stops"]) or not clean(r["junction"]) or ncbi_region_wrong(r["junction"], jc)
    gs, gl, _ = m.intron(k, ex)
    why = (f"clade consensus placement 5' {d5:+d} / 3' {d3:+d} nt: score {gain:+.1f}, better for {wins}/{len(homs)}, "
           f"junction {junction_text(jc)}")
    if not ncbi_wrong and not (wins >= STRONG_WINS and gain >= STRONG_GAIN):
        g.update(gold_s=r["gap_s"], gold_len=r["gap_len"])     # weak evidence: the gold keeps NCBI's placement
        row["gold"] = coords(r["gap_s"], r["gap_len"], n)
        set_verdict(row, "consensus")
        row["reason"] = row["reason"] + (f"; NCBI product and junction clean; the {why} is better but the evidence is "
                                         "not strong, so the gold keeps NCBI's placement")
        row["gold_handling"] = "keep"
        return f"consensus placement differs (kept NCBI's): {why}"
    g.update(gold_s=gs, gold_len=gl)
    row["gold"] = coords(gs, gl, n)
    if ncbi_wrong:
        set_verdict(row, "D-corrected")
        row["reason"] = row["reason"].replace("; no canonical placement within 60 nt removes it (lineage-specific indel?)", "") + \
            (f"; the family lacks omega-G, so its consensus overrides the generic rule: corrected to the {why}" if override
             else f"; corrected to the {why}")
    else:
        set_verdict(row, "consensus")
        row["reason"] = row["reason"] + f"; NCBI product and junction clean; strong evidence for the {why} (gold moved)"
    row["gold_handling"] = "keep"
    why = ("corrected to the " if ncbi_wrong else "moved to the ") + why
    return why


_LIN = {}


def V_lineage(V, r):
    if r["accession"] not in _LIN:
        _LIN[r["accession"]] = SeqIO.read(str(GOLD.parent / "ncbi_gb" / f"{r['accession']}.gb"),
                                          "genbank").annotations.get("taxonomy", [])
    return _LIN[r["accession"]]


def V_ref(V, gene):
    return next(p for x in V.HOM if x["acc"] == V.REF_ACC for g, (p, _) in x["prot"].items() if g == gene)


if __name__ == "__main__":
    sys.exit(main())
