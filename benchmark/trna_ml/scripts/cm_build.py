#!/usr/bin/env python3
"""Build clade-specific mt-tRNA covariance models from training tRNAs.

1. label: every training tRNA is aligned to the MITOS2 pan-metazoan CM of its amino acid; the anticodon is
   read from the anticodon loop; tRNAs whose anticodon does not fit the annotated amino acid are dropped
   (annotation errors), Leu/Ser are split into L1/L2, S1/S2 by anticodon.
2. clades: lineage tokens (phylum and below) with enough records and genera; a token whose record set is
   (almost) the same as its parent's is skipped.
3. per clade x tRNA family with >= MIN_SEQ distinct sequences: cmalign to the pan CM -> cmbuild ->
   cmcalibrate -L 0.1.

Usage: cm_build.py <dataset_dir> <out_dir> [--folds K]   (--folds: also build K cross-fit model sets,
       each without one fold of training families, for unbiased classifier features)
"""
import argparse
import hashlib
import json
import random
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import FAMS, INF, PAN, anticodon_from_alignment, parse_sto, read_tsv, sh, write_tsv  # noqa: E402

MIN_REC, MIN_GEN, MIN_SEQ, MAX_SEQ = 40, 10, 15, 200
# anticodon -> amino acid, invertebrate mito code 5 / 9 / 2 / 4 / 13 / 14 / 24 / 33 differences handled loosely:
# accept the standard assignment plus the known mito reassignments
STD = {}
_b = "tcag"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, a in enumerate(_b):
    for j, b in enumerate(_b):
        for k, c in enumerate(_b):
            STD[a + b + c] = _aa[16 * i + 4 * j + k]
MITO_EXTRA = {"tga": {"W"}, "aga": {"S", "G", "R"}, "agg": {"S", "G", "R", "K"}, "ata": {"M", "I"},
              "aaa": {"N", "K"}, "tag": {"L"}, "taa": {"L", "Y"}}


def rc(s):
    return s[::-1].translate(str.maketrans("acgt", "tgca"))


def fits(aa, ac):
    if len(ac) != 3 or set(ac) - set("acgt"):
        return False
    codon = rc(ac)
    # wobble: anticodon position 1 pairs with codon position 3; check codon family of the 2 first bases
    return aa == STD[codon] or aa in MITO_EXTRA.get(codon, set()) or any(
        STD[codon[:2] + x] == aa for x in "tcag")


def label(ds, out):
    """Anticodon + family for every tRNA of the data set (train and test: test labels are used for scoring)."""
    out.mkdir(parents=True, exist_ok=True)
    tr = read_tsv(ds / "trnas.tsv")
    by_aa = defaultdict(list)
    for i, t in enumerate(tr):
        if t["aa"] and 40 <= int(t["length"]) <= 110 and set(t["seq"]) <= set("ACGT"):
            by_aa[t["aa"]].append(i)
    jobs, chunks = [], defaultdict(list)
    for aa, idx in by_aa.items():
        cm = PAN / f"trn{ {'L': 'L1', 'S': 'S1'}.get(aa, aa)}.cm"
        for c in range(0, len(idx), 500):              # small chunks: cmalign time grows faster than linear
            fa = out / f"lab_{aa}_{c}.fa"
            sto = out / f"lab_{aa}_{c}.sto"
            chunks[aa].append(sto)
            if sto.exists() and sto.stat().st_size:
                continue
            fa.write_text("".join(f">t{i}\n{tr[i]['seq']}\n" for i in idx[c:c + 500]))
            jobs.append((cm, fa, sto))
    with ProcessPoolExecutor(16) as ex:
        list(ex.map(_align, jobs))
    for aa, idx in by_aa.items():
        seqs = {}
        for sto in chunks[aa]:
            s_, ss, _ = parse_sto(sto.read_text())
            seqs.update((k, (v, ss)) for k, v in s_.items())
        for i in idx:
            ac = anticodon_from_alignment(*seqs[f"t{i}"]) if f"t{i}" in seqs else ""
            tr[i]["ac_aln"] = ac
    for t in tr:
        ac = t.get("ac_aln") or t["anticodon"]
        t["ac_final"] = ac if fits(t["aa"], ac) else (t["anticodon"] if fits(t["aa"], t["anticodon"]) else "")
        if t["aa"] in ("L", "S"):
            fam = {"L": {"tag": "L1", "taa": "L2"}, "S": {"tct": "S1", "gct": "S1", "tga": "S2", "cct": "S1",
                                                          "tcc": "S1", "act": "S1"}}[t["aa"]].get(t["ac_final"], "")
        else:
            fam = t["aa"] if t["ac_final"] else ""
        t["family"] = fam
        t["clean"] = int(bool(fam))
    write_tsv(ds / "trnas_labeled.tsv", tr)
    return tr


def _align(job):
    cm, fa, sto = job
    sh([str(INF / "cmalign"), "--cpu", "1", "--outformat", "Stockholm", "-o", str(sto), str(cm), str(fa)])


def clades(recs):
    """Lineage tokens (phylum and below) -> set of record accessions; nested near-duplicates removed."""
    tok = defaultdict(set)
    parent = {}
    for r in recs:
        lin = r["lineage"].split(";")
        if r["phylum"] not in lin:
            continue
        sub = lin[lin.index(r["phylum"]):]
        for d, t in enumerate(sub):
            tok[t].add(r["acc"])
            if d:
                parent.setdefault(t, sub[d - 1])
    gen = {r["acc"]: r["genus"] for r in recs}
    # data-rich insects need more records per clade; small, hard clades (Enoplea, mites, ...) get models too
    insect_tok = {t for r in recs if ";Insecta;" in r["lineage"] + ";"
                  for t in r["lineage"].split(";")[r["lineage"].split(";").index("Insecta"):]}
    keep = {}
    for t, s in tok.items():
        insect = t in insect_tok
        if len(s) < (MIN_REC if insect else 12) or len({gen[a] for a in s}) < (MIN_GEN if insect else 5):
            continue
        keep[t] = s
    out = {}
    for t, s in keep.items():
        p = parent.get(t)
        dup = False
        while p:
            if p in keep and len(s) >= 0.8 * len(keep[p]):
                dup = True
                break
            p = parent.get(p)
        if not dup:
            out[t] = s
    return out, parent


def build_one(job):
    name, fam, seqs, odir = job
    odir = Path(odir)
    fa, sto, cm = odir / f"{name}.fa", odir / f"{name}.sto", odir / f"{name}.cm"
    if cm.exists() and cm.stat().st_size and "ECMLC" in cm.read_text()[:4000]:
        return name, "cached"
    fa.write_text("".join(f">s{i}\n{s}\n" for i, s in enumerate(seqs)))
    pan = PAN / f"trn{fam}.cm"
    sh([str(INF / "cmalign"), "--cpu", "1", "--outformat", "Stockholm", "-o", str(sto), str(pan), str(fa)])
    r = sh([str(INF / "cmbuild"), "-F", "-n", name, str(cm), str(sto)])
    if r.returncode:
        return name, "cmbuild failed: " + r.stderr[-200:]
    r = sh([str(INF / "cmcalibrate"), "-L", "0.02", "--cpu", "1", str(cm)])
    fa.unlink(missing_ok=True)
    return name, "ok" if r.returncode == 0 else "calibrate failed"


def build_set(recs_train, tr, odir, jobs=16, reuse=None):
    odir.mkdir(parents=True, exist_ok=True)
    cl, parent = clades(recs_train)
    if reuse and (reuse / "models.json").exists():   # same clade record set -> same models: copy them
        old = json.load(open(reuse / "models.json"))["clades"]
        for c, s in cl.items():
            if sorted(s) == old.get(c):
                for f in reuse.glob(f"{c}__*.cm"):
                    dst = odir / f.name
                    if not dst.exists():
                        dst.write_bytes(f.read_bytes())
    acc_ok = {r["acc"] for r in recs_train}
    by = defaultdict(lambda: defaultdict(set))          # clade -> fam -> seqs
    rec_cl = defaultdict(list)
    for t, s in cl.items():
        for a in s:
            rec_cl[a].append(t)
    for t in tr:
        if t["acc"] in acc_ok and str(t["clean"]) == "1":
            for c in rec_cl[t["acc"]]:
                by[c][t["family"]].add(t["seq"])
    work, meta = [], {}
    for c, fams in by.items():
        for fam in FAMS:
            seqs = sorted(fams.get(fam, ()))
            if len(seqs) > MAX_SEQ:                     # cmalign time grows fast; a seeded sample keeps diversity
                seqs = sorted(random.Random(f"{c}{fam}").sample(seqs, MAX_SEQ))
            if len(seqs) >= MIN_SEQ:
                name = f"{c}__{fam}"
                work.append((name, fam, seqs, str(odir)))
                meta[name] = {"clade": c, "family": fam, "n_seq": len(seqs), "n_rec": len(cl[c])}
    with ProcessPoolExecutor(jobs) as ex:
        res = dict(ex.map(build_one, work))
    for k, v in res.items():
        meta[k]["status"] = v
    json.dump({"models": meta, "clades": {c: sorted(s) for c, s in cl.items()},
               "parent": {c: parent.get(c, "") for c in cl}}, open(odir / "models.json", "w"), indent=1)
    bad = [k for k, v in res.items() if v not in ("ok", "cached")]
    print(f"{odir}: {len(cl)} clades, {len(work)} models, {len(bad)} failed", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset")
    ap.add_argument("out")
    ap.add_argument("--folds", type=int, default=0)
    ap.add_argument("--exclude_families", default="", help="comma list: build only out/<--name> without these families")
    ap.add_argument("--name", default="hard")
    ap.add_argument("--jobs", type=int, default=16)
    a = ap.parse_args()
    ds, out = Path(a.dataset), Path(a.out)
    lab = ds / "trnas_labeled.tsv"
    tr = read_tsv(lab) if lab.exists() else label(ds, out / "label")
    recs = [r for r in read_tsv(ds / "records.tsv") if r["split"] == "train"]
    if a.exclude_families:
        ex = set(a.exclude_families.split(","))
        build_set([r for r in recs if r["family"] not in ex], tr, out / a.name, a.jobs, reuse=out / "full")
        return
    build_set(recs, tr, out / "full", a.jobs)
    for k in range(a.folds):
        sub = [r for r in recs if int(hashlib.md5((r["family"] or r["genus"]).encode()).hexdigest()[:8], 16) % a.folds != k]
        build_set(sub, tr, out / f"fold{k}", a.jobs)


if __name__ == "__main__":
    main()
