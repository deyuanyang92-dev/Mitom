#!/usr/bin/env python3
"""Recall of mito_orf on make_truth.py: each truth CDS -> the reported ORF on its strand covering >= 50 % of it
(written / warned:<class> / missed). Left out: NC_137528.1 (its 13 genes carry names gene_key does not know) and
products "ATPase subunit 6 / 8" (standard genes). Usage: score_recall.py <work>/recall <out folder name>
-> <work>/recall/<out>_scored.tsv"""
import csv,sys,collections,re
W=sys.argv[1]; out=sys.argv[2] if len(sys.argv)>2 else "out"
truth=list(csv.DictReader(open(f"{W}/truth.tsv"),delimiter="\t"))
res=collections.Counter(); bycat=collections.defaultdict(collections.Counter); rows=[]
def cat(p):
    p=p.lower()
    if "muts" in p or "mismatch" in p: return "mutS"
    if "subunit 9" in p or "atp9" in p: return "atp9"
    if "atpase subunit" in p: return "std-named"
    if "hypothetical" in p or "orf" in p: return "hypothetical/ORF"
    if "endonuclease" in p or "laglidadg" in p or "giy" in p: return "homing endonuclease"
    return "other"
for t in truth:
    c=cat(t["product"])
    if c=="std-named" or t["record"]=="NC_137528.1": continue
    s,e=int(t["start"]),int(t["end"]); st="+" if t["strand"]=="1" else "-"
    try: rep=list(csv.DictReader(open(f"{W}/{out}/{t['record']}/orf_report.tsv"),delimiter="\t"))
    except FileNotFoundError: rep=[]
    best=None
    for r in rep:
        if r["strand"]!=st: continue
        a,b=int(r["start"]),int(r["end"])
        ov=max(0,min(b,e)-max(a,s)+1)
        if ov>=0.5*(e-s+1) and (best is None or ov>best[0]): best=(ov,r)
    if best is None: k="missed"
    elif best[1]["written"]=="yes": k="written"
    else: k="warned:"+best[1]["class"]
    if t["parts"]!="1": k+=" (truth has introns)"
    res[k]+=1; bycat[c][k.split(" ")[0].split(":")[0]]+=1
    rows.append([t["record"],t["phylum"],t["product"][:30],t["codons"],t["parts"],k,best[1]["class"] if best else "",best[1]["evidence"][:100] if best else ""])
print(sum(res.values()),dict(res))
for c,v in bycat.items(): print(" ",c,dict(v))
csv.writer(open(f"{W}/{out}_scored.tsv","w"),delimiter="\t").writerows(rows)
