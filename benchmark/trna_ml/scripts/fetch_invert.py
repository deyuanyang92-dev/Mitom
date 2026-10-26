#!/usr/bin/env python3
"""Download all invertebrate RefSeq (NC_) complete mitogenomes as GenBank, in batches, resumable."""
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
Q = ('Metazoa[Organism] NOT Vertebrata[Organism] AND mitochondrion[Filter] AND "complete genome"[Title] '
     'AND refseq[Filter] AND 10000:60000[SLEN]')
OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "fetch")
OUT.mkdir(parents=True, exist_ok=True)


def get(url, data=None, tries=6):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, data=data, timeout=180) as r:
                return r.read().decode()
        except Exception as e:  # noqa: BLE001
            print("retry", k, e, flush=True)
            time.sleep(5 * (k + 1))
    raise RuntimeError(url)


ids_file = OUT / "ids.txt"
if not ids_file.exists():
    r = get(E + "esearch.fcgi?" + urllib.parse.urlencode({"db": "nuccore", "term": Q, "retmax": 100000, "idtype": "acc"}))
    ids = [x.split("</Id>")[0] for x in r.split("<Id>")[1:]]
    ids = [i for i in ids if i.startswith("NC_")]
    ids_file.write_text("\n".join(ids) + "\n")
ids = ids_file.read_text().split()
print(len(ids), "records", flush=True)
B = 100
for s in range(0, len(ids), B):
    f = OUT / f"batch_{s // B:04d}.gb"
    if f.exists() and f.stat().st_size > 0:
        continue
    txt = get(E + "efetch.fcgi", urllib.parse.urlencode(
        {"db": "nuccore", "id": ",".join(ids[s:s + B]), "rettype": "gbwithparts", "retmode": "text"}).encode())
    if txt.count("\n//") < len(ids[s:s + B]) * 0.9:
        print("short batch", s, flush=True)
        continue
    f.write_text(txt)
    print("batch", s // B, flush=True)
    time.sleep(0.4)
print("done", flush=True)
