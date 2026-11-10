"""All NCBI animal mitogenomes assembled by genome projects (title 'genome assembly, organelle: mitochondrion':
DToL / ERGA / VGP ..., mostly PacBio HiFi + MitoHiFi). esearch with history, efetch gbwithparts 200 per request
(tool=mitom, no e-mail)."""
import re, subprocess, time, sys
from pathlib import Path
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TERM = "Metazoa[Organism] AND mitochondrion[filter] AND genome assembly[Title] AND organelle[Title]"
def curl(*a):
    for _ in range(5):
        p = subprocess.run(["curl", "-s", "--max-time", "900", *a], capture_output=True, text=True)
        if p.returncode == 0 and p.stdout.strip():
            return p.stdout
        time.sleep(10)
    raise SystemExit("NCBI request failed")
x = curl("-G", f"{E}/esearch.fcgi", "--data-urlencode", "db=nuccore", "--data-urlencode", "tool=mitom",
         "--data-urlencode", "usehistory=y", "--data-urlencode", "retmax=0", "--data-urlencode", f"term={TERM}")
n = int(re.search(r"<Count>(\d+)", x).group(1)); web = re.search(r"<WebEnv>(\S+?)</WebEnv>", x).group(1)
key = re.search(r"<QueryKey>(\d+)", x).group(1)
print(n, "records", file=sys.stderr)
for i in range(0, n, 200):
    f = Path("gb") / f"chunk_{i // 200:03d}.gb"
    if f.exists() and f.read_text()[-200:].rstrip().endswith("//"):
        continue
    f.write_text(curl(f"{E}/efetch.fcgi?db=nuccore&tool=mitom&rettype=gbwithparts&retmode=text&WebEnv={web}"
                      f"&query_key={key}&retstart={i}&retmax=200"))
    time.sleep(0.5)
print("done", file=sys.stderr)
