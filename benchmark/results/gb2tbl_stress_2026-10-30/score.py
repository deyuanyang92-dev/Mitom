"""table2asn messages (ERROR / WARNING by code) of both converters per case; submitter / publication messages left out."""
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, "/media/deyuan/18T_B/yzz/MitoM")
import ncbi_validate as nv

R = Path(__file__).resolve().parent
for case in sorted(p.stem for p in (R / "in").glob("*.gb")):
    line = [f"{case:30s}"]
    for tool in ("mitom", "ncbi"):
        vals = list((R / tool / case).glob("*.val"))
        if not vals:
            line.append(f"{tool}: no output")
            continue
        c = Counter()
        for v in vals:
            for ln in v.read_text(errors="ignore").splitlines():
                m = nv.MSG_RE.match(ln.strip())
                if m and not nv.IGNORE_CODES.search(m.group(2)) and m.group(1).upper() in ("ERROR", "REJECT", "FATAL", "WARNING"):
                    c[(m.group(1).upper()[0], m.group(2).split(".")[-1])] += 1
        errs = sum(v for k, v in c.items() if k[0] != "W")
        warns = sum(v for k, v in c.items() if k[0] == "W")
        top = ", ".join(f"{k[0]}:{k[1]}x{v}" for k, v in c.most_common(4))
        line.append(f"{tool}: E{errs:3d} W{warns:3d} [{top}]")
    print("  ".join(line))
