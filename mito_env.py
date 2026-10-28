#!/usr/bin/env python3
"""mito_env.py - find the external programs of the mitogenome annotation scripts, one place for all.

Every script (batch_mitoz.py, mito_intron.py, mito_trna.py, mito_curate.py, ncbi_validate.py) asks
this module where a program is, so no --bin / --intron_bin / --mitos_bin / --mitos_refdir has to be
typed. Order: command-line value > config file > PATH > automatic scan of conda / mamba environments.

Config file (JSON, all keys optional), first found of:
  $YZZ_MITO_CONFIG, <yzz root>/config/mitoz.json (inside the yzz workflow), <this dir>/mito_config.json,
  ~/.config/yzz/mito.json
  {"mitos_bin": ".../envs/mitos/bin",
   "mitos_refdir": "/opt/mitos_refs", "tools_bin": ".../envs/mitotools/bin",
   "intron_bin": ".../envs/mitointron/bin", "trna_bin": ".../envs/mitotrna/bin"}

YZZ_MITO_NO_SCAN=1 turns the conda scan off (only command line, config file and PATH are used).

Usage:
  python3 mito_env.py doctor                  # what is found, what is missing, how to install it
  python3 mito_env.py doctor --write-config   # also save the found locations to the config file
"""
from __future__ import annotations

import json
import os
import shutil
import sys
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional

HERE = Path(__file__).resolve().parent
MITOS_REF = "refseq89m"
# what each program is needed for, and the conda environment (env/*.yml) that installs it
TOOLS = {
    "runmitos":    ("optional: MITOS2 rescue of missed genes (mito_rescue.py, --rescue yes), mito_trna.py", "mitos"),
    "exonerate":   ("protein gene boundaries (mito_core.py); intron models (mito_intron.py)", "mitotools"),
    "cmsearch":    ("tRNA / rRNA covariance models (mito_core.py); intron typing (mito_intron.py)", "mitotools"),
    "blastx":      ("protein gene search (mito_core.py)", "mitotools"),
    "java":        ("MiTFi tRNA selection (mito_core.py, the default core)", "mitotools"),
    "blastp":      ("intron homolog search (mito_intron.py)", "mitotools"),
    "table2asn":   ("NCBI validation (mito_curate.py, ncbi_validate.py)", "mitotools"),
    "nhmmer":      ("per-phylum rRNA profile HMMs (mito_core.py, mito_rrna_hmm.py)", "mitotools"),
    "mafft":       ("optional: rebuild the rRNA HMM library (mito_rrna_hmm.py build)", "mitotools"),
    "arwen":       ("optional extra tRNA method (mito_trna.py)", "mitotools"),
    "aragorn":     ("optional extra tRNA method (mito_trna.py)", "mitotools"),
    "tRNAscan-SE": ("optional extra tRNA method (mito_trna.py)", "mitotools"),
    "blastn":      ("rDNA 18S-28S annotation (mito_nuc.py, mitom nuc)", "mitotools"),
    "tblastn":     ("histone annotation (mito_nuc.py, mitom nuc)", "mitotools"),
    "bwa":         ("optional: read depth ring of the genome maps (mito_plot.py, --infq)", "mitotools"),
    "minimap2":    ("optional: used instead of bwa for the depth ring when found (about 8x faster; on Linux "
                    "mamba install -n mitotools -c bioconda minimap2 adds 0.5 GB)", "mitotools"),
}
REQUIRED = ("blastx", "exonerate", "cmsearch", "nhmmer", "java", "blastp", "table2asn")   # the native core (default)


def config_path() -> Path:
    env = os.environ.get("YZZ_MITO_CONFIG")
    yzz = HERE.parents[1] / "config" / "mitoz.json" if (HERE.parents[1] / "yzz.py").exists() else None
    if env:
        return Path(env)                                       # explicit: used (and written) even before it exists
    cands = [c for c in (yzz, HERE / "mito_config.json", Path.home() / ".config" / "yzz" / "mito.json") if c]
    return next((c for c in cands if c.exists()), cands[0])   # new file: yzz config, or here


@lru_cache(maxsize=None)
def config() -> Dict[str, str]:
    p = config_path()
    try:
        return json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}
    except (OSError, ValueError) as exc:
        print(f"[mito_env] cannot read {p}: {exc}", file=sys.stderr)
        return {}


@lru_cache(maxsize=None)
def conda_roots() -> List[Path]:
    """Base directories of conda / mamba / micromamba installations (none when YZZ_MITO_NO_SCAN=1)."""
    if os.environ.get("YZZ_MITO_NO_SCAN") == "1":
        return []
    roots = []
    for var in ("CONDA_EXE", "MAMBA_EXE"):
        if os.environ.get(var):
            roots.append(Path(os.environ[var]).resolve().parents[1])
    for var in ("CONDA_PREFIX", "MAMBA_ROOT_PREFIX", "CONDA_ROOT"):
        if os.environ.get(var):
            p = Path(os.environ[var])
            roots.append(p.parents[1] if p.parent.name == "envs" else p)
    home = Path.home()
    roots += [home / d for d in ("miniforge3", "mambaforge", "miniconda3", "anaconda3", "micromamba", ".local/share/mamba")]
    roots += [Path("/opt/conda"), Path("/opt/miniforge3"), Path("/opt/miniconda3")]
    exe = shutil.which("conda") or shutil.which("mamba")
    if exe:
        roots.append(Path(exe).resolve().parents[1])
    out = []
    for r in roots:
        if r.is_dir() and r not in out:
            out.append(r)
    return out


@lru_cache(maxsize=None)
def env_bins() -> List[Path]:
    """bin/ of every conda environment (named environments first, then the base)."""
    out = []
    for r in conda_roots():
        envs = r / "envs"
        if envs.is_dir():
            out += sorted(e / "bin" for e in envs.iterdir() if (e / "bin").is_dir())
        if (r / "bin").is_dir():
            out.append(r / "bin")
    return out


def _exe(p: Path) -> bool:
    return p.is_file() and os.access(p, os.X_OK)


def find(name: str, explicit: Optional[str] = None, *config_keys: str) -> Optional[str]:
    """Path of program `name`: explicit dir/file (only there: a wrong explicit path is an error, never
    silently replaced) > config bins > PATH > conda environments."""
    cfg = config()
    dirs = [explicit] if explicit else \
        [cfg.get(k, "") for k in (*config_keys, "tools_bin", "intron_bin", "trna_bin", "mitos_bin")]
    for d in dirs:
        if not d:
            continue
        p = Path(d).expanduser()
        if p.is_file() and p.name == name:
            return str(p)
        if _exe(p / name):
            return str(p / name)
    if explicit:
        return None
    w = shutil.which(name)
    if w:
        return w
    for b in env_bins():
        if _exe(b / name):
            return str(b / name)
    return None


def bin_with(*names: str, explicit: Optional[str] = None) -> Optional[str]:
    """One bin directory that holds all `names` (e.g. exonerate + cmsearch + blastp for mito_intron)."""
    if explicit:
        return explicit if all(_exe(Path(explicit).expanduser() / n) for n in names) else None
    cfg = config()
    cands = [cfg.get(k, "") for k in ("tools_bin", "intron_bin", "trna_bin")]
    first = find(names[0])
    if first:
        cands.append(str(Path(first).parent))
    cands += [str(b) for b in env_bins()]
    for d in cands:
        if d and all(_exe(Path(d).expanduser() / n) for n in names):
            return str(Path(d).expanduser())
    return None


def mitos_bin(explicit: Optional[str] = None) -> Optional[str]:
    exe = find("runmitos", explicit, "mitos_bin")
    return str(Path(exe).parent) if exe else None


def mitos_refdir(explicit: Optional[str] = None, ref: str = MITOS_REF) -> Optional[str]:
    """Directory that contains the MITOS2 reference data folder `ref` (e.g. refseq89m)."""
    if explicit:                                  # used as given; a wrong one makes MITOS2 fail visibly
        return explicit
    cands = [config().get("mitos_refdir"), os.environ.get("MITOS_REFDIR"),
             "/opt/mitos_refs", str(Path.home() / "mitos_refs"), str(Path.home() / ".local/share/mitos_refs"),
             str(HERE / "mitos_refs")]
    mb = mitos_bin()
    if mb:
        cands += [str(Path(mb).parent / "share" / "mitos_refs"), str(Path(mb).parent / "mitos_refs")]
    for c in cands:
        if c and (Path(c).expanduser() / ref).is_dir():
            return str(Path(c).expanduser())
    return None


MITOS_MODES = ("fast", "sensitive", "full")


def mitos_mode(value: Optional[str], default: str = "sensitive") -> str:
    """YZZ_MITOS_MODE: case-insensitive; an unknown value is reported on stderr and replaced by `default`
    (it used to behave as 'sensitive' silently)."""
    v = (value or "").strip().lower()
    if not v:
        return default
    if v not in MITOS_MODES:
        print(f"[mito_env] YZZ_MITOS_MODE='{value}' is not one of {', '.join(MITOS_MODES)}; using {default}",
              file=sys.stderr)
        return default
    return v


def biopython_version() -> Optional[str]:
    """Version of Biopython importable by this interpreter, None when it cannot be imported (every MitoM
    script needs it)."""
    try:
        import Bio
        return str(Bio.__version__)
    except ImportError:
        return None


def matplotlib_version() -> Optional[str]:
    """Version of matplotlib importable by this interpreter (the genome maps of mito_plot.py), None when missing."""
    try:
        import matplotlib
        return str(matplotlib.__version__)
    except ImportError:
        return None


def doctor(write: bool = False) -> int:
    print(f"config: {config_path()}{'' if config_path().exists() else ' (not present)'}")
    print(f"conda/mamba: {', '.join(map(str, conda_roots())) or 'not found'}")
    missing = []
    for name, (use, env_name) in TOOLS.items():
        exe = find(name)
        flag = "OK " if exe else ("-- " if name not in REQUIRED else "!! ")
        print(f"  {flag}{name:12} {exe or 'not found':60} {use}")
        if not exe:
            missing.append((name, env_name))
    ref = mitos_refdir()
    print(f"  {'OK ' if ref else '!! '}{'MITOS2 refs':12} {ref or 'not found (' + MITOS_REF + ')':60} MITOS2 reference data")
    bio = biopython_version()          # required: all MitoM scripts read GenBank files with it
    print(f"  {'OK ' if bio else '!! '}{'Biopython':12} "
          f"{('version ' + bio) if bio else 'not found for ' + sys.executable:60} sequence input / output (all scripts)")
    mpl = matplotlib_version()         # optional: without it run with --plot no
    print(f"  {'OK ' if mpl else '-- '}{'matplotlib':12} "
          f"{('version ' + mpl) if mpl else 'not found for ' + sys.executable:60} genome maps (mito_plot.py, --plot)")
    req = [n for n, _ in missing if n in REQUIRED]
    if req or not ref or not bio:
        envs = sorted({e for n, e in missing if n in REQUIRED})
        print("\nMissing required parts. Install everything with one command:\n"
              f"  bash {HERE / 'install.sh'}" + (f"   (environments: {', '.join(envs)})" if envs else "")
              + ("\n  MITOS2 reference data is downloaded by the same script." if not ref else "")
              + ("\n  Biopython is missing: pip install biopython (install.sh does this too)." if not bio else ""))
    else:
        print("\nAll required programs found.")
    if write:
        cfg_file = config_path()
        cfg = dict(config())
        mb, tb = mitos_bin(), bin_with("exonerate", "cmsearch", "blastp", "makeblastdb")
        for k, v in (("mitos_bin", mb), ("mitos_refdir", ref), ("tools_bin", tb)):
            if v and not cfg.get(k):
                cfg[k] = v
        cfg_file.parent.mkdir(parents=True, exist_ok=True)
        cfg_file.write_text(json.dumps(cfg, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"config written: {cfg_file}")
    return 1 if req or not ref or not bio else 0


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "doctor":
        sys.exit(doctor(write="--write-config" in sys.argv))
    print(__doc__)
