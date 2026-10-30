#!/usr/bin/env python3
"""mito_env.py - find the external programs of the mitogenome annotation scripts, one place for all.

Every script (mitom_pipeline.py, mito_intron.py, mito_trna.py, mito_curate.py, ncbi_validate.py) asks
this module where a program is, so no --bin / --intron_bin / --mitos_bin / --mitos_refdir has to be
typed. Order: command-line value > config file > the conda environment of the running Python (MitoM's own,
env/mitotools.yml) > PATH > automatic scan of conda / mamba environments.

Config file (JSON, all keys optional), first found of:
  $YZZ_MITO_CONFIG, <yzz root>/config/mitoz.json (inside the yzz workflow), <this dir>/mito_config.json,
  ~/.config/yzz/mito.json
  {"mitos_bin": ".../envs/mitos/bin",
   "mitos_refdir": "/opt/mitos_refs", "tools_bin": ".../envs/mitotools/bin",
   "intron_bin": ".../envs/mitointron/bin", "trna_bin": ".../envs/mitotrna/bin"}

YZZ_MITO_NO_SCAN=1 turns the conda scan off (only command line, config file and PATH are used).

Usage:
  python3 mito_env.py setup [--refdir DIR]    # download the MITOS2 reference data refseq89m, write the config
  python3 mito_env.py doctor                  # what is found, what is missing, how to install it
  python3 mito_env.py doctor --write-config   # also save the found locations to the config file
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys
import tarfile
import urllib.request
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional

# Numerical libraries (NumPy / OpenBLAS / OpenMP of scikit-learn, used by the tRNA classifier) start one thread per
# machine core in every process unless told otherwise. `mitom run` annotates several records at once, each process
# with --threads for the external tools: on a 128-core machine 4 records ran 99 s instead of 26 s (OMP_NUM_THREADS=1),
# and a whole `mitom run` of 4 mitogenomes took 47 min. MitoM's numerical work is small; the external programs use
# --threads. Set before NumPy / scikit-learn are imported (every module imports mito_env first); values the user set
# are kept.
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_var, "1")

HERE = Path(__file__).resolve().parent
MITOS_REF = "refseq89m"
REF_URL = "https://zenodo.org/records/4284483/files/refseq89m.tar.bz2"      # MITOS2 reference data, 19.6 MB
REF_MD5 = "bb6325e27612e61a6995e88e8ffcecf8"
# what each program is needed for, and the conda environment (env/*.yml) that installs it
TOOLS = {
    "runmitos":    ("optional: MITOS2 rescue of missed genes (mito_rescue.py, --rescue yes), mito_trna.py", "mitos"),
    "exonerate":   ("protein gene boundaries (mito_core.py); intron models (mito_intron.py)", "mitotools"),
    "cmsearch":    ("tRNA / rRNA covariance models (mito_core.py); intron typing (mito_intron.py)", "mitotools"),
    "cmalign":     ("MitoM-ML tRNAs, the default tRNA caller (trna_ml.py; without it MiTFi only)", "mitotools"),
    "blastx":      ("protein gene search (mito_core.py)", "mitotools"),
    "java":        ("MiTFi tRNA selection (mito_core.py, the default core); any Java >= 8", "mitotools"),
    "blastp":      ("intron homolog search (mito_intron.py)", "mitotools"),
    "makeblastdb": ("intron homolog library (mito_intron.py); nuclear marker references (mitom nuc)", "mitotools"),
    "table2asn":   ("NCBI validation (mito_curate.py, ncbi_validate.py)", "mitotools"),
    "nhmmer":      ("per-phylum rRNA profile HMMs (mito_core.py, mito_rrna_hmm.py)", "mitotools"),
    "mafft":       ("optional: rebuild the rRNA HMM library (mito_rrna_hmm.py build)", "mitotools"),
    "hmmbuild":    ("optional: rebuild the rRNA HMM library (mito_rrna_hmm.py build)", "mitotools"),
    "arwen":       ("optional extra tRNA method (mito_trna.py; extra evidence of MitoM-ML)", "mitotools"),
    "aragorn":     ("optional extra tRNA method (mito_trna.py; extra evidence of MitoM-ML)", "mitotools"),
    "tRNAscan-SE": ("optional extra tRNA method (mito_trna.py)", "mitotools"),
    "blastn":      ("intron site probes (mito_intron_ref.py); rDNA 18S-28S (mitom nuc)", "mitotools"),
    "tblastn":     ("intron site anchors (mito_intron_ref.py); histones (mitom nuc)", "mitotools"),
    "bwa":         ("optional: read depth ring of the genome maps (mito_plot.py, --infq)", "mitotools"),
    "minimap2":    ("optional: used instead of bwa for the depth ring when found (about 8x faster; on Linux "
                    "mamba install -n mitotools -c bioconda minimap2 adds 0.5 GB)", "mitotools"),
}
# the programs each step of `mitom run` needs: doctor checks every one, the preflight of a run the steps it will do
REQUIRES = {
    "core": ("blastx", "exonerate", "cmsearch", "nhmmer", "java"),            # java: MiTFi
    "trna_ml": ("cmsearch", "cmalign"),                                       # MitoM-ML tRNAs (else MiTFi only)
    "introns": ("exonerate", "blastp", "makeblastdb", "cmsearch", "tblastn", "blastn"),
    "submission": ("table2asn",),
    "rrna_hmm": ("nhmmer",),
}
STEP_NAMES = {"core": "the annotation core", "trna_ml": "MitoM-ML tRNAs", "introns": "the intron step (--intron_fix)",
              "submission": "the NCBI submission files", "rrna_hmm": "the rRNA HMM step (--rrna_hmm)"}
REQUIRED = tuple(dict.fromkeys(p for progs in REQUIRES.values() for p in progs))   # the default `mitom run`


def missing_programs(steps) -> List[str]:
    """'<program> (for <step>)' for every program the given REQUIRES steps need and find() does not find."""
    out = []
    for step in steps:
        for prog in REQUIRES[step]:
            if not find(prog):
                out.append(f"{prog} (for {STEP_NAMES[step]})")
    return list(dict.fromkeys(out))


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


def own_env_bin() -> Optional[str]:
    """bin/ of the conda environment this Python runs in (MitoM's own when started by `mitom`), else None."""
    return str(Path(sys.executable).parent) if (Path(sys.prefix) / "conda-meta").is_dir() else None


def find(name: str, explicit: Optional[str] = None, *config_keys: str) -> Optional[str]:
    """Path of program `name`: explicit dir/file (only there: a wrong explicit path is an error, never
    silently replaced) > config bins > the running Python's conda environment > PATH > conda environments."""
    cfg = config()
    dirs = [explicit] if explicit else \
        [cfg.get(k, "") for k in (*config_keys, "tools_bin", "intron_bin", "trna_bin", "mitos_bin")] + [own_env_bin()]
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


def data_dir() -> Path:
    """Folder for the data MitoM makes on first use (MitoM-ML: unpacked clade models, trained classifiers), so that
    the install folder may be read-only: $MITOM_DATA, else $XDG_DATA_HOME/mitom, else ~/.local/share/mitom."""
    if os.environ.get("MITOM_DATA"):
        return Path(os.environ["MITOM_DATA"]).expanduser()
    xdg = os.environ.get("XDG_DATA_HOME")
    return (Path(xdg).expanduser() if xdg else Path.home() / ".local" / "share") / "mitom"


def refseq_dir(refdir: Optional[str] = None, ref: str = MITOS_REF) -> Optional[Path]:
    """The MITOS2 reference data folder itself (refseq89m: featureProt/, ncRNA/), given that folder or the folder
    that holds it (as mitos_refdir returns it); found automatically when refdir is None; None when not found."""
    d = refdir or mitos_refdir(None, ref)
    if not d:
        return None
    p = Path(d).expanduser()
    return p / ref if (p / ref).is_dir() else p


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


# conda packages of the programs (the environment of install.sh, env/mitotools.yml): the fix `mitom debug` names
PACKAGES = {"exonerate": "exonerate", "cmsearch": "infernal", "cmalign": "infernal", "blastx": "blast", "blastp": "blast",
            "makeblastdb": "blast", "blastn": "blast", "tblastn": "blast", "table2asn": "table2asn", "nhmmer": "hmmer",
            "hmmbuild": "hmmer", "mafft": "mafft", "aragorn": "aragorn", "tRNAscan-SE": "trnascan-se", "bwa": "bwa",
            "minimap2": "minimap2", "java": "openjdk=8"}
VERSION_ARG = {"exonerate": "--version", "cmsearch": "-h", "cmalign": "-h", "nhmmer": "-h", "hmmbuild": "-h",
               "java": "-version", "mafft": "--version", "aragorn": "-h", "arwen": "-h", "tRNAscan-SE": "-h", "bwa": "",
               "minimap2": "--version", "runmitos": "--version"}          # others (BLAST, table2asn): -version


def install_hint(name: str) -> str:
    """How to install program `name` into MitoM's environment."""
    if name == "runmitos":
        return "optional: MITOS2 in its own conda environment (mamba create -n mitos -c bioconda mitos=2.1.10)"
    if name in PACKAGES:
        return (f"bash {HERE / 'install.sh'}, or: mamba install -n mitotools -c conda-forge -c bioconda {PACKAGES[name]}"
                + ("  (or any system Java >= 8)" if name == "java" else ""))
    return f"bash {HERE / 'install.sh'}"


def program_version(exe: str, name: Optional[str] = None) -> str:
    """The line with the version number that program `exe` prints (VERSION_ARG; BLAST style -version by default):
    the first one naming the program or a version, else the first with a number; '' when there is none or it does not
    answer within 10 s."""
    import re
    import subprocess
    arg = VERSION_ARG.get(name or Path(exe).name, "-version")
    try:
        p = subprocess.run([exe] + ([arg] if arg else []), stdin=subprocess.DEVNULL, capture_output=True, text=True,
                           errors="replace", timeout=10)
    except (OSError, subprocess.SubprocessError):
        return ""
    lines = [ln.strip() for ln in (p.stdout + "\n" + p.stderr).splitlines() if re.search(r"\d+\.\d+", ln)]
    named = [ln for ln in lines if re.search(rf"version|infernal|hmmer|{re.escape(name or Path(exe).name)}", ln, re.I)]
    return (named or lines or [""])[0][:100]


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
              + ("\n  MITOS2 reference data: mitom setup (install.sh runs it)." if not ref else "")
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


def setup(refdir: Optional[str] = None) -> int:
    """Download the MITOS2 reference data refseq89m (Zenodo 4284483, md5 checked) into refdir (default
    $MITOS_REFDIR or ~/mitos_refs) unless it is found, then save the found locations (doctor --write-config)."""
    have = refdir if refdir and (Path(refdir).expanduser() / MITOS_REF).is_dir() else None if refdir else mitos_refdir()
    if have:
        print(f"MITOS2 reference data present: {Path(have).expanduser() / MITOS_REF}")
    else:
        target = Path(refdir or os.environ.get("MITOS_REFDIR") or Path.home() / "mitos_refs").expanduser()
        target.mkdir(parents=True, exist_ok=True)
        tgz = target / f"{MITOS_REF}.tar.bz2"
        print(f"downloading the MITOS2 reference data {MITOS_REF} (19.6 MB) to {target}", flush=True)
        try:
            urllib.request.urlretrieve(REF_URL, tgz)            # honours https_proxy
        except OSError as exc:
            print(f"download failed: {exc}\n  get {REF_URL} yourself and run: mitom setup --refdir <its folder>",
                  file=sys.stderr)
            return 1
        if hashlib.md5(tgz.read_bytes()).hexdigest() != REF_MD5:
            tgz.unlink()
            print(f"checksum mismatch for {tgz.name}: run mitom setup again", file=sys.stderr)
            return 1
        with tarfile.open(tgz) as tf:
            tf.extractall(target, **({"filter": "data"} if hasattr(tarfile, "data_filter") else {}))
        tgz.unlink()
        have = str(target)
    os.environ["MITOS_REFDIR"] = str(Path(have).expanduser())  # doctor --write-config records it
    config.cache_clear()
    prepare_trna_ml()
    return doctor(write=True)


def prepare_trna_ml() -> None:
    """MitoM-ML's clade models unpacked and classifiers trained now (into data_dir()), not by the first run's
    parallel jobs; without scikit-learn the core uses MiTFi and nothing is prepared."""
    try:
        sys.path.insert(0, str(HERE))
        import trna_ml
        print(f"MitoM-ML tRNA models ready: {trna_ml.prepare()}", flush=True)
    except Exception as exc:                       # noqa: BLE001  (setup goes on; the first run tries again)
        print(f"MitoM-ML tRNA models not prepared ({exc!r}); the first run prepares them", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "doctor":
        sys.exit(doctor(write="--write-config" in sys.argv))
    if len(sys.argv) > 1 and sys.argv[1] == "setup":
        import argparse
        a = argparse.ArgumentParser(prog="mitom setup", description=setup.__doc__)
        a.add_argument("--refdir", default=None, help="folder for refseq89m (default: $MITOS_REFDIR or ~/mitos_refs)")
        sys.exit(setup(a.parse_args(sys.argv[2:]).refdir))
    print(__doc__)
