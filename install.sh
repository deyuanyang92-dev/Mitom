#!/usr/bin/env bash
# install.sh - one command installs everything MitoM needs into ONE conda environment, mitotools (env/mitotools.yml):
#   Python with Biopython, scikit-learn (MitoM-ML tRNAs) and matplotlib (genome maps) - the `mitom` command runs on
#   this Python, nothing is installed into the system Python - and the programs (BLAST, exonerate, Infernal, HMMER,
#   table2asn, ARAGORN, tRNAscan-SE, MAFFT, bwa). Java for MiTFi: the system Java when there is one, else openjdk 8
#   is added to the environment. Then the `mitom` command (in the environment and
#   in ~/.local/bin) and `mitom setup`: the MITOS2 reference data refseq89m (Zenodo 4284483, md5 checked) and the
#   config file, so no tool path has to be typed.
#   `--envs "mitotools mitos"` also installs MITOS2 2.1.10 (optional: --engine mitos2, --rescue yes, mito_trna.py).
#   MitoZ is never installed (optional --engine mitoz: install it yourself, see mito_engines.py).
# Existing environments are kept (--force recreates them). Safe to run again.
#
# Platforms (checked with mamba solves, 2026-10-20):
#   Linux x86_64 / Windows WSL2   everything
#   macOS Intel                   mitotools; the optional MITOS2 has no macOS build (viennarna<2)
#   macOS Apple Silicon           as macOS Intel, through Rosetta 2 (osx-64 packages)
#   Windows (native)              not supported by the bioinformatics packages: use WSL2 or Docker
#
# Usage: bash install.sh [--force] [--refdir DIR] [--envs "mitotools mitos"]
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REFDIR="${MITOS_REFDIR:-$HOME/mitos_refs}"
ENVS="mitotools"            # MitoM needs only this one; --envs "mitotools mitos" adds the optional MITOS2
FORCE=0

while [ $# -gt 0 ]; do
  case "$1" in
    --force) FORCE=1 ;;
    --refdir) REFDIR="$2"; shift ;;
    --envs) ENVS="$2"; shift ;;
    -h|--help) sed -n '2,16p' "$0"; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
  shift
done

say() { printf '\n[install] %s\n' "$*"; }

# ── platform ─────────────────────────────────────────────────────────────────
OS="$(uname -s)"; ARCH="$(uname -m)"
case "$OS" in
  Linux) ;;
  Darwin)
    if [ "$ARCH" = "arm64" ]; then
      if ! arch -x86_64 /usr/bin/true 2>/dev/null; then
        echo "Apple Silicon: Rosetta 2 is needed. Run: softwareupdate --install-rosetta --agree-to-license" >&2
        exit 1
      fi
      export CONDA_SUBDIR=osx-64          # bioconda has no arm64 build of table2asn
      say "Apple Silicon: installing Intel (osx-64) packages, run through Rosetta 2"
    fi
    keep=""; for e in $ENVS; do [ "$e" != mitos ] && keep="$keep $e"; done; ENVS="${keep# }"
    say "macOS: the optional MITOS2 cannot be installed (no macOS build of viennarna<2)"
    ;;
  MINGW*|MSYS*|CYGWIN*)
    echo "Windows: please use WSL2 (wsl --install -d Ubuntu, then run this script inside Ubuntu) or Docker." >&2
    exit 1 ;;
  *) echo "unsupported system: $OS" >&2; exit 1 ;;
esac

# ── conda / mamba (Miniforge installed into ~/miniforge3 when none is present) ─
SOLVER="$(command -v mamba || command -v micromamba || command -v conda || true)"
if [ -z "$SOLVER" ]; then
  for d in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3" /opt/conda; do
    for x in mamba conda; do [ -z "$SOLVER" ] && [ -x "$d/bin/$x" ] && SOLVER="$d/bin/$x"; done
  done
fi
if [ -z "$SOLVER" ]; then
  say "no conda/mamba found: installing Miniforge into $HOME/miniforge3"
  MF_OS="$OS"; MF_ARCH="$ARCH"; [ -n "${CONDA_SUBDIR:-}" ] && MF_ARCH="x86_64"
  curl -fsSL -o /tmp/miniforge.sh \
    "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-${MF_OS}-${MF_ARCH}.sh"
  bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"
  rm -f /tmp/miniforge.sh
  SOLVER="$HOME/miniforge3/bin/mamba"
fi
say "using $SOLVER"

env_exists() { "$SOLVER" env list 2>/dev/null | awk '{print $1}' | grep -qx "$1"; }

for e in $ENVS; do
  if env_exists "$e" && [ "$FORCE" = 0 ]; then
    say "environment '$e' exists - updating it from env/$e.yml (use --force to recreate)"
    "$SOLVER" env update -n "$e" -f "$HERE/env/$e.yml"
    continue
  fi
  say "creating environment '$e' from env/$e.yml"
  if env_exists "$e"; then "$SOLVER" env remove -y -n "$e"; fi
  "$SOLVER" env create -y -n "$e" -f "$HERE/env/$e.yml" 2>/dev/null || "$SOLVER" env create -n "$e" -f "$HERE/env/$e.yml"
done

# ── Java for MiTFi (its classes are Java 5: any Java runs it): the system one, else openjdk 8 in the environment ──
if command -v java >/dev/null 2>&1 && java -version >/dev/null 2>&1; then
  say "Java found ($(command -v java)): used for MiTFi"
else
  say "no Java found: adding openjdk 8 to mitotools (MiTFi)"
  "$SOLVER" install -y -n mitotools -c conda-forge "openjdk=8"
fi

# ── the mitom command: runs on the environment's Python (in the environment's bin/ and in ~/.local/bin) ──
ENV_PREFIX="$("$SOLVER" run -n mitotools python -c 'import sys; print(sys.prefix)' | tail -n 1)"
PY="$ENV_PREFIX/bin/python"
chmod +x "$HERE/mitom"
for d in "$ENV_PREFIX/bin" "$HOME/.local/bin"; do
  mkdir -p "$d"
  rm -f "$d/mitom"                       # an older install left a symlink to $HERE/mitom: never write through it
  printf '#!/usr/bin/env bash\nexec "%s" "%s/mitom" "$@"\n' "$PY" "$HERE" > "$d/mitom"
  chmod +x "$d/mitom"
done
case ":$PATH:" in
  *":$HOME/.local/bin:"*) say "command installed: mitom" ;;
  *) say "command installed: $HOME/.local/bin/mitom - add it to PATH: echo 'export PATH=\$HOME/.local/bin:\$PATH' >> ~/.bashrc   (zsh, the macOS default shell: ~/.zshrc instead of ~/.bashrc)" ;;
esac

# ── MITOS2 reference data and the config file (what was found): mitom setup ──
say "mitom setup: MITOS2 reference data in $REFDIR, config file"
"$PY" "$HERE/mito_env.py" setup --refdir "$REFDIR"
