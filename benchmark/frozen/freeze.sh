#!/usr/bin/env bash
# freeze.sh - write the manifest that fixes everything the frozen evaluation depends on (PLAN.md section 3).
# Run it once, after the test set is built and before any tool runs; commit the manifest folder.
#
# Usage: bash benchmark/frozen/freeze.sh <set.tsv> <manifest_dir> [env ...]
#        env = conda environments to record (default: mitotools mitos mitoz mitofinder; missing ones are noted)
set -euo pipefail
SET=$1; OUT=$2; shift 2
ENVS=${*:-"mitotools mitos mitoz mitofinder"}
REPO=$(cd "$(dirname "$0")/../.." && pwd)
mkdir -p "$OUT"
cd "$REPO"
PY=${PYTHON:-python3}                  # the mitotools Python (conda activate mitotools)

if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
    echo "ERROR: uncommitted changes in $REPO - commit (and tag) the version to evaluate first" >&2; exit 1
fi
CONDA=$(command -v mamba || command -v conda || echo "$HOME/miniforge3/bin/mamba")
REFDIR=$("$PY" -c "import sys; sys.path.insert(0, '$REPO'); import mito_env; print(mito_env.refseq_dir() or '')")

{
  echo "date	$(date -Iseconds)"
  echo "host	$(hostname)"
  echo "disk	$(df -P "$REPO" | awk 'NR==2{print $6}')"
  echo "commit	$(git rev-parse HEAD)"
  echo "tag	$(git describe --tags --exact-match 2>/dev/null || echo 'none (tag the commit, e.g. v0.0.<n>-frozen)')"
  echo "set	$SET"
  echo "set_md5	$(md5sum "$SET" | cut -d' ' -f1)"
  echo "set_records	$(($(wc -l < "$SET") - 1))"
  echo "refseq89m	$REFDIR"
  echo "models_md5	$(find models -type f -print0 | sort -z | xargs -0 md5sum | md5sum | cut -d' ' -f1)"
  echo "uname	$(uname -srm)"
  echo "cpus	$(nproc)"
} > "$OUT/manifest.tsv"

# reference data: md5 of every file (refseq89m: MitoM's and MITOS2's protein / RNA references)
if [ -n "$REFDIR" ]; then
    (cd "$REFDIR" && find . -type f -print0 | sort -z | xargs -0 md5sum) > "$OUT/refseq89m.md5"
    echo "refseq89m_md5	$(md5sum < "$OUT/refseq89m.md5" | cut -d' ' -f1)" >> "$OUT/manifest.tsv"
fi
find models -type f -print0 | sort -z | xargs -0 md5sum > "$OUT/models.md5"
cp "$SET" "$OUT/set.tsv"

# programs: mitom doctor (paths and versions) and the exact conda environments
"$PY" "$REPO/mitom_pipeline.py" doctor > "$OUT/doctor.txt" 2>&1 || true
for e in $ENVS; do
    if "$CONDA" env export -n "$e" > "$OUT/env_$e.yml" 2>/dev/null; then
        echo "env_$e	$(md5sum < "$OUT/env_$e.yml" | cut -d' ' -f1)" >> "$OUT/manifest.tsv"
    else
        rm -f "$OUT/env_$e.yml"; echo "env_$e	not installed" >> "$OUT/manifest.tsv"
    fi
done
for tool in runmitos runmitos.py mitoz MitoFinder; do
    exe=$("$PY" -c "import sys; sys.path.insert(0, '$REPO'); import mito_env; print(mito_env.find('$tool') or '')")
    [ -n "$exe" ] && echo "$tool	$exe	$("$exe" --version 2>&1 | head -1)" >> "$OUT/programs.tsv" || true
done
echo "[freeze] manifest -> $OUT (commit it before running: git add $OUT && git commit)"
