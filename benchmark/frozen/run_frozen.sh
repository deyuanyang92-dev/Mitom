#!/usr/bin/env bash
# run_frozen.sh - the one run of the frozen evaluation (PLAN.md section 4): every tool on every test record, scored once.
# Refuses to run when the code is not the manifest's commit or has uncommitted changes. Resumable (records with a
# status.json are skipped), so an interrupted run is continued with the same command, never restarted with other code.
#
# Usage: bash benchmark/frozen/run_frozen.sh <manifest_dir> <data_root> [jobs] [threads]
#   <data_root> holds gb/ and bare/ from make_large_set.py (PLAN.md step 2); results go to <data_root>/runs, eval/all
#   MitoFinder: MITOFINDER=<script> PY2=<python2.7> REF_POOL=<RefSeq pool .gb> (left out when MITOFINDER is unset)
set -euo pipefail
MAN=$1; ROOT=$2; JOBS=${3:-8}; THREADS=${4:-2}
REPO=$(cd "$(dirname "$0")/../.." && pwd)
SET="$MAN/set.tsv"
cd "$REPO"
PY=${PYTHON:-python3}                  # the mitotools Python (conda activate mitotools)

want=$(awk -F'\t' '$1=="commit"{print $2}' "$MAN/manifest.tsv")
[ "$(git rev-parse HEAD)" = "$want" ] || { echo "ERROR: HEAD $(git rev-parse --short HEAD) is not the frozen commit ${want:0:7}" >&2; exit 1; }
[ -z "$(git status --porcelain --untracked-files=no)" ] || { echo "ERROR: uncommitted changes" >&2; exit 1; }
[ "$(md5sum < "$SET" | cut -d' ' -f1)" = "$(awk -F'\t' '$1=="set_md5"{print $2}' "$MAN/manifest.tsv")" ] || { echo "ERROR: set.tsv changed" >&2; exit 1; }

# primary: mitom_native vs mitos2 (both given the NCBI record's genetic code); mitom_blind = sequence only
# (TOOLS=... overrides the list, e.g. without MitoZ on a machine that lacks it; say so in the report)
TOOLS=${TOOLS:-"mitom_native,mitom_blind,mitos2,mitos2_cur,mitoz,mitoz_cur"}
MF=()
if [ -n "${MITOFINDER:-}" ]; then
    TOOLS="$TOOLS,mitofinder"; MF=(--mitofinder "$MITOFINDER" ${PY2:+--py2 "$PY2"} ${REF_POOL:+--ref_pool "$REF_POOL"})
fi
FIRST=$(tr , '\n' <<< "$TOOLS" | grep -v '_cur$' | paste -sd,)
SECOND=$(tr , '\n' <<< "$TOOLS" | grep '_cur$' | paste -sd, || true)
LOG="$MAN/commands.log"
run() { echo "$(date -Iseconds)	$*" >> "$LOG"; "$@"; }

# mitos2 and mitoz run before their *_cur (which import their output)
run "$PY" benchmark/mito_bench_run.py --set "$SET" --root "$ROOT" --split all --jobs "$JOBS" --threads "$THREADS" \
    --timeout 7200 --tools "$FIRST" "${MF[@]}"
[ -z "$SECOND" ] || run "$PY" benchmark/mito_bench_run.py --set "$SET" --root "$ROOT" --split all --jobs "$JOBS" \
    --threads "$THREADS" --timeout 7200 --tools "$SECOND"
run "$PY" benchmark/mito_bench_eval.py --set "$SET" --root "$ROOT" --split all --tools "$TOOLS" --table2asn
run "$PY" benchmark/bench_report.py --eval "$ROOT/eval/all" --tools "$TOOLS" --set "$SET" --out "$MAN/report.md"
run "$PY" benchmark/frozen/frozen_stats.py --eval "$ROOT/eval/all" --set "$SET" --tools "$TOOLS" --ref mitos2 \
    -o "$MAN/stats.tsv"
cp "$ROOT"/eval/all/{summary,errors,timing}.tsv "$MAN/"
gzip -c "$ROOT/eval/all/features.tsv" > "$MAN/features.tsv.gz"
echo "[frozen] done: $MAN/report.md, $MAN/stats.tsv (commit $MAN; add one row to knowledge/TEST_REGISTRY.md)"
