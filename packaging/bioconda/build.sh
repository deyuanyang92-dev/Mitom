#!/bin/bash
# The scripts and models go to $PREFIX/share/mitom; $PREFIX/bin/mitom runs them with the environment's Python.
set -euo pipefail
dest="$PREFIX/share/mitom"
mkdir -p "$dest/scripts" "$PREFIX/bin"
cp mitom *.py "$dest/"
cp scripts/runmitos_fast.py "$dest/scripts/"
cp -r models "$dest/"
rm -rf "$dest/models/trna_ml/clades" "$dest/models/trna_ml/cache"      # unpacked / trained on first use
cat > "$PREFIX/bin/mitom" <<'WRAP'
#!/usr/bin/env bash
d="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
exec "$d/bin/python" "$d/share/mitom/mitom" "$@"
WRAP
chmod +x "$PREFIX/bin/mitom" "$dest/mitom"
