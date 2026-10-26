#!/bin/bash
# Final comparison after the --max pan rescue (rescue_max.py) and the divergent-clade fold (hardx).
# A = classifier on excluded-family genomes + shuffled negatives (w 0.1); B = A + hardx genomes.
set -e
cd "$(dirname "$0")/.."
L=scripts/learn.py; B=scripts/bounds.py; E=scripts/evaluate.py
[ -d scan/mx_hardx ] || python3 scripts/rescue_max.py data scan/hardx scan/mx_hardx --jobs 12
f() { [ -s feat/$2.tsv ] || python3 $L features scan/$1 $3 feat/$2.tsv --ext $4; }
f mx_excl mx_excl data baselines/fast_excl & f mx_hardx mx_hardx data baselines/fast_hardx &
f mx_shufx mx_shufx data_shufx baselines/fast_shufx & f mx_test mx_test data baselines/fast &
f mx_shuf mx_shuf data_shuf baselines/fast_shuf & f mx_hard mx_hard data baselines/fast_hard & wait
head -1 feat/mx_excl.tsv > feat/mx_exclhx.tsv; tail -n +2 -q feat/mx_excl.tsv feat/mx_hardx.tsv >> feat/mx_exclhx.tsv
python3 $L train feat/mx_excl.tsv ml/mxA.pkl --neg feat/mx_shufx.tsv
python3 $L train feat/mx_exclhx.tsv ml/mxB.pkl --neg feat/mx_shufx.tsv
python3 $B pairs feat/mx_exclhx.tsv data feat/mx_pairs_train.tsv --positives_only
python3 $B train feat/mx_pairs_train.tsv ml/mx_bounds.pkl
for m in mxA mxB; do
  for s in test hard; do
    python3 $L apply ml/$m.pkl feat/mx_$s.tsv preds/${m}_${s}_call.tsv
    python3 $B pairs feat/mx_$s.tsv data feat/${m}_pairs_$s.tsv --keys preds/${m}_${s}_call.tsv
    python3 $B apply ml/mx_bounds.pkl feat/${m}_pairs_$s.tsv preds/${m}_${s}_call.tsv preds/${m}_$s.tsv
  done
  python3 $L apply ml/$m.pkl feat/mx_shuf.tsv preds/${m}_shuf.tsv
  echo "$m shuffled calls: $(( $(wc -l < preds/${m}_shuf.tsv) - 1 )) / 150"
  python3 $E data preds/${m}_test.tsv --accs test_accs.txt --groups
  python3 $E data preds/${m}_hard.tsv --accs hard_accs.txt
done
