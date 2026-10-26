#!/bin/bash
# run_ml.sh <tag> <train_hits_dir> <test_hits_dir> [threshold]
# features -> call classifier -> boundary model -> evaluation (train = excluded families, test = test split)
set -e
cd "$(dirname "$0")/.."
tag=$1; tr=$2; te=$3; thr=${4:-0.5}
python3 scripts/learn.py features "$tr" data feat/${tag}_excl.tsv --ext baselines/fast_excl
python3 scripts/learn.py features "$te" data feat/${tag}_test.tsv --ext baselines/fast
python3 scripts/learn.py train feat/${tag}_excl.tsv ml/${tag}.pkl
python3 scripts/learn.py apply ml/${tag}.pkl feat/${tag}_test.tsv preds/${tag}_call.tsv --threshold "$thr"
python3 scripts/bounds.py pairs feat/${tag}_excl.tsv data feat/${tag}_pairs_excl.tsv --positives_only
python3 scripts/bounds.py pairs feat/${tag}_test.tsv data feat/${tag}_pairs_test.tsv --keys preds/${tag}_call.tsv
python3 scripts/bounds.py train feat/${tag}_pairs_excl.tsv ml/${tag}_bounds.pkl
python3 scripts/bounds.py apply ml/${tag}_bounds.pkl feat/${tag}_pairs_test.tsv preds/${tag}_call.tsv preds/${tag}.tsv
python3 scripts/evaluate.py data preds/${tag}.tsv --accs test_accs.txt --groups
