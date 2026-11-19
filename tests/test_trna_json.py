"""The JSON trees of the tRNA classifiers (mito_gbdt, models/trna_ml/*_model.json) give scikit-learn's probabilities."""
import csv
import gzip
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]

import trna_ml  # noqa: E402


def test_shipped_json_matches_tables():
    for name in trna_ml.PARAMS:
        assert trna_ml._json_model(name, trna_ml.MODELS) is not None, f"{name}_model.json stale: python3 trna_ml.py --train"


def test_json_inference_needs_no_sklearn(monkeypatch):
    monkeypatch.setitem(sys.modules, "sklearn", None)
    clf, bclf = trna_ml.classifiers()
    assert clf.predict_proba([[0.0] * len(trna_ml.FEATS)]).shape == (1, 2)


@pytest.mark.parametrize("name", ["caller", "bounds"])
def test_json_equals_sklearn(name):
    pytest.importorskip("sklearn")
    import numpy as np
    feats, kw = trna_ml.PARAMS[name]
    table = trna_ml.MODELS / f"{name}_train.tsv.gz"
    with gzip.open(table, "rt", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    x = np.array([[float(r[k]) for k in feats] for r in rows])
    sk = trna_ml._fit(table, feats, **kw)
    trees = trna_ml._json_model(name, trna_ml.MODELS)
    assert np.abs(trees.raw(x) - sk.decision_function(x)).max() < 1e-9
    assert np.abs(trees.predict_proba(x) - sk.predict_proba(x)).max() < 1e-9
    # missing values go where scikit-learn sends them
    xm = x[:2000].copy()
    xm[::3, ::5] = np.nan
    assert np.abs(trees.predict_proba(xm) - sk.predict_proba(xm)).max() < 1e-9
