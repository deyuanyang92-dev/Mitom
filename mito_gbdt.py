"""mito_gbdt.py - gradient-boosted trees as plain JSON, evaluated with numpy only (no scikit-learn, no pickle).

Shared by the tRNA classifiers (trna_ml.py) and the start-codon model (mito_start.py). A scikit-learn
HistGradientBoostingClassifier is read once, at training (export_hgb), into per-tree arrays: left / right child
(-1 at a leaf), split feature, threshold, leaf value, and whether a missing value goes left. Evaluation walks the same
trees: raw score = baseline + sum of the leaf values; the probability of a binary model is the logistic function of the
raw score (scikit-learn's own link). Numerical features only (what both models use).
"""
from __future__ import annotations

from typing import Dict, List, Sequence


def export_hgb(model) -> Dict:
    """{baseline, trees: [{left, right, feature, threshold, value, missing_left}]} of a fitted
    HistGradientBoostingClassifier (binary: one tree per iteration). Training time only: reads scikit-learn internals."""
    import numpy as np
    base = np.asarray(model._baseline_prediction, dtype=float).ravel()
    if len(base) != 1 or any(len(p) != 1 for p in model._predictors):
        raise ValueError("only binary HistGradientBoostingClassifier models are exported")
    trees = []
    for (pred,) in model._predictors:
        nodes = pred.nodes
        if pred.is_categorical_.any() if hasattr(pred, "is_categorical_") else False:
            raise ValueError("categorical splits are not exported")
        leaf = nodes["is_leaf"].astype(bool)
        trees.append({"left": np.where(leaf, -1, nodes["left"].astype(np.int64)).tolist(),
                      "right": np.where(leaf, -1, nodes["right"].astype(np.int64)).tolist(),
                      "feature": nodes["feature_idx"].astype(int).tolist(),
                      "threshold": [float(x) for x in nodes["num_threshold"]],
                      "value": [float(x) for x in nodes["value"]],
                      "missing_left": nodes["missing_go_to_left"].astype(int).tolist()})
    return {"baseline": float(base[0]), "trees": trees}


class Trees:
    """The exported trees of a model; raw() and predict_proba() give what scikit-learn's decision_function /
    predict_proba give (to rounding: the sums are added in the same order, tree by tree)."""

    def __init__(self, meta: Dict):
        import numpy as np
        self.baseline = float(meta.get("baseline", 0.0))
        self.trees = []
        for t in meta["trees"]:
            ml = t.get("missing_left")
            self.trees.append((np.asarray(t["left"]), np.asarray(t["right"]), np.asarray(t["feature"]),
                               np.asarray(t["threshold"], dtype=float), np.asarray(t["value"], dtype=float),
                               None if ml is None else np.asarray(ml, dtype=bool)))

    def raw(self, x):
        """Raw scores (log-odds) of the rows of matrix x (n x features)."""
        import numpy as np
        x = np.asarray(x, dtype=float)
        n = len(x)
        total = np.full(n, self.baseline)
        rows = np.arange(n)
        for left, right, feature, threshold, value, miss_left in self.trees:
            node = np.zeros(n, dtype=int)
            while True:
                inner = left[node] >= 0
                if not inner.any():
                    break
                v = x[rows, np.where(inner, feature[node], 0)]
                go_left = v <= threshold[node]
                if miss_left is not None:
                    go_left = np.where(np.isnan(v), miss_left[node], go_left)
                node = np.where(inner, np.where(go_left, left[node], right[node]), node)
            total += value[node]
        return total

    def predict_proba(self, rows: Sequence[Sequence[float]]):
        """n x 2 array [P(0), P(1)] like scikit-learn's predict_proba."""
        import numpy as np
        if len(rows) == 0:
            return np.zeros((0, 2))
        p = 1.0 / (1.0 + np.exp(-self.raw(rows)))
        return np.column_stack([1.0 - p, p])
