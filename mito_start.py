"""mito_start.py - the start codon of a protein gene as NCBI annotates it, learned from NCBI records (mito_curate P11).

Mitochondrial genes start with ATG and with ATA / ATT / ATC / GTG / TTG ..., often several in frame within a few
codons, and NCBI submitters choose among them by conventions no single rule captures (start nearest the preceding
gene, ATG preferred, the homologue's start ...). P11 scores every candidate - an in-frame start codon of the genetic
code within MAX_SHIFT codons of the start the boundary rules chose (P1), reading frame open in between - with
gradient-boosted trees and takes the best one.

Only genes found by MitoM's core (with its alignment note): the starts of other annotators (MitoZ, MITOS2, mitom
import) are left as they are.

Features of a candidate: its codon; its offset from P1's start and from the annotation core's start; the residues of
the reference protein missing before the core's alignment and the share it covers (from the core's note: a start
the alignment reaches is reliable, one before a divergent N-terminus is a guess); the gap to the end of the preceding
feature, that feature's type and strand, whether the candidate is the first start codon after it; the gene, the
phylum and the genetic code.

The model is plain JSON (models/start_model/start_model.json: columns and the trees' split features, thresholds and
leaf values), evaluated here with numpy: no pickle, nothing executed from the file, no scikit-learn needed to run it.
Training: benchmark/train_start_model.py (scikit-learn GradientBoostingClassifier on mito_curate.py --dump_starts of
the leave-out benchmark records, each candidate labelled by the NCBI record's start).
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

MAX_SHIFT = 30                                # codons searched on each side of P1's start
MODEL = Path(__file__).resolve().parent / "models" / "start_model" / "start_model.json"
NUMERIC = ["cand", "rank", "is_p1", "core_off", "aln_missing5", "aln_cover", "up_gap", "up_same", "first_after_up",
           "code"]
CATEGORICAL = ["codon", "gene", "phylum", "up_type"]
MISSING = -999.0                              # a numeric feature that does not apply (no core note, no neighbour)


def candidates(s: str, a: int, b: int, starts: set, stops: set, core_a: Optional[int],
               aln: Optional[Tuple[int, int, int]], up: Optional[Tuple[int, str, bool]], gene: str, phylum: str,
               code: int) -> List[Dict]:
    """Candidate starts of the CDS a..b (1-based, transcript orientation, s the transcript-oriented sequence) whose
    start P1 put at a: rows with the features of every in-frame start codon within MAX_SHIFT codons of a, reading
    frame open between it and a (a itself is always a candidate). core_a: the annotation core's start; aln: (first
    residue, last residue, reference length) of its alignment; up: (end, type, same strand) of the preceding feature."""
    cands = []
    for c in range(a - 3 * MAX_SHIFT, a + 3 * MAX_SHIFT + 1, 3):
        if c < 1 or c > b - 30 or s[c - 1:c + 2] not in starts and c != a:
            continue
        lo, hi = min(c, a), max(c, a)
        if any(s[k - 1:k + 2] in stops for k in range(lo, hi, 3)):
            continue
        cands.append(c)
    first_after = min((c for c in cands if up is not None and c > up[0]), default=None)
    rows = []
    for c in cands:
        rows.append({"pos": c, "codon": s[c - 1:c + 2], "gene": gene, "phylum": phylum or "?",
                     "up_type": up[1] if up else "none",
                     "cand": c - a, "rank": abs(c - a) // 3, "is_p1": int(c == a),
                     "core_off": MISSING if core_a is None else c - core_a,
                     "aln_missing5": MISSING if not aln else aln[0] - 1,
                     "aln_cover": MISSING if not aln else round((aln[1] - aln[0] + 1) / aln[2], 3),
                     "up_gap": MISSING if up is None else c - up[0] - 1,
                     "up_same": int(bool(up and up[2])), "first_after_up": int(c == first_after), "code": code})
    return rows


def columns(categories: Dict[str, List[str]]) -> List[str]:
    """Model columns: the numeric features, then one 0 / 1 column per category value (codon=ATG ...)."""
    return NUMERIC + [f"{c}={v}" for c in CATEGORICAL for v in categories[c]]


def matrix(rows: Sequence[Dict], categories: Dict[str, List[str]]):
    """Feature matrix in columns() order."""
    import numpy as np
    cols = {k: i for i, k in enumerate(columns(categories))}
    x = np.zeros((len(rows), len(cols)))
    for i, r in enumerate(rows):
        for k in NUMERIC:
            x[i, cols[k]] = float(r[k])
        for c in CATEGORICAL:
            j = cols.get(f"{c}={r[c]}")
            if j is not None:                          # (a value the training never saw: all its columns 0)
                x[i, j] = 1.0
    return x


class StartModel:
    """Gradient-boosted trees from start_model.json: score = sum of the leaf values each tree gives a candidate
    (the ranking needs no link function)."""

    def __init__(self, meta: Dict):
        import numpy as np
        self.meta = meta
        self.trees = [tuple(np.asarray(t[k]) for k in ("left", "right", "feature", "threshold", "value"))
                      for t in meta["trees"]]

    def score(self, rows: Sequence[Dict]) -> List[float]:
        import numpy as np
        x = matrix(rows, self.meta["categories"])
        total = np.zeros(len(rows))
        for left, right, feature, threshold, value in self.trees:
            node = np.zeros(len(rows), dtype=int)
            while True:
                inner = left[node] >= 0
                if not inner.any():
                    break
                go_left = x[np.arange(len(rows)), np.where(inner, feature[node], 0)] <= threshold[node]
                node = np.where(inner, np.where(go_left, left[node], right[node]), node)
            total += value[node]
        return [float(v) for v in total]

    def best(self, rows: Sequence[Dict]) -> Tuple[int, List[float]]:
        """Index of the best-scoring candidate (P1's start on a tie) and every candidate's score."""
        scores = self.score(rows)
        return max(range(len(rows)), key=lambda i: (scores[i], rows[i]["is_p1"])), scores


def load(path: Optional[str] = None) -> Optional[StartModel]:
    """The model in `path` (default models/start_model/start_model.json); None when there is none."""
    p = Path(path) if path else MODEL
    if not p.is_file():
        return None
    return StartModel(json.loads(p.read_text()))
