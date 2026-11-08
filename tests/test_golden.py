"""Whole-workflow regression (scripts/golden.py): the same inputs give the same outputs as tests/data/golden.
About 1.5 min with every program installed; run with MITOM_GOLDEN=1 (skipped otherwise, and on another machine than
tests/data/golden/MACHINE.txt, where the reference outputs are first rebuilt with scripts/golden.py --update)."""
import os
import platform
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import golden  # noqa: E402

MACHINE = (golden.GOLDEN / "MACHINE.txt").read_text().strip()


@pytest.mark.skipif(not os.environ.get("MITOM_GOLDEN"), reason="set MITOM_GOLDEN=1 (about 1.5 min)")
@pytest.mark.skipif(platform.node() != MACHINE, reason=f"reference outputs made on {MACHINE}")
def test_outputs_unchanged():
    assert golden.main([]) == 0
