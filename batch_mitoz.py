#!/usr/bin/env python3
"""Old name of mitom_pipeline.py (the MitoM workflow; use the `mitom` command). Kept so that old scripts keep working:
`python3 batch_mitoz.py run ...` runs mitom_pipeline.py, `import batch_mitoz` gives the mitom_pipeline module."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mitom_pipeline  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(mitom_pipeline.main())
sys.modules[__name__] = mitom_pipeline
