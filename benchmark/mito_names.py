"""Moved to the top-level mito_names.py (shipped with MitoM); this name keeps the benchmark scripts that import
mito_names from benchmark/ working: importing it gives the top-level module."""
import importlib.util
import sys
from pathlib import Path

_spec = importlib.util.spec_from_file_location("mito_names", Path(__file__).resolve().parents[1] / "mito_names.py")
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
sys.modules[__name__] = _mod
