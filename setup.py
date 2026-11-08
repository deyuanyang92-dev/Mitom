"""Installs every top-level module (pyproject.toml has the rest of the metadata): a new module needs no entry."""
from pathlib import Path

from setuptools import setup

HERE = Path(__file__).resolve().parent
if __name__ == "__main__":                      # (pip runs it so; importing it - install tests - does nothing)
    setup(py_modules=sorted(p.stem for p in HERE.glob("*.py") if p.stem != "setup"))
