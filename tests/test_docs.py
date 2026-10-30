"""The documentation keeps up with the command lines: every option of every mitom command is in docs/OPTIONS.md
(regenerate it with python3 scripts/make_options_doc.py), and the manual has a section for every command.

Run: python3 -m pytest -q tests/test_docs.py
"""
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
import make_options_doc as mod  # noqa: E402

OPTION = re.compile(r"(?<![\w-])(--[a-z][a-z0-9_-]+)")


def sections(text: str) -> dict:
    return {m.group(1): m.group(2) for m in re.finditer(r"^## mitom (\S+)\n(.*?)(?=^## mitom |\Z)", text, re.M | re.S)}


def test_options_doc_has_every_option():
    cmds = mod.commands()
    with ThreadPoolExecutor(8) as ex:
        helps = dict(zip(cmds, ex.map(mod.help_text, cmds)))
    doc = sections((ROOT / "docs" / "OPTIONS.md").read_text())
    for c, text in helps.items():
        assert c in doc, f"mitom {c} not in docs/OPTIONS.md: python3 scripts/make_options_doc.py"
        missing = set(OPTION.findall(text)) - set(OPTION.findall(doc[c]))
        assert not missing, f"mitom {c}: {sorted(missing)} not in docs/OPTIONS.md: python3 scripts/make_options_doc.py"


def test_manual_has_every_command():
    usage = (ROOT / "docs" / "USAGE_zh.md").read_text()
    for c in mod.commands():
        assert f"`mitom {c}`" in usage, f"mitom {c} not in docs/USAGE_zh.md"
