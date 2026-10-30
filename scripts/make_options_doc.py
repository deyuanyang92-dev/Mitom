#!/usr/bin/env python3
"""make_options_doc.py - docs/OPTIONS.md: every MitoM command with all its options, as `mitom <command> -h` prints
them. Run it after changing a command line (python3 scripts/make_options_doc.py); tests/test_docs.py fails while an
option is missing from docs/OPTIONS.md."""
import importlib.machinery
import importlib.util
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs" / "OPTIONS.md"
# `mitom doctor` takes no -h (it runs the check): its one option, by hand
DOCTOR = ("usage: mitom doctor [--write-config]\n\n"
          "Which external programs MitoM finds (config, PATH, conda environments): OK / -- (optional, missing) /\n"
          "!! (required, missing), with what each one is used for.\n\n"
          "options:\n"
          "  --write-config  save the locations found to mito_config.json (mitom setup does it too)\n")


def commands():
    """{command: (script, args)} of the mitom launcher, in its order."""
    loader = importlib.machinery.SourceFileLoader("mitom_launcher", str(ROOT / "mitom"))
    spec = importlib.util.spec_from_loader("mitom_launcher", loader)
    mod = importlib.util.module_from_spec(spec)
    loader.exec_module(mod)
    return mod.COMMANDS


def help_text(cmd: str) -> str:
    """`mitom <cmd> -h`, wrapped at 110 columns whatever the terminal."""
    if cmd == "doctor":
        return DOCTOR
    env = {**os.environ, "COLUMNS": "110"}
    p = subprocess.run([sys.executable, str(ROOT / "mitom"), cmd, "-h"], capture_output=True, text=True, env=env,
                       timeout=300)
    return (p.stdout or p.stderr).rstrip() + "\n"


def main() -> int:
    cmds = commands()
    with ThreadPoolExecutor(8) as ex:
        helps = dict(zip(cmds, ex.map(help_text, cmds)))
    lines = ["# MitoM 全部命令与参数 / All commands and options", "",
             "由 `python3 scripts/make_options_doc.py` 从每个命令的 `-h` 生成（命令行改动后重新生成；"
             "`tests/test_docs.py` 检查这里有没有漏掉参数）。各命令的用法、示例和输出见 [USAGE_zh.md](USAGE_zh.md)。",
             "",
             "Generated from each command's `-h` by `scripts/make_options_doc.py`; usage, examples and outputs: "
             "[USAGE_zh.md](USAGE_zh.md).", "",
             "| 命令 / command | 程序 / script |", "|---|---|"]
    lines += [f"| [`mitom {c}`](#mitom-{c}) | `{s}{' ' + ' '.join(a) if a else ''}` |" for c, (s, a) in cmds.items()]
    for c in cmds:
        lines += ["", f"## mitom {c}", "", "```text", helps[c].rstrip(), "```"]
    OUT.write_text("\n".join(lines) + "\n")
    print(f"{OUT}: {len(cmds)} commands")
    return 0


if __name__ == "__main__":
    sys.exit(main())
