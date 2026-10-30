"""mitom_log.py - what MitoM shows and logs while it runs: the run log (setup_root_logger: everything goes to the log
file, the console gets the numbered steps), numbered steps and progress lines on the terminal (Steps, Progress, say),
one log file per sample of the annotation stage.

Lines of the run log that `mitom debug` reads (their form stays as it is):
  [STEP 6/11 add_missing] start: <title>            a step starts (Steps.start)
  [STEP 6/11 add_missing] introns done 5.4s: ...    a part of a step (Steps.part)
  [STEP 6/11 add_missing] done 9.2s: <result>       ... ends ('failed' instead of 'done': it could not do its work)
  [CMD rc=1 12.3s] <command line>  (in <folder>)    every external program (mito_common.run_tool); when it
      | <a line of its stderr>                      failed, the end of its stderr below (debug mode: always)
  [RECORD <name>] <what happened to this record>    a record that failed or was left out (record_event)
  [ENV] ...                                         the environment (programs, versions with --debug)
  [ERROR] <exception> + traceback                   the run stopped on an error
  [DONE] exit=<code> <time>                         the run ended (mitom_pipeline.main)
"""
from __future__ import annotations

import logging
import threading
import time
from pathlib import Path
from typing import Callable, Optional


_LOG_LOCK = threading.Lock()


class _LockedHandler(logging.Handler):
    def __init__(self, inner: logging.Handler):
        super().__init__()
        self._inner = inner
        self.setFormatter(inner.formatter)
        self.setLevel(inner.level)

    def emit(self, record: logging.LogRecord) -> None:
        with _LOG_LOCK:
            self._inner.emit(record)


# libraries whose INFO lines are noise in run.log (fontTools: ~40 lines on the font subset of every PDF map)
QUIET_LOGGERS = ("fontTools", "matplotlib", "PIL")


def setup_root_logger(log_file: Path, console_level: int = logging.INFO, mode: str = "a", debug: bool = False) -> None:
    """Everything goes to log_file (mode "w": a new log for this run), every external program too ([CMD] lines of
    'mitom.tools'); the console gets console_level and above (run: WARNING, its steps are printed by Steps / say());
    QUIET_LOGGERS only from WARNING. debug: DEBUG lines of every module in log_file."""
    log_file.parent.mkdir(parents=True, exist_ok=True)
    root = logging.getLogger()
    root.setLevel(logging.DEBUG if debug else logging.INFO)
    tools = logging.getLogger("mitom.tools")    # the programs run: in the log file, not on the console
    for h in list(tools.handlers):              # (mito_common's stderr handler of a step run in debug mode)
        tools.removeHandler(h)
    tools.setLevel(logging.DEBUG)
    tools.propagate = True
    for h in list(root.handlers):
        root.removeHandler(h)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    fh = logging.FileHandler(str(log_file), mode=mode, encoding="utf-8")
    fh.setFormatter(fmt)
    sh = logging.StreamHandler()
    sh.setFormatter(logging.Formatter("    ! %(message)s") if console_level > logging.INFO else fmt)
    sh.setLevel(console_level)
    root.addHandler(_LockedHandler(fh))
    root.addHandler(_LockedHandler(sh))
    for name in QUIET_LOGGERS:
        logging.getLogger(name).setLevel(logging.WARNING)


# ── what the user sees on the terminal during `run` (GetOrganelle style: numbered steps, counts, time, files) ──
def say(text: str = "") -> None:
    """One line on the terminal and in the run log."""
    print(text, flush=True)
    if text.strip():
        logging.info(text.strip())


def duration(sec: float) -> str:
    sec = int(round(sec))
    if sec < 60:
        return f"{sec} s"
    if sec < 3600:
        return f"{sec // 60} min {sec % 60:02d} s"
    return f"{sec // 3600} h {sec % 3600 // 60:02d} min"


class Steps:
    """Numbered steps: start() prints '[Step n/N] title', info() an indented line, done() the result and the time; the
    run log gets '[STEP n/N key] start: title' and '[STEP n/N key] done|failed <seconds>s: result' (see above)."""
    current: Optional["Steps"] = None    # the steps of the running workflow: main() logs the step an error stopped

    def __init__(self, total: int):
        self.total, self.n, self.t0, self.run_t0 = total, 0, time.time(), time.time()
        self.key = ""                    # the running step ('' between steps)
        Steps.current = self

    def tag(self) -> str:
        return f"[STEP {self.n}/{self.total} {self.key}]"

    def start(self, key: str, title: str) -> None:
        self.n += 1
        self.key, self.t0 = key, time.time()
        print("", flush=True)
        print(f"[Step {self.n}/{self.total}] {title}", flush=True)
        logging.info(f"{self.tag()} start: {title}")

    def info(self, text: str) -> None:
        say(f"    {text}")

    def part(self, name: str, fn: Callable, *args):
        """Run one part of the step (fn returns a StepResult): its time and whether it did its work in the log."""
        t0 = time.time()
        res = fn(*args)
        logging.info(f"{self.tag()} {name} {'done' if getattr(res, 'ok', True) else 'failed'} {time.time() - t0:.1f}s: "
                     f"{getattr(res, 'summary', res)}")
        return res

    def done(self, text: str = "", ok: bool = True) -> None:
        dt = time.time() - self.t0
        print(f"    -> {text}  ({duration(dt)})" if text else f"    -> done ({duration(dt)})", flush=True)
        logging.info(f"{self.tag()} {'done' if ok else 'failed'} {dt:.1f}s" + (f": {text}" if text else ""))
        self.key = ""

    def failed(self, exc: BaseException) -> None:
        """The running step stopped on exc (logged by main() before the traceback)."""
        if self.key:
            logging.error(f"{self.tag()} failed {time.time() - self.t0:.1f}s: {type(exc).__name__}: {exc}")
            self.key = ""


def record_event(name: str, text: str, level: int = logging.WARNING) -> None:
    """One line about one record, as `mitom debug` reads it: '[RECORD <name>] <text>' (a record that failed or was
    left out of a step)."""
    logging.log(level, f"[RECORD {name}] {text}")


class Progress:
    """'    label: 120/636 (19%), 5 min 02 s so far, about 18 min left' - at most every `every` seconds, and at the end."""

    def __init__(self, label: str, total: int, every: float = 30.0):
        self.label, self.total, self.every = label, max(1, total), every
        self.t0 = self.last = time.time()
        self.shown = -1
        self.lock = threading.Lock()

    def update(self, done: int) -> None:
        with self.lock:
            now = time.time()
            if done == self.shown or (done < self.total and now - self.last < self.every):
                return
            self.last, self.shown = now, done
        el = now - self.t0
        left = f", about {duration(el / done * (self.total - done))} left" if 0 < done < self.total else ""
        say(f"    {self.label}: {done}/{self.total} ({100 * done // self.total}%), {duration(el)} so far{left}")


def watch_files(pattern_dir: Path, pattern: str, progress: Progress, stop: threading.Event) -> threading.Thread:
    """Background thread: progress.update(number of files matching pattern) every 10 s until stop is set."""
    def loop() -> None:
        while not stop.wait(10):
            progress.update(len(list(pattern_dir.glob(pattern))))
    t = threading.Thread(target=loop, daemon=True)
    t.start()
    return t


def get_sample_logger(name: str, log_file: Path) -> logging.Logger:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    lg = logging.getLogger(name)
    lg.setLevel(logging.INFO)
    lg.propagate = False
    for h in list(lg.handlers):
        try:
            h.close()
        except Exception:
            pass
        lg.removeHandler(h)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    fh = logging.FileHandler(str(log_file), encoding="utf-8")
    fh.setFormatter(fmt)
    lg.addHandler(fh)
    return lg


def close_sample_logger(name: str) -> None:
    lg = logging.getLogger(name)
    for h in list(lg.handlers):
        try:
            h.close()
        except Exception:
            pass
        lg.removeHandler(h)


def tail_text(s: str, n_lines: int = 120) -> str:
    return "\n".join((s or "").splitlines()[-n_lines:])
