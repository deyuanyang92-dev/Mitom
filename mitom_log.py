"""mitom_log.py - what MitoM shows and logs while it runs: the run log (setup_root_logger: everything goes to the log
file, the console gets the numbered steps), numbered steps and progress lines on the terminal (Steps, Progress, say),
one log file per sample of the annotation stage."""
from __future__ import annotations

import logging
import threading
import time
from pathlib import Path


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


def setup_root_logger(log_file: Path, console_level: int = logging.INFO, mode: str = "a") -> None:
    """Everything goes to log_file (mode "w": a new log for this run); the console gets console_level and above
    (run: WARNING, its steps are printed by Steps / say())."""
    log_file.parent.mkdir(parents=True, exist_ok=True)
    root = logging.getLogger()
    root.setLevel(logging.INFO)
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
    """Numbered steps: start() prints '[Step n/N] title', info() an indented line, done() the result and the time."""

    def __init__(self, total: int):
        self.total, self.n, self.t0, self.run_t0 = total, 0, time.time(), time.time()

    def start(self, title: str) -> None:
        self.n += 1
        self.t0 = time.time()
        say("")
        say(f"[Step {self.n}/{self.total}] {title}")

    def info(self, text: str) -> None:
        say(f"    {text}")

    def done(self, text: str = "") -> None:
        say(f"    -> {text}  ({duration(time.time() - self.t0)})" if text else f"    -> done ({duration(time.time() - self.t0)})")


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
