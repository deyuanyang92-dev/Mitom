"""Tests for mito_env.py (reference data download).

Run: python3 -m pytest -q tests/test_mito_env.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import mito_env  # noqa: E402


def test_refseq89m_download_falls_back_to_the_next_source(tmp_path, monkeypatch):
    """Issue #11: Zenodo unreachable -> the mirror; a source with a wrong md5 is skipped too."""
    import hashlib
    import io
    import urllib.request
    good = b"refseq89m stand-in"
    monkeypatch.setattr(mito_env, "REF_MD5", hashlib.md5(good).hexdigest())
    served = {"http://a/zenodo": OSError("unreachable"), "http://b/wrong": b"another file", "http://c/mirror": good}

    def fake(url, timeout=0):
        v = served[url]
        if isinstance(v, Exception):
            raise v
        return io.BytesIO(v)
    monkeypatch.setattr(urllib.request, "urlopen", fake)
    monkeypatch.setattr(mito_env.shutil, "which", lambda name: None)          # (no curl fallback here)
    tgz = tmp_path / "r.tar.bz2"
    assert mito_env.download_ref(tgz, ("http://a/zenodo", "http://b/wrong", "http://c/mirror"))
    assert tgz.read_bytes() == good
    tgz.unlink()
    assert not mito_env.download_ref(tgz, ("http://a/zenodo", "http://b/wrong")) and not tgz.exists()


def test_refseq89m_download_uses_curl_when_urllib_fails(tmp_path, monkeypatch):
    """Issue #11: urllib through the proxy fails the TLS handshake where curl works."""
    import hashlib
    import subprocess
    import urllib.request
    good = b"refseq89m stand-in"
    monkeypatch.setattr(mito_env, "REF_MD5", hashlib.md5(good).hexdigest())
    monkeypatch.setattr(urllib.request, "urlopen", lambda url, timeout=0: (_ for _ in ()).throw(OSError("SSL EOF")))
    monkeypatch.setattr(mito_env.shutil, "which", lambda name: "/usr/bin/curl")

    def fake_run(cmd, **kw):
        Path(cmd[cmd.index("-o") + 1]).write_bytes(good)
        return subprocess.CompletedProcess(cmd, 0)
    monkeypatch.setattr(mito_env.subprocess, "run", fake_run)
    tgz = tmp_path / "r.tar.bz2"
    assert mito_env.download_ref(tgz, ("http://a/x",)) and tgz.read_bytes() == good
