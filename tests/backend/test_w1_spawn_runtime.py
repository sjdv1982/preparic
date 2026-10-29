"""W1 execution acceptance using real workers, with TCP services unavailable.

Run in a fresh interpreter so earlier tests cannot supply an existing pool,
configured remote clients, or in-memory transformation results.
"""

import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


@pytest.mark.integration
def test_backend_runs_real_eight_worker_pool_without_tcp_services(tmp_path):
    repo = Path(__file__).resolve().parents[2]
    denied = tmp_path / "denied-tcp.jsonl"
    # sitecustomize is loaded by the fresh interpreter and spawned workers.
    # AF_UNIX multiprocessing transport remains available; hashserver/database
    # TCP connections fail immediately and are recorded even if caught upstream.
    (tmp_path / "sitecustomize.py").write_text('''
import json
import os
import socket

_connect = socket.socket.connect
_connect_ex = socket.socket.connect_ex

def _deny(sock, address):
    if sock.family in (socket.AF_INET, socket.AF_INET6):
        with open(os.environ["W1_DENIED_TCP_LOG"], "a", encoding="utf-8") as stream:
            stream.write(json.dumps({"pid": os.getpid(), "address": address}) + "\\n")
        raise OSError("TCP services deliberately unavailable in W1 acceptance")

def connect(sock, address):
    _deny(sock, address)
    return _connect(sock, address)

def connect_ex(sock, address):
    _deny(sock, address)
    return _connect_ex(sock, address)

socket.socket.connect = connect
socket.socket.connect_ex = connect_ex
''', encoding="utf-8")
    env = dict(os.environ)
    for name in ("SEAMLESS_CACHE", "SEAMLESS_REMOTE_CLIENTS", "SEAMLESS_DASK_SCHEDULER"):
        env.pop(name, None)
    env["W1_DENIED_TCP_LOG"] = str(denied)
    env["PYTHONPATH"] = os.pathsep.join((str(tmp_path), str(repo), env.get("PYTHONPATH", "")))
    completed = subprocess.run(
        [sys.executable, str(Path(__file__).with_name("w1_spawn_runtime_probe.py"))],
        cwd=tmp_path, env=env, capture_output=True, text=True, timeout=120,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    payloads = [line.removeprefix("W1_SPAWN_RESULT=") for line in completed.stdout.splitlines()
                if line.startswith("W1_SPAWN_RESULT=")]
    assert len(payloads) == 1, completed.stdout + completed.stderr
    result = json.loads(payloads[0])
    assert result["configured_workers"] == result["worker_handles"] == 8
    assert sorted(item["index"] for item in result["results"]) == list(range(16))
    assert all(item["pid"] != result["parent_pid"] for item in result["results"])
    assert all(item["tcp_disabled"] for item in result["results"])
    assert result["clients"] == {"database": [], "buffer": []}
    assert result["persistent_cache"] is None
    assert not denied.exists() or not denied.read_text(encoding="utf-8"), "execution attempted a TCP service connection"
