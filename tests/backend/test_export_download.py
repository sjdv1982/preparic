"""E0: real Context/HTTP checks for the export download integration contract.

These probes intentionally do not import Preparic's export implementation.
They establish the Seamless assumptions needed before wiring that graph.
The actual Vite proxy and preparic-start/preparic-stop restart checks belong
to the live E5/E7 checks; a direct share-server URL does not prove the proxy.
"""

from concurrent.futures import ThreadPoolExecutor
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import threading
import time
import uuid
import zipfile

import pytest
import requests

from seamless import Cell
from seamless_workflow import Context, shareserver


pytestmark = pytest.mark.integration


def _probe_tree(decision, validated_pairs, validated_classes, validated_energies):
    if any(value is not None for value in
           (validated_pairs, validated_classes, validated_energies)):
        raise ValueError("validator gates must deliver successful null")
    if decision < 0:
        raise ValueError("E0 deliberate producer failure")
    return {
        "DATA/reference.pdb": b"ATOM exact bytes\r\n\x00\xff\n",
        "DATA/export-manifest.json": b'{"schema":1}\n',
        "nested/deeper/value.txt": str(decision).encode("ascii") + b"\n",
    }


def _copy_folder(value):
    return value


def _slow_probe_tree(decision, validated_pairs, validated_classes, validated_energies):
    # Gives the HTTP probe a deterministic window to observe an in-flight run.
    import time as _time
    _time.sleep(1.5)
    if any(value is not None for value in
           (validated_pairs, validated_classes, validated_energies)):
        raise ValueError("validator gates must deliver successful null")
    if decision < 0:
        raise ValueError("E0 deliberate producer failure")
    return {
        "DATA/reference.pdb": b"ATOM exact bytes\r\n\x00\xff\n",
        "DATA/export-manifest.json": b'{"schema":1}\n',
        "nested/deeper/value.txt": str(decision).encode("ascii") + b"\n",
    }


def _expected(decision):
    return {"alaric_project/" + name: value for name, value in
            _probe_tree(decision, None, None, None).items()}


def _members(response):
    with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
        return {name: archive.read(name) for name in archive.namelist()}


def _get(url):
    return requests.get(url, timeout=20, headers={"Cache-Control": "no-cache"})


def _build_context(namespace, decision=1):
    context = Context()
    context.shares.namespace = namespace
    context.decision = Cell(celltype="int")
    context.decision.set(decision)
    context.decision.share(readonly=False)
    for name in ("validated_pairs", "validated_classes", "validated_energies"):
        gate = Cell(celltype="plain")
        gate.set(None)
        setattr(context, name, gate)
    context.render = _probe_tree
    context.render.pins.decision = context.decision
    context.render.celltypes.decision = "int"
    for name in ("validated_pairs", "validated_classes", "validated_energies"):
        setattr(context.render.pins, name, getattr(context, name))
        setattr(context.render.celltypes, name, "plain")
    context.render.celltypes.result = "folder"
    context.export_tree = context.render.result
    url = context.shares.add_download(
        {"alaric_project": "export_tree"}, path="export.zip", format="zip")
    return context, url


@pytest.fixture(scope="module")
def download_server(seamless_runtime):
    shareserver.close()
    shareserver.configure(host="127.0.0.1", port=0)
    yield
    shareserver.close()


@pytest.fixture
def probe(download_server):
    context, url = _build_context("e0-" + uuid.uuid4().hex)
    try:
        yield context, url
    finally:
        context.close()


def test_e0_folder_zip_exact_bytes_and_required_null_plain_pins(probe):
    context, url = probe
    context.compute(timeout=30)
    assert context.render.state == "complete", context.render.exception
    for name in ("validated_pairs", "validated_classes", "validated_energies"):
        assert getattr(context.render.pins, name).state == "complete"
        assert getattr(context, name).value is None
    before = context.decision.checksum, context.export_tree.checksum
    response = _get(url)
    assert response.status_code == 200, response.text
    assert response.headers["Content-Type"] == "application/zip"
    assert response.headers["Content-Disposition"] == 'attachment; filename="export.zip"'
    assert response.headers["Cache-Control"] == "no-cache"
    assert _members(response) == _expected(1)
    assert _get(url).content == response.content
    assert (context.decision.checksum, context.export_tree.checksum) == before


def test_e0_changed_put_cannot_serve_old_tree_and_equal_put_keeps_tree(probe):
    context, url = probe
    context.compute(timeout=30)
    old = _get(url)
    assert old.status_code == 200
    put = requests.put(context.decision.share.url, data=b"2", timeout=20)
    assert put.status_code == 200, put.text
    # No compute/capture call is allowed between the accepted write and GET.
    # A very fast producer may already finish; then only the new cut is legal.
    immediate = _get(url)
    if immediate.status_code == 200:
        assert _members(immediate) == _expected(2)
    else:
        assert immediate.status_code == 409, immediate.text
        assert immediate.json()["cells"]["export_tree"] in {
            "waiting", "computing", "blocked"}
    context.compute(timeout=30)
    current = _get(url)
    assert current.status_code == 200, current.text
    assert _members(current) == _expected(2)
    equal = requests.put(context.decision.share.url, data=b"2", timeout=20)
    assert equal.status_code == 200, equal.text
    assert equal.json()["marker"] == put.json()["marker"]
    assert _get(url).content == current.content


def test_e0_incomplete_download_refuses_instead_of_empty_zip(probe):
    context, url = probe
    # Make the required pin state deterministic. A graph with all inputs
    # available may finish before the first GET, even without an explicit
    # compute call, so an unscheduled transformation is not a reliable way to
    # exercise the non-complete download response.
    context.unwired_gate = Cell(celltype="plain")
    context.render.pins.validated_pairs = context.unwired_gate
    context.compute(timeout=30)
    assert context.export_tree.state == "blocked"
    response = _get(url)
    assert response.status_code == 409, response.text
    assert response.json()["cells"]["export_tree"] == "blocked"


def test_e0_computing_download_refuses_until_tree_completes(probe):
    context, url = probe
    context.compute(timeout=30)
    context.render.code = _slow_probe_tree
    with ThreadPoolExecutor(max_workers=1) as pool:
        computation = pool.submit(context.compute, timeout=30)
        deadline = time.monotonic() + 10
        while context.export_tree.state not in {"waiting", "computing"}:
            assert not computation.done(), "slow producer completed before becoming pending"
            assert time.monotonic() < deadline, context.export_tree.state
            time.sleep(0.01)
        response = _get(url)
        assert response.status_code == 409, response.text
        assert response.json()["cells"]["export_tree"] in {"waiting", "computing"}
        computation.result(timeout=30)
    assert context.render.state == "complete", context.render.exception
    complete = _get(url)
    assert complete.status_code == 200, complete.text
    assert _members(complete) == _expected(1)


def test_e0_failed_and_blocked_producers_refuse_download(download_server):
    context, url = _build_context("e0-failure-" + uuid.uuid4().hex, decision=-1)
    try:
        context.copy = _copy_folder
        context.copy.pins.value = context.export_tree
        context.copy.celltypes.value = "folder"
        context.copy.celltypes.result = "folder"
        context.blocked_tree = context.copy.result
        blocked_url = context.shares.add_download(
            {"alaric_project": "blocked_tree"}, path="blocked.zip", format="zip")
        context.compute(timeout=30)
        assert context.render.state == "failed", context.render.exception
        assert context.copy.state == "blocked"
        for target, address in ((url, "export_tree"), (blocked_url, "blocked_tree")):
            response = _get(target)
            assert response.status_code == 409, response.text
            assert response.json()["cells"][address] in {"failed", "blocked"}
    finally:
        context.close()


def test_e0_edit_after_cut_keeps_old_archive_consistent(probe, monkeypatch):
    from seamless_workflow.attachments.share import server

    context, url = probe
    context.compute(timeout=30)
    entered, release = threading.Event(), threading.Event()
    original_pack = server.pack_archive

    def paused_pack(format, members):
        # The server resolves every leaf from the held cut before packing.
        entered.set()
        if not release.wait(15):
            raise RuntimeError("E0 test did not release archive packing")
        return original_pack(format, members)

    monkeypatch.setattr(server, "pack_archive", paused_pack)
    with ThreadPoolExecutor(max_workers=1) as pool:
        pending = pool.submit(_get, url)
        try:
            assert entered.wait(15), "download never reached packing after its cut"
            context.decision.set(2)
            context.compute(timeout=30)
        finally:
            release.set()
        old = pending.result(timeout=20)
    assert old.status_code == 200, old.text
    assert _members(old) == _expected(1)
    assert _members(_get(url)) == _expected(2)


def test_e0_missing_folder_leaf_answers_503(download_server):
    with Context() as context:
        context.shares.namespace = "e0-missing-" + uuid.uuid4().hex
        context.export_tree = Cell(celltype="folder")
        absent = hashlib.sha256(b"E0 deliberately absent buffer " + uuid.uuid4().bytes).hexdigest()
        context.export_tree.set({"DATA/reference.pdb": absent})
        context.compute(timeout=30)
        url = context.shares.add_download(
            {"alaric_project": "export_tree"}, path="export.zip", format="zip")
        assert context.export_tree.state == "complete"
        response = _get(url)
        assert response.status_code == 503, response.text


def test_e0_folder_results_resolve_after_fresh_interpreter_restart(tmp_path):
    # Each process has new worker/cache/server state and the same durable
    # decision. This must compute and resolve the folder without old owners.
    decision_file = tmp_path / "decision.json"
    decision_file.write_text("1", encoding="utf-8")
    results = []
    for _ in range(2):
        completed = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--restart-probe", str(decision_file)],
            capture_output=True, text=True, timeout=90,
        )
        assert completed.returncode == 0, completed.stdout + completed.stderr
        rows = [line.removeprefix("E0_RESTART_RESULT=")
                for line in completed.stdout.splitlines()
                if line.startswith("E0_RESTART_RESULT=")]
        assert len(rows) == 1, completed.stdout + completed.stderr
        results.append(json.loads(rows[0]))
    expected = {name: content.hex() for name, content in _expected(1).items()}
    assert results == [expected, expected]


def _restart_probe(decision_file):
    import seamless
    from seamless_config.select import select_execution, select_persistent

    select_execution("process")
    select_persistent(False)
    shareserver.configure(host="127.0.0.1", port=0)
    context, url = _build_context("e0-restart", json.loads(Path(decision_file).read_text()))
    try:
        context.compute(timeout=30)
        response = _get(url)
        assert response.status_code == 200, response.text
        print("E0_RESTART_RESULT=" + json.dumps(
            {name: content.hex() for name, content in _members(response).items()}, sort_keys=True))
    finally:
        context.close()
        shareserver.close()
        seamless.close()


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] != "--restart-probe":
        raise SystemExit("usage: test_export_download.py --restart-probe DECISION.json")
    _restart_probe(sys.argv[2])
