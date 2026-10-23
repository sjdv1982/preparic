"""Integration checks for the live backend graph and its stacking attachments."""

import json
from urllib.request import urlopen

import pytest

import backend
from detect_stacking import _detect_stacking


pytestmark = [pytest.mark.integration, pytest.mark.usefixtures("seamless_runtime")]


def _redirect_mounts(monkeypatch, mount_directory):
    monkeypatch.setattr(backend, "MOUNT_DIRECTORY", mount_directory)
    monkeypatch.setattr(backend, "PDBCODE_PATH", mount_directory / "pdbcode.txt")
    monkeypatch.setattr(backend, "BIOUNIT_PATH", mount_directory / "biounit.json")
    monkeypatch.setattr(backend, "RANGES_PATH", mount_directory / "ranges.yaml")
    monkeypatch.setattr(backend, "PDBS_PATH", mount_directory / "pdbs")
    monkeypatch.setattr(backend, "PDBS_AA_PATH", mount_directory / "pdbs-aa")
    monkeypatch.setenv("SEAMLESS_SHARE_HOST", "127.0.0.1")
    monkeypatch.setenv("SEAMLESS_SHARE_PORT", "0")


@pytest.mark.parametrize(
    ("serialized", "expected"),
    [
        ("[]", []),
        ("[true, false]", [True, False]),
        ("[true, 1]", []),
        ('{"selected": true}', []),
        ("not-json", []),
    ],
)
def test_initial_stacking_selection_accepts_only_boolean_lists(
    tmp_path, serialized, expected
):
    selection_path = tmp_path / "stacking-selection.json"
    selection_path.write_text(serialized, encoding="utf-8")

    assert backend._initial_stacking_selection(selection_path) == expected


def test_backend_graph_mounts_and_shares_stacking_cells(
    monkeypatch, tmp_path, cif_text, range_specs, expected_pdbs
):
    _redirect_mounts(monkeypatch, tmp_path)
    selection_path = tmp_path / "stacking-selection.json"
    selection_path.write_text("[true, 1]", encoding="utf-8")

    context = None
    try:
        context = backend.build_context()

        graph = context.get_graph()
        connections = graph["connections"]
        assert {
            "type": "connection",
            "source": ["range_pdbs"],
            "target": ["detect_stacking", "range_pdbs"],
        } in connections
        assert not any(
            edge["source"] == ["range_pdb_files_aa"]
            and edge["target"] == ["detect_stacking", "range_pdbs"]
            for edge in connections
        )
        assert any(
            node["path"] == ["detect_stacking"] for node in graph["nodes"]
        )

        assert context.stackings.celltype == "plain"
        assert context.stackings.mount.spec.mode == "w"
        assert context.stackings.mount.spec.authority == "cell"
        assert context.stackings.mount.spec.path == str(tmp_path / "stackings.json")
        assert context.stacking_selection.celltype == "plain"
        assert context.stacking_selection.mount.spec.mode == "rw"
        assert context.stacking_selection.mount.spec.authority == "file"
        assert context.stacking_selection.share.spec.readonly is False
        assert context.stackings.share.spec.readonly is True

        assert context.stacking_selection.value == []
        assert selection_path.read_bytes() == b"[]\n"

        # Empty input is a valid initialized pipeline value and mount payload.
        assert context.stackings.run() == []
        context.mounts.sync(timeout=10)
        stackings_path = tmp_path / "stackings.json"
        assert json.loads(stackings_path.read_text(encoding="utf-8")) == []
        assert stackings_path.read_bytes() == b"[]\n"

        # Exercise the actual extraction -> detector dependency and compare it
        # to the same isolated transformation over its extracted PDB inputs.
        context.cif.set(cif_text)
        context.range_specs.set(range_specs)
        expected_stackings = _detect_stacking(expected_pdbs)
        assert context.stackings.run() == expected_stackings
        context.mounts.sync(timeout=10)
        assert json.loads(stackings_path.read_text(encoding="utf-8")) == expected_stackings

        with urlopen(f"{context.shares.url}/state-graph", timeout=5) as response:
            state_graph = json.load(response)
        transformer_paths = {
            tuple(node["path"])
            for node in state_graph["nodes"]
            if node["type"] == "transformer"
        }
        assert ("ranges_to_pdb",) in transformer_paths
        assert ("detect_stacking",) in transformer_paths
    finally:
        if context is not None:
            context.close()
            backend.ctx = None
