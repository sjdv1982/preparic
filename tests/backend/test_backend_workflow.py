"""Integration checks for the live backend graph and its stacking attachments."""

import json
from urllib.request import urlopen

import pytest
import yaml

from alaric_dir import resolve_alaric_repo_root
import backend
from detect_stacking import _detect_stacking


pytestmark = [pytest.mark.integration, pytest.mark.usefixtures("seamless_runtime")]


def _redirect_mounts(monkeypatch, mount_directory):
    monkeypatch.setattr(backend, "MOUNT_DIRECTORY", mount_directory)
    monkeypatch.setattr(backend, "PDBCODE_PATH", mount_directory / "pdbcode.txt")
    monkeypatch.setattr(backend, "BIOUNIT_PATH", mount_directory / "biounit.json")
    monkeypatch.setattr(backend, "CIF_PATH", mount_directory / "cif.cif")
    monkeypatch.setattr(backend, "RANGES_PATH", mount_directory / "ranges.yaml")
    monkeypatch.setattr(backend, "PDBS_PATH", mount_directory / "pdbs")
    monkeypatch.setattr(backend, "PDBS_AA_PATH", mount_directory / "pdbs-aa")
    marker_files = {
        "STACKING_SELECTION_PATH": "stacking-selection.json",
        "THRESHOLD_EDITS_PATH": "best-fit-threshold-edits.json",
        "BEST_FIT_CLIP_PATH": "best-fit-clip.json",
        "INTERACTION_MATRIX_PATH": "interaction-matrix.json",
        "GAUSSIAN_CLASSES_PATH": "gaussian-classes.json",
        "ENERGY_THRESHOLD_EDITS_PATH": "energy-threshold-edits.json",
    }
    for name, filename in marker_files.items():
        monkeypatch.setattr(backend, name, mount_directory / filename)
    monkeypatch.setenv("SEAMLESS_SHARE_HOST", "127.0.0.1")
    monkeypatch.setenv("SEAMLESS_SHARE_PORT", "0")


@pytest.mark.parametrize(
    ("serialized", "expected"),
    [
        ('{"schema":1,"basis":{"pdbcode":"1B7F","biounit":1,"rows":[{"chain":"A","start":101,"end":104}]},"active":true,"selected":{"anchor-A":true,"anchor-B":false}}',
         {"schema": 1, "basis": {"pdbcode": "1B7F", "biounit": 1, "rows": [{"chain": "A", "start": 101, "end": 104}]}, "active": True, "selected": {"anchor-A": True, "anchor-B": False}}),
        ("[]", None),
        ('{"schema":1,"basis":{"pdbcode":"1B7F","biounit":1,"rows":[{"chain":"A","start":101,"end":104}]},"active":true,"selected":{"anchor-A":1}}', None),
        ("not-json", None),
    ],
)
def test_initial_stacking_selection_uses_canonical_keyed_markers(
    tmp_path, serialized, expected
):
    selection_path = tmp_path / "stacking-selection.json"
    selection_path.write_text(serialized, encoding="utf-8")

    payload, diagnostics = backend.load_marker_payload(selection_path, "stacking_selection")
    assert payload == expected
    assert bool(diagnostics) is (expected is None)


def test_backend_graph_mounts_and_shares_stacking_cells(
    monkeypatch, tmp_path, cif_text, expected_pdbs
):
    _redirect_mounts(monkeypatch, tmp_path)
    range_rows = [
        {"chain": "A", "start": 1, "end": 4},
        {"chain": "C", "start": 2, "end": 8},
        {"chain": "C", "start": 90, "end": 96},
    ]
    selection_path = tmp_path / "stacking-selection.json"
    selection = {
        "schema": 1,
        "basis": {"pdbcode": "1B7F", "biounit": 1, "rows": range_rows},
        "active": True,
        "selected": {"anchor-A": True},
    }
    selection_path.write_text(json.dumps(selection), encoding="utf-8")

    context = None
    try:
        context = backend.build_context()
        assert context.cif.mount.spec.path == str(tmp_path / "cif.cif")

        # T4 must pin the installed files directly, so its inputs follow the
        # mounted paths and remain available to spawned workers.
        from pathlib import Path
        installed_references = resolve_alaric_repo_root() / "ring-refe"
        for name in ("A", "C", "PHE", "HIS", "ARG", "TRP"):
            filename = f"refe-{name}.pdb"
            mounted_file = Path(getattr(context, f"ring_refe_{name.lower()}").mount.spec.path)
            assert mounted_file.samefile(installed_references / filename)

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

        assert context.stacking_selection.value == selection
        assert json.loads(selection_path.read_text(encoding="utf-8")) == selection

        context.cif.set(cif_text)
        # Empty input still carries the S7 identity echo, so consumers can
        # distinguish a fresh empty result from stale stacking rows.
        empty_stackings = context.stackings.run()
        empty_ranges = context.range_specs.value
        assert empty_stackings == {
            "schema": 1,
            "echo": {
                "structure_id": context.structure_id_result.value,
                "pdbcode": empty_ranges["pdbcode"],
                "biounit": empty_ranges["biounit"],
                "rows": [],
            },
            "rows": [],
        }
        context.mounts.sync(timeout=10)
        stackings_path = tmp_path / "stackings.json"
        assert json.loads(stackings_path.read_text(encoding="utf-8")) == empty_stackings

        # Exercise the actual extraction -> detector dependency and compare it
        # to the same isolated transformation over its extracted PDB inputs.
        decision = {
            "schema": 1,
            "basis": {"pdbcode": "1B7F", "biounit": 1},
            "active": True,
            "rows": range_rows,
        }
        context.pdbcode.set("1B7F")
        context.biounit.set(1)
        context.ranges.set(yaml.safe_dump(decision, sort_keys=False))
        # T3's prepared RNA must reach both B4 and T24 in the live graph.
        # This catches a nonempty S23 whose separate B4 pin is empty or stale.
        prepared_folder = context.range_pdb_files_aa.run()
        assert "range-1-aa.pdb" in prepared_folder
        prepared_rna = context.rna_pdb.run()
        assert isinstance(prepared_rna, bytes) and prepared_rna.strip()
        assert prepared_folder["range-1-aa.pdb"] == backend.seamless.Buffer(
            prepared_rna, "bytes"
        ).get_checksum().hex()
        shared_prepared = context.range_pdbs_aa.run()
        assert shared_prepared["pdbs"][0]["pdb"].encode() == prepared_rna
        fragments = context.fragments.run()
        assert fragments["fragments"]
        assert fragments["residues"]
        # Run the isolated oracle over T1's exact emitted text and the same
        # checksum-pinned ring references consumed by the graph transformer.
        ring_root = resolve_alaric_repo_root() / "ring-refe"
        ring_names = ("A", "C", "PHE", "HIS", "ARG", "TRP")
        ring_references = {
            f"refe-{name}.pdb": getattr(context, f"ring_refe_{name.lower()}").value
            for name in ring_names
        }
        assert {
            filename: backend.seamless.Buffer(data, "bytes").get_checksum()
            for filename, data in ring_references.items()
        } == {
            f"refe-{name}.pdb": backend.seamless.Buffer(
                (ring_root / f"refe-{name}.pdb").read_bytes(), "bytes"
            ).get_checksum()
            for name in ring_names
        }
        t1_range_pdbs = context.range_pdbs.run()
        expected_stackings = {
            "schema": 1,
            "echo": {
                "structure_id": backend._structure_id(cif_text),
                "pdbcode": "1B7F",
                "biounit": 1,
                "rows": range_rows,
            },
            "rows": _detect_stacking(t1_range_pdbs, ring_references=ring_references),
        }
        actual_stackings = context.stackings.run()
        assert actual_stackings == expected_stackings
        context.mounts.sync(timeout=60)
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
