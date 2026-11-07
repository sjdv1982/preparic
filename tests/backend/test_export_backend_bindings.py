"""E5 checks for the export transformer bindings and registered download."""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest
import requests

import backend
from seamless_workflow import shareserver


EXPECTED_INPUT_CELLTYPES = {
    "pdbcode": "str",
    "biounit": "plain",
    "cif": "str",
    "ranges": "plain",
    "stacking_selection": "plain",
    "threshold_edits": "plain",
    "best_fit_clip": "plain",
    "interaction_matrix": "plain",
    "gaussian_classes": "plain",
    "energy_threshold_edits": "plain",
    "range_specs": "plain",
    "structure_id": "plain",
    "range_pdb_files": "folder",
    "range_pdb_files_aa": "folder",
    "rna_pdb": "bytes",
    "fragments": "plain",
    "best_fit_thresholds": "plain",
    "best_fits": "plain",
    "clipped_fragments": "folder",
    "stackings": "plain",
    "interactions": "plain",
    "energetics": "plain",
    "bestfit_stackings": "plain",
    "gaussian_scores": "plain",
    "gaussian_validation": "plain",
    "gaussian_check": "plain",
    "gaussian_model_json": "bytes",
    "fraglib_identities": "str",
    "exclude_pdb_code": "str",
    "ring_refe_a": "bytes",
    "ring_refe_c": "bytes",
    "ring_refe_phe": "bytes",
    "ring_refe_his": "bytes",
    "ring_refe_arg": "bytes",
    "ring_refe_trp": "bytes",
    "range_pdbs_aa": "plain",
}
VALIDATOR_PINS = ("validated_pairs", "validated_classes", "validated_energies")


def _close_context_without_closing_runtime(context):
    """Close this test Context without invalidating Seamless for later tests."""
    stop = backend._STAGE5_RANGE_SYNC_STOP
    thread = backend._STAGE5_RANGE_SYNC_THREAD
    if stop is not None:
        stop.set()
    if thread is not None:
        thread.join(timeout=2)
    backend._STAGE5_RANGE_SYNC_STOP = None
    backend._STAGE5_RANGE_SYNC_THREAD = None
    context.close()
    if backend.ctx is context:
        backend.ctx = None


def test_registry_and_both_export_wrappers_use_the_explicit_plan_inputs():
    assert backend.EXPORT_INPUT_CELLTYPES == EXPECTED_INPUT_CELLTYPES

    admission_parameters = tuple(inspect.signature(backend._export_admission).parameters)
    tree_parameters = tuple(inspect.signature(backend._export_tree).parameters)
    registry_parameters = tuple(EXPECTED_INPUT_CELLTYPES)
    assert admission_parameters == registry_parameters + VALIDATOR_PINS
    assert tree_parameters == registry_parameters + ("admission",) + VALIDATOR_PINS
    for wrapper in (backend._export_admission, backend._export_tree):
        parameters = inspect.signature(wrapper).parameters
        assert all(parameter.kind is inspect.Parameter.POSITIONAL_OR_KEYWORD
                   for parameter in parameters.values())


@pytest.fixture
def wired_export_context(monkeypatch, tmp_path, seamless_runtime):
    """Build the real context with isolated mount paths and share server."""
    old_mount = backend.MOUNT_DIRECTORY
    for name, value in tuple(vars(backend).items()):
        if isinstance(value, Path) and value.parent == old_mount:
            monkeypatch.setattr(backend, name, tmp_path / value.name)
    monkeypatch.setattr(backend, "MOUNT_DIRECTORY", tmp_path)
    # The wiring inspection itself should not start the configured worker pool.
    monkeypatch.setattr(backend, "configure_local_execution", lambda: 8)
    monkeypatch.setenv("SEAMLESS_SHARE_HOST", "127.0.0.1")
    monkeypatch.setenv("SEAMLESS_SHARE_PORT", "0")
    shareserver.close()

    context = backend.build_context()
    try:
        yield context
    finally:
        _close_context_without_closing_runtime(context)
        shareserver.close()


@pytest.mark.integration
def test_build_context_wires_matching_pins_modules_gates_and_live_download(
    wired_export_context,
):
    context = wired_export_context
    admission = context.export_admission_transformer
    tree = context.export_tree_transformer

    assert admission.celltypes.result == "plain"
    assert tree.celltypes.result == "folder"
    for pin, celltype in EXPECTED_INPUT_CELLTYPES.items():
        assert admission.celltypes[pin] == celltype, pin
        assert tree.celltypes[pin] == celltype, pin
    for pin in VALIDATOR_PINS:
        assert admission.celltypes[pin] == "plain", pin
        assert tree.celltypes[pin] == "plain", pin
    assert tree.celltypes.admission == "plain"

    graph_connections = {
        tuple(connection["target"]): tuple(connection["source"])
        for connection in context.get_graph()["connections"]
        if connection.get("type") == "connection"
    }
    for pin in tuple(EXPECTED_INPUT_CELLTYPES) + VALIDATOR_PINS:
        source = graph_connections[("export_admission_transformer", pin)]
        assert graph_connections[("export_tree_transformer", pin)] == source, pin
    assert graph_connections[("export_tree_transformer", "admission")] == (
        "export_admission",
    )

    assert context.export_admission.share.spec.readonly is True
    assert context.export_tree.share.spec is None
    download = context.shares.downloads["export.zip"]
    assert download.format == "zip"
    assert download.cells == (("alaric_project", "export_tree"),)

    # The incomplete initial graph still has a registered route. GET observes
    # the live tree state and refuses it instead of returning an empty ZIP.
    response = requests.get(
        context.shares.url.rstrip("/") + "/export.zip",
        timeout=20,
        headers={"Cache-Control": "no-cache"},
    )
    assert response.status_code == 409, response.text
    assert response.json()["cells"]["export_tree"] in {
        "waiting", "computing", "failed", "blocked",
    }

    expected_admission_modules = {
        "export_snapshot": backend.PROJECT_ROOT / "export_snapshot.py",
        "preparic_stage4": backend.PROJECT_ROOT / "stage4.py",
        "preparic_stage5": backend.PROJECT_ROOT / "stage5.py",
        "preparic_stage6": backend.PROJECT_ROOT / "stage6.py",
        "alaric_gaussian_filter": backend.resolve_alaric_dir() / "gaussian_filter.py",
    }
    for module_name, source_path in expected_admission_modules.items():
        for wrapper in (admission, tree):
            module = getattr(wrapper.modules, module_name)
            assert module["type"] == "interpreted"
            assert module["language"] == "python"
            assert module["code"] == source_path.read_text(encoding="utf-8")

    for module_name, source_name in (
        ("export_gens", "export_gens.py"),
        ("export_tree", "export_tree.py"),
    ):
        module = getattr(tree.modules, module_name)
        assert module["code"] == (backend.PROJECT_ROOT / source_name).read_text(
            encoding="utf-8"
        )

