"""W1 acceptance: durable decision loading and local execution policy.

These tests cover the shared-project contract, independently of scientific fits.
The frontend companion set covers pure stage/freshness/error semantics.
"""

import importlib.util
import json
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import urlopen

import pytest
import seamless.config
import seamless_config
import yaml

import backend


RNA = {"chain": "A", "start": 101, "end": 104}
PROTEIN = {"chain": "B", "start": 10, "end": 30}
ROWS = [RNA, PROTEIN]


def marker(name):
    values = {
        "ranges": {"schema": 1, "basis": {"pdbcode": "1B7F", "biounit": 1}, "active": True, "rows": ROWS},
        "stacking_selection": {"schema": 1, "basis": {"pdbcode": "1B7F", "biounit": 1, "rows": ROWS}, "active": False, "selected": {"anchor-A": True}},
        "threshold_edits": {"schema": 1, "basis": {"pdbcode": "1B7F", "rna": RNA}, "active": True, "pairs": []},
        "best_fit_clip": {"schema": 1, "basis": {"pdbcode": "1B7F", "rna": RNA}, "clip1": 0, "clip2": 0},
        "interaction_matrix": {"schema": 1, "basis": ["anchor-A"], "active": False, "rows": ["frag-1-AC"], "ranges": ["range-2"], "grid": [[True]]},
        "gaussian_classes": {"schema": 1, "basis": ["anchor-A"], "active": True, "rows": [{"key": "anchor-A", "rrm": 50}]},
        "energy_threshold_edits": {"schema": 1, "basis": {"pdbcode": "1B7F", "rna": RNA}, "ranges": [{"row": PROTEIN, "values": {"frag-1-AC": None}}]},
    }
    return values[name]


@pytest.mark.parametrize("name", ["ranges", "stacking_selection", "threshold_edits", "best_fit_clip", "interaction_matrix", "gaussian_classes", "energy_threshold_edits"])
def test_restart_load_preserves_each_canonical_decision_including_inactive_and_blank(tmp_path, name):
    """OC P2/P4/P8: restart retains mounted edits; blanks are not undecided."""
    path = tmp_path / ("ranges.yaml" if name == "ranges" else name + ".json")
    value = marker(name)
    path.write_text(yaml.safe_dump(value) if name == "ranges" else json.dumps(value), encoding="utf-8")
    loaded, diagnostics = backend.load_marker_payload(path, name)
    assert loaded == value
    assert diagnostics == []


@pytest.mark.parametrize("name", ["ranges", "stacking_selection", "threshold_edits", "interaction_matrix", "gaussian_classes", "energy_threshold_edits"])
def test_explicit_undecided_survives_load(tmp_path, name):
    path = tmp_path / ("ranges.yaml" if name == "ranges" else name + ".json")
    path.write_text("null\n", encoding="utf-8")
    assert backend.load_marker_payload(path, name) == (None, [])


@pytest.mark.parametrize("raw", ['{not-json', '[true,false]', '{"schema":1,"active":true,"selected":{"anchor-A":true}}', '{"schema":2,"basis":{},"active":true,"selected":{}}'])
def test_malformed_legacy_and_unsupported_decisions_reset_with_preserved_diagnostic(tmp_path, raw):
    """Author I7: do not silently assign a basis or swallow parse failures."""
    path = tmp_path / "stacking-selection.json"
    path.write_text(raw, encoding="utf-8")
    payload, diagnostics = backend.load_marker_payload(path, "stacking_selection")
    assert payload is None
    assert diagnostics
    assert all(isinstance(message, str) and message.strip() for message in diagnostics)
    # Normalization must not destroy the only copy before reporting diagnostics.
    assert path.read_text(encoding="utf-8") == raw


def test_missing_decision_has_no_diagnostic_and_clip_has_explicit_startup_fallback(tmp_path):
    assert backend.load_marker_payload(tmp_path / "absent.json", "stacking_selection") == (None, [])
    fallback = {"schema": 1, "basis": None, "clip1": 0, "clip2": 0}
    assert backend.load_marker_payload(tmp_path / "clip.json", "best_fit_clip", fallback) == (fallback, [])


def test_t0_s5_rows_envelope_is_consumed_by_t1(cif_text, range_specs, expected_pdbs):
    decision = {
        "schema": 1,
        "basis": {"pdbcode": "1B7F", "biounit": 1},
        "active": True,
        "rows": range_specs["ranges"],
    }
    s5 = backend._parse_ranges(decision, "1B7F", 1)
    assert s5 == {"schema": 1, "pdbcode": "1B7F", "biounit": 1, "rows": range_specs["ranges"]}
    assert backend._ranges_to_pdb(cif_text, s5) == expected_pdbs


def test_execution_selects_eight_local_spawn_workers_once_without_service_initialization(monkeypatch):
    """D14: injected spawn observes the actual startup scheduling decision."""
    # A fresh module has not previously configured workers in another test.
    spec = importlib.util.spec_from_file_location("_w1_backend_execution", backend.__file__)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    calls = []
    monkeypatch.setattr(module.seamless.config, "init", lambda *args, **kwargs: pytest.fail("startup contacted configured execution/cache services"), raising=False)
    monkeypatch.setattr(seamless_config, "init", lambda *args, **kwargs: pytest.fail("startup initialized configured services"))
    assert module.configure_local_execution(spawn_fn=lambda count: calls.append(count)) == 8
    assert module.configure_local_execution(spawn_fn=lambda count: calls.append(count)) == 8
    assert calls == [8]


def test_all_decisions_and_cif_are_file_authoritative_mounts(monkeypatch, tmp_path, seamless_runtime):
    """P2/CB S3: restoration must use CIF bytes on disk, without browser repair."""
    old_directory = backend.MOUNT_DIRECTORY
    for name, value in tuple(vars(backend).items()):
        if isinstance(value, Path) and value.parent == old_directory:
            monkeypatch.setattr(backend, name, tmp_path / value.name)
    monkeypatch.setattr(backend, "MOUNT_DIRECTORY", tmp_path)
    monkeypatch.setattr(backend, "configure_local_execution", lambda *args, **kwargs: 8)
    monkeypatch.setenv("SEAMLESS_SHARE_HOST", "127.0.0.1")
    monkeypatch.setenv("SEAMLESS_SHARE_PORT", "0")
    context = backend.build_context()
    try:
        for name in ["pdbcode", "biounit", "cif", "ranges", "stacking_selection", "threshold_edits", "best_fit_clip", "interaction_matrix", "gaussian_classes", "energy_threshold_edits"]:
            cell = getattr(context, name)
            assert cell.mount.spec.authority == "file", name
            assert cell.mount.spec.mode == "rw", name
            assert Path(cell.mount.spec.path).parent == tmp_path, name
        for name in ["ranges", "stacking_selection", "threshold_edits", "best_fit_clip", "interaction_matrix", "gaussian_classes", "energy_threshold_edits"]:
            assert getattr(context, name).share.spec.readonly is False, name
        assert context.range_specs.share.spec.readonly is True
        assert context.range_pdbs.share.spec is None
        # S26 is shared view state, deliberately not persisted (CB §8.2).
        assert context.stage4_active_nucleotide.share.spec.readonly is False
        assert context.stage4_active_nucleotide.mount.spec is None

        cell_files = {
            "best_fit_thresholds": "best-fit-thresholds.json",
            "energy_thresholds": "computed-energy-thresholds.json",
            "interactions": "interactions.json",
            "gaussian_validation": "gaussian-validation.json",
        }
        for name, filename in cell_files.items():
            mount = getattr(context, name).mount.spec
            assert mount.authority == "cell", name
            assert mount.mode == "w", name
            assert Path(mount.path) == tmp_path / filename, name

        published = {
            "best_fit_thresholds": context.best_fit_thresholds,
            "best_fits": context.best_fits,
            "bestfit_stackings": context.bestfit_stackings,
            "energy_thresholds": context.energy_thresholds,
            "energetics": context.energetics,
            "interactions": context.interactions,
            "gaussian_validation": context.gaussian_validation,
            "gaussian_scores": context.gaussian_scores,
            "reference_fragment_pdbs": context.reference_fragment_pdbs,
        }
        server_paths = context.shares.openapi()["paths"]
        for name, cell in published.items():
            share = cell.share
            assert share.spec.readonly is True, name
            assert share.spec.path == name, name
            route = f"/preparic/{name}"
            assert route in server_paths, name
            assert share.url.endswith(route), name

        for name in ("best_fits", "bestfit_stackings", "gaussian_scores"):
            assert getattr(context, name).mount.spec is None, name

        # Exercise the actual T8 URL as a browser does. 404 means the route is
        # registered but its computed value is not available yet; OpenAPI
        # above independently proves the route itself is live.
        try:
            with urlopen(context.best_fit_thresholds.share.url, timeout=5) as response:
                assert response.status in (200, 204)
        except HTTPError as error:
            assert error.code == 404
    finally:
        context.close()
        backend.ctx = None


def test_legacy_threshold_bases_load_as_undecided_with_preserved_diagnostics(tmp_path):
    path = tmp_path / "best-fit-threshold-edits.json"
    raw = json.dumps({**marker("threshold_edits"), "bases": []})
    path.write_text(raw)
    value, diagnostics = backend.load_marker_payload(path, "threshold_edits")
    assert value is None
    assert diagnostics
    assert path.read_text() == raw
