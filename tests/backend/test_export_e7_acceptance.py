"""Transferred W8 seam obligations against real named export input pins.

The E2 fixture models producer payloads rather than the retired bundle API.
These checks bind admission to the rendered archive payload and manifest.
"""

from copy import deepcopy
import hashlib
import json

import pytest
import yaml

import export_snapshot
import export_tree
from test_export_snapshot import _attached_helper_modules, _inputs  # noqa: F401


def _build(inputs):
    admission = export_snapshot.build_admission(inputs)
    assert admission["ready"] is True, admission["blocks"]
    tree = export_tree.build_tree(inputs, admission)
    manifest = json.loads(tree["DATA/export-manifest.json"])
    return admission, tree, manifest


def _actions(tree):
    return {path: yaml.safe_load(value) for path, value in tree.items()
            if path.endswith("/alaric.yaml")}


def test_e7_x1_x2_real_pin_cut_preserves_exact_payload_and_identity():
    inputs = _inputs()
    original = deepcopy(inputs)
    admission, tree, manifest = _build(inputs)
    assert inputs == original
    assert tree["DATA/reference.pdb"] == inputs["rna_pdb"]
    assert tree["DATA/reference.pdb"] == inputs["range_pdb_files_aa"]["range-1-aa.pdb"]
    assert b"\r\n" in tree["DATA/reference.pdb"]
    assert manifest["reference_sha256"] == hashlib.sha256(inputs["rna_pdb"]).hexdigest()
    assert manifest["snapshot_sha256"] == admission["digest"]
    assert tree["DATA/pdbcode.txt"] == b"1b7f\n"
    for pin in ("range_pdb_files", "range_pdb_files_aa"):
        for name, value in inputs[pin].items():
            assert tree[f"DATA/{name}"] == value
    for name, value in inputs["clipped_fragments"].items():
        if name.startswith(("reference/", "bestfit/")):
            assert tree[f"DATA/clipped-fragments/{name}"] == value
    assert all(b"ALARIC_ANCHOR_MIN_CELL_FRACTION" not in value for value in tree.values())
    assert _build(deepcopy(inputs))[1] == tree


@pytest.mark.parametrize("change", ["exclude", "fraction", "rna-bytes", "displayed-bytes"])
def test_e7_x1_x2_conflicting_preparation_evidence_is_non_overridable(change):
    inputs = _inputs()
    if change == "exclude":
        inputs["exclude_pdb_code"] = "3sxl"
    elif change == "fraction":
        inputs["gaussian_check"]["min_cell_fraction"] = 0.5
    elif change == "rna-bytes":
        inputs["rna_pdb"] += b"REMARK changed\n"
    else:
        inputs["range_pdbs_aa"]["pdbs"][0]["pdb"] += "REMARK changed\n"
    admission = export_snapshot.build_admission(inputs)
    assert admission["ready"] is False
    assert admission["blocks"] and admission["warnings"] == []
    # Consent is local to a browser and cannot make a blocked folder build.
    admission["proceed_anyway"] = True
    with pytest.raises(export_snapshot.ExportNotReady):
        export_tree.build_tree(inputs, admission)


@pytest.mark.parametrize("rounded,raw,sentinel", [
    ("1.73", 1.72, False), ("1.74", 1.731, True), ("2.25", 2.24, True),
])
def test_e7_x3_computed_reference_thresholds_bind_yaml_and_manifest(rounded, raw, sentinel):
    inputs = _inputs()
    base = inputs["best_fit_thresholds"]["bases"][0]
    base.update(nuc2_ovrmsd=rounded, nuc2_ovrmsd_raw=raw)
    admission, tree, manifest = _build(inputs)
    actions = _actions(tree)
    matches = {path: action for path, action in actions.items()
               if action["action"] == "anchor-refe" and action["fragment"] == 1
               and action["nucleotide"] == "second"}
    assert matches
    for path, action in matches.items():
        assert float(action["ovrmsd"]) == (-1 if sentinel else float(rounded))
        text = tree[path].decode()
        assert "comment:" not in text
        comments = " ".join(line[1:].strip() for line in text.splitlines() if line.startswith("#"))
        assert bool(comments) is sentinel
        if sentinel:
            assert str(raw) in comments and rounded in comments
            assert "disabled positive control" in comments.lower()
            assert "without misses" in comments.lower()
            if raw < 3 ** 0.5:
                assert "best-fit pose itself lies beyond" not in comments.lower()
    disabled = [row for row in manifest["disabled_reference_anchors"]
                if row["candidate"] == {"fragment": "frag-1-AC", "nucleotide": "second"}]
    assert len(disabled) == int(sentinel)
    if sentinel:
        assert set(disabled[0]["action_paths"]) == set(matches)
    # Changing only a reference candidate does not change the Gaussian graph
    # or remove any grow/identity/filter descendants from either alternative.
    baseline = _actions(_build(_inputs())[1])
    assert {path: action["action"] for path, action in actions.items()} == {
        path: action["action"] for path, action in baseline.items()}
    assert admission["warnings"] == _build(_inputs())[0]["warnings"]
    for gen in manifest["gens"]:
        raw = [actions[f"{gen['path']}/{name}/alaric.yaml"] for name in gen["actions"]]
        refs = [(row["fragment"], row["nucleotide"]) for row in raw
                if row["action"] == "anchor-refe"]
        assert len(refs) == len(set(refs)), gen["path"]


@pytest.mark.parametrize("pin", [
    "threshold_edits", "best_fit_clip", "interaction_matrix", "interactions",
    "energy_threshold_edits", "gaussian_classes", "gaussian_validation",
    "best_fit_thresholds", "best_fits", "bestfit_stackings", "gaussian_scores",
])
@pytest.mark.parametrize("replacement", [None, {}], ids=["null", "empty"])
def test_e7_x6_missing_real_evidence_refuses_tree_with_addressed_blocks(pin, replacement):
    inputs = _inputs()
    inputs[pin] = replacement
    admission = export_snapshot.build_admission(inputs)
    assert admission["ready"] is False
    assert admission["blocks"] and admission["warnings"] == []
    for block in admission["blocks"]:
        assert set(block) == {"code", "address", "message", "remedy"}
        assert block["message"] and block["remedy"]
    with pytest.raises(export_snapshot.ExportNotReady):
        export_tree.build_tree(inputs, admission)


@pytest.mark.parametrize("field", ["crmsd", "ovrmsd"])
@pytest.mark.parametrize("value", [None, "auto", True, 0, -1, "NaN", "inf"])
def test_e7_x6_kept_pair_thresholds_refuse_malformed_values(field, value):
    inputs = _inputs()
    inputs["threshold_edits"]["pairs"][0][field] = value
    admission = export_snapshot.build_admission(inputs)
    assert admission["ready"] is False and admission["blocks"]
    with pytest.raises(export_snapshot.ExportNotReady):
        export_tree.build_tree(inputs, admission)


def test_e7_x6_kept_pair_limit_and_missing_true_cell_threshold_cannot_be_overridden():
    for change in ("pair-limit", "missing-energy"):
        inputs = _inputs()
        if change == "pair-limit":
            inputs["threshold_edits"]["pairs"][0]["ovrmsd"] = "1.74"
        else:
            del inputs["energy_threshold_edits"]["ranges"][0]["values"]["frag-1-AC"]
        admission = export_snapshot.build_admission(inputs)
        assert admission["ready"] is False and admission["blocks"]
        with pytest.raises(export_snapshot.ExportNotReady):
            export_tree.build_tree(inputs, {**admission, "ready": True, "blocks": []})


@pytest.mark.parametrize("field", ["nuc2_ovrmsd", "nuc2_ovrmsd_raw"])
@pytest.mark.parametrize("value", [None, "auto", True, "NaN", float("inf")])
def test_e7_x3_missing_or_malformed_computed_candidate_evidence_blocks(field, value):
    inputs = _inputs()
    inputs["best_fit_thresholds"]["bases"][0][field] = value
    admission = export_snapshot.build_admission(inputs)
    assert admission["ready"] is False and admission["blocks"]
    with pytest.raises(export_snapshot.ExportNotReady):
        export_tree.build_tree(inputs, admission)


@pytest.mark.parametrize("change", ["missing-base", "duplicate-base", "missing-pair", "bad-address", "bad-sequence"])
def test_e7_x3_x6_retained_candidate_and_adjacent_pair_coverage_is_complete(change):
    inputs = _inputs()
    if change == "missing-base":
        inputs["best_fit_thresholds"]["bases"].pop(0)
    elif change == "duplicate-base":
        duplicate = deepcopy(inputs["best_fit_thresholds"]["bases"][0])
        duplicate["nuc2_ovrmsd"] = "2.26"
        inputs["best_fit_thresholds"]["bases"].append(duplicate)
    elif change == "missing-pair":
        inputs["threshold_edits"]["pairs"] = []
    elif change == "bad-address":
        inputs["threshold_edits"]["pairs"][0]["up"] = 3
    else:
        inputs["threshold_edits"]["pairs"][0]["sequence"] = "UUU"
    admission = export_snapshot.build_admission(inputs)
    assert admission["ready"] is False and admission["blocks"]
    with pytest.raises(export_snapshot.ExportNotReady):
        export_tree.build_tree(inputs, admission)


def test_e7_x7_changed_decision_invalidates_admission_and_old_tree_stays_immutable():
    inputs = _inputs()
    admission, tree, manifest = _build(inputs)
    frozen = deepcopy(tree)
    inputs["energy_threshold_edits"]["ranges"][0]["values"]["frag-1-AC"] = "-3"
    current, changed_tree, changed_manifest = _build(inputs)
    assert current["digest"] != admission["digest"]
    assert changed_manifest["snapshot_sha256"] != manifest["snapshot_sha256"]
    assert changed_tree["DATA/energy-thresholds.json"] != tree["DATA/energy-thresholds.json"]
    with pytest.raises(export_snapshot.ExportNotReady, match="digest"):
        export_tree.build_tree(inputs, admission)
    assert tree == frozen


@pytest.mark.parametrize("change", ["pair", "range-file", "ring", "library", "model"])
def test_e7_x7_output_and_validation_inputs_change_manifest_identity(change):
    inputs = _inputs()
    admission, tree, _ = _build(inputs)
    if change == "pair":
        inputs["threshold_edits"]["pairs"][0]["crmsd"] = "0.15"
        inputs["clipped_fragments"]["thresholds.json"] = json.dumps({
            "schema": 1, "pairs": inputs["threshold_edits"]["pairs"],
        }, sort_keys=True).encode() + b"\n"
    elif change == "range-file":
        inputs["range_pdb_files"]["range-2.pdb"] += b"REMARK identity edit\n"
    elif change == "ring":
        inputs["ring_refe_a"] += b"changed"
    elif change == "library":
        inputs["fraglib_identities"] += "-changed"
    else:
        model = json.loads(inputs["gaussian_model_json"])
        model["thresholds"]["rrm90"]["cutoff"] = 0.02
        inputs["gaussian_model_json"] = json.dumps(model, sort_keys=True).encode()
    current, changed_tree, manifest = _build(inputs)
    assert current["digest"] != admission["digest"]
    assert manifest["snapshot_sha256"] == current["digest"]
    assert changed_tree["DATA/export-manifest.json"] != tree["DATA/export-manifest.json"]
    with pytest.raises(export_snapshot.ExportNotReady, match="digest"):
        export_tree.build_tree(inputs, admission)


def test_e7_x7_view_and_local_consent_do_not_change_archive_identity():
    inputs = _inputs()
    admission, tree, _ = _build(inputs)
    inputs.update(stage4_active_nucleotide=102, proceed_anyway=True, acknowledgedDigest="local")
    current, current_tree, _ = _build(inputs)
    assert current["digest"] == admission["digest"]
    assert current_tree == tree
