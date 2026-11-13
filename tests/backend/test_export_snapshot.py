"""Pure E2 export snapshot acceptance fixtures.

The input values below mirror the named Preparic pins and producer payloads.
They deliberately avoid Context/Seamless construction and filesystem reads.
"""

from copy import deepcopy
import json
import sys
from types import ModuleType

import pytest

import export_snapshot


RNA = {"chain": "A", "start": 101, "end": 103}
PROTEIN = {"chain": "B", "start": 1, "end": 20}
RANGE_SPECS = {"schema": 1, "pdbcode": "1B7F", "biounit": 1,
               "rows": [RNA, PROTEIN]}
RESIDUES = [
    {"chain": "A", "resid": 101, "insertion_code": "", "base": "A"},
    {"chain": "A", "resid": 102, "insertion_code": "", "base": "C"},
    {"chain": "A", "resid": 103, "insertion_code": "", "base": "G"},
]
FIT_ECHO = {"pdb_code": "1b7f", "rna": RNA, "resids": RESIDUES}
FRAGMENT_ECHO = {"pdbcode": "1B7F", "rna": RNA}
FRAGMENT_KEYS = ["frag-1-AC", "frag-2-CG"]
STACKING_KEY = json.dumps(["C", 102, 2, 214], separators=(",", ":"))
CLIP = {"schema": 1, "basis": {"pdbcode": "1B7F", "rna": RNA},
        "clip1": 0, "clip2": 0}
CLASS_ROWS = [{"key": STACKING_KEY, "rrm": 90, "tilt30": False, "disp45": False}]
GAUSSIAN_NAMES = [
    "rrm50", "rrm60", "rrm70", "rrm80", "rrm90", "rrm90-tilt30",
    "rrm90-disp45", "rrm90-tilt30-disp45", "rrm95", "rrm95-tilt30",
    "rrm95-disp45", "rrm95-tilt30-disp45",
]
DIGEST_GOLDEN_EXPORTER_VERSION = 2
DIGEST_GOLDEN_SHA256 = "508508a364893bd691cb1051924c6a004457b18c9b532a4118ea65ed940a2ab7"


@pytest.fixture(autouse=True)
def _attached_helper_modules(monkeypatch):
    """Provide transformer-attached helper names without Seamless or Alaric reads."""
    import stage4 as preparic_stage4
    import stage5 as preparic_stage5

    stage6 = ModuleType("preparic_stage6")

    def validate_gaussian_classes(marker, model_json, selection, gaussian_filter):
        model = json.loads(model_json.decode("utf-8") if isinstance(model_json, bytes) else model_json)
        selected = sorted(key for key, value in selection["selected"].items() if value is True)
        assert model["schema"] == 1
        assert set(model["thresholds"]) == set(GAUSSIAN_NAMES)
        assert marker["schema"] == 1 and marker["active"] is True
        assert marker["basis"] == selected

    stage6.validate_gaussian_classes = validate_gaussian_classes
    gaussian_filter = ModuleType("alaric_gaussian_filter")
    monkeypatch.setitem(sys.modules, "preparic_stage4", preparic_stage4)
    monkeypatch.setitem(sys.modules, "preparic_stage5", preparic_stage5)
    monkeypatch.setitem(sys.modules, "preparic_stage6", stage6)
    monkeypatch.setitem(sys.modules, "alaric_gaussian_filter", gaussian_filter)


def _pdb(label):
    return (f"REMARK E2 fixture {label}\r\n"
            "ATOM      1  P     A A 101       0.000   0.000   0.000  1.00  0.00           P\r\n"
            "END\r\n").encode("ascii")


def _rna_pdb(label):
    """A prepared-RNA fixture with complete base atoms for every residue in RESIDUES."""
    lines = [f"REMARK E2 fixture {label}"]
    serial = 0
    for residue in RESIDUES:
        names = ("P", "C1'", *export_snapshot.ANCHOR_REFE_BASE_ATOMS[residue["base"]])
        for name in names:
            serial += 1
            element = name[0]
            atom_name = f" {name:<3}" if len(name) < 4 else name
            lines.append(
                f"ATOM  {serial:5d} {atom_name} {'R' + residue['base']:>3} A{residue['resid']:4d}    "
                f"{serial * 0.5:8.3f}{serial * -0.25:8.3f}{serial * 0.125:8.3f}  1.00  0.00"
                f"           {element}"
            )
    lines.append("END")
    return ("\r\n".join(lines) + "\r\n").encode("ascii")


def _inputs():
    """A consistent stage-6 pin cut with one selected, two-candidate anchor."""
    rna_pdb = _pdb("exact prepared RNA")
    protein_pdb = _pdb("prepared range 2")
    fit_rows = [
        {"key": "frag-1-AC", "fragment": 1, "sequence": "AC", "conformer": 1,
         "rotamer": 1, "grid": [0.0, 0.0, 0.0], "rmsd": 0.1, "pdb": "END\n"},
        {"key": "frag-2-CG", "fragment": 2, "sequence": "CG", "conformer": 1,
         "rotamer": 1, "grid": [0.0, 0.0, 0.0], "rmsd": 0.2, "pdb": "END\n"},
    ]
    best_fits = {"schema": 1, "echo": deepcopy(FIT_ECHO), "rows": fit_rows}
    fragments = {
        "schema": 1, "echo": deepcopy(FRAGMENT_ECHO), "residues": deepcopy(RESIDUES),
        "fragments": [
            {"number": 1, "key": "frag-1-AC", "sequence": "AC",
             "residues": deepcopy(RESIDUES[:2])},
            {"number": 2, "key": "frag-2-CG", "sequence": "CG",
             "residues": deepcopy(RESIDUES[1:])},
        ],
        "adjacent_pairs": [{"down": 1, "up": 2, "sequence": "ACG"}],
    }
    pair = {"down": 1, "up": 2, "sequence": "ACG", "ovrmsd": "0.14",
            "crmsd": "0.13", "ovrmsd_raw": 0.130001, "crmsd_raw": 0.129995}
    bases = [
        {"fragment": 1, "sequence": "AC", "nuc1_ovrmsd": "0.20",
         "nuc2_ovrmsd": "2.25", "nuc1_ovrmsd_raw": 0.1999,
         "nuc2_ovrmsd_raw": 2.249999},
        {"fragment": 2, "sequence": "CG", "nuc1_ovrmsd": "0.20",
         "nuc2_ovrmsd": "0.30", "nuc1_ovrmsd_raw": 0.1999,
         "nuc2_ovrmsd_raw": 0.2999},
    ]
    best_fit_thresholds = {"schema": 1, "echo": deepcopy(FIT_ECHO),
                           "pairs": [pair], "bases": bases}
    selections = {
        "schema": 1,
        "basis": {"pdbcode": "1B7F", "biounit": 1, "rows": [RNA, PROTEIN]},
        "active": True,
        "selected": {STACKING_KEY: True},
    }
    stacking_row = {
        "rna_resname": "C", "rna_resid": 102, "protein_range_pdb_id": 2,
        "protein_resid": 214, "distance": 4.0, "tilt_angle": 10.0,
        "slip_angle": 20.0, "rna_centroid": [0.0, 0.0, 0.0],
        "protein_centroid": [1.0, 0.0, 0.0],
    }
    stackings = {
        "schema": 1,
        "echo": {"structure_id": "1B7F", "pdbcode": "1B7F", "biounit": 1,
                 "rows": [RNA, PROTEIN]},
        "rows": [stacking_row],
    }
    matrix = {"schema": 1, "basis": [STACKING_KEY], "active": True,
              "rows": list(FRAGMENT_KEYS), "ranges": ["range-2"],
              "grid": [[True], [True]]}
    interactions = {
        "schema": 1, "clip1": 0, "clip2": 0,
        "kept_residues": deepcopy(RESIDUES),
        "ranges": [{"index": 2, **deepcopy(PROTEIN)}],
        "runs": [{"range": "range-2", "first_fragment": 1, "last_fragment": 2,
                  "anchors": [{"fragment": "frag-1-AC", "nucleotide": 2},
                              {"fragment": "frag-2-CG", "nucleotide": 1}]}],
        "interactions": {key: ["range-2"] for key in FRAGMENT_KEYS},
    }
    protein_range_echo = [deepcopy(PROTEIN)]
    energetics = {
        "schema": 1,
        "echo": {"best_fits": deepcopy(FIT_ECHO), "protein_ranges": protein_range_echo},
        "ranges": [{
            "range": "range-2", "row": deepcopy(PROTEIN),
            "receptor": {"beads": 1},
            "rows": [
                {"key": key, "reference": {"energy": -2.0, "error": None},
                 "bestfit": {"energy": energy, "error": None}}
                for key, energy in zip(FRAGMENT_KEYS, (-1.0, -2.0))
            ],
        }],
    }
    energy_threshold_edits = {
        "schema": 1, "basis": deepcopy(FRAGMENT_ECHO),
        "ranges": [{"row": deepcopy(PROTEIN),
                    "values": {"frag-1-AC": "-1", "frag-2-CG": "-1"}}],
    }
    bestfit_stackings_echo = {"best_fits": deepcopy(FIT_ECHO),
                              "protein_ranges": protein_range_echo}
    bestfit_stackings = {
        "schema": 1, "echo": bestfit_stackings_echo,
        "rows": [
            {"key": STACKING_KEY, "fragment": 1, "nucleotide": "second",
             "range": 2, "protein_resid": 214, "distance": 4.0,
             "phi_deg": 20.0, "theta_deg": 10.0},
            {"key": STACKING_KEY, "fragment": 2, "nucleotide": "first",
             "range": 2, "protein_resid": 214, "distance": 4.0,
             "phi_deg": 20.0, "theta_deg": 10.0},
        ],
    }
    gaussian_model_json = json.dumps({
        "schema": 1, "filter_set_version": "e2-model-v1",
        "thresholds": {name: {"accept": ">=", "cutoff": 0.01, "filters": []}
                       for name in GAUSSIAN_NAMES},
    }, sort_keys=True).encode("utf-8")
    gaussian_scores = {
        "schema": 1,
        "echo": {"bestfit_stackings": deepcopy(bestfit_stackings_echo),
                 "selected_keys": [STACKING_KEY], "clip": deepcopy(CLIP)},
        "filter_set_version": "e2-model-v1", "filter_names": list(GAUSSIAN_NAMES),
        "classes": [{"key": STACKING_KEY, "filter_name": "rrm90",
                     "class_pass": True, "below_rrm95": False}],
        "candidate_scores": [
            {"key": STACKING_KEY, "fragment": 1, "nucleotide": "second", "score": 0.3},
            {"key": STACKING_KEY, "fragment": 2, "nucleotide": "first", "score": 0.3},
        ],
        "reference_scores": [{"key": STACKING_KEY, "score": 0.3}],
    }
    gaussian_check = {
        "schema": 1, "filter_set_version": "e2-model-v1", "min_cell_fraction": 1 / 3,
        "jobs": [
            {"key": STACKING_KEY, "fragment": 1, "nucleotide": "second",
             "filters": {"rrm90": {"fill": 0.2, "emitted": True, "pose_pass": True}}},
            {"key": STACKING_KEY, "fragment": 2, "nucleotide": "first",
             "filters": {"rrm90": {"fill": 0.2, "emitted": True, "pose_pass": True}}},
        ],
    }
    gaussian_validation = {
        "schema": 1,
        "echo": {"classes": deepcopy(CLASS_ROWS), "clip": deepcopy(CLIP),
                 "selected_keys": [STACKING_KEY], "best_fits": deepcopy(FIT_ECHO)},
        "filter_set_version": "e2-model-v1", "filter_names": list(GAUSSIAN_NAMES),
        "rows": [{"key": STACKING_KEY, "filter_name": "rrm90", "worst_state": "amber",
                  "candidates": [
                      {"fragment": 1, "nucleotide": "second", "state": "amber",
                       "fill": 0.2, "emitted": True, "pose_pass": True, "score": 0.3},
                      {"fragment": 2, "nucleotide": "first", "state": "amber",
                       "fill": 0.2, "emitted": True, "pose_pass": True, "score": 0.3},
                  ]}],
    }
    threshold_edits = {
        "schema": 1, "basis": {"pdbcode": "1B7F", "rna": RNA}, "active": True,
        "pairs": [{"down": 1, "up": 2, "sequence": "ACG",
                   "ovrmsd": "0.13", "crmsd": "0.13"}],
    }
    reference_bytes = _rna_pdb("reference")
    protein_bytes = _pdb("protein")
    range_pdbs_aa = {
        "schema": 1, "pdb_code": "1b7f", "range_specs": deepcopy(RANGE_SPECS),
        "pdbs": [
            {"range": "range-1", "name": "range-1-aa.pdb", "pdb": reference_bytes.decode("ascii")},
            {"range": "range-2", "name": "range-2-aa.pdb", "pdb": protein_bytes.decode("ascii")},
        ],
    }
    clipped_fragments = {
        "reference/frag-1-AC-refe.pdb": _pdb("reference frag1"),
        "bestfit/frag-1-AC-bestfit.pdb": b"END\n",
        "reference/frag-2-CG-refe.pdb": _pdb("reference frag2"),
        "bestfit/frag-2-CG-bestfit.pdb": b"END\n",
        "thresholds.json": json.dumps({"schema": 1, "pairs": threshold_edits["pairs"]},
                                       sort_keys=True).encode("utf-8") + b"\n",
        "clip.json": json.dumps({"schema": 1, "clip1": 0, "clip2": 0,
                                 "kept_residues": RESIDUES}, sort_keys=True).encode("utf-8") + b"\n",
    }
    return {
        "pdbcode": "1B7F", "biounit": 1, "cif": "data_1B7F\n_entry.id 1B7F\n",
        "ranges": {"schema": 1, "basis": {"pdbcode": "1B7F", "biounit": 1},
                   "active": True, "rows": [RNA, PROTEIN]},
        "range_specs": deepcopy(RANGE_SPECS), "structure_id": "1B7F",
        "exclude_pdb_code": "1b7f",
        "range_pdb_files": {"range-1.pdb": _pdb("raw RNA"),
                             "range-2.pdb": _pdb("raw protein")},
        "range_pdb_files_aa": {"range-1-aa.pdb": reference_bytes,
                               "range-2-aa.pdb": protein_bytes},
        "range_pdbs_aa": range_pdbs_aa, "rna_pdb": reference_bytes,
        "fragments": fragments, "best_fit_thresholds": best_fit_thresholds,
        "best_fits": best_fits, "clipped_fragments": clipped_fragments,
        "stackings": stackings,
        "stacking_selection": selections, "threshold_edits": threshold_edits,
        "best_fit_clip": deepcopy(CLIP), "interaction_matrix": matrix, "interactions": interactions,
        "energetics": energetics, "energy_threshold_edits": energy_threshold_edits,
        "bestfit_stackings": bestfit_stackings, "gaussian_scores": gaussian_scores,
        "gaussian_validation": gaussian_validation, "gaussian_check": gaussian_check,
        "gaussian_classes": {"schema": 1, "basis": [STACKING_KEY], "active": True,
                              "rows": deepcopy(CLASS_ROWS)},
        "gaussian_model_json": gaussian_model_json,
        "fraglib_identities": "fixture-fraglib-identity-v1",
        "ring_refe_a": b"A", "ring_refe_c": b"C", "ring_refe_phe": b"PHE",
        "ring_refe_his": b"HIS", "ring_refe_arg": b"ARG", "ring_refe_trp": b"TRP",
        "validated_pairs": None, "validated_classes": None, "validated_energies": None,
    }


def _inputs_with_two_stacking_classes():
    inputs = _inputs()
    second_key = json.dumps(["A", 101, 2, 215], separators=(",", ":"))
    inputs["stackings"]["rows"].insert(0, {
        "rna_resname": "A", "rna_resid": 101, "protein_range_pdb_id": 2,
        "protein_resid": 215, "distance": 4.1, "tilt_angle": 12.0,
        "slip_angle": 18.0, "rna_centroid": [0.0, 0.0, 0.0],
        "protein_centroid": [1.0, 0.0, 0.0],
    })
    selected = sorted([STACKING_KEY, second_key])
    inputs["stacking_selection"]["selected"][second_key] = True
    inputs["interaction_matrix"]["basis"] = selected
    inputs["interactions"]["runs"][0]["anchors"].insert(
        0, {"fragment": "frag-1-AC", "nucleotide": 1})

    class_rows = [
        {"key": second_key, "rrm": 80},
        deepcopy(CLASS_ROWS[0]),
    ]
    inputs["gaussian_classes"]["basis"] = selected
    inputs["gaussian_classes"]["rows"] = deepcopy(class_rows)
    inputs["gaussian_validation"]["echo"]["classes"] = deepcopy(class_rows)
    inputs["gaussian_validation"]["echo"]["selected_keys"] = selected
    inputs["gaussian_validation"]["rows"].insert(0, {
        "key": second_key, "filter_name": "rrm80", "worst_state": "green",
        "candidates": [{"fragment": 1, "nucleotide": "first", "state": "green",
                        "fill": 0.5, "emitted": True, "pose_pass": True, "score": 0.4}],
    })
    inputs["gaussian_scores"]["echo"]["selected_keys"] = selected
    inputs["gaussian_scores"]["classes"].insert(0, {
        "key": second_key, "rrm": 80, "filter_name": "rrm80",
        "class_pass": True, "below_rrm95": False,
    })
    inputs["gaussian_scores"]["candidate_scores"].insert(0, {
        "key": second_key, "fragment": 1, "nucleotide": "first", "score": 0.4,
    })
    inputs["bestfit_stackings"]["rows"].insert(0, {
        "key": second_key, "fragment": 1, "nucleotide": "first", "range": 2,
        "protein_resid": 215, "distance": 4.1, "phi_deg": 18.0, "theta_deg": 12.0,
    })
    return inputs


def _normalised(inputs=None):
    snapshot, blocks = export_snapshot.normalise_export_inputs(_inputs() if inputs is None else inputs)
    assert blocks == []
    assert snapshot is not None
    return snapshot


def _reverse_mappings(value):
    if isinstance(value, dict):
        return {key: _reverse_mappings(item) for key, item in reversed(list(value.items()))}
    if isinstance(value, list):
        return [_reverse_mappings(item) for item in value]
    return value


def _assert_blocked(inputs):
    snapshot, blocks = export_snapshot.normalise_export_inputs(inputs)
    assert blocks
    assert snapshot is None or isinstance(snapshot, dict)
    for block in blocks:
        assert set(block) == {"code", "address", "message", "remedy"}
        assert isinstance(block["code"], str) and block["code"]
        assert block["address"] is None or isinstance(block["address"], str)
        assert isinstance(block["message"], str) and block["message"]
        assert isinstance(block["remedy"], str) and block["remedy"]


def test_e2_real_named_pin_shapes_normalise_and_emit_schema_one_admission():
    inputs = _inputs()
    snapshot = _normalised(inputs)
    assert snapshot["schema"] == 1

    admission = export_snapshot.build_admission(inputs)
    assert set(admission) == {"schema", "ready", "blocks", "warnings", "digest"}
    assert admission["schema"] == 1
    assert admission["ready"] is True  # X5 warnings never change readiness.
    assert admission["blocks"] == []
    assert admission["digest"] == export_snapshot.snapshot_digest(snapshot)


def test_e2_stage_six_gate_reports_not_ready_for_undecided_inputs():
    cases = []
    no_pair_decision = _inputs()
    no_pair_decision["threshold_edits"] = None
    cases.append(no_pair_decision)
    no_classes = _inputs()
    no_classes["gaussian_classes"] = None
    cases.append(no_classes)
    inactive_classes = _inputs()
    inactive_classes["gaussian_classes"] = {
        "schema": 1, "basis": [STACKING_KEY], "active": False, "rows": [],
    }
    cases.append(inactive_classes)

    for inputs in cases:
        admission = export_snapshot.build_admission(inputs)
        assert admission["schema"] == 1
        assert admission["ready"] is False
        assert admission["digest"] is None
        assert admission["blocks"]
        assert any(block["code"] == "stage-gate" for block in admission["blocks"])
        assert all("malformed" not in block["message"].lower()
                   for block in admission["blocks"])


def test_e2_x5_warning_order_content_and_unrounded_values():
    snapshot = _normalised()
    warnings = export_snapshot.export_warnings(snapshot)
    assert warnings
    assert all(set(item) == {"kind", "address", "message"} for item in warnings)
    assert all(item["kind"] in {"gaussian", "pair-boundary", "energy"} for item in warnings)
    kinds = [item["kind"] for item in warnings]
    assert kinds[0] == "gaussian"
    assert "pair-boundary" in kinds
    assert kinds[-1] == "energy"
    assert kinds == sorted(kinds, key={"gaussian": 0, "pair-boundary": 1, "energy": 2}.get)
    messages = "\n".join(item["message"] for item in warnings)
    assert "0.130001" in messages and "0.129995" in messages
    assert "amber" in messages
    assert any(item["address"] and "range-2" in item["address"] for item in warnings
               if item["kind"] == "energy")
    # The bases sentinel above L_ov is deliberately not an X5 reason.
    assert all("2.249999" not in item["message"] for item in warnings)

    admission = export_snapshot.build_admission(_inputs())
    assert admission["ready"] is True
    assert admission["blocks"] == []
    assert admission["warnings"] == warnings


def test_e2_x5_implied_gaussian_validation_is_a_warning():
    inputs = _inputs()
    validation = inputs["gaussian_validation"]["rows"][0]
    validation["worst_state"] = "implied"
    validation["candidates"][0]["state"] = "implied"
    validation["candidates"][1]["state"] = "red"

    admission = export_snapshot.build_admission(inputs)
    assert admission["ready"] is True
    assert any(item["kind"] == "gaussian" and "implied" in item["message"]
               for item in admission["warnings"])


def test_e2_x5_computed_class_below_rrm95_warns_even_when_validation_is_green():
    inputs = _inputs()
    inputs["gaussian_validation"]["rows"][0]["worst_state"] = "green"
    for candidate in inputs["gaussian_validation"]["rows"][0]["candidates"]:
        candidate["state"] = "green"
    inputs["gaussian_scores"]["classes"][0]["below_rrm95"] = True

    admission = export_snapshot.build_admission(inputs)
    assert admission["ready"] is True
    assert any(item["kind"] == "gaussian" and "RRM-95" in item["message"]
               for item in admission["warnings"])


@pytest.mark.parametrize("pin", [
    "ranges", "range_specs", "range_pdbs_aa", "fragments", "best_fit_thresholds", "best_fits",
    "stackings", "threshold_edits", "best_fit_clip", "interaction_matrix",
    "interactions", "energetics", "energy_threshold_edits", "bestfit_stackings",
    "gaussian_classes", "gaussian_scores", "gaussian_validation", "gaussian_check",
])
def test_e2_strict_schema_one_boundary_jsons(pin):
    inputs = _inputs()
    inputs[pin]["schema"] = 2
    _assert_blocked(inputs)
    admission = export_snapshot.build_admission(inputs)
    assert admission["schema"] == 1
    assert admission["ready"] is False
    assert admission["blocks"]
    assert admission["digest"] is None
    assert admission["warnings"] == []


def test_e2_gaussian_model_boundary_requires_schema_one():
    inputs = _inputs()
    model = json.loads(inputs["gaussian_model_json"])
    model["schema"] = 2
    inputs["gaussian_model_json"] = json.dumps(model).encode("utf-8")
    _assert_blocked(inputs)


def test_e2_addressed_hard_blocks_for_lost_anchor_lock_and_exclude_mismatch():
    inputs = _inputs()
    inputs["interaction_matrix"]["grid"][0][0] = False
    _assert_blocked(inputs)
    result = export_snapshot.build_admission(inputs)
    assert result["ready"] is False
    assert any(block["address"] and "frag-1-AC" in block["address"]
               for block in result["blocks"])

    inputs = _inputs()
    inputs["exclude_pdb_code"] = "1B7F"
    _assert_blocked(inputs)


def test_e2_freshness_echoes_and_exact_prepared_reference_bytes_are_hard_gates():
    inputs = _inputs()
    inputs["gaussian_scores"]["echo"]["bestfit_stackings"]["best_fits"]["pdb_code"] = "stale"
    _assert_blocked(inputs)

    inputs = _inputs()
    inputs["range_pdb_files_aa"]["range-1-aa.pdb"] += b"changed"
    _assert_blocked(inputs)


def test_e2_b11_clip_residue_list_must_match_current_kept_fragments():
    inputs = _inputs()
    clip_record = json.loads(inputs["clipped_fragments"]["clip.json"])
    clip_record["kept_residues"].pop()
    inputs["clipped_fragments"]["clip.json"] = json.dumps(
        clip_record, sort_keys=True).encode("utf-8") + b"\n"
    _assert_blocked(inputs)


def test_e2_selected_stacking_base_must_match_fragment_residue_identity():
    inputs = _inputs()
    for fragment in inputs["fragments"]["fragments"]:
        for residue in fragment["residues"]:
            if residue["resid"] == 102:
                residue["base"] = "A"
    snapshot, blocks = export_snapshot.normalise_export_inputs(inputs)
    assert blocks
    assert snapshot is None or isinstance(snapshot, dict)
    assert any(block["address"] == STACKING_KEY for block in blocks)


@pytest.mark.parametrize(("field", "value"), [
    ("structure_id", "OTHER"),
    ("pdbcode", "OTHER"),
    ("biounit", 2),
    ("rows", [RNA, {"chain": "C", "start": 1, "end": 20}]),
])
def test_e2_stacking_echo_must_match_current_structure_and_ranges(field, value):
    inputs = _inputs()
    inputs["stackings"]["echo"][field] = value
    snapshot, blocks = export_snapshot.normalise_export_inputs(inputs)
    assert blocks
    assert snapshot is None or isinstance(snapshot, dict)
    assert any(isinstance(block.get("address"), str) and block["address"]
               for block in blocks)


@pytest.mark.parametrize("change", ["run-endpoint", "run-anchor", "kept-residue"])
def test_e2_interaction_handoff_must_match_matrix_and_kept_fragments(change):
    inputs = _inputs()
    original_matrix = deepcopy(inputs["interaction_matrix"])
    handoff = inputs["interactions"]
    if change == "run-endpoint":
        handoff["runs"][0]["last_fragment"] = 1
    elif change == "run-anchor":
        handoff["runs"][0]["anchors"].pop()
    else:
        handoff["kept_residues"].pop()

    assert inputs["interaction_matrix"] == original_matrix
    snapshot, blocks = export_snapshot.normalise_export_inputs(inputs)
    assert blocks
    assert snapshot is None or isinstance(snapshot, dict)
    assert any(isinstance(block.get("address"), str) and block["address"]
               for block in blocks)


@pytest.mark.parametrize("change", ["missing", "duplicate"])
def test_e2_gaussian_validation_has_exactly_one_current_row_per_class(change):
    inputs = _inputs()
    current_echo = deepcopy(inputs["gaussian_validation"]["echo"])
    rows = inputs["gaussian_validation"]["rows"]
    if change == "missing":
        rows.clear()
    else:
        rows.append(deepcopy(rows[0]))

    assert inputs["gaussian_validation"]["echo"] == current_echo
    snapshot, blocks = export_snapshot.normalise_export_inputs(inputs)
    assert blocks
    assert snapshot is None or isinstance(snapshot, dict)
    assert any(isinstance(block.get("address"), str) and block["address"]
               for block in blocks)


@pytest.mark.parametrize("mutation", [
    lambda value: value["best_fit_thresholds"]["pairs"][0].__setitem__("ovrmsd_raw", True),
    lambda value: value["energetics"]["ranges"][0]["rows"][0]["bestfit"].__setitem__("energy", True),
    lambda value: value["gaussian_check"].__setitem__("min_cell_fraction", True),
])
def test_e2_boolean_values_are_not_numeric_measurements(mutation):
    inputs = _inputs()
    mutation(inputs)
    _assert_blocked(inputs)


def test_e2_zero_raw_rmsds_are_valid_measurements():
    # Alaric's strictly_above_cent rounds a raw 0 to "0.01". 1CVJ has cRMSD
    # 0.000000 for three adjacent best-fit pairs.
    inputs = _inputs()
    pair = inputs["best_fit_thresholds"]["pairs"][0]
    pair["crmsd_raw"], pair["crmsd"] = 0.0, "0.01"
    base = inputs["best_fit_thresholds"]["bases"][0]
    base["nuc1_ovrmsd_raw"], base["nuc1_ovrmsd"] = 0.0, "0.01"
    admission = export_snapshot.build_admission(inputs)
    assert admission["blocks"] == []
    assert admission["ready"] is True


@pytest.mark.parametrize("mutation", [
    lambda value: value["best_fit_thresholds"]["pairs"][0].__setitem__("crmsd_raw", -0.001),
    lambda value: value["best_fit_thresholds"]["pairs"][0].__setitem__("ovrmsd_raw", -0.001),
    lambda value: value["best_fit_thresholds"]["bases"][0].__setitem__("nuc1_ovrmsd_raw", -0.001),
])
def test_e2_negative_raw_rmsds_are_blocked(mutation):
    inputs = _inputs()
    mutation(inputs)
    _assert_blocked(inputs)


@pytest.mark.parametrize("bad_value", [None, {"schema": 2}])
def test_e2_null_or_bad_schema_computed_thresholds_become_admission_blocks(bad_value):
    inputs = _inputs()
    inputs["best_fit_thresholds"] = bad_value
    admission = export_snapshot.build_admission(inputs)
    assert admission["ready"] is False
    assert admission["digest"] is None
    assert admission["blocks"]


@pytest.mark.parametrize("mutation", [
    lambda value: value["energetics"]["ranges"][0]["rows"][0].__setitem__("key", []),
    lambda value: value["energetics"]["ranges"][0]["rows"][0]["bestfit"].__setitem__("energy", 10 ** 400),
])
def test_e2_malformed_or_overflowing_energy_rows_become_addressed_blocks(mutation):
    inputs = _inputs()
    mutation(inputs)
    snapshot, blocks = export_snapshot.normalise_export_inputs(inputs)
    assert blocks
    assert snapshot is None or isinstance(snapshot, dict)
    assert any(isinstance(block["address"], str) and "frag-1-AC" in block["address"]
               for block in blocks)


@pytest.mark.parametrize("mutation", [
    lambda value: value["gaussian_check"].__setitem__("min_cell_fraction", 10 ** 400),
    lambda value: value["gaussian_validation"]["filter_names"].__setitem__(0, ["rrm50"]),
    lambda value: value["gaussian_scores"]["filter_names"].__setitem__(0, ["rrm50"]),
])
def test_e2_malformed_gaussian_check_and_filter_names_become_admission_blocks(mutation):
    inputs = _inputs()
    mutation(inputs)

    snapshot, blocks = export_snapshot.normalise_export_inputs(inputs)
    assert blocks
    assert snapshot is None or isinstance(snapshot, dict)

    admission = export_snapshot.build_admission(inputs)
    assert admission["ready"] is False
    assert admission["digest"] is None
    assert admission["blocks"]


def test_e2_snapshot_digest_canonicalizes_mapping_order_and_tracks_bytes_and_version():
    snapshot = _normalised()
    digest = export_snapshot.snapshot_digest(snapshot)
    assert isinstance(digest, str) and len(digest) == 64
    int(digest, 16)
    assert export_snapshot.snapshot_digest(_reverse_mappings(snapshot)) == digest

    changed = _inputs()
    changed["rna_pdb"] += b"\x00"
    changed["range_pdb_files_aa"]["range-1-aa.pdb"] += b"\x00"
    changed["range_pdbs_aa"]["pdbs"][0]["pdb"] += "\x00"
    changed_digest = export_snapshot.snapshot_digest(_normalised(changed))
    assert changed_digest != digest

    old_version = export_snapshot.EXPORTER_VERSION
    try:
        export_snapshot.EXPORTER_VERSION = old_version + 1
        assert export_snapshot.snapshot_digest(snapshot) != digest
    finally:
        export_snapshot.EXPORTER_VERSION = old_version


def test_e2_snapshot_digest_golden_is_tied_to_exporter_version_and_byte_fingerprint():
    assert export_snapshot.EXPORTER_VERSION == DIGEST_GOLDEN_EXPORTER_VERSION
    assert export_snapshot.snapshot_digest({"schema": 1, "reference_pdb": b"E2"}) \
        == DIGEST_GOLDEN_SHA256


def test_e2_snapshot_digest_ignores_order_of_semantically_unordered_gaussian_rows():
    inputs = _inputs_with_two_stacking_classes()
    first = _normalised(inputs)
    reordered = deepcopy(inputs)
    reordered["gaussian_classes"]["rows"].reverse()
    reordered["gaussian_validation"]["echo"]["classes"].reverse()
    reordered["gaussian_validation"]["rows"].reverse()
    reordered["gaussian_scores"]["classes"].reverse()
    reordered["gaussian_scores"]["candidate_scores"].reverse()
    reordered["bestfit_stackings"]["rows"].reverse()
    for validation_row in reordered["gaussian_validation"]["rows"]:
        validation_row["candidates"].reverse()
    second = _normalised(reordered)
    assert export_snapshot.snapshot_digest(first) == export_snapshot.snapshot_digest(second)


def test_e2_validator_null_pins_and_view_state_are_excluded_from_snapshot_identity():
    inputs = _inputs()
    snapshot = _normalised(inputs)
    for name in ("validated_pairs", "validated_classes", "validated_energies",
                 "stage4_active_nucleotide"):
        assert name not in snapshot
    inputs["stage4_active_nucleotide"] = 3
    changed_snapshot = _normalised(inputs)
    assert export_snapshot.snapshot_digest(changed_snapshot) == export_snapshot.snapshot_digest(snapshot)


def test_best_fit_statistics_collect_kept_rmsd_energies_and_pair_values():
    inputs = _inputs()
    inputs["energetics"]["ranges"][0]["rows"][0]["bestfit"]["energy"] = -1.5
    snapshot = _normalised(inputs)
    expected_energies = [{"range-2": {"energy": -1.5, "error": None}},
                         {"range-2": {"energy": -2.0, "error": None}}]
    statistics = snapshot["best_fit_statistics"]
    assert statistics == {
        "schema": 1,
        "fragments": [
            {"fragment": "frag-1-AC", "rmsd": 0.1,
             "nuc1_ovrmsd_raw": 0.1999, "nuc1_ovrmsd": "0.20",
             "nuc2_ovrmsd_raw": 2.249999, "nuc2_ovrmsd": "2.25",
             "energies": expected_energies[0]},
            {"fragment": "frag-2-CG", "rmsd": 0.2,
             "nuc1_ovrmsd_raw": 0.1999, "nuc1_ovrmsd": "0.20",
             "nuc2_ovrmsd_raw": 0.2999, "nuc2_ovrmsd": "0.30",
             "energies": expected_energies[1]},
        ],
        "pairs": [
            {"down": 1, "up": 2, "sequence": "ACG",
             "crmsd_raw": 0.129995, "crmsd_computed": "0.13", "crmsd": "0.13",
             "ovrmsd_raw": 0.130001, "ovrmsd_computed": "0.14", "ovrmsd": "0.13"},
        ],
    }


def test_best_fit_statistics_keep_a_scoring_error_of_a_non_interacting_cell():
    # A false matrix cell may carry a scoring error instead of an energy.
    inputs = _inputs()
    inputs["energetics"]["ranges"][0]["rows"][1]["bestfit"] = {"energy": None, "error": "minimisation failed"}
    statistics = export_snapshot._best_fit_statistics(
        inputs["best_fits"], inputs["energetics"], inputs["best_fit_thresholds"],
        [row for row in inputs["fragments"]["fragments"]], [],
    )
    assert [row["energies"] for row in statistics["fragments"]] == [
        {"range-2": {"energy": -1.0, "error": None}},
        {"range-2": {"energy": None, "error": "minimisation failed"}},
    ]
    assert statistics["pairs"] == []


def test_best_fit_statistics_block_a_kept_fragment_without_finite_rmsd():
    inputs = _inputs()
    inputs["best_fits"]["rows"][1]["rmsd"] = float("nan")
    snapshot, blocks = export_snapshot.normalise_export_inputs(inputs)
    assert snapshot is None
    assert [(block["code"], block["address"]) for block in blocks] == [("best-fit-rmsd", "frag-2-CG")]


def _replace_reference(inputs, reference):
    inputs["rna_pdb"] = reference
    inputs["range_pdb_files_aa"]["range-1-aa.pdb"] = reference
    inputs["range_pdbs_aa"]["pdbs"][0]["pdb"] = reference.decode("ascii")


@pytest.mark.parametrize("change", ["missing-atom", "duplicate-atom", "shifted-residue", "wrong-base"])
def test_reference_nucleobases_block_an_unusable_anchor_refe_reference(change):
    inputs = _inputs()
    lines = inputs["rna_pdb"].decode("ascii").split("\r\n")
    if change == "missing-atom":
        lines = [line for line in lines if not (line[22:26] == " 102" and line[12:16] == " N4 ")]
        message = "lacks base atom N4"
    elif change == "duplicate-atom":
        lines.insert(lines.index(next(line for line in lines if line[12:16] == " N4 ")),
                     next(line for line in lines if line[12:16] == " N4 "))
        message = "repeats base atom N4"
    elif change == "shifted-residue":
        lines = [line for line in lines if line[22:26] != " 101"]
        message = "not A A 101"
    else:
        lines = [line.replace(" RC A 102", " RU A 102") for line in lines]
        message = "not C A 102"
    _replace_reference(inputs, "\r\n".join(lines).encode("ascii"))
    snapshot, blocks = export_snapshot.normalise_export_inputs(inputs)
    assert snapshot is None
    assert len(blocks) == 1 and blocks[0]["code"] == "reference-nucleobase"
    assert message in blocks[0]["message"]


def test_reference_nucleobases_return_verbatim_records_in_template_order():
    reference = _rna_pdb("reference")
    fragments = _inputs()["fragments"]["fragments"]
    bases = export_snapshot.reference_nucleobases(reference, fragments, [3, 2])
    assert list(bases) == [2, 3]
    assert bases[3]["residue"] == RESIDUES[2]
    source = reference.decode("ascii").split("\r\n")
    assert [line[12:16].strip() for line in bases[3]["lines"]] \
        == list(export_snapshot.ANCHOR_REFE_BASE_ATOMS["G"])
    assert all(line in source for line in bases[3]["lines"])


def test_anchor_refe_base_atoms_mirror_alaric_only_base_template_selection():
    np = pytest.importorskip("numpy")
    from alaric_dir import resolve_alaric_repo_root

    try:
        templates = resolve_alaric_repo_root() / "fraglib" / "templates"
    except FileNotFoundError:
        pytest.skip("The Alaric checkout is not available.")
    if not templates.is_dir():
        pytest.skip("The Alaric fragment-library templates are not available.")
    for base, names in export_snapshot.ANCHOR_REFE_BASE_ATOMS.items():
        template = np.load(templates / f"{base}-ppdb.npy", allow_pickle=False)
        # Alaric's Library.create(only_base=True): template atoms with two-character names.
        expected = tuple(name.decode() for name in template["name"] if len(name.decode()) == 2)
        assert names == expected, base
