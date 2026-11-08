"""W5 pure acceptance vectors; authored separately from implementation, not run.

Authority: OC P5.3–P5.8, D3/D12/D13; Alaric tests/W45_ACCEPTANCE.md M/E.
Browser-only caps/gates/write-order checks remain in that acceptance protocol.
"""
from copy import deepcopy
import json
from types import SimpleNamespace

import pytest

import backend
from stage5 import (
    adapt_interaction_matrix, candidate_anchors, columns_to_matrix,
    interaction_handoff, make_interaction_matrix, matrix_to_columns,
    recompute_interaction_matrix, reconcile_energy_threshold_edits,
    validate_energy_threshold_edits, validate_interaction_matrix,
)


KEYS = [f"frag-{i}-GU" for i in range(1, 10)]
RANGES = ["range-2", "range-3", "range-4"]
RNA = {"chain": "A", "start": 101, "end": 119}
SPECS = {"schema": 1, "rows": [RNA, *[
    {"chain": chain, "start": 1, "end": 20} for chain in "BCD"
]]}


def fragments():
    # Deliberately non-1 residue labels; no lock may use resid-1 arithmetic.
    return {"schema": 1, "echo": {"pdbcode": "1b7f", "rna": deepcopy(RNA)},
            "fragments": [{"key": key, "residues": [
                {"chain": "A", "resid": 101 + 2 * i, "base": "G"},
                {"chain": "A", "resid": 102 + 2 * i, "base": "U"},
            ]} for i, key in enumerate(KEYS)]}


def anchor(n, range_id="range-2", locks=None):
    return {"resid": 99 + 2 * n, "range": range_id, "stackings": [f"a{n}"],
            "locks": [{"fragment": KEYS[i - 1], "nucleotide": 1}
                      for i in (locks or [n])]}


def selection_and_stackings(*addresses):
    rows = [{"rna_resname": "G", "rna_resid": resid,
             "protein_range_pdb_id": rng, "protein_resid": protein}
            for resid, rng, protein in addresses]
    selected = {json.dumps(["G", resid, rng, protein], separators=(",", ":")): True
                for resid, rng, protein in addresses}
    return {"active": True, "selected": selected}, {"rows": rows}


def matrix(rows=None, columns=None):
    return {"schema": 1, "basis": ["a4"], "active": True,
            "rows": list(KEYS if rows is None else rows), "ranges": list(RANGES if columns is None else columns),
            "grid": [[False] * len(RANGES if columns is None else columns) for _ in (KEYS if rows is None else rows)]}


def energies(high=False):
    values = [(2.999, 0), (2.999, -2.999), (3, 0), (0, 0), (0, -3),
              (25, 0), None, (-4, 0), (0, 0)]
    result = {"schema": 1, "ranges": []}
    for j, range_id in enumerate(RANGES):
        rows = []
        for i, key in enumerate(KEYS):
            value = (3, 3) if high else values[i] if j == 0 else (20, 0) if j == 1 and i == 4 else (0, 0)
            rows.append({"key": key, **{
                label: {"energy": None if value is None else value[k],
                        "error": "bad ligand" if value is None else None}
                for k, label in enumerate(("reference", "bestfit"))}})
        result["ranges"].append({"range": range_id, "row": SPECS["rows"][j + 1], "rows": rows})
    return result


def true_rows(value, range_id):
    columns = [i for i, name in enumerate(value["ranges"]) if name == range_id]
    return [key for key, row in zip(value["rows"], value["grid"]) if any(row[i] for i in columns)]


def test_m1_exact_energy_boundary_stop_row_and_unanchored_ranges():
    actual = recompute_interaction_matrix(matrix(), [anchor(4)], energies())
    assert true_rows(actual, "range-2") == KEYS[2:6]
    assert true_rows(actual, "range-3") == []
    assert true_rows(actual, "range-4") == []
    assert validate_interaction_matrix(actual, KEYS, RANGES, [anchor(4)], energies()) == actual


def test_t13_range_scorers_pin_receptors_from_the_protein_folder_result(monkeypatch, tmp_path):
    class Folder:
        def __init__(self, files=None):
            self.files = files or {}

        def __getitem__(self, key):
            return self.files[key]

        def __setitem__(self, key, value):
            self.files[key] = value

    class Scorer:
        def __init__(self):
            self.pins = SimpleNamespace()
            self.celltypes = SimpleNamespace()
            self.result = {"energetics.json": object()}

    receptors = {f"range-{index}-aa.pdb": object() for index in (2, 3)}
    context = SimpleNamespace(
        range_specs=SimpleNamespace(value={"schema": 1, "rows": [RNA, *SPECS["rows"][1:3]]}),
        energetics_ligands=SimpleNamespace(result=object()),
        protein_pdbs=SimpleNamespace(result=Folder(receptors)),
        energetics_range_results=Folder(),
    )
    scorers = []

    def make_scorer(_alaric_dir):
        scorer = Scorer()
        scorers.append(scorer)
        return scorer

    class OnePassStop:
        calls = 0

        def wait(self, _timeout):
            self.calls += 1
            return self.calls > 1

    monkeypatch.setattr(backend, "energetics_range_transformer", make_scorer)
    backend._reconcile_stage5_range_scorers(context, tmp_path, OnePassStop())

    assert [scorer.pins.receptor for scorer in scorers] == [receptors["range-2-aa.pdb"], receptors["range-3-aa.pdb"]]
    assert all(scorer.pins.ligands is context.energetics_ligands.result for scorer in scorers)
    assert set(context.energetics_range_results.files) == {"range-2.json", "range-3.json"}


def test_t14_merge_imports_json_in_an_isolated_transformer_namespace(monkeypatch):
    import inspect
    import sys
    import textwrap
    import stage5 as preparic_stage5

    monkeypatch.setitem(sys.modules, "preparic_stage5", preparic_stage5)
    namespace = {}
    source = textwrap.dedent(inspect.getsource(backend._stage5_merge_energetics))
    exec(source, namespace)
    merge = namespace["_stage5_merge_energetics"]
    assert "json" not in merge.__globals__

    energy_result = {
        "schema": 1,
        "receptor": {"beads": 4},
        "pairs": [
            {"ligand": "frag-1-GU-refe.pdb", "energy": -1.0, "error": None},
            {"ligand": "frag-1-GU-bestfit.pdb", "energy": -2.0, "error": None},
        ],
    }
    result = merge(
        {"range-2.json": json.dumps(energy_result).encode("utf-8")},
        {"schema": 1, "rows": [RNA, SPECS["rows"][1]]},
        {"echo": {"fit": "current"}},
        {"schema": 1, "fragments": [{"key": "frag-1-GU"}]},
    )

    assert result["ranges"][0]["rows"] == [{
        "key": "frag-1-GU",
        "reference": {"energy": -1.0, "error": None},
        "bestfit": {"energy": -2.0, "error": None},
    }]


def test_t15_is_a_sentinel_until_s14_is_active():
    inactive = matrix()
    inactive["active"] = False
    assert interaction_handoff(None, None, None, None, None, None) is None
    assert interaction_handoff(inactive, None, None, None, None, None) is None


def test_m2_independent_ranges_can_overlap():
    anchors = [anchor(4), anchor(5, "range-3")]
    actual = recompute_interaction_matrix(matrix(), anchors, energies())
    assert actual["ranges"] == RANGES
    assert true_rows(actual, "range-2") == KEYS[2:6]
    assert true_rows(actual, "range-3") == [KEYS[4]]


def test_m2_candidate_collapse_keeps_original_stacking_keys_and_mapped_nucleotides():
    frags = fragments()
    # A residue can be the second nucleotide of one fragment and first of next.
    frags["fragments"][2]["residues"][1]["resid"] = 107
    selected, stacks = selection_and_stackings((107, 2, 8), (107, 2, 9), (107, 3, 8))
    actual = candidate_anchors(frags, selected, stacks, KEYS)
    assert len(actual) == 2
    assert actual[0]["locks"] == [{"fragment": KEYS[2], "nucleotide": 2},
                                   {"fragment": KEYS[3], "nucleotide": 1}]
    assert len(actual[0]["stackings"]) == 2
    assert set(actual[0]["stackings"] + actual[1]["stackings"]) == set(selected["selected"])
    # Equal residue/range anchors collapse, distinct residue/range anchors do not.
    selected_before = deepcopy(selected)
    empty = make_interaction_matrix(frags, {"clip1": 0, "clip2": 0}, SPECS, selected, stacks)
    assert empty["ranges"] == RANGES
    assert empty["basis"] == sorted(selected["selected"])
    assert selected == selected_before  # Original Gaussian job identities survive.
    assert candidate_anchors(frags, selected, stacks, KEYS[3:])[0]["locks"] == [
        {"fragment": KEYS[3], "nucleotide": 1}]
    assert candidate_anchors(frags, selected, stacks, KEYS[:3])[0]["locks"] == [
        {"fragment": KEYS[2], "nucleotide": 2}]
    with pytest.raises(ValueError, match="clipped|unmapped"):
        candidate_anchors(frags, selected, stacks, KEYS[4:])


def test_s14_multi_anchor_per_range_uses_one_column_and_keeps_each_lock_in_runs():
    anchors = [anchor(3, "range-2"), anchor(8, "range-2"),
               anchor(4, "range-3"), anchor(8, "range-3")]
    value = matrix()
    value["basis"] = [f"stack-{index}" for index in range(4)]
    # Two separated runs for each range, with one selected anchor lock in each.
    value["grid"][2][0] = value["grid"][7][0] = True
    value["grid"][3][1] = value["grid"][7][1] = True
    assert value["ranges"] == RANGES
    assert all(len(row) == len(RANGES) for row in value["grid"])
    assert validate_interaction_matrix(value, KEYS, RANGES, anchors, energies(high=True)) == value
    columns = matrix_to_columns(value, anchors)
    assert [column["range"] for column in columns] == RANGES
    assert len(columns[0]["runs"]) == 2
    assert len(columns[1]["runs"]) == 2
    assert columns[2]["runs"] == []
    assert columns_to_matrix(columns, KEYS, RANGES, value["basis"]) == value


def test_anchor_locks_keep_unique_protein_range_columns_and_repeat_anchor_ranges():
    frags = fragments()
    selected, stacks = selection_and_stackings((115, 2, 8), (107, 3, 9), (105, 2, 10), (105, 2, 11))
    actual = make_interaction_matrix(frags, {"clip1": 0, "clip2": 0}, SPECS, selected, stacks)
    assert actual["ranges"] == RANGES
    assert len(actual["grid"]) == 9
    assert all(row == [False, False, False] for row in actual["grid"])
    assert len(actual["basis"]) == 4  # four original Gaussian anchors, three matrix anchors
    assert len(actual["ranges"]) == len(set(actual["ranges"]))
    anchors = candidate_anchors(frags, selected, stacks, KEYS)
    assert [item["resid"] for item in anchors] == [105, 107, 115]
    assert [len(item["stackings"]) for item in anchors] == [2, 1, 1]


@pytest.mark.parametrize("left,right", [([3, 4], [5, 6]), ([3, 4], [4, 5])])
def test_m3_touching_or_overlapping_locks_form_one_run(left, right):
    anchors = [anchor(3, locks=left), anchor(5, locks=right)]
    result = recompute_interaction_matrix(matrix(), anchors, energies(high=True))
    assert result["ranges"] == RANGES
    assert all(len(row) == len(RANGES) for row in result["grid"])
    columns = matrix_to_columns(result, anchors)
    assert len(columns) == len(RANGES)
    assert [column["range"] for column in columns] == RANGES
    assert {(lock["fragment"], lock["nucleotide"]) for lock in columns[0]["runs"][0]["anchors"]} \
        == {(KEYS[index - 1], 1) for index in set(left + right)}
    assert true_rows(result, "range-2") == KEYS
    assert columns_to_matrix(columns, KEYS, result["ranges"], result["basis"]) == result


def test_m4_roundtrip_separated_joined_and_cross_range_runs():
    # Multiple runs in range-2 share one column, alongside the range-3 column.
    anchors = [anchor(3), anchor(4, "range-3"), anchor(8)]
    column_ranges = RANGES
    value = matrix(columns=column_ranges)
    for n in (2, 3, 4):
        value["grid"][n - 1][0] = True
    for n in (3, 4, 5):
        value["grid"][n - 1][1] = True
    for n in (7, 8):
        value["grid"][n - 1][2] = True
    columns = matrix_to_columns(value, anchors)
    assert len(columns) == len(RANGES)
    assert [column["range"] for column in columns] == RANGES
    assert value["ranges"] == column_ranges
    assert columns_to_matrix(columns, KEYS, column_ranges, value["basis"]) == value


@pytest.mark.parametrize("case", ["row-count", "ragged", "int-zero", "int-one", "null",
    "string", "missing-row", "duplicate-row", "reverse-rows", "range-order", "schema",
    "lock", "unanchored-run", "true-error"])
def test_m5_each_corrupted_read_is_refused_without_mutation(case):
    anchors, energy = [anchor(4), anchor(5, "range-3")], energies()
    value = recompute_interaction_matrix(matrix(), anchors, energy)
    if case == "row-count": value["grid"].pop()
    elif case == "ragged": value["grid"][0].pop()
    elif case in ("int-zero", "int-one", "null", "string"):
        value["grid"][0][0] = {"int-zero": 0, "int-one": 1, "null": None, "string": "true"}[case]
    elif case == "missing-row": value["rows"].pop()
    elif case == "duplicate-row": value["rows"][0] = value["rows"][1]
    elif case == "reverse-rows": value["rows"].reverse()
    elif case == "range-order": value["ranges"].reverse()
    elif case == "schema": value["schema"] = 2
    elif case == "lock": value["grid"][3][0] = False
    elif case == "unanchored-run": value["grid"][0][1] = True
    elif case == "true-error": value["grid"][6][0] = True
    before = deepcopy(value)
    with pytest.raises(ValueError):
        validate_interaction_matrix(value, KEYS, RANGES, anchors, energy)
    assert value == before


def test_m5_error_at_anchor_lock_names_anchor_and_never_drops_it():
    anchors = [anchor(7)]
    before = deepcopy(anchors)
    with pytest.raises(ValueError, match=str(anchors[0]["resid"])):
        recompute_interaction_matrix(matrix(), anchors, energies())
    assert anchors == before


def test_m5_opposite_end_clip_with_same_count_is_refused():
    value = matrix(KEYS[1:])
    with pytest.raises(ValueError, match="extent|order"):
        validate_interaction_matrix(value, KEYS[:-1], RANGES, [], energies())


@pytest.mark.parametrize("stop", [None, "low", "error"])
def test_m6_only_boundary_touching_outer_ends_regrow(stop):
    old = matrix(KEYS[2:7])
    for row in old["grid"]:
        row[0] = True
    for row in old["grid"][1:-1]:
        row[1] = True  # user-capped range never touched either old boundary
    anchors = [anchor(4), anchor(5, "range-3")]
    energy = energies(high=True)
    if stop:
        cell = energy["ranges"][0]["rows"][1]
        for label in ("reference", "bestfit"):
            cell[label] = {"energy": None if stop == "error" else 2.999,
                           "error": "bad" if stop == "error" else None}
    actual = adapt_interaction_matrix(old, KEYS, anchors, energy)
    assert true_rows(actual, "range-2") == (KEYS if stop is None else KEYS[2:])
    assert true_rows(actual, "range-3") == KEYS[3:6]
    assert actual["grid"][2:7] == old["grid"]
    validate_interaction_matrix(actual, KEYS, RANGES, anchors, energy)
    assert adapt_interaction_matrix(actual, KEYS[2:7], anchors, energy) == old
    assert adapt_interaction_matrix(old, old["rows"], anchors, energy) == old
    with pytest.raises(ValueError, match="lock|anchor"):
        adapt_interaction_matrix(old, KEYS[5:], anchors, energy)


def test_m6_adapts_shifted_extent_and_keeps_the_inner_caps():
    old = matrix(KEYS[2:7])
    for index, row in enumerate(old["grid"]):
        row[0] = index < 4  # the 3' cap stopped short of the old boundary
        row[1] = index > 0  # the 5' cap stopped short; the 3' end touched the boundary
    anchors = [anchor(4), anchor(5, "range-3")]
    energy = energies(high=True)
    inner = deepcopy(old["grid"][1:4])
    shifted_rows = KEYS[1:6]  # adds frag-2 and drops frag-7
    shifted = adapt_interaction_matrix(old, shifted_rows, anchors, energy)
    assert shifted["rows"] == shifted_rows
    assert [row[0] for row in shifted["grid"]] == [True, True, True, True, True]
    assert [row[1] for row in shifted["grid"]] == [False, False, True, True, True]
    assert shifted["grid"][2:5] == inner
    validate_interaction_matrix(shifted, shifted_rows, RANGES, anchors, energy)


def test_m6_refuses_legacy_duplicate_anchor_columns_without_normalizing():
    old = matrix(KEYS[2:7], columns=["range-2", "range-2"])
    anchors = [anchor(4), anchor(5)]
    before = deepcopy(old)
    with pytest.raises(ValueError, match="grid|range"):
        adapt_interaction_matrix(old, old["rows"], anchors, energies(high=True))
    assert old == before


def test_handoff_contains_only_kept_rows_and_actual_candidate_anchors():
    frags = fragments()
    selected, stacks = selection_and_stackings((107, 2, 8))
    clip = {"clip1": 1, "clip2": 1}
    value = make_interaction_matrix(frags, clip, SPECS, selected, stacks)
    anchors = candidate_anchors(frags, selected, stacks, value["rows"])
    assert value["ranges"] == RANGES
    value = recompute_interaction_matrix(value, anchors, energies())
    actual = interaction_handoff(value, frags, selected, stacks, clip, SPECS)
    assert actual == {
        "schema": 1, "clip1": 1, "clip2": 1,
        "kept_residues": [res for row in frags["fragments"][1:-1] for res in row["residues"]],
        "ranges": [{"index": i, **row} for i, row in enumerate(SPECS["rows"][1:], 2)],
        "runs": [{"range": "range-2", "first_fragment": 3, "last_fragment": 6,
                  "anchors": [{"fragment": KEYS[3], "nucleotide": 1}]}],
        "interactions": {key: ["range-2"] for key in KEYS[2:6]},
    }


def test_t15_handoff_defers_anchor_admission_to_export():
    """T15 accepts structurally valid S14 and leaves lock admission for export."""
    frags = fragments()
    selected, stacks = selection_and_stackings((107, 2, 8))
    clip = {"clip1": 0, "clip2": 0}
    value = make_interaction_matrix(frags, clip, SPECS, selected, stacks)
    anchors = candidate_anchors(frags, selected, stacks, value["rows"])
    # This has valid schema, addressed extent, ranges, and boolean cells, but
    # the selected row is not the anchor lock (fragment 4).
    value["grid"][0][0] = True
    with pytest.raises(ValueError, match="locked|anchor"):
        validate_interaction_matrix(value, KEYS, value["ranges"], anchors, energies())

    handoff = interaction_handoff(value, frags, selected, stacks, clip, SPECS)
    assert handoff["runs"] == [{
        "range": "range-2", "first_fragment": 1, "last_fragment": 1, "anchors": [],
    }]
    assert handoff["interactions"] == {KEYS[0]: ["range-2"]}


def test_e1_empty_skeleton_and_row_based_reconciliation():
    frags = fragments()
    edits = reconcile_energy_threshold_edits(None, frags, SPECS)
    assert edits["basis"] == frags["echo"]
    assert len(edits["ranges"]) == 3
    assert all(entry["values"] == dict.fromkeys(KEYS) for entry in edits["ranges"])
    edits["ranges"][0]["values"][KEYS[0]] = "-16.60"
    edits["ranges"][1]["values"][KEYS[1]] = "21.00"
    reordered = {"schema": 1, "rows": [RNA, SPECS["rows"][2], SPECS["rows"][1],
        {"chain": "Z", "start": 5, "end": 10}]}
    result = reconcile_energy_threshold_edits(edits, frags, reordered)
    assert result["ranges"][0]["values"][KEYS[1]] == "21.00"
    assert result["ranges"][1]["values"][KEYS[0]] == "-16.60"
    assert result["ranges"][2]["values"] == dict.fromkeys(KEYS)
    changed = deepcopy(frags)
    changed["echo"]["pdbcode"] = "other"
    assert all(entry["values"] == dict.fromkeys(KEYS)
               for entry in reconcile_energy_threshold_edits(edits, changed, SPECS)["ranges"])


@pytest.mark.parametrize("value", [None, "-16.60", "0.00", "+21", "0.1"])
def test_e3_t27_accepts_partial_current_values_without_scoring(value):
    frags = fragments()
    edits = reconcile_energy_threshold_edits(None, frags, SPECS)
    edits["ranges"][0]["values"][KEYS[0]] = value
    assert validate_energy_threshold_edits(edits, frags, SPECS) is None
    assert validate_energy_threshold_edits(None, frags, SPECS) is None


@pytest.mark.parametrize("value", ["1.234", "NaN", "Infinity", 1.0, True, "bad"])
def test_e3_t27_rejects_invalid_current_values(value):
    frags = fragments()
    edits = reconcile_energy_threshold_edits(None, frags, SPECS)
    edits["ranges"][0]["values"][KEYS[0]] = value
    with pytest.raises(ValueError):
        validate_energy_threshold_edits(edits, frags, SPECS)


def test_e3_t27_unknown_fragment_fails_but_stale_basis_and_ranges_defer():
    frags = fragments()
    edits = reconcile_energy_threshold_edits(None, frags, SPECS)
    edits["ranges"][0]["values"]["frag-999-GU"] = "1.00"
    with pytest.raises(ValueError, match="fragment"):
        validate_energy_threshold_edits(edits, frags, SPECS)
    stale = deepcopy(edits)
    stale["basis"]["pdbcode"] = "old"
    assert validate_energy_threshold_edits(stale, frags, SPECS) is None
    edits["ranges"][0]["row"] = {"chain": "gone", "start": 1, "end": 2}
    assert validate_energy_threshold_edits(edits, frags, SPECS) is None


def test_w8a_matrix_reuse_restores_uncovered_lock_despite_low_energy_and_preserves_caps():
    frags = fragments()
    frags["fragments"][2]["residues"][1]["resid"] = 107
    selected, stacks = selection_and_stackings((107, 2, 8), (111, 3, 9))
    clipped_rows = KEYS[3:]
    old_anchors = candidate_anchors(frags, selected, stacks, clipped_rows)
    new_anchors = candidate_anchors(frags, selected, stacks, KEYS)
    energy = energies()
    for column in energy["ranges"]:
        for row in column["rows"]:
            row.update(reference={"energy": 0, "error": None}, bestfit={"energy": 0, "error": None})
    old = make_interaction_matrix(frags, {"clip1": 3, "clip2": 0}, SPECS, selected, stacks)
    old = recompute_interaction_matrix(old, old_anchors, energy)
    old["grid"][1][1] = True  # a user cap beside the interior range-3 lock
    before = deepcopy(old)
    restored = adapt_interaction_matrix(old, KEYS, new_anchors, energy)
    assert restored["basis"] == old["basis"] == sorted(selected["selected"])
    assert restored["grid"][2][0] is True  # restored candidate lock ignores low energy
    assert restored["grid"][3:] == old["grid"]
    assert adapt_interaction_matrix(restored, clipped_rows, old_anchors, energy) == old
    validate_interaction_matrix(restored, KEYS, RANGES, new_anchors, energy)
    assert old == before
    handoff = interaction_handoff(old, frags, selected, stacks, {"clip1": 3, "clip2": 0}, SPECS)
    assert {(item["fragment"], item["nucleotide"]) for run in handoff["runs"] for item in run["anchors"]} == {
        (KEYS[3], 1), (KEYS[5], 1)}
    uncovered_handoff = interaction_handoff(restored, frags, selected, stacks, {"clip1": 0, "clip2": 0}, SPECS)
    assert {(item["fragment"], item["nucleotide"]) for run in uncovered_handoff["runs"] for item in run["anchors"]} == {
        (KEYS[2], 2), (KEYS[3], 1), (KEYS[5], 1)}
    energy["ranges"][0]["rows"][2]["bestfit"] = {"energy": None, "error": "bad candidate"}
    restored_with_error = adapt_interaction_matrix(old, KEYS, new_anchors, energy)
    assert restored_with_error["grid"][2][0] is True
    assert restored_with_error["grid"][3:] == old["grid"]
    with pytest.raises(ValueError, match="107"):
        validate_interaction_matrix(restored_with_error, KEYS, RANGES, new_anchors, energy)
    with pytest.raises(ValueError, match="107"):
        validate_interaction_matrix(restored, KEYS, RANGES, new_anchors, energy)
    # Explicit clipping removes only the failing candidate, retaining its stacking.
    validate_interaction_matrix(old, clipped_rows, RANGES, old_anchors, energy)
