"""W3 early-fragment, editable-threshold, and byte-exact handoff oracles."""

import copy
import json

import pytest

from stage4 import (
    _candidate_fragments, _clamp_clip, _clip_fragments, _list_fragments,
    _selected_candidate_numbers, _selected_candidate_groups, _share_range_pdbs_aa, _valid_clip,
    _validate_threshold_edits,
)


RNA = {"chain": "A", "start": 101, "end": 104}
BASIS = {"pdbcode": "1B7F", "rna": RNA}
SEQUENCES = ["AC", "CG", "GU"]


def _pdb(residues, chain="A"):
    return "".join(
        f"ATOM  {i:5d}  P   {base:>3s} {chain}{resid:4d}    "
        f"{float(i):8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{20.0:6.2f}           P\n"
        for i, (resid, base) in enumerate(residues, 1)
    ) + "END\n"


def _identity(resid, base):
    return {"chain": "A", "resid": resid, "insertion_code": "", "base": base}


def fragment_map():
    residues = [_identity(101 + i, base) for i, base in enumerate("ACGU")]
    return {
        "schema": 1, "echo": copy.deepcopy(BASIS), "residues": residues,
        "fragments": [{"number": i + 1, "key": f"frag-{i + 1}-{seq}", "sequence": seq,
                       "residues": residues[i:i + 2]} for i, seq in enumerate(SEQUENCES)],
        "adjacent_pairs": [{"down": 1, "up": 2, "sequence": "ACG"}, {"down": 2, "up": 3, "sequence": "CGU"}],
    }


def edits(blank=False):
    value = None if blank else "0.25"
    return {
        "schema": 1, "basis": copy.deepcopy(BASIS), "active": True,
        "pairs": [{**row, "ovrmsd": value, "crmsd": value} for row in fragment_map()["adjacent_pairs"]],
    }


def test_t24_uses_t6_keys_but_original_rna_residue_identities_without_a_fit():
    refe = {f"frag-{i + 1}-{seq}-refe.pdb": _pdb([(1, seq[0]), (2, seq[1])], "B").encode()
            for i, seq in reversed(list(enumerate(SEQUENCES)))}
    stacking = {"rna_resname": "C", "rna_resid": 102, "protein_range_pdb_id": 2, "protein_resid": 30}
    result = _list_fragments(refe, _pdb(list(zip(range(101, 105), "ACGU"))), "1B7F", RNA)
    expected = fragment_map()
    assert result == expected
    assert _candidate_fragments(result, [stacking]) == [
        {"key": '["C",102,2,30]', "fragment": "frag-1-AC", "nucleotide": 2},
        {"key": '["C",102,2,30]', "fragment": "frag-2-CG", "nucleotide": 1},
    ]


def test_candidate_membership_tracks_fresh_protein_stacking_inputs_not_t24():
    fragments = fragment_map()
    selection = {"schema": 1, "basis": BASIS, "active": True,
                 "selected": {'["C",102,2,30]': True}}
    original = [{"rna_resname": "C", "rna_resid": 102,
                 "protein_range_pdb_id": 2, "protein_resid": 30}]
    assert _candidate_fragments(fragments, original) == [
        {"key": '["C",102,2,30]', "fragment": "frag-1-AC", "nucleotide": 2},
        {"key": '["C",102,2,30]', "fragment": "frag-2-CG", "nucleotide": 1},
    ]
    assert _selected_candidate_numbers(fragments, selection, original) == [1, 2]
    edited_protein_rows = [{**original[0], "protein_resid": 31}]
    with pytest.raises(ValueError, match="outside the fresh stacking table"):
        _selected_candidate_numbers(fragments, selection, edited_protein_rows)


def test_t25_keeps_all_atom_pdb_bytes_and_range_order_with_explicit_echo():
    # Identical internal chain IDs cannot be used as range identities.
    rna = b"ATOM      1  P     A A 101      1.000   2.000   3.000\r\nEND  \r\n"
    protein = b"ATOM      1  CA  PHE A  30      7.000   8.000   9.000\nEND\n"
    specs = {"schema": 1, "pdbcode": "1B7F", "biounit": 2, "rows": [RNA, {"chain": "C-2", "start": 30, "end": 40}]}
    result = _share_range_pdbs_aa({"range-2-aa.pdb": protein, "range-1-aa.pdb": rna}, specs, "1b7f")
    assert result == {
        "schema": 1, "pdb_code": "1b7f", "range_specs": specs,
        "pdbs": [{"range": "range-1", "name": "range-1-aa.pdb", "pdb": rna.decode()},
                 {"range": "range-2", "name": "range-2-aa.pdb", "pdb": protein.decode()}],
    }
    assert result["pdbs"][0]["pdb"].encode() == rna


def test_t10_accepts_undecided_and_blank_skeleton_before_fit_completes():
    assert _validate_threshold_edits(None, fragment_map()) is None
    assert _validate_threshold_edits(edits(blank=True), fragment_map()) is None


@pytest.mark.parametrize("value", ["0", "-0.1", "0.251", "NaN", "Infinity", "1e-2", 0.25, True])
def test_t10_rejects_invalid_value_without_rounding(value):
    marker = edits()
    marker["pairs"][0]["ovrmsd"] = value
    with pytest.raises(ValueError):
        _validate_threshold_edits(marker, fragment_map())


@pytest.mark.parametrize("schema", [None, 2, True])
def test_t10_rejects_boundary_schema_without_partial_use(schema):
    marker = edits()
    if schema is None:
        del marker["schema"]
    else:
        marker["schema"] = schema
    with pytest.raises(ValueError, match=r"(?i)schema"):
        _validate_threshold_edits(marker, fragment_map())


def test_t10_refuses_missing_duplicate_or_extra_pair_rows():
    for mutation in (lambda m: m["pairs"].pop(),
                     lambda m: m["pairs"].append(copy.deepcopy(m["pairs"][0])),
                     lambda m: m["pairs"].append({"down": 3, "up": 4, "sequence": "GUA", "ovrmsd": "0.2", "crmsd": "0.2"})):
        marker = edits()
        mutation(marker)
        with pytest.raises(ValueError):
            _validate_threshold_edits(marker, fragment_map())


def test_clip_validation_and_entering_stage4_clamp_keep_one_candidate_per_anchor():
    clip = {"schema": 1, "basis": copy.deepcopy(BASIS), "clip1": 2, "clip2": 0}
    assert not _valid_clip(clip, 3, [[1, 2]])
    clamped = _clamp_clip(clip, 3, [[1, 2]])
    assert clamped == {**clip, "clip1": 1}
    assert clip["clip1"] == 2
    assert _valid_clip(clamped, 3, [[1, 2]])
    assert _clamp_clip(clamped, 3, [[1, 2]]) == clamped
    assert _clamp_clip({**clip, "clip1": 0, "clip2": 2}, 3, [[2, 3]]) == {**clip, "clip1": 0, "clip2": 1}
    assert not _valid_clip({**clamped, "clip1": True}, 3, [[1, 2]])
    assert not _valid_clip({**clamped, "clip2": -1}, 3, [[1, 2]])
    assert not _valid_clip({**clamped, "clip1": 0.5}, 3, [[1, 2]])
    assert not _valid_clip({**clamped, "clip1": 2, "clip2": 1}, 3, [[1, 2]])
    assert _clamp_clip({**clip, "clip1": 1, "clip2": 1}, 3, [[1], [3]]) == {**clip, "clip1": 0, "clip2": 0}


def test_selected_candidate_groups_preserve_stacking_identity():
    stackings = [{"rna_resname": "C", "rna_resid": 102, "protein_range_pdb_id": 2, "protein_resid": 30},
                 {"rna_resname": "G", "rna_resid": 103, "protein_range_pdb_id": 2, "protein_resid": 31}]
    selection = {"schema": 1, "basis": BASIS, "active": True, "selected": {json.dumps([row["rna_resname"], row["rna_resid"], 2, row["protein_resid"]], separators=(",", ":")): True for row in stackings}}
    assert _selected_candidate_groups(fragment_map(), selection, stackings) == [[1, 2], [2, 3]]


def test_t11_filters_edited_rows_and_keeps_pdb_bytes_and_actual_residues():
    fragments = fragment_map()
    marker = edits()
    marker["pairs"][0]["crmsd"] = None
    reference = {f"{f['key']}-refe.pdb": f"REFERENCE {f['key']}\r\n".encode() for f in fragments["fragments"]}
    fit_echo = {"pdb_code": "1b7f", "rna": RNA, "resids": fragments["residues"]}
    best_fits = {"schema": 1, "echo": fit_echo,
                 "rows": [{"key": f["key"], "fragment": f["number"], "sequence": f["sequence"],
                           "conformer": 1, "rotamer": 1, "grid": [0, 0, 0], "pdb": f"BEST {f['key']}\n"}
                          for f in fragments["fragments"]]}
    clip = {"schema": 1, "basis": copy.deepcopy(BASIS), "clip1": 1, "clip2": 0}
    before = copy.deepcopy((best_fits, reference, marker, clip, fragments))
    folder = _clip_fragments(best_fits, reference, marker, clip, fragments)
    assert set(folder) == {"reference/frag-2-CG-refe.pdb", "reference/frag-3-GU-refe.pdb",
                           "bestfit/frag-2-CG-bestfit.pdb", "bestfit/frag-3-GU-bestfit.pdb",
                           "thresholds.json", "clip.json"}
    assert folder["reference/frag-2-CG-refe.pdb"] == reference["frag-2-CG-refe.pdb"]
    assert folder["bestfit/frag-2-CG-bestfit.pdb"] == b"BEST frag-2-CG\n"
    assert json.loads(folder["thresholds.json"]) == {"schema": 1, "pairs": marker["pairs"][1:]}
    assert json.loads(folder["clip.json"]) == {"schema": 1, "clip1": 1, "clip2": 0,
                                               "kept_residues": fragments["residues"][1:]}
    assert (best_fits, reference, marker, clip, fragments) == before


def test_t10_refuses_legacy_bases_instead_of_migrating():
    marker = {**edits(), "bases": []}
    with pytest.raises(ValueError):
        _validate_threshold_edits(marker, fragment_map())
