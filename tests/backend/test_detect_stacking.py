"""Plain-Python tests for the deterministic stacking detector."""

import math

import numpy as np
import pytest

from conftest import expected_pdbs
from detect_stacking import _detect_stacking


_RNA_ATOMS = ("N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9")
_PHE_ATOMS = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")


def _ring_coordinates(atom_names, center=(0.0, 0.0, 0.0), tilt_degrees=0.0):
    tilt = math.radians(tilt_degrees)
    cos_tilt, sin_tilt = math.cos(tilt), math.sin(tilt)
    center = np.asarray(center, dtype=float)
    coordinates = {}
    for index, atom_name in enumerate(atom_names):
        angle = 2 * math.pi * index / len(atom_names)
        x, y = math.cos(angle), math.sin(angle)
        coordinates[atom_name] = center + np.array(
            [x, y * cos_tilt, y * sin_tilt]
        )
    return coordinates


def _pdb(*residues):
    lines = []
    serial = 1
    for resname, resid, coordinates, chain in residues:
        for atom_name, (x, y, z) in coordinates.items():
            element = next((char for char in atom_name if char.isalpha()), "C")
            lines.append(
                f"ATOM  {serial:5d} {atom_name:>4s} {resname:>3s} {chain}{resid:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{20.0:6.2f}          {element:>2s}\n"
            )
            serial += 1
    lines.append("END\n")
    return "".join(lines)


def _pair_pdbs(distance=5.0, protein_resid=93, rna_resid=4, tilt=0.0, slip=0.0,
               rna_resname="A", protein_resname="PHE"):
    # slip is the angle between the centroid displacement and the protein ring
    # normal (z): 0 deg = zenith, 90 deg = in the ring plane.
    slip_radians = math.radians(slip)
    rna_center = (
        distance * math.sin(slip_radians),
        0.0,
        distance * math.cos(slip_radians),
    )
    return [
        _pdb((rna_resname, rna_resid, _ring_coordinates(_RNA_ATOMS, rna_center, tilt), "A")),
        _pdb((protein_resname, protein_resid, _ring_coordinates(_PHE_ATOMS), "C")),
    ]


def _plain_types(value):
    if isinstance(value, dict):
        return all(_plain_types(key) and _plain_types(item) for key, item in value.items())
    if isinstance(value, list):
        return all(_plain_types(item) for item in value)
    return value is None or type(value) in {str, int, float, bool}


def test_positive_pair_reports_independently_known_geometry():
    range_pdbs = _pair_pdbs(distance=4.8, tilt=20.0, slip=30.0)

    result = _detect_stacking(range_pdbs)

    assert len(result) == 1
    row = result[0]
    assert row["rna_resname"] == "A"
    assert row["rna_resid"] == 4
    assert row["protein_range_pdb_id"] == 2
    assert row["protein_resid"] == 93
    assert row["distance"] == pytest.approx(4.8, abs=0.002)
    assert row["tilt_angle"] == pytest.approx(20.0, abs=0.1)
    assert row["slip_angle"] == pytest.approx(30.0, abs=0.1)
    assert row["rna_centroid"] == pytest.approx([2.4, 0.0, 4.8 * math.cos(math.radians(30))], abs=0.002)
    assert row["protein_centroid"] == pytest.approx([0.0, 0.0, 0.0], abs=0.002)
    assert _plain_types(result)


def test_slip_angle_is_measured_from_the_protein_ring_normal():
    zenith = _detect_stacking(_pair_pdbs(distance=4.0, slip=0.0))
    horizon = _detect_stacking(_pair_pdbs(distance=4.0, slip=90.0))

    assert zenith[0]["slip_angle"] == pytest.approx(0.0, abs=0.1)
    assert horizon[0]["slip_angle"] == pytest.approx(90.0, abs=0.1)


def test_distance_cutoff_is_inclusive_and_just_outside_is_excluded():
    at_cutoff = _detect_stacking(_pair_pdbs(distance=6.0))
    outside = _detect_stacking(_pair_pdbs(distance=6.01))

    assert len(at_cutoff) == 1
    assert at_cutoff[0]["distance"] == pytest.approx(6.0, abs=0.002)
    assert outside == []


def test_tilt_cutoff_is_inclusive_and_just_outside_is_excluded():
    at_cutoff = _detect_stacking(_pair_pdbs(distance=5.0, tilt=45.0))
    outside = _detect_stacking(_pair_pdbs(distance=5.0, tilt=45.1))

    assert len(at_cutoff) == 1
    assert at_cutoff[0]["tilt_angle"] == pytest.approx(45.0, abs=0.1)
    assert outside == []


def test_rna_residue_with_missing_ring_atom_is_ignored():
    incomplete_rna = dict(_ring_coordinates(_RNA_ATOMS))
    incomplete_rna.pop("N9")
    range_pdbs = [
        _pdb(("A", 4, incomplete_rna, "A")),
        _pdb(("PHE", 93, _ring_coordinates(_PHE_ATOMS), "C")),
    ]

    assert _detect_stacking(range_pdbs) == []


def test_protein_residue_with_missing_ring_atom_is_ignored():
    incomplete_protein = dict(_ring_coordinates(_PHE_ATOMS))
    incomplete_protein.pop("CZ")
    range_pdbs = [
        _pdb(("A", 4, _ring_coordinates(_RNA_ATOMS, (0.0, 0.0, 5.0)), "A")),
        _pdb(("PHE", 93, incomplete_protein, "C")),
    ]

    assert _detect_stacking(range_pdbs) == []


def test_non_ring_residues_are_ignored():
    alanine = {"N": (-1.0, 0.0, 0.0), "CA": (0.0, 0.0, 0.0), "C": (1.0, 0.0, 0.0)}
    range_pdbs = [
        _pdb(("A", 4, _ring_coordinates(_RNA_ATOMS, (0.0, 0.0, 5.0)), "A")),
        _pdb(("ALA", 93, alanine, "C")),
    ]

    assert _detect_stacking(range_pdbs) == []


@pytest.mark.parametrize("alias,canonical", [("RA", "A"), ("RC", "C"), ("RG", "G"), ("RU", "U")])
def test_rna_aliases_are_canonicalized(alias, canonical):
    result = _detect_stacking(_pair_pdbs(rna_resname=alias))

    assert len(result) == 1
    assert result[0]["rna_resname"] == canonical


def test_multiple_protein_files_keep_one_based_original_range_ids():
    rna_pdb = _pair_pdbs()[0]
    protein_pdb = _pair_pdbs()[1]

    result = _detect_stacking([rna_pdb, protein_pdb, protein_pdb])

    assert [row["protein_range_pdb_id"] for row in result] == [2, 3]


def test_results_are_sorted_by_protein_range_then_residue_ids():
    rna = _pdb(
        ("A", 10, _ring_coordinates(_RNA_ATOMS, (0.0, 0.0, 5.0)), "A"),
        ("G", 2, _ring_coordinates(_RNA_ATOMS, (0.0, 0.0, 5.0)), "A"),
    )
    protein = _pdb(
        ("PHE", 9, _ring_coordinates(_PHE_ATOMS), "C"),
        ("TYR", 3, _ring_coordinates(_PHE_ATOMS), "C"),
    )

    result = _detect_stacking([rna, protein])

    assert [(row["protein_resid"], row["rna_resid"]) for row in result] == [
        (3, 2), (3, 10), (9, 2), (9, 10),
    ]
    assert result == _detect_stacking([rna, protein])


def test_only_the_first_model_is_used():
    def with_models(*model_pdbs):
        parts = []
        for index, model_pdb in enumerate(model_pdbs, start=1):
            parts.extend((f"MODEL     {index:4d}\n", model_pdb.removesuffix("END\n"), "ENDMDL\n"))
        return "".join(parts) + "END\n"

    empty_rna_model = _pdb(("ALA", 4, {"CA": (0.0, 0.0, 0.0)}, "A"))
    ring_rna_model = _pdb(("A", 4, _ring_coordinates(_RNA_ATOMS, (0.0, 0.0, 5.0)), "A"))
    protein_model = _pdb(("PHE", 93, _ring_coordinates(_PHE_ATOMS), "C"))

    assert _detect_stacking([
        with_models(empty_rna_model, ring_rna_model),
        with_models(protein_model, protein_model),
    ]) == []


def test_empty_inputs_and_too_few_ranges_return_empty():
    assert _detect_stacking([]) == []
    assert _detect_stacking(["valid PDB text\n"]) == []
    assert _detect_stacking(["", "not parsed"]) == []
    assert _detect_stacking([" \n", " "]) == []


@pytest.mark.parametrize("value", [None, (), "PDB", {"pdb": "text"}, ["pdb", 1]])
def test_rejects_invalid_boundary_types(value):
    with pytest.raises(ValueError, match="range_pdbs|PDB string"):
        _detect_stacking(value)


@pytest.mark.parametrize("malformed", ["not a PDB", "END\n", "ATOM      1  C   ALA A   1      nope\n"])
def test_rejects_malformed_nonempty_pdbs(malformed):
    with pytest.raises(ValueError, match=r"range_pdbs\[0\]|range_pdbs\[1\]"):
        _detect_stacking([malformed, _pdb(("PHE", 1, _ring_coordinates(_PHE_ATOMS), "C"))])


def test_repository_1b7f_fixture_detects_rna_range_against_protein_ranges(expected_pdbs):
    result = _detect_stacking(expected_pdbs)

    assert [
        (row["rna_resname"], row["rna_resid"], row["protein_range_pdb_id"], row["protein_resid"])
        for row in result
    ] == [("G", 4, 3, 93)]
    assert result[0]["distance"] == pytest.approx(3.699, abs=0.001)
    assert result[0]["tilt_angle"] == pytest.approx(15.719, abs=0.001)
    assert result[0]["slip_angle"] == pytest.approx(21.949, abs=0.001)
    assert result == sorted(
        result,
        key=lambda row: (row["protein_range_pdb_id"], row["protein_resid"], row["rna_resid"]),
    )
