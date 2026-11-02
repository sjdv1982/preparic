"""D7 / OC P2.1–P2.2: composite copies must match the Mol* picked frame.

Oracle: Mol* assembly.js expandOperators enumerates the leftmost group fastest;
getAssemblyOperators right-multiplies matrices. Thus (1,2)(3,4) gives
M1 M3, M2 M3, M1 M4, M2 M4 (rightmost operation acts first).
src/lib/chains.ts names these in that order, reserving bare A for identity.
Expected coordinates below are literal, independent of backend expansion/math.
The 4→5 UI gate belongs to W5 after W4 and is outside this W3 correction.
"""

from io import StringIO

import numpy as np
import pytest
from Bio.PDB import MMCIFIO, PDBParser

from range_pdbs import _ranges_to_pdb


def _composite_cif(translations, rotations=None, expression="(1,2)(3,4)"):
    """One RNA residue, three non-collinear atoms, four assembly operators."""
    data = {
        "data_": "SYNTH",
        "_struct_asym.id": ["A"],
        "_struct_asym.entity_id": ["1"],
        "_entity_poly.entity_id": ["1"],
        "_entity_poly.type": ["polyribonucleotide"],
        "_atom_site.group_PDB": ["ATOM"] * 3,
        "_atom_site.id": ["1", "2", "3"],
        "_atom_site.type_symbol": ["P", "C", "O"],
        "_atom_site.label_atom_id": ["P", "C4'", "O3'"],
        "_atom_site.label_alt_id": ["."] * 3,
        "_atom_site.label_comp_id": ["A"] * 3,
        "_atom_site.label_asym_id": ["A"] * 3,
        "_atom_site.label_entity_id": ["1"] * 3,
        "_atom_site.label_seq_id": ["7"] * 3,
        "_atom_site.auth_seq_id": ["7"] * 3,
        "_atom_site.auth_asym_id": ["A"] * 3,
        "_atom_site.pdbx_PDB_ins_code": ["?"] * 3,
        "_atom_site.Cartn_x": ["1", "2", "1"],
        "_atom_site.Cartn_y": ["2", "2", "3"],
        "_atom_site.Cartn_z": ["3"] * 3,
        "_atom_site.occupancy": ["1"] * 3,
        "_atom_site.B_iso_or_equiv": ["10"] * 3,
        "_atom_site.pdbx_PDB_model_num": ["1"] * 3,
        "_pdbx_struct_assembly.id": ["1"],
        "_pdbx_struct_assembly_gen.assembly_id": ["1"],
        "_pdbx_struct_assembly_gen.oper_expression": [expression],
        "_pdbx_struct_assembly_gen.asym_id_list": ["A"],
        "_pdbx_struct_oper_list.id": ["1", "2", "3", "4"],
    }
    if rotations is None:
        rotations = [np.eye(3)] * 4
    for i in range(3):
        data[f"_pdbx_struct_oper_list.vector[{i + 1}]"] = [str(t[i]) for t in translations]
        for j in range(3):
            data[f"_pdbx_struct_oper_list.matrix[{i + 1}][{j + 1}]"] = [str(r[i][j]) for r in rotations]
    writer = MMCIFIO()
    writer.set_dict(data)
    stream = StringIO()
    writer.save(stream)
    return stream.getvalue()


def _assert_copy(cif, chain, expected):
    pdbs = _ranges_to_pdb(cif, {
        "schema": 1, "pdbcode": "SYNTH", "biounit": 1,
        "rows": [{"chain": chain, "start": 7, "end": 7}],
    })
    assert len(pdbs) == 1
    atoms = list(PDBParser(QUIET=True).get_structure("copy", StringIO(pdbs[0])).get_atoms())
    assert [(a.name, a.parent.id[1]) for a in atoms] == [("P", 7), ("C4'", 7), ("O3'", 7)]
    assert np.asarray([a.coord for a in atoms]) == pytest.approx(np.asarray(expected), abs=0.001)


@pytest.mark.parametrize("expression", ["(1,2)(3,4)", "(1-2)(3-4)"])
@pytest.mark.parametrize("chain,offset", [
    ("A-1", (10, 30, 0)),  # (1,3)
    ("A-2", (20, 30, 0)),  # (2,3), not (1,4)
    ("A-3", (10, 40, 0)),  # (1,4), not (2,3)
    ("A-4", (20, 40, 0)),  # (2,4)
])
def test_composite_copy_names_follow_molstar_cartesian_order(expression, chain, offset):
    # Commuting translations isolate enumeration from composition order.
    cif = _composite_cif([(10, 0, 0), (20, 0, 0), (0, 30, 0), (0, 40, 0)], expression=expression)
    _assert_copy(cif, chain, np.array([(1, 2, 3), (2, 2, 3), (1, 3, 3)]) + offset)


@pytest.mark.parametrize("chain,expected", [
    ("A-1", [(5, 1, 3), (5, 2, 3), (4, 1, 3)]),
    ("A-2", [(25, -1, 3), (25, -2, 3), (26, -1, 3)]),
    ("A-3", [(4, 1, 3), (4, 2, 3), (3, 1, 3)]),
    ("A-4", [(26, -1, 3), (26, -2, 3), (27, -1, 3)]),
])
def test_composite_coordinates_apply_rightmost_operator_first(chain, expected):
    # Rotation and translation do not commute, so reversing multiplication fails
    # even for A-1/A-4, whose positions in the Cartesian enumeration do not change.
    cif = _composite_cif(
        [(10, 0, 0), (20, 0, 0), (0, 3, 0), (0, 4, 0)],
        [[[0, -1, 0], [1, 0, 0], [0, 0, 1]],
         [[0, 1, 0], [-1, 0, 0], [0, 0, 1]], np.eye(3), np.eye(3)],
    )
    _assert_copy(cif, chain, expected)


@pytest.mark.parametrize("chain,offset", [
    ("A-1", (10, -20, 0)),  # (1,3)
    ("A", (0, 0, 0)),      # (2,3): identity is second, not first
    ("A-2", (10, 0, 30)),  # (1,4): suffix skips the identity
    ("A-3", (0, 20, 30)),  # (2,4)
])
def test_composite_identity_copy_reserves_bare_chain_name(chain, offset):
    cif = _composite_cif([(10, 0, 0), (0, 20, 0), (0, -20, 0), (0, 0, 30)])
    _assert_copy(cif, chain, np.array([(1, 2, 3), (2, 2, 3), (1, 3, 3)]) + offset)
