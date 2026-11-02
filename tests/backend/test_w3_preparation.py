"""W3 range fidelity oracles (OC G8, P2.1, D7/D15), never RCSB-dependent."""

from io import StringIO
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
from Bio.PDB import MMCIFIO, PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from range_pdbs import _ranges_to_pdb
from alaric_dir import resolve_alaric_dir, resolve_alaric_repo_root

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _cif(data):
    writer = MMCIFIO()
    writer.set_dict(data)
    stream = StringIO()
    writer.save(stream)
    return stream.getvalue()


def _specs(*rows, biounit=1):
    return {"schema": 1, "pdbcode": "1B7F", "biounit": biounit, "rows": list(rows)}


def _atoms(text):
    return list(PDBParser(QUIET=True).get_structure("range", StringIO(text)).get_atoms())


def _assembly_cif(cif_text):
    data = MMCIF2Dict(StringIO(cif_text))
    rotations = [np.eye(3), np.eye(3), np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])]
    translations = [np.zeros(3), np.array([30, 0, 0]), np.array([11, -7, 3])]
    data["_pdbx_struct_oper_list.id"] = ["1", "2", "3"]
    data["_pdbx_struct_oper_list.type"] = ["identity operation", "point symmetry operation", "point symmetry operation"]
    data["_pdbx_struct_oper_list.name"] = ["1_555", "copy2", "copy3"]
    data["_pdbx_struct_oper_list.symmetry_operation"] = ["x,y,z", "x+30,y,z", "-y+11,x-7,z+3"]
    for i in range(3):
        data[f"_pdbx_struct_oper_list.vector[{i + 1}]"] = [str(t[i]) for t in translations]
        for j in range(3):
            data[f"_pdbx_struct_oper_list.matrix[{i + 1}][{j + 1}]"] = [str(r[i, j]) for r in rotations]
    data["_pdbx_struct_assembly_gen.assembly_id"] = ["1", "2"]
    data["_pdbx_struct_assembly_gen.oper_expression"] = ["(1,2,3)", "1"]
    data["_pdbx_struct_assembly_gen.asym_id_list"] = ["A,C", "A"]
    return _cif(data), rotations[2], translations[2]


def _remove_residue(cif_text, chain, residue):
    data = MMCIF2Dict(StringIO(cif_text))
    keep = [not (asym == chain and seq == str(residue)) for asym, seq in zip(
        data["_atom_site.label_asym_id"], data["_atom_site.label_seq_id"], strict=True
    )]
    for key in list(data):
        if key.startswith("_atom_site."):
            data[key] = [value for value, retain in zip(data[key], keep, strict=True) if retain]
    return _cif(data)


def test_assembly_copy_is_the_browser_named_copy_in_its_transformed_frame(cif_text):
    synthetic, rotation, translation = _assembly_cif(cif_text)
    original = _ranges_to_pdb(synthetic, _specs({"chain": "A", "start": 1, "end": 4}))[0]
    copied = _ranges_to_pdb(synthetic, _specs({"chain": "A-2", "start": 1, "end": 4}))[0]
    original_atoms, copied_atoms = _atoms(original), _atoms(copied)
    assert [atom.full_id[3:] for atom in copied_atoms] == [atom.full_id[3:] for atom in original_atoms]
    assert np.asarray([atom.coord for atom in copied_atoms]) == pytest.approx(
        np.asarray([atom.coord for atom in original_atoms]) @ rotation.T + translation, abs=0.0015
    )


def test_row_validity_uses_selected_assembly_not_asymmetric_unit(cif_text):
    synthetic, _, _ = _assembly_cif(cif_text)
    with pytest.raises(ValueError, match=r"(?i)(assembly|biounit|not present|chain)"):
        _ranges_to_pdb(synthetic, _specs(
            {"chain": "A", "start": 1, "end": 4},
            {"chain": "C", "start": 2, "end": 8}, biounit=2,
        ))
    with pytest.raises(ValueError, match=r"(?i)(assembly|biounit|not present|chain)"):
        _ranges_to_pdb(synthetic, _specs({"chain": "A-2", "start": 1, "end": 4}, biounit=2))


def test_missing_internal_rna_residue_is_named_and_refused(cif_text):
    missing = _remove_residue(cif_text, "A", 2)
    with pytest.raises(ValueError, match=r"(?i)(unobserved|missing).*2|2.*(unobserved|missing)"):
        _ranges_to_pdb(missing, _specs({"chain": "A", "start": 1, "end": 4}))


def test_missing_protein_residue_does_not_apply_the_rna_gap_rule(cif_text):
    missing = _remove_residue(cif_text, "C", 4)
    pdbs = _ranges_to_pdb(missing, _specs(
        {"chain": "A", "start": 1, "end": 4}, {"chain": "C", "start": 2, "end": 8},
    ))
    assert {atom.parent.id[1] for atom in _atoms(pdbs[1])} == {2, 3, 5, 6, 7, 8}


def test_consecutive_labels_do_not_hide_a_chemical_rna_chain_break(cif_text):
    data = MMCIF2Dict(StringIO(cif_text))
    for i, (chain, seq) in enumerate(zip(data["_atom_site.label_asym_id"], data["_atom_site.label_seq_id"], strict=True)):
        if chain == "A" and seq == "3":
            data["_atom_site.Cartn_x"][i] = str(float(data["_atom_site.Cartn_x"][i]) + 100)
    with pytest.raises(ValueError, match=r"(?i)(break|continu|connect|bond)"):
        _ranges_to_pdb(_cif(data), _specs({"chain": "A", "start": 1, "end": 4}))


@pytest.mark.parametrize("schema", [None, 2, True])
def test_canonical_range_envelope_rejects_unknown_or_untyped_schema(cif_text, schema):
    specs = _specs({"chain": "A", "start": 1, "end": 4})
    if schema is None:
        del specs["schema"]
    else:
        specs["schema"] = schema
    with pytest.raises(ValueError, match=r"(?i)schema"):
        _ranges_to_pdb(cif_text, specs)


def _ring_pdb(reference, coordinates, resname, chain, resid):
    lines = []
    for serial, (atom, coord) in enumerate(zip(_atoms(reference), coordinates, strict=True), 1):
        x, y, z = coord
        lines.append(
            f"ATOM  {serial:5d} {atom.name:>4s} {resname:>3s} {chain}{resid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{20.0:6.2f}          {atom.element:>2s}\n"
        )
    return "".join(lines) + "END\n"


def _fitted_reference_normal(reference, observed):
    reference = reference - reference.mean(axis=0)
    observed = observed - observed.mean(axis=0)
    left, _, right = np.linalg.svd(reference.T @ observed)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(left @ right)
    row_rotation = left @ correction @ right
    normal = np.linalg.svd(reference, full_matrices=False)[2][-1] @ row_rotation
    return normal / np.linalg.norm(normal)


def test_stage3_descriptors_use_pinned_kabsch_reference_normals():
    from detect_stacking import _detect_stacking

    alaric = resolve_alaric_repo_root()
    references = {name: (alaric / "ring-refe" / f"refe-{name}.pdb").read_bytes()
                  for name in ("A", "C", "PHE", "HIS", "ARG", "TRP")}
    a_text, arg_text = references["A"].decode(), references["ARG"].decode()
    a_reference = np.asarray([atom.coord for atom in _atoms(a_text)], dtype=float)
    arg_reference = np.asarray([atom.coord for atom in _atoms(arg_text)], dtype=float)
    a_observed = a_reference - a_reference.mean(axis=0) + [2, 0, 4]
    arg_observed = arg_reference - arg_reference.mean(axis=0)
    arg_observed[:, 0] *= 1.4
    arg_observed[:, 1] += 0.2 * arg_observed[:, 0]
    arg_observed[0, 2] += 0.9
    rna = _ring_pdb(a_text, a_observed, "A", "A", 101)
    protein = _ring_pdb(arg_text, arg_observed, "ARG", "C", 30)
    # Derive the oracle from the delivered PDB precision, rather than a looser
    # comparison against the pre-format coordinates.
    a_observed = np.asarray([atom.coord for atom in _atoms(rna)], dtype=float)
    arg_observed = np.asarray([atom.coord for atom in _atoms(protein)], dtype=float)
    rna_normal = _fitted_reference_normal(a_reference, a_observed)
    protein_normal = _fitted_reference_normal(arg_reference, arg_observed)
    displacement = a_observed.mean(axis=0) - arg_observed.mean(axis=0)
    expected_phi = math.degrees(math.acos(np.clip(abs(displacement @ protein_normal) / np.linalg.norm(displacement), 0, 1)))
    expected_tilt = math.degrees(math.acos(np.clip(abs(rna_normal @ protein_normal), 0, 1)))
    rows = _detect_stacking([rna, protein], ring_references=references)
    assert len(rows) == 1
    assert rows[0]["slip_angle"] == pytest.approx(expected_phi, abs=1e-5)
    assert rows[0]["tilt_angle"] == pytest.approx(expected_tilt, abs=1e-5)


def test_guanine_uses_the_pinned_adenine_ring_reference():
    from detect_stacking import _detect_stacking

    alaric = resolve_alaric_repo_root()
    references = {name: (alaric / "ring-refe" / f"refe-{name}.pdb").read_bytes()
                  for name in ("A", "C", "PHE", "HIS", "ARG", "TRP")}
    # The contract pins six templates and explicitly maps G to A. No separate
    # G file is needed because both rings use the same nine atom names.
    assert "G" not in references
    a_text, arg_text = references["A"].decode(), references["ARG"].decode()
    a_reference = np.asarray([atom.coord for atom in _atoms(a_text)], dtype=float)
    arg_reference = np.asarray([atom.coord for atom in _atoms(arg_text)], dtype=float)
    g_observed = a_reference - a_reference.mean(axis=0) + [2, 0, 4]
    g_observed[:, 0] *= 1.4
    g_observed[:, 1] += 0.2 * g_observed[:, 0]
    g_observed[0, 2] += 0.9
    arg_observed = arg_reference - arg_reference.mean(axis=0)
    arg_observed[:, 0] *= 1.4
    arg_observed[:, 1] += 0.2 * arg_observed[:, 0]
    arg_observed[0, 2] += 0.9
    rna = _ring_pdb(a_text, g_observed, "G", "A", 101)
    protein = _ring_pdb(arg_text, arg_observed, "ARG", "C", 30)
    g_observed = np.asarray([atom.coord for atom in _atoms(rna)], dtype=float)
    arg_observed = np.asarray([atom.coord for atom in _atoms(protein)], dtype=float)
    g_normal = _fitted_reference_normal(a_reference, g_observed)
    protein_normal = _fitted_reference_normal(arg_reference, arg_observed)
    displacement = g_observed.mean(axis=0) - arg_observed.mean(axis=0)
    expected_phi = math.degrees(math.acos(np.clip(
        abs(displacement @ protein_normal) / np.linalg.norm(displacement), 0, 1
    )))
    expected_tilt = math.degrees(math.acos(np.clip(abs(g_normal @ protein_normal), 0, 1)))

    rows = _detect_stacking([rna, protein], ring_references=references)
    assert len(rows) == 1
    assert rows[0]["slip_angle"] == pytest.approx(expected_phi, abs=1e-5)
    assert rows[0]["tilt_angle"] == pytest.approx(expected_tilt, abs=1e-5)


def test_one_alaric_setting_resolves_override_env_symlink_and_repo_root(tmp_path, monkeypatch):
    first, second = tmp_path / "first" / "alaric", tmp_path / "second" / "alaric"
    first.mkdir(parents=True)
    second.mkdir(parents=True)
    link = tmp_path / "chosen"
    link.symlink_to(first, target_is_directory=True)
    monkeypatch.setenv("ALARIC_DIR", str(link))
    assert resolve_alaric_dir() == first.resolve()
    assert resolve_alaric_repo_root() == first.parent.resolve()
    assert resolve_alaric_dir(second) == second.resolve()
    assert resolve_alaric_repo_root(second) == second.parent.resolve()
    with pytest.raises(FileNotFoundError, match="ALARIC_DIR"):
        resolve_alaric_dir(tmp_path / "missing")


def test_backend_import_does_not_import_installed_alaric_code():
    program = """
import importlib.abc
import sys
class NoInstalledAlaric(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'alaric' or fullname.startswith('alaric.'):
            raise AssertionError('backend imported installed Alaric: ' + fullname)
sys.meta_path.insert(0, NoInstalledAlaric())
import backend
"""
    result = subprocess.run([sys.executable, "-c", program], cwd=PROJECT_ROOT, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr


def test_real_preparation_preserves_rna_ids_and_is_accepted_without_lenient_flags(tmp_path, expected_pdbs):
    pdbs = tmp_path / "pdbs"
    pdbs.mkdir()
    source = "".join(
        line[:22] + f"{int(line[22:26]) + 100:4d}" + line[26:] if line.startswith("ATOM") else line
        for line in expected_pdbs[0].splitlines(keepends=True)
    )
    (pdbs / "range-1.pdb").write_text(source)
    (tmp_path / "aareduce").symlink_to(PROJECT_ROOT / "aareduce", target_is_directory=True)
    prepared = subprocess.run(["bash", str(PROJECT_ROOT / "range-pdbs-aareduce.sh")],
                              cwd=tmp_path, capture_output=True, text=True)
    assert prepared.returncode == 0, prepared.stdout + prepared.stderr
    rna_path = tmp_path / "RESULT" / "range-1-aa.pdb"
    atoms = _atoms(rna_path.read_text())
    assert sorted({atom.parent.id[1] for atom in atoms}) == [101, 102, 103, 104]
    assert all(atom.element not in {"H", "D"} for atom in atoms)
    repo = resolve_alaric_repo_root()
    strict = subprocess.run([
        sys.executable, str(repo / "util" / "refe-fragments.py"), str(rna_path),
        "--templates", str(repo / "fraglib" / "templates"), "--output", str(tmp_path / "fragments"),
    ], cwd=tmp_path, capture_output=True, text=True)
    assert strict.returncode == 0, strict.stdout + strict.stderr
    assert len(list((tmp_path / "fragments").glob("frag-*-refe.pdb"))) == 3
