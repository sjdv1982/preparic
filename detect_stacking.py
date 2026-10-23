"""Detect RNA--protein aromatic stacking pairs in range PDB structures."""

from __future__ import annotations

def _detect_stacking(range_pdbs: list[str]) -> list[dict]:
    """Return deterministic RNA--protein ring pairs from range PDB text.

    ``range_pdbs[0]`` is the RNA structure and subsequent entries are protein
    structures. The result contains only ordinary Python values suitable for a
    Seamless ``plain`` cell.
    """
    import math
    from io import StringIO

    import numpy as np
    from Bio.PDB import PDBParser

    # Keep all implementation dependencies and configuration within this
    # function. Seamless can then bind its source as a self-contained,
    # content-addressed transformation without implicit module globals.
    rings = {
        "T": ("C6", "C5", "C7", "N3", "O2", "O4"),
        "U": ("N1", "C2", "N3", "C4", "C5", "C6"),
        "C": ("N1", "C2", "N3", "C4", "C5", "C6"),
        "G": ("N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"),
        "A": ("N1", "C2", "N3", "C4", "C5", "C6", "N7", "C8", "N9"),
        "PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
        "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
        "HIS": ("CG", "ND1", "CD2", "CE1", "NE2"),
        "ARG": ("CD", "NE", "CZ", "NH1", "NH2"),
        "TRP": ("CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
    }
    rna_rings = frozenset(("A", "C", "G", "U", "T"))
    protein_rings = frozenset(("PHE", "TYR", "HIS", "ARG", "TRP"))
    max_centroid_distance = 6.0
    max_plane_tilt = 45.0
    distance_boundary_tolerance = 1e-12
    angle_boundary_tolerance = 1e-5

    def canonical_resname(resname: str) -> str:
        normalized = resname.strip().upper()
        if (
            len(normalized) == 2
            and normalized[0] in {"R", "D"}
            and normalized[1] in rna_rings
        ):
            return normalized[1]
        return normalized

    def fit_plane_normal(coordinates: np.ndarray) -> np.ndarray:
        centered = coordinates - coordinates.mean(axis=0)
        _, singular_values, right_vectors = np.linalg.svd(
            centered, full_matrices=False
        )
        if (
            len(singular_values) < 2
            or singular_values[0] == 0
            or singular_values[1]
            <= np.finfo(float).eps * max(centered.shape) * singular_values[0]
        ):
            raise ValueError("ring atoms do not define a plane")
        normal = right_vectors[-1]
        norm = float(np.linalg.norm(normal))
        if not math.isfinite(norm) or norm == 0:
            raise ValueError("ring atoms do not define a finite plane")
        return normal / norm

    def complete_rings(model, allowed_resnames: frozenset[str]):
        rings_found = []
        for residue in model.get_residues():
            resname = canonical_resname(residue.get_resname())
            if resname not in allowed_resnames:
                continue
            atom_names = rings[resname]
            if not all(atom_name in residue for atom_name in atom_names):
                continue

            coordinates = np.asarray(
                [residue[atom_name].coord for atom_name in atom_names], dtype=float
            )
            if not np.isfinite(coordinates).all():
                raise ValueError(
                    f"residue {resname} {residue.id[1]} has non-finite coordinates"
                )
            centroid = coordinates.mean(axis=0)
            normal = fit_plane_normal(coordinates)
            rings_found.append((resname, int(residue.id[1]), centroid, normal))
        return rings_found

    def clamp_unit(value: float) -> float:
        return min(1.0, max(-1.0, value))

    def acute_angle_degrees(dot_product: float) -> float:
        return math.degrees(math.acos(abs(clamp_unit(dot_product))))

    def slip_angle_degrees(
        displacement: np.ndarray, protein_normal: np.ndarray
    ) -> float:
        normal_component = abs(float(np.dot(displacement, protein_normal)))
        distance = float(np.linalg.norm(displacement))
        return math.degrees(math.asin(clamp_unit(normal_component / distance)))

    if not isinstance(range_pdbs, list):
        raise ValueError("range_pdbs must be provided as a list of PDB strings")
    if any(not isinstance(pdb_text, str) for pdb_text in range_pdbs):
        raise ValueError("every range_pdbs entry must be a PDB string")
    if len(range_pdbs) < 2 or any(not pdb_text.strip() for pdb_text in range_pdbs):
        return []

    models = []
    for index, pdb_text in enumerate(range_pdbs, start=1):
        try:
            structure = PDBParser(QUIET=True).get_structure(
                f"range-{index}", StringIO(pdb_text)
            )
        except Exception as exc:
            raise ValueError(
                f"range_pdbs[{index - 1}] is not valid PDB text: {exc}"
            ) from exc
        model = next(structure.get_models(), None)
        if model is None or not any(model.get_atoms()):
            raise ValueError(
                f"range_pdbs[{index - 1}] contains no atoms in its first model"
            )
        models.append(model)

    rna_ring_values = complete_rings(models[0], rna_rings)
    if not rna_ring_values:
        return []

    results = []
    for protein_range_pdb_id, model in enumerate(models[1:], start=2):
        protein_ring_values = complete_rings(model, protein_rings)
        for (
            protein_resname,
            protein_resid,
            protein_centroid,
            protein_normal,
        ) in protein_ring_values:
            for rna_resname, rna_resid, rna_centroid, rna_normal in rna_ring_values:
                displacement = rna_centroid - protein_centroid
                distance = float(np.linalg.norm(displacement))
                if distance == 0:
                    continue
                tilt_angle = acute_angle_degrees(
                    float(np.dot(protein_normal, rna_normal))
                )
                if distance > max_centroid_distance + distance_boundary_tolerance:
                    continue
                if tilt_angle > max_plane_tilt + angle_boundary_tolerance:
                    continue
                results.append(
                    {
                        "rna_resname": rna_resname,
                        "rna_resid": rna_resid,
                        "protein_range_pdb_id": protein_range_pdb_id,
                        "protein_resid": protein_resid,
                        "distance": distance,
                        "tilt_angle": tilt_angle,
                        "slip_angle": slip_angle_degrees(
                            displacement, protein_normal
                        ),
                        "rna_centroid": [float(value) for value in rna_centroid],
                        "protein_centroid": [
                            float(value) for value in protein_centroid
                        ],
                    }
                )

    results.sort(
        key=lambda row: (
            row["protein_range_pdb_id"],
            row["protein_resid"],
            row["rna_resid"],
        )
    )
    return results
