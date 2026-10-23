"""Backend-style standalone transformations, without workflow or mounts."""

from io import StringIO

from Bio.PDB import PDBParser
import pytest
from seamless import Checksum
from seamless_transformer import Transformer, delayed
from seamless_transformer.transformation_class import Transformation, TransformationError

from detect_stacking import _detect_stacking
from range_pdbs import _pdbs_to_folder, _ranges_to_pdb
from conftest import PROJECT_ROOT


pytestmark = [pytest.mark.integration, pytest.mark.usefixtures("seamless_runtime")]


def folder_transformation(pdb_texts):
    builder = delayed(_pdbs_to_folder)
    builder.celltypes.pdb_texts = "plain"
    builder.celltypes.result = "folder"
    return builder.build(pdb_texts=pdb_texts)


def reducer(pdbs, aareduce):
    builder = Transformer("bash")
    builder.code = (PROJECT_ROOT / "range-pdbs-aareduce.sh").read_text()
    builder.environment.set_conda_env("preparic")
    builder.pins.pdbs = pdbs
    builder.celltypes.pdbs = "folder"
    builder.pins.aareduce = aareduce
    builder.celltypes.aareduce = "folder"
    builder.celltypes.result = "folder"
    return builder.build()


def folder_contents(transformation):
    index = transformation.run()
    assert all(isinstance(value, Checksum) for value in index.values())
    # Folder results are indexes; resolve contents explicitly for assertions.
    return {name: checksum.resolve("bytes").content for name, checksum in index.items()}


def assert_reduced(result, expected_pdbs):
    assert set(result) == {f"range-{i}-aa.pdb" for i in (1, 2, 3)}
    parser = PDBParser(QUIET=True)
    for i, original in enumerate(expected_pdbs, 1):
        source = parser.get_structure("source", StringIO(original))
        reduced = parser.get_structure(
            "reduced", StringIO(result[f"range-{i}-aa.pdb"].decode()),
        )
        assert [(r.parent.id, r.id) for r in reduced.get_residues()] == [
            (r.parent.id, r.id) for r in source.get_residues()
        ]
        atoms = list(reduced.get_atoms())
        assert atoms
        assert all(atom.element not in {"H", "D"} for atom in atoms)
        if i > 1:
            assert all({"N", "CA", "C", "O"} <= set(r.child_dict)
                       for r in reduced.get_residues())


def test_python_transformations(cif_text, range_specs, expected_pdbs):
    extraction = delayed(_ranges_to_pdb).build(cif=cif_text, ranges=range_specs)
    assert isinstance(extraction, Transformation)
    assert extraction.run() == expected_pdbs
    folder = folder_transformation(extraction)
    assert folder_contents(folder) == {
        f"range-{i}.pdb": pdb.encode() for i, pdb in enumerate(expected_pdbs, 1)
    }


@pytest.mark.xfail(
    strict=True,
    raises=AttributeError,
    reason="Seamless Checksum.resolution creates a trivial empty-list Buffer without calculating its checksum",
)
def test_python_empty_transformation():
    extraction = delayed(_ranges_to_pdb).build(cif="", ranges={"ranges": []})
    assert extraction.run() == []
    assert folder_contents(folder_transformation(extraction)) == {}


def test_python_validation_failure(cif_text):
    transformation = delayed(_ranges_to_pdb).build(cif=cif_text, ranges={"ranges": [
        {"chain": "C", "start": 2, "end": 8},
    ]})
    with pytest.raises(TransformationError, match="expected RNA"):
        transformation.run()


def test_bash_transformation(aareduce_folder, expected_pdbs):
    # Independent of CIF extraction: use the checked-in PDB inputs.
    pdbs = folder_transformation(expected_pdbs)
    transformation = reducer(pdbs, aareduce_folder)
    assert isinstance(transformation, Transformation)
    assert_reduced(folder_contents(transformation), expected_pdbs)


def test_bash_empty_transformation(aareduce_folder):
    result = folder_contents(reducer(folder_transformation([]), aareduce_folder))
    assert result == {"no-ranges.txt": b"No range PDB files were provided.\n"}


def test_bash_failure_propagates(aareduce_folder):
    transformation = reducer(folder_transformation(["invalid PDB\n"]), aareduce_folder)
    with pytest.raises(TransformationError, match="Bash transformer exception"):
        transformation.run()


def test_extraction_to_reduction_dependencies(cif_text, range_specs, expected_pdbs, aareduce_folder):
    extraction = delayed(_ranges_to_pdb).build(cif=cif_text, ranges=range_specs)
    pdbs = folder_transformation(extraction)
    # Pass handles all the way through; the final run executes dependencies.
    assert_reduced(folder_contents(reducer(pdbs, aareduce_folder)), expected_pdbs)


def test_extraction_to_stacking_dependency_is_content_bound(
    cif_text, range_specs, expected_pdbs
):
    extraction = delayed(_ranges_to_pdb).build(cif=cif_text, ranges=range_specs)
    detector = delayed(_detect_stacking)

    stacking_transformation = detector.build(range_pdbs=extraction)
    assert stacking_transformation.run() == _detect_stacking(expected_pdbs)


def test_wrapped_stacking_validation_failure():
    transformation = delayed(_detect_stacking).build(range_pdbs="not a list")
    with pytest.raises(TransformationError, match="range_pdbs must be provided"):
        transformation.run()


def test_wrapped_stacking_identity_and_result_are_deterministic(expected_pdbs):
    first = delayed(_detect_stacking).build(range_pdbs=expected_pdbs)
    second = delayed(_detect_stacking).build(range_pdbs=expected_pdbs)
    first.construct()
    second.construct()

    assert first.transformation_checksum == second.transformation_checksum
    assert first.run() == second.run() == _detect_stacking(expected_pdbs)


@pytest.mark.xfail(
    strict=True,
    raises=AttributeError,
    reason="Seamless Checksum.resolution creates a trivial empty-list Buffer without calculating its checksum",
)
def test_wrapped_stacking_empty_pipeline_is_stable():
    assert delayed(_detect_stacking).build(range_pdbs=[]).run() == []
