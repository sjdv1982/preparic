"""Direct Python tests for selection, validation and folder serialization."""

from io import StringIO

from Bio.PDB import MMCIFIO, PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import pytest
import yaml

from range_pdbs import _pdbs_to_folder, _ranges_to_pdb
from conftest import DATA


def residues(text):
    structure = PDBParser(QUIET=True).get_structure("test", StringIO(text))
    return [(r.parent.id, r.id[1]) for r in structure.get_residues()]


def test_extracts_rna_then_proteins_in_range_order(cif_text, range_specs, expected_pdbs):
    result = _ranges_to_pdb(cif_text, range_specs)
    assert result == expected_pdbs
    assert [residues(pdb) for pdb in result] == [
        [("A", n) for n in range(1, 5)],
        [("C", n) for n in range(2, 9)],
        [("C", n) for n in range(90, 97)],
    ]


def test_original_deployment_selections(cif_text):
    specs = yaml.safe_load((DATA / "deployment-ranges.yaml").read_text())
    result = _ranges_to_pdb(cif_text, specs)
    assert [residues(pdb) for pdb in result] == [
        [("A", n) for n in range(1, 13)],
        [("C", n) for n in range(2, 81)],
        [("C", n) for n in range(90, 169)],
    ]


@pytest.mark.parametrize("cif,specs", [
    ("not parsed", {"ranges": []}),
    ("", {"ranges": [{"chain": "A", "start": 1, "end": 4}]}),
    ("  \n", {"ranges": [{"chain": "A", "start": 1, "end": 4}]}),
])
def test_empty_inputs(cif, specs):
    assert _ranges_to_pdb(cif, specs) == []


@pytest.mark.parametrize("specs", [None, [], {}, {"ranges": {}}, {"ranges": None}])
def test_requires_ranges_list(specs):
    with pytest.raises(ValueError, match="must contain a ranges list"):
        _ranges_to_pdb("", specs)


@pytest.mark.parametrize("spec,message", [
    (None, "Range 1 must be a mapping"),
    ({}, "Range 1 is missing its chain"),
    ({"chain": 1}, "Range 1 is missing its chain"),
    ({"chain": "A", "end": 4}, "missing its start residue"),
    ({"chain": "A", "start": 1}, "missing its end residue"),
    ({"chain": "A", "start": 4, "end": 1}, "starts after it ends"),
    ({"chain": "absent", "start": 1, "end": 4}, "not present in the CIF"),
    ({"chain": "A", "start": 1000, "end": 1001}, "contains no atoms"),
    ({"chain": "C", "start": 2, "end": 4}, "expected RNA"),
])
def test_invalid_range(cif_text, spec, message):
    with pytest.raises(ValueError, match=message):
        _ranges_to_pdb(cif_text, {"ranges": [spec]})


@pytest.mark.parametrize("value", [True, False, 1.5, "1.5", "1e0", "", None, []])
@pytest.mark.parametrize("boundary", ["start", "end"])
def test_rejects_noninteger_boundaries(cif_text, value, boundary):
    spec = {"chain": "A", "start": 1, "end": 4, boundary: value}
    with pytest.raises(ValueError, match="must be integer residue numbers"):
        _ranges_to_pdb(cif_text, {"ranges": [spec]})


@pytest.mark.parametrize("start,end", [(1.0, 4.0), (" +1 ", "04")])
def test_accepts_integer_valued_boundaries(cif_text, expected_pdbs, start, end):
    assert _ranges_to_pdb(cif_text, {"ranges": [
        {"chain": "A", "start": start, "end": end},
    ]}) == expected_pdbs[:1]


def test_assembly_copy_suffix(cif_text, range_specs, expected_pdbs):
    range_specs["ranges"][0]["chain"] = "A-2"
    assert _ranges_to_pdb(cif_text, range_specs) == expected_pdbs


def test_later_ranges_must_be_protein(cif_text):
    with pytest.raises(ValueError, match="Range 2.*expected protein"):
        _ranges_to_pdb(cif_text, {"ranges": [
            {"chain": "A", "start": 1, "end": 4},
            {"chain": "A", "start": 5, "end": 6},
        ]})


@pytest.mark.parametrize("annotation,message", [
    (None, "no polymer annotation"),
    ("polydeoxyribonucleotide", "expected RNA"),
])
def test_polymer_annotations(cif_text, annotation, message):
    data = MMCIF2Dict(StringIO(cif_text))
    if annotation is None:
        del data["_entity_poly.type"]
    else:
        data["_entity_poly.type"] = [
            annotation if value == "polyribonucleotide" else value
            for value in data["_entity_poly.type"]
        ]
    writer = MMCIFIO()
    writer.set_dict(data)
    output = StringIO()
    writer.save(output)
    with pytest.raises(ValueError, match=message):
        _ranges_to_pdb(output.getvalue(), {"ranges": [
            {"chain": "A", "start": 1, "end": 4},
        ]})


def test_multicharacter_chain_is_written_as_a(cif_text, range_specs, expected_pdbs):
    data = MMCIF2Dict(StringIO(cif_text))
    for field in ("_struct_asym.id", "_atom_site.label_asym_id"):
        data[field] = ["RNA" if value == "A" else value for value in data[field]]
    writer = MMCIFIO()
    writer.set_dict(data)
    output = StringIO()
    writer.save(output)
    range_specs["ranges"][0]["chain"] = "RNA"
    assert _ranges_to_pdb(output.getvalue(), range_specs) == expected_pdbs


def test_numbered_folder_contains_raw_pdb_bytes(expected_pdbs):
    assert _pdbs_to_folder(expected_pdbs) == {
        f"range-{i}.pdb": text.encode("utf-8")
        for i, text in enumerate(expected_pdbs, 1)
    }
    assert _pdbs_to_folder([]) == {}


@pytest.mark.parametrize("value", [None, {}, "PDB", ("PDB",)])
def test_folder_requires_list(value):
    with pytest.raises(ValueError, match="must be provided as a list"):
        _pdbs_to_folder(value)
