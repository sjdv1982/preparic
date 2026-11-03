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


def test_in_process_range_extraction_does_not_require_module_globals(cif_text, monkeypatch):
    import range_pdbs

    # Delayed/in-process execution may retain the function body without its
    # original module namespace. The backbone cutoff must remain available.
    monkeypatch.delattr(range_pdbs, "MAX_RNA_BACKBONE_BOND_DISTANCE")
    result = _ranges_to_pdb(cif_text, {"ranges": [
        {"chain": "A", "start": 1, "end": 2},
    ]})
    assert len(result) == 1
    assert residues(result[0]) == [("A", 1), ("A", 2)]


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
    ({"chain": "absent", "start": 1, "end": 4}, "not present in the selected biological assembly"),
    ({"chain": "A", "start": 1000, "end": 1001}, "unobserved residues: 1000, 1001"),
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
    data = MMCIF2Dict(StringIO(cif_text))
    data["_pdbx_struct_oper_list.id"] = ["1", "2", "3"]
    data["_pdbx_struct_oper_list.type"] = [
        "identity operation", "translation operation", "translation operation",
    ]
    data["_pdbx_struct_oper_list.name"] = ["1_555", "copy2", "copy3"]
    data["_pdbx_struct_oper_list.symmetry_operation"] = ["x,y,z", "x+10,y,z", "x+20,y,z"]
    translations = [(0, 0, 0), (10, 0, 0), (20, 0, 0)]
    for row in range(1, 4):
        data[f"_pdbx_struct_oper_list.vector[{row}]"] = [str(value[row - 1]) for value in translations]
        for column in range(1, 4):
            data[f"_pdbx_struct_oper_list.matrix[{row}][{column}]"] = [
                "1" if row == column else "0" for _ in translations
            ]
    data["_pdbx_struct_assembly_gen.assembly_id"] = ["1"]
    data["_pdbx_struct_assembly_gen.oper_expression"] = ["(1,2,3)"]
    data["_pdbx_struct_assembly_gen.asym_id_list"] = ["A,C"]
    writer = MMCIFIO()
    writer.set_dict(data)
    synthetic = StringIO()
    writer.save(synthetic)

    range_specs["ranges"][0]["chain"] = "A-2"
    extracted = _ranges_to_pdb(synthetic.getvalue(), range_specs)
    assert [residues(pdb) for pdb in extracted] == [residues(pdb) for pdb in expected_pdbs]


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
    data["_pdbx_struct_assembly_gen.asym_id_list"] = [
        ",".join("RNA" if asym == "A" else asym for asym in value.split(","))
        for value in data["_pdbx_struct_assembly_gen.asym_id_list"]
    ]
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


@pytest.mark.parametrize("bond_distance,accepted", [(2.1, True), (2.3, False)],
                         ids=["below-shared-2.2-cutoff", "above-shared-2.2-cutoff"])
def test_consecutive_rna_labels_respect_shared_backbone_distance(cif_text, bond_distance, accepted):
    """Preparic must reject the 2.3 Å link that Alaric's 2.2 Å rule rejects."""
    data = MMCIF2Dict(StringIO(cif_text))

    def atom_index(residue, name):
        matches = [index for index, (chain, seq, atom) in enumerate(zip(
            data["_atom_site.label_asym_id"], data["_atom_site.label_seq_id"],
            data["_atom_site.label_atom_id"], strict=True,
        )) if chain == "A" and seq == str(residue) and atom.replace("*", "'") == name]
        assert len(matches) == 1
        return matches[0]

    oxygen = atom_index(1, "O3'")
    phosphorus = atom_index(2, "P")
    # Keep labels 1→2 consecutive and alter only the connecting phosphorus.
    # An axis-aligned displacement fixes the intended distance independently
    # of the original crystallographic bond direction or length.
    for axis in "xyz":
        field = f"_atom_site.Cartn_{axis}"
        coordinate = float(data[field][oxygen])
        data[field][phosphorus] = str(coordinate + (bond_distance if axis == "x" else 0.0))
    writer = MMCIFIO()
    writer.set_dict(data)
    output = StringIO()
    writer.save(output)
    specs = {"schema": 1, "pdbcode": "1B7F", "biounit": 1,
             "rows": [{"chain": "A", "start": 1, "end": 2}]}
    if accepted:
        extracted = _ranges_to_pdb(output.getvalue(), specs)
        assert len(extracted) == 1
        assert residues(extracted[0]) == [("A", 1), ("A", 2)]
    else:
        with pytest.raises(ValueError, match=r"(?i)(break|continu|connect|bond)"):
            _ranges_to_pdb(output.getvalue(), specs)
