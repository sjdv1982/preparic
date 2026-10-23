"""Range extraction and PDB-file materialization transformers."""


def _ranges_to_pdb(cif: str, ranges: dict) -> list[str]:
    """Extract selected CIF ranges as PDB text, validating polymer annotations."""
    import copy
    import re
    from io import StringIO

    from Bio.PDB import MMCIFParser, PDBIO

    range_specs = ranges.get("ranges") if isinstance(ranges, dict) else None
    if not isinstance(range_specs, list):
        raise ValueError("The ranges cell must contain a ranges list.")
    if not range_specs:
        return []
    if not isinstance(cif, str) or not cif.strip():
        return []

    parser = MMCIFParser(auth_chains=False, auth_residues=False, QUIET=True)
    structure = parser.get_structure("cif", StringIO(cif))
    cif_dictionary = parser._mmcif_dict
    # MMCIFParser keeps the parsed categories here; reuse them instead of
    # tokenizing a potentially large CIF a second time.
    asym_ids = cif_dictionary.get("_struct_asym.id", [])
    asym_entity_ids = cif_dictionary.get("_struct_asym.entity_id", [])
    entity_ids = cif_dictionary.get("_entity_poly.entity_id", [])
    entity_types = cif_dictionary.get("_entity_poly.type", [])
    chain_entities = dict(zip(asym_ids, asym_entity_ids))
    polymer_types = dict(zip(entity_ids, entity_types))
    parsed_chain_ids = {
        chain.id for model in structure for chain in model
    }

    pdb_texts = []
    for index, range_spec in enumerate(range_specs):
        if not isinstance(range_spec, dict):
            raise ValueError(f"Range {index + 1} must be a mapping.")
        requested_chain = range_spec.get("chain")
        if not isinstance(requested_chain, str) or not requested_chain:
            raise ValueError(f"Range {index + 1} is missing its chain.")
        try:
            raw_start = range_spec["start"]
            raw_end = range_spec["end"]
        except KeyError as exc:
            raise ValueError(
                f"Range {index + 1} is missing its {exc.args[0]} residue."
            ) from exc

        def residue_number(value):
            if isinstance(value, bool):
                raise ValueError
            if isinstance(value, int):
                return value
            if isinstance(value, float) and value.is_integer():
                return int(value)
            if isinstance(value, str) and re.fullmatch(r"[+-]?[0-9]+", value.strip()):
                return int(value)
            raise ValueError

        try:
            start = residue_number(raw_start)
            end = residue_number(raw_end)
        except ValueError as exc:
            raise ValueError(
                f"Range {index + 1} start and end must be integer residue numbers."
            ) from exc
        if start > end:
            raise ValueError(
                f"Range {index + 1} starts after it ends ({start} > {end})."
            )

        chain_id = requested_chain
        if chain_id not in parsed_chain_ids:
            suffix_match = re.fullmatch(r"(.+)-[0-9]+", chain_id)
            if suffix_match and suffix_match.group(1) in parsed_chain_ids:
                chain_id = suffix_match.group(1)
            else:
                raise ValueError(
                    f"Range {index + 1} refers to chain {requested_chain!r}, "
                    "which is not present in the CIF."
                )

        pdb_chain_id = chain_id if len(chain_id) == 1 else "A"
        entity_id = chain_entities.get(chain_id)
        annotation = polymer_types.get(entity_id)
        if annotation is None:
            raise ValueError(
                f"Chain {requested_chain!r} has no polymer annotation in the CIF."
            )
        normalized_annotation = annotation.lower()
        if normalized_annotation == "polyribonucleotide":
            polymer_kind = "RNA"
        elif normalized_annotation.startswith("polypeptide"):
            polymer_kind = "protein"
        else:
            polymer_kind = annotation
        expected_kind = "RNA" if index == 0 else "protein"
        if polymer_kind != expected_kind:
            raise ValueError(
                f"Range {index + 1} chain {requested_chain!r} is annotated as "
                f"{annotation!r}; expected {expected_kind}."
            )

        selected_structure = copy.deepcopy(structure)
        selected_atom_count = 0
        for model in list(selected_structure):
            selected_chain = next(
                (chain for chain in model if chain.id == chain_id), None
            )
            for chain in list(model):
                if chain is not selected_chain:
                    model.detach_child(chain.id)
            if selected_chain is None:
                selected_structure.detach_child(model.id)
                continue

            for residue in list(selected_chain):
                residue_id = residue.id[1]
                if not isinstance(residue_id, int) or not start <= residue_id <= end:
                    selected_chain.detach_child(residue.id)
                else:
                    selected_atom_count += sum(1 for _ in residue.get_atoms())

            if not any(selected_chain.get_atoms()):
                selected_structure.detach_child(model.id)
                continue
            model.detach_child(chain_id)
            selected_chain.id = pdb_chain_id
            model.add(selected_chain)

        if selected_atom_count == 0:
            raise ValueError(
                f"Range {index + 1} chain {requested_chain!r} residues "
                f"{start}-{end} contains no atoms."
            )

        output = StringIO()
        writer = PDBIO()
        writer.set_structure(selected_structure)
        writer.save(output)
        pdb_texts.append(output.getvalue())

    return pdb_texts


def _pdbs_to_folder(pdb_texts: list[str]) -> dict[str, bytes]:
    """Store each range PDB as a numbered file in a folder cell."""
    if not isinstance(pdb_texts, list):
        raise ValueError("Range PDBs must be provided as a list.")
    return {
        f"range-{index + 1}.pdb": pdb_text.encode("utf-8")
        for index, pdb_text in enumerate(pdb_texts)
    }
