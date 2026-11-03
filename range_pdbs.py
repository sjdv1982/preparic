"""Range extraction and PDB-file materialization transformers."""

MAX_RNA_BACKBONE_BOND_DISTANCE = 2.2


def _ranges_to_pdb(cif: str, ranges: dict) -> list[str]:
    """Extract selected CIF ranges as PDB text, validating polymer annotations."""
    import copy
    from collections import defaultdict
    import re
    from io import StringIO

    from Bio.PDB import MMCIFParser, PDBIO
    import numpy as np

    # Keep this value local: Seamless may execute this function from its body
    # without the module globals that normally accompany an in-process call.
    max_rna_backbone_bond_distance = 2.2

    # T0 now publishes the canonical S5 envelope. Keep the legacy field as a
    # read-only compatibility path for callers outside the shared project.
    range_specs = None
    if isinstance(ranges, dict):
        if "schema" in ranges or "rows" in ranges:
            if type(ranges.get("schema")) is not int or ranges.get("schema") != 1:
                raise ValueError("The ranges cell must use schema: 1.")
            range_specs = ranges.get("rows")
        else:
            range_specs = ranges.get("ranges")
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
    parsed_chain_ids = {chain.id for model in structure for chain in model}

    def category_values(key):
        value = cif_dictionary.get(key, [])
        if isinstance(value, list):
            return value
        return [value] if value not in (None, "") else []

    def operator_table():
        identifiers = category_values("_pdbx_struct_oper_list.id")
        result = {}
        for index, identifier in enumerate(identifiers):
            matrix = np.asarray([
                [float(category_values(f"_pdbx_struct_oper_list.matrix[{row}][{col}]")[index]) for col in (1, 2, 3)]
                for row in (1, 2, 3)
            ], dtype=float)
            vector = np.asarray([
                float(category_values(f"_pdbx_struct_oper_list.vector[{row}]")[index])
                for row in (1, 2, 3)
            ], dtype=float)
            result[str(identifier)] = (matrix, vector)
        return result

    def expand_operator_expression(expression):
        """Expand mmCIF operator lists, ranges, and parenthesized products."""
        expression = str(expression).strip()
        groups = re.findall(r"\(([^()]*)\)", expression)
        if not groups:
            groups = [expression]

        def expand_group(group):
            values = []
            for item in group.split(","):
                item = item.strip()
                match = re.fullmatch(r"([^,()]+?)\s*-\s*([^,()]+)", item)
                if match and match.group(1).isdigit() and match.group(2).isdigit():
                    start, end = int(match.group(1)), int(match.group(2))
                    step = 1 if end >= start else -1
                    values.extend(str(number) for number in range(start, end + step, step))
                elif item:
                    values.append(item)
            return values

        choices = [expand_group(group) for group in groups]
        if not choices or any(not group for group in choices):
            raise ValueError(f"Invalid assembly operator expression {expression!r}.")
        combinations = [[]]
        for group in choices:
            # Mol* expandOperators varies the leftmost group fastest.
            combinations = [prefix + [identifier] for identifier in group for prefix in combinations]
        return combinations

    def assembly_chain_copies(biounit):
        """Map Mol*'s displayed chain-copy names to source asym ids and transforms."""
        generators = category_values("_pdbx_struct_assembly_gen.assembly_id")
        assembly_ids = [str(value) for value in category_values("_pdbx_struct_assembly.id")]
        if not assembly_ids and generators:
            # Small synthetic mmCIFs sometimes omit the descriptive assembly
            # category while retaining the generator rows used for expansion.
            assembly_ids = list(dict.fromkeys(str(value) for value in generators))
        asym_lists = category_values("_pdbx_struct_assembly_gen.asym_id_list")
        expressions = category_values("_pdbx_struct_assembly_gen.oper_expression")
        operators = operator_table()
        if not assembly_ids:
            return {
                chain_id: (chain_id, np.eye(3), np.zeros(3))
                for chain_id in parsed_chain_ids
            }
        if type(biounit) is not int or not 1 <= biounit <= len(assembly_ids):
            raise ValueError(f"Biological assembly {biounit!r} is not available in the CIF.")
        selected_id = assembly_ids[biounit - 1]
        copies = defaultdict(list)
        for assembly_id, asym_list, expression in zip(generators, asym_lists, expressions):
            if str(assembly_id) != selected_id:
                continue
            asym_ids = [item.strip() for item in str(asym_list).split(",") if item.strip() not in (".", "?")]
            for operation_ids in expand_operator_expression(expression):
                rotation = np.eye(3)
                translation = np.zeros(3)
                for operation_id in operation_ids:
                    if operation_id not in operators:
                        raise ValueError(f"Assembly {selected_id} refers to missing operator {operation_id!r}.")
                    next_rotation, next_translation = operators[operation_id]
                    # Mol* getAssemblyOperators right-multiplies each matrix.
                    translation = rotation @ next_translation + translation
                    rotation = rotation @ next_rotation
                for asym_id in asym_ids:
                    copies[asym_id].append((rotation, translation))

        result = {}
        for asym_id, transforms in copies.items():
            identity_index = next((index for index, (rotation, translation) in enumerate(transforms)
                                   if np.allclose(rotation, np.eye(3), atol=1e-8, rtol=0)
                                   and np.allclose(translation, 0, atol=1e-8, rtol=0)), -1)
            for index, (rotation, translation) in enumerate(transforms):
                if len(transforms) == 1 or index == identity_index:
                    display_name = asym_id
                else:
                    display_index = index + (0 if identity_index >= 0 and index > identity_index else 1)
                    display_name = f"{asym_id}-{display_index}"
                result[display_name] = (asym_id, rotation, translation)
        return result

    pdb_texts = []
    selected_biounit = ranges.get("biounit", 1) if isinstance(ranges, dict) else 1
    chain_copies = assembly_chain_copies(selected_biounit)

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

        copy_entry = chain_copies.get(requested_chain)
        if copy_entry is None:
            raise ValueError(
                f"Range {index + 1} refers to chain {requested_chain!r}, "
                "which is not present in the selected biological assembly."
            )
        chain_id, assembly_rotation, assembly_translation = copy_entry

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
        observed_residues = set()
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
                    observed_residues.add(residue_id)
                    selected_atom_count += sum(1 for _ in residue.get_atoms())

            if index == 0:
                ordered = sorted(
                    (residue for residue in selected_chain.get_residues()
                     if isinstance(residue.id[1], int) and start <= residue.id[1] <= end),
                    key=lambda residue: residue.id[1],
                )
                for before, after in zip(ordered, ordered[1:]):
                    if after.id[1] != before.id[1] + 1:
                        continue
                    before_atoms = {atom.get_name().strip().replace("*", "'"): atom for atom in before.get_atoms()}
                    after_atoms = {atom.get_name().strip().replace("*", "'"): atom for atom in after.get_atoms()}
                    oxygen = before_atoms.get("O3'")
                    phosphorus = after_atoms.get("P")
                    if oxygen is not None and phosphorus is not None:
                        bond_distance = float(oxygen - phosphorus)
                        if (
                            not np.isfinite(bond_distance)
                            or bond_distance > max_rna_backbone_bond_distance
                        ):
                            raise ValueError(
                                f"RNA range {index + 1} chain {requested_chain!r} has a chain break "
                                f"between residues {before.id[1]} and {after.id[1]}."
                            )

            if not any(selected_chain.get_atoms()):
                selected_structure.detach_child(model.id)
                continue
            model.detach_child(chain_id)
            selected_chain.id = pdb_chain_id
            if not (np.allclose(assembly_rotation, np.eye(3), atol=1e-8, rtol=0)
                    and np.allclose(assembly_translation, 0, atol=1e-8, rtol=0)):
                for atom in selected_chain.get_atoms():
                    atom.set_coord(assembly_rotation @ atom.get_coord() + assembly_translation)
            model.add(selected_chain)

        if index == 0:
            missing_residues = [resid for resid in range(start, end + 1) if resid not in observed_residues]
            if missing_residues:
                missing_text = ", ".join(str(resid) for resid in missing_residues)
                raise ValueError(
                    f"RNA range {index + 1} chain {requested_chain!r} has unobserved residues: {missing_text}."
                )

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
