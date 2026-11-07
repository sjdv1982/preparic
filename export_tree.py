"""Pure renderer for the export archive's folder-cell mapping.

This module does not read files, construct a Context, or write a ZIP.  It
receives the already-normalised snapshot and checked generator catalogue and
returns the exact relative paths and bytes held by the ``export_tree`` cell.
"""

from __future__ import annotations

import hashlib
import json
import re
from typing import Any

import yaml

import export_gens
import export_snapshot


_GEN_RECIPE_KEYS = {
    "path", "modes", "range", "direction", "final_pools", "external_inputs", "actions",
}
_RANGE_FILE = re.compile(r"range-([1-9][0-9]*)\.pdb\Z")
_RANGE_AA_FILE = re.compile(r"range-([1-9][0-9]*)-aa\.pdb\Z")
_FRAGMENT_KEY = re.compile(r"frag-([1-9][0-9]*)-([A-Z]{2})\Z")


def _json_bytes(value: Any) -> bytes:
    """Serialize a JSON document deterministically as UTF-8 with one newline."""
    return (json.dumps(value, ensure_ascii=False, allow_nan=False, indent=2) + "\n").encode("utf-8")


def _yaml_bytes(value: dict, comments: list[str] | None = None) -> bytes:
    """Serialize an action mapping, keeping generated comments as YAML comments."""
    if type(value) is not dict:
        raise ValueError("An action mapping must be a dictionary.")
    comment_lines: list[str] = []
    if comments is not None:
        if type(comments) is not list or any(not isinstance(comment, str) for comment in comments):
            raise ValueError("Action comments must be a list of strings.")
        for comment in comments:
            lines = comment.splitlines()
            if not lines:
                lines = [""]
            comment_lines.extend("#" + (f" {line}" if line else "") for line in lines)
    body = yaml.safe_dump(value, sort_keys=False, allow_unicode=True,
                          default_flow_style=False)
    if not body.endswith("\n"):
        body += "\n"
    prefix = "\n".join(comment_lines)
    if prefix:
        prefix += "\n"
    return (prefix + body).encode("utf-8")


def _checked_bytes(value: Any, address: str) -> bytes:
    if not isinstance(value, bytes):
        raise ValueError(f"Expected exact byte content for {address}.")
    return value


def _add_member(tree: dict[str, bytes], path: str, content: bytes) -> None:
    if not isinstance(path, str) or not path or path.startswith("/") or "\\" in path or "\x00" in path:
        raise ValueError(f"Unsafe generated archive path: {path!r}.")
    parts = path.split("/")
    if any(not part or part in (".", "..")
           or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", part) is None for part in parts):
        raise ValueError(f"Unsafe generated archive path: {path!r}.")
    if type(content) is not bytes:
        raise ValueError(f"Generated archive member {path} is not bytes.")
    if path in tree:
        raise ValueError(f"Duplicate generated archive path: {path}.")
    tree[path] = content


def _validate_member_tree(tree: dict[str, bytes]) -> None:
    keys = set(tree)
    for path in keys:
        parts = path.split("/")
        for index in range(1, len(parts)):
            prefix = "/".join(parts[:index])
            if prefix in keys:
                raise ValueError(f"Generated archive path {prefix} is a directory prefix of {path}.")


def _validate_recipes(recipes: Any) -> list[dict]:
    if type(recipes) is not list:
        raise ValueError("A generator catalogue must be a list.")
    seen_paths: set[str] = set()
    for recipe in recipes:
        if type(recipe) is not dict or set(recipe) != _GEN_RECIPE_KEYS:
            raise ValueError("A generator recipe has an unexpected shape.")
        export_gens.validate_recipe(recipe)
        path = recipe["path"]
        if path in seen_paths:
            raise ValueError(f"Duplicate generator path: {path}.")
        for previous in seen_paths:
            if path.startswith(previous + "/") or previous.startswith(path + "/"):
                raise ValueError(f"Generator path {path} collides by directory prefix with {previous}.")
        seen_paths.add(path)

    checked: list[dict] = []
    for recipe in recipes:
        path = recipe["path"]

        modes = recipe["modes"]
        if modes not in (["fixed"], ["unfixed"], ["fixed", "unfixed"]):
            raise ValueError(f"Generator {path} has invalid mode metadata.")
        if recipe["range"] is not None and not isinstance(recipe["range"], str):
            raise ValueError(f"Generator {path} has invalid range metadata.")
        if recipe["direction"] not in (None, "fwd", "bwd"):
            raise ValueError(f"Generator {path} has invalid direction metadata.")
        if not isinstance(recipe["final_pools"], list) or not recipe["final_pools"] \
                or any(not isinstance(pool, str) or not pool for pool in recipe["final_pools"]):
            raise ValueError(f"Generator {path} has invalid final-pool metadata.")

        parts = path.split("/")
        root = parts[0]
        if root == "common-conformer":
            if len(parts) != 2:
                raise ValueError(f"Common-conformer generator path has extra levels: {path!r}.")
            match = re.fullmatch(
                r"gen-frag([1-9][0-9]*)-range-([1-9][0-9]*)-([1-9][0-9]*)", parts[1]
            )
            if match is None or int(match.group(2)) >= int(match.group(3)):
                raise ValueError(f"Common-conformer generator path is malformed: {path!r}.")
            if modes != ["unfixed"] or recipe["range"] is not None \
                    or recipe["direction"] is not None:
                raise ValueError(f"Common-conformer generator {path} has inconsistent metadata.")
            checked.append(recipe)
            continue

        if root not in {"anchor", "anchor-refe", "grow", "merge"}:
            raise ValueError(f"Generator path is outside the export catalogue: {path!r}.")

        index = 1
        direction = None
        if root in {"anchor", "anchor-refe"}:
            if index < len(parts) - 1 and parts[index] in {"fwd", "bwd"}:
                direction = parts[index]
                index += 1
        else:
            if index >= len(parts) - 1 or parts[index] not in {"fwd", "bwd"}:
                raise ValueError(f"{root} generator path must specify a direction: {path!r}.")
            direction = parts[index]
            index += 1

        mode = None
        range_level = None
        if index < len(parts) - 1 and parts[index] in {"fixed", "unfixed"}:
            mode = parts[index]
            index += 1
            if mode == "unfixed" and index < len(parts) - 1:
                range_match = re.fullmatch(r"range-[1-9][0-9]*", parts[index])
                if range_match is not None:
                    range_level = parts[index]
                    index += 1

        if root == "merge" and index < len(parts) - 1 \
                and parts[index] in {"anchor", "anchor-refe"}:
            index += 1

        if index != len(parts) - 1:
            raise ValueError(f"Generator path has misplaced or unsupported levels: {path!r}.")

        terminal = parts[index]
        if root in {"anchor", "anchor-refe"} and direction is None:
            valid_terminal = re.fullmatch(r"gen-frag[1-9][0-9]*", terminal) is not None
        elif root in {"anchor", "anchor-refe"}:
            terminal_match = re.fullmatch(r"gen-frag([1-9][0-9]*)-([1-9][0-9]*)", terminal)
            valid_terminal = terminal_match is not None \
                and int(terminal_match.group(1)) < int(terminal_match.group(2))
        else:
            valid_terminal = re.fullmatch(r"gen-frag[1-9][0-9]*", terminal) is not None
        if not valid_terminal:
            raise ValueError(f"Generator path has an invalid terminal gen folder: {path!r}.")

        if recipe["direction"] != direction:
            raise ValueError(f"Generator {path} has direction metadata inconsistent with its path.")
        if mode is None:
            expected_modes = ["fixed", "unfixed"]
            expected_range = None
        elif mode == "fixed":
            expected_modes = ["fixed"]
            expected_range = None
        else:
            expected_modes = ["unfixed"]
            expected_range = recipe["range"]
            if not isinstance(expected_range, str) \
                    or re.fullmatch(r"range-[1-9][0-9]*", expected_range) is None:
                raise ValueError(f"Unfixed generator {path} has invalid range metadata.")
            if range_level is not None and range_level != expected_range:
                raise ValueError(f"Generator {path} has a range level inconsistent with its metadata.")
        if modes != expected_modes or recipe["range"] != expected_range:
            raise ValueError(f"Generator {path} has mode/range metadata inconsistent with its path.")

        checked.append(recipe)
    return checked


def _range_sources(snapshot: dict) -> list[tuple[int, bytes, bytes]]:
    range_specs = snapshot.get("range_specs")
    range_rows = range_specs.get("rows") if type(range_specs) is dict else None
    plain = snapshot.get("range_pdb_files")
    all_atom = snapshot.get("range_pdb_files_aa")
    if not isinstance(range_rows, list) or not range_rows:
        raise ValueError("The snapshot has no normalised range rows.")
    if type(plain) is not dict or type(all_atom) is not dict:
        raise ValueError("The snapshot has no prepared range-file mappings.")
    expected_plain = {f"range-{index}.pdb" for index in range(1, len(range_rows) + 1)}
    expected_all_atom = {f"range-{index}-aa.pdb" for index in range(1, len(range_rows) + 1)}
    if set(plain) != expected_plain or set(all_atom) != expected_all_atom:
        raise ValueError("Prepared range files do not match the range specifications.")
    result: list[tuple[int, bytes, bytes]] = []
    for index in range(1, len(range_rows) + 1):
        plain_name = f"range-{index}.pdb"
        all_atom_name = f"range-{index}-aa.pdb"
        result.append((index,
                       _checked_bytes(plain[plain_name], plain_name),
                       _checked_bytes(all_atom[all_atom_name], all_atom_name)))
    return result


def _constraints(snapshot: dict) -> dict:
    kept = snapshot.get("kept_fragments")
    pairs = snapshot.get("kept_pairs")
    if not isinstance(kept, list) or not isinstance(pairs, list):
        raise ValueError("The snapshot has no kept fragment or pair table.")
    result: dict[str, Any] = {"schema": 1}
    kept_numbers: set[int] = set()
    sequence_by_number: dict[int, str] = {}
    ordered_kept = sorted(kept, key=lambda row: row.get("number", -1) if type(row) is dict else -1)
    for fragment in ordered_kept:
        if type(fragment) is not dict:
            raise ValueError("A kept fragment is malformed.")
        key = fragment.get("key")
        sequence = fragment.get("sequence")
        number = fragment.get("number")
        match = _FRAGMENT_KEY.fullmatch(key) if isinstance(key, str) else None
        if match is None or type(number) is not int or int(match.group(1)) != number \
                or not isinstance(sequence, str) or sequence != match.group(2):
            raise ValueError(f"Invalid kept fragment identity: {key!r}.")
        if number in kept_numbers:
            raise ValueError(f"Duplicate kept fragment number: {number}.")
        kept_numbers.add(number)
        sequence_by_number[number] = sequence
        result[f"frag{number}"] = {"sequence": sequence}
    pair_rows = []
    seen_pairs: set[tuple[int, int]] = set()
    for pair in pairs:
        if type(pair) is not dict:
            raise ValueError("A kept grow pair is malformed.")
        down, up = pair.get("down"), pair.get("up")
        sequence = pair.get("sequence")
        crmsd, ovrmsd = pair.get("crmsd"), pair.get("ovrmsd")
        if type(down) is not int or type(up) is not int or up != down + 1 \
                or down not in kept_numbers or up not in kept_numbers \
                or not isinstance(sequence, str) or len(sequence) != 3 \
                or not isinstance(crmsd, str) or not isinstance(ovrmsd, str):
            raise ValueError(f"Invalid kept grow-pair row: {pair!r}.")
        address = (down, up)
        if address in seen_pairs:
            raise ValueError(f"Duplicate kept grow pair: {down}->{up}.")
        seen_pairs.add(address)
        down_sequence = sequence_by_number[down]
        up_sequence = sequence_by_number[up]
        if down_sequence[1] != up_sequence[0] \
                or sequence != down_sequence + up_sequence[1]:
            raise ValueError(f"Grow-pair sequence does not match the kept fragments {down}->{up}.")
        pair_rows.append({"down": down, "up": up, "sequence": sequence,
                          "crmsd": crmsd, "ovrmsd": ovrmsd})
    pair_rows.sort(key=lambda row: row["down"])
    expected_pairs = {(number, number + 1) for number in kept_numbers if number + 1 in kept_numbers}
    if seen_pairs != expected_pairs:
        raise ValueError("The kept grow-pair table does not contain exactly the adjacent kept pairs.")
    result["pairs"] = pair_rows
    return result


def _energy_thresholds(snapshot: dict) -> dict:
    rows = snapshot.get("energy_cells")
    kept = snapshot.get("kept_fragments")
    range_specs = snapshot.get("range_specs")
    ranges = range_specs.get("rows") if type(range_specs) is dict else None
    if not isinstance(rows, list) or not isinstance(kept, list) or not isinstance(ranges, list):
        raise ValueError("The snapshot has no validated energy or range table.")
    kept_keys = {row.get("key") for row in kept if type(row) is dict}
    fragment_numbers = {row.get("key"): row.get("number") for row in kept if type(row) is dict}
    range_names = {f"range-{index}" for index in range(1, len(ranges) + 1)}
    interaction = snapshot.get("interaction_matrix")
    if type(interaction) is not dict or not isinstance(interaction.get("rows"), list) \
            or not isinstance(interaction.get("ranges"), list) \
            or not isinstance(interaction.get("grid"), list):
        raise ValueError("The snapshot has no validated interaction matrix.")
    matrix_rows, matrix_ranges, grid = interaction["rows"], interaction["ranges"], interaction["grid"]
    if len(grid) != len(matrix_rows) or any(not isinstance(row, list) or len(row) != len(matrix_ranges)
                                             for row in grid):
        raise ValueError("The interaction matrix has an invalid shape.")
    expected_cells: set[tuple[str, str]] = set()
    for row_index, fragment in enumerate(matrix_rows):
        if not isinstance(fragment, str):
            raise ValueError("The interaction matrix has an invalid fragment address.")
        if fragment not in kept_keys:
            continue
        for range_index, range_name in enumerate(matrix_ranges):
            if not isinstance(range_name, str):
                raise ValueError("The interaction matrix has an invalid range address.")
            if type(grid[row_index][range_index]) is not bool:
                raise ValueError("The interaction matrix contains a non-boolean cell.")
            if grid[row_index][range_index]:
                expected_cells.add((fragment, range_name))
    result_rows = []
    seen: set[tuple[str, str]] = set()
    for row in rows:
        if type(row) is not dict:
            raise ValueError("An energy threshold row is malformed.")
        fragment, range_name, threshold = row.get("fragment"), row.get("range"), row.get("threshold")
        if not isinstance(fragment, str) or fragment not in kept_keys \
                or not isinstance(range_name, str) or range_name not in range_names \
                or not isinstance(threshold, str):
            raise ValueError("An energy threshold is not addressed to a kept true cell.")
        address = (fragment, range_name)
        if address not in expected_cells:
            raise ValueError("An energy threshold is not addressed to a true interaction cell.")
        if address in seen:
            raise ValueError(f"Duplicate energy threshold for {fragment}/{range_name}.")
        seen.add(address)
        result_rows.append((fragment, range_name, threshold))
    if seen != expected_cells:
        raise ValueError("The energy threshold table does not cover exactly the true kept cells.")
    result: dict[str, Any] = {"schema": 1}
    result_rows.sort(key=lambda row: (fragment_numbers[row[0]], int(row[1].split("-")[1])))
    for fragment, range_name, threshold in result_rows:
        result.setdefault(fragment, {})[range_name] = threshold
    return result


def _disabled_reference_anchors(snapshot: dict, recipes: list[dict]) -> list[dict]:
    threshold_rows = snapshot.get("threshold_bases")
    fragment_rows = snapshot.get("fragments")
    if not isinstance(threshold_rows, list) or not isinstance(fragment_rows, list):
        raise ValueError("The snapshot has no reference-threshold evidence.")
    fragments_by_number = {row.get("number"): row for row in fragment_rows if type(row) is dict}
    thresholds_by_address = {
        (row["fragment"], row["sequence"]): row
        for row in threshold_rows
        if type(row) is dict and type(row.get("fragment")) is int
        and isinstance(row.get("sequence"), str)
    }
    candidates: dict[tuple[str, str], dict] = {}
    for recipe in recipes:
        for action in recipe["actions"]:
            parameters = action.get("parameters")
            if type(parameters) is not dict or parameters.get("action") != "anchor-refe" \
                    or type(parameters.get("ovrmsd")) is not int or parameters.get("ovrmsd") != -1:
                continue
            number = parameters.get("fragment")
            nucleotide = parameters.get("nucleotide")
            fragment = fragments_by_number.get(number)
            if fragment is None or nucleotide not in ("first", "second"):
                raise ValueError(f"Malformed disabled reference candidate in {recipe['path']}.")
            key = fragment.get("key")
            if not isinstance(key, str):
                raise ValueError(f"Malformed disabled reference candidate in {recipe['path']}.")
            threshold = thresholds_by_address.get((number, fragment.get("sequence")))
            suffix = "nuc1" if nucleotide == "first" else "nuc2"
            raw = threshold.get(f"{suffix}_ovrmsd_raw") if threshold else None
            rounded = threshold.get(f"{suffix}_ovrmsd") if threshold else None
            if threshold is None or raw is None or rounded is None:
                raise ValueError(f"Missing base threshold for disabled candidate {key}/{nucleotide}.")
            candidate_key = (key, nucleotide)
            record = candidates.get(candidate_key)
            if record is None:
                record = {"candidate": {"fragment": key, "nucleotide": nucleotide},
                          "raw": raw, "rounded": rounded,
                          "limit": export_snapshot.L_OV, "action_paths": []}
                candidates[candidate_key] = record
            elif record["raw"] != raw or record["rounded"] != rounded:
                raise ValueError(f"Conflicting disabled-reference evidence for {key}/{nucleotide}.")
            record["action_paths"].append(f"{recipe['path']}/{action['name']}/alaric.yaml")
    result = []
    def candidate_order(item: tuple[str, str]) -> tuple[int, str]:
        match = _FRAGMENT_KEY.fullmatch(item[0])
        if match is None:
            raise ValueError(f"Malformed disabled reference candidate address: {item[0]!r}.")
        return int(match.group(1)), item[1]

    for key in sorted(candidates, key=candidate_order):
        record = candidates[key]
        record["action_paths"] = sorted(set(record["action_paths"]))
        result.append(record)
    return result


def _manifest(snapshot: dict, recipes: list[dict], digest: str) -> dict:
    if not isinstance(snapshot.get("range_specs"), dict) or not isinstance(snapshot.get("handoff"), dict) \
            or not isinstance(snapshot.get("class_table"), list):
        raise ValueError("The snapshot is missing manifest inputs.")
    gens = []
    for recipe in recipes:
        gens.append({
            "path": recipe["path"],
            "modes": list(recipe["modes"]),
            "range": recipe["range"],
            "direction": recipe["direction"],
            "final_pools": list(recipe["final_pools"]),
            "external_inputs": [dict(row) for row in recipe["external_inputs"]],
            "actions": [action["name"] for action in recipe["actions"]],
        })
    return {
        "schema": 1,
        "snapshot_sha256": digest,
        "reference_sha256": snapshot["reference_sha256"],
        "range_specs": snapshot["range_specs"],
        "handoff": snapshot["handoff"],
        "class_table": snapshot["class_table"],
        "gens": gens,
        "disabled_reference_anchors": _disabled_reference_anchors(snapshot, recipes),
    }


def render_tree(snapshot: Any, recipes: Any, digest: Any) -> dict[str, bytes]:
    """Render a validated snapshot/catalogue to the folder-cell mapping.

    The result uses paths relative to ``alaric_project/``.  It holds source PDB
    buffers unchanged and creates only the JSON/YAML configuration members.
    """
    if type(snapshot) is not dict or type(snapshot.get("schema")) is not int or snapshot.get("schema") != 1:
        raise ValueError("A schema-1 export snapshot is required.")
    if not isinstance(digest, str) or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
        raise ValueError("A lowercase SHA-256 snapshot digest is required.")
    recomputed_digest = export_snapshot.snapshot_digest(snapshot)
    if digest != recomputed_digest:
        raise ValueError("The supplied snapshot digest does not match the normalized snapshot.")
    checked_recipes = _validate_recipes(recipes)

    reference = _checked_bytes(snapshot.get("reference_pdb"), "reference.pdb")
    reference_sha256 = hashlib.sha256(reference).hexdigest()
    if snapshot.get("reference_sha256") != reference_sha256:
        raise ValueError("The reference PDB checksum does not match the normalized snapshot.")

    tree: dict[str, bytes] = {}
    _add_member(tree, "DATA/reference.pdb", reference)
    exclude = snapshot.get("exclude_pdb_code")
    structure = snapshot.get("structure")
    if not isinstance(exclude, str) or re.fullmatch(r"[a-z0-9]{4}", exclude) is None \
            or type(structure) is not dict or structure.get("pdbcode") != exclude:
        raise ValueError("The lower-case export code is malformed or stale.")
    _add_member(tree, "DATA/pdbcode.txt", f"{exclude}\n".encode("utf-8"))

    for number, plain_bytes, all_atom_bytes in _range_sources(snapshot):
        _add_member(tree, f"DATA/range-{number}.pdb", plain_bytes)
        _add_member(tree, f"DATA/range-{number}-aa.pdb", all_atom_bytes)
    all_atom_range_one = snapshot["range_pdb_files_aa"]["range-1-aa.pdb"]
    if reference != all_atom_range_one:
        raise ValueError("DATA/reference.pdb must be the exact prepared range-1 all-atom buffer.")

    _add_member(tree, "DATA/constraints.json", _json_bytes(_constraints(snapshot)))
    _add_member(tree, "DATA/energy-thresholds.json", _json_bytes(_energy_thresholds(snapshot)))
    score_filter = {
        "min_cascade": 1000,
        "cascade": [
            {"plateau_distance": 5, "neighbour_distance": 7, "max_error": 1.0},
            {"plateau_distance": 10, "neighbour_distance": 12, "max_error": 0.05},
        ],
    }
    _add_member(tree, "DATA/score-filter.yaml", _yaml_bytes(score_filter))

    clipped = snapshot.get("clipped_fragments")
    kept_fragments = snapshot.get("kept_fragments")
    if type(clipped) is not dict or not isinstance(kept_fragments, list):
        raise ValueError("The snapshot has no kept clipped-fragment files.")
    expected_clipped: set[str] = {"clip.json", "thresholds.json"}
    for fragment in kept_fragments:
        if type(fragment) is not dict:
            raise ValueError("A kept fragment is malformed.")
        key, number, sequence = fragment.get("key"), fragment.get("number"), fragment.get("sequence")
        match = _FRAGMENT_KEY.fullmatch(key) if isinstance(key, str) else None
        if match is None or type(number) is not int or int(match.group(1)) != number \
                or sequence != match.group(2):
            raise ValueError(f"Invalid kept fragment identity: {key!r}.")
        ref_name = f"reference/frag-{number}-{sequence}-refe.pdb"
        fit_name = f"bestfit/frag-{number}-{sequence}-bestfit.pdb"
        expected_clipped.update((ref_name, fit_name))
        for name in (ref_name, fit_name):
            _add_member(tree, f"DATA/clipped-fragments/{name}", _checked_bytes(clipped.get(name), name))
    if set(clipped) != expected_clipped:
        raise ValueError("The clipped-fragment folder contains missing or unexpected members.")
    _checked_bytes(clipped["clip.json"], "clip.json")
    _checked_bytes(clipped["thresholds.json"], "thresholds.json")

    for recipe in checked_recipes:
        # Recipes are pure configuration; their paths are checked by the
        # planner and again here before becoming archive keys.
        for action in recipe["actions"]:
            action_path = f"{recipe['path']}/{action['name']}/alaric.yaml"
            parameters = action.get("parameters")
            comments = action.get("comments")
            if type(parameters) is not dict or not isinstance(parameters.get("action"), str):
                raise ValueError(f"Malformed action mapping in {recipe['path']}/{action['name']}.")
            _add_member(tree, action_path, _yaml_bytes(parameters, comments))

    manifest = _manifest(snapshot, checked_recipes, digest)
    _add_member(tree, "DATA/export-manifest.json", _json_bytes(manifest))
    _validate_member_tree(tree)
    return {key: tree[key] for key in sorted(tree)}
