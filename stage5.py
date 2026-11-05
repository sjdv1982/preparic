"""Pure helpers for stage-5 energetics, interaction selection and thresholds."""

from __future__ import annotations

import json
import math
import re
from typing import Any


_FRAGMENT_KEY = re.compile(r"^frag-([1-9][0-9]*)-([ACGU]{2})$")
_ENERGY_VALUE = re.compile(r"^[+-]?(?:0|[1-9][0-9]*)(?:\.[0-9]{1,2})?$")


def _fragment_rows(fragments: object) -> list[dict]:
    if type(fragments) is not dict or fragments.get("schema") != 1 \
        or not isinstance(fragments.get("fragments"), list):
        raise ValueError("The canonical fragment list is not ready.")
    rows = fragments["fragments"]
    if not rows or any(type(row) is not dict or not isinstance(row.get("key"), str)
                       or _FRAGMENT_KEY.fullmatch(row["key"]) is None for row in rows):
        raise ValueError("The canonical fragment list contains an invalid row.")
    keys = [row["key"] for row in rows]
    if len(set(keys)) != len(keys):
        raise ValueError("The canonical fragment list contains duplicate keys.")
    return rows


def _protein_rows(range_specs: object) -> list[dict]:
    if type(range_specs) is not dict or range_specs.get("schema") != 1 \
        or not isinstance(range_specs.get("rows"), list) or len(range_specs["rows"]) < 2:
        raise ValueError("At least one RNA and one protein range are required.")
    rows = range_specs["rows"]
    if any(type(row) is not dict or set(row) != {"chain", "start", "end"}
           or not isinstance(row["chain"], str) or not row["chain"]
           or type(row["start"]) is not int or type(row["end"]) is not int
           or row["end"] < row["start"] for row in rows):
        raise ValueError("The current range specification is invalid.")
    return rows[1:]


def energetics_ligands(
    refe_fragments: object,
    best_fit: object,
    fragments: object,
) -> dict[str, bytes]:
    """T12: stage all reference and best-fit PDBs, including clipped rows."""
    rows = _fragment_rows(fragments)
    if type(refe_fragments) is not dict or type(best_fit) is not dict:
        raise ValueError("The reference-fragment and best-fit folders are required.")
    output: dict[str, bytes] = {}
    for row in rows:
        key = row["key"]
        reference_name = f"{key}-refe.pdb"
        bestfit_name = f"pdbs/{key}-bestfit.pdb"
        reference = refe_fragments.get(reference_name)
        bestfit = best_fit.get(bestfit_name)
        if not isinstance(reference, bytes):
            raise ValueError(f"The reference-fragment folder is missing {reference_name}.")
        if not isinstance(bestfit, bytes):
            raise ValueError(f"The best-fit folder is missing {bestfit_name}.")
        output[reference_name] = reference
        output[f"{key}-bestfit.pdb"] = bestfit
    return output


def merge_energetics(
    range_results: object,
    range_specs: object,
    best_fits: object,
    fragments: object,
) -> dict:
    """T14: join range-local T13 outputs and publish the S13 freshness echo."""
    protein_rows = _protein_rows(range_specs)
    fragment_rows = _fragment_rows(fragments)
    if type(best_fits) is not dict or not isinstance(best_fits.get("echo"), dict):
        raise ValueError("The fresh best-fit result is not ready.")
    if type(range_results) is not dict:
        raise ValueError("Per-range energetics results are not ready.")

    ranges = []
    for index, row in enumerate(protein_rows, start=2):
        range_id = f"range-{index}"
        result = range_results.get(range_id)
        if type(result) is not dict or result.get("schema") != 1 \
            or type(result.get("receptor")) is not dict \
            or type(result.get("receptor", {}).get("beads")) is not int \
            or not isinstance(result.get("pairs"), list):
            raise ValueError(f"Energetics for {range_id} is not ready.")
        by_name: dict[str, dict] = {}
        for pair in result["pairs"]:
            if type(pair) is not dict or not isinstance(pair.get("ligand"), str):
                raise ValueError(f"Energetics for {range_id} contains an invalid ligand row.")
            if pair["ligand"] in by_name:
                raise ValueError(f"Energetics for {range_id} contains duplicate {pair['ligand']}.")
            energy = pair.get("energy")
            error = pair.get("error")
            if energy is not None and (type(energy) not in (int, float) or not math.isfinite(energy)):
                raise ValueError(f"Energetics for {range_id} contains a non-finite energy.")
            if error is not None and not isinstance(error, str):
                raise ValueError(f"Energetics for {range_id} contains an invalid error.")
            by_name[pair["ligand"]] = {"energy": energy, "error": error}

        ligand_rows = []
        for fragment in fragment_rows:
            key = fragment["key"]
            reference_name = f"{key}-refe.pdb"
            bestfit_name = f"{key}-bestfit.pdb"
            reference = by_name.get(reference_name)
            bestfit = by_name.get(bestfit_name)
            if reference is None or bestfit is None:
                raise ValueError(f"Energetics for {range_id} is missing one or more rows for {key}.")
            ligand_rows.append({"key": key, "reference": reference, "bestfit": bestfit})
        ranges.append({
            "range": range_id,
            "row": dict(row),
            "receptor": dict(result["receptor"]),
            "rows": ligand_rows,
        })

    return {
        "schema": 1,
        "echo": {
            "best_fits": best_fits["echo"],
            "protein_ranges": [dict(row) for row in protein_rows],
        },
        "ranges": ranges,
    }


def energy_thresholds(energetics: object, best_fits: object, energy_threshold_module: object) -> dict:
    """T16: calculate S15 only from each best-fit float64 energy."""
    if type(energetics) is not dict or energetics.get("schema") != 1 \
        or not isinstance(energetics.get("echo"), dict) or not isinstance(energetics.get("ranges"), list):
        raise ValueError("Fresh energetics are not ready.")
    if type(best_fits) is not dict or not isinstance(best_fits.get("echo"), dict) \
        or best_fits["echo"] != energetics["echo"].get("best_fits"):
        raise ValueError("Energetics do not match the current best-fit geometry.")
    round_energy = getattr(energy_threshold_module, "energy_threshold", None)
    if not callable(round_energy):
        raise ValueError("The pinned energy-threshold module is unavailable.")
    values: dict[str, dict[str, str]] = {}
    for range_row in energetics["ranges"]:
        if type(range_row) is not dict or not isinstance(range_row.get("range"), str) \
            or not isinstance(range_row.get("rows"), list):
            raise ValueError("An energetics range is malformed.")
        range_id = range_row["range"]
        for row in range_row["rows"]:
            if type(row) is not dict or not isinstance(row.get("key"), str):
                raise ValueError("An energetics row is malformed.")
            bestfit = row.get("bestfit")
            energy = bestfit.get("energy") if type(bestfit) is dict else None
            error = bestfit.get("error") if type(bestfit) is dict else "missing best-fit result"
            if error is not None or energy is None:
                continue
            if type(energy) not in (int, float) or not math.isfinite(energy):
                raise ValueError(f"Best-fit energy for {row['key']} / {range_id} is invalid.")
            values.setdefault(row["key"], {})[range_id] = round_energy(float(energy))
    return {"schema": 1, "echo": energetics["echo"], "values": values}


def validate_energy_threshold_edits(
    edits: object,
    fragments: object,
    range_specs: object,
) -> None:
    """T27: validate current S24 rows while tolerating stale basis/range rows."""
    if edits is None:
        return None
    if type(edits) is not dict or set(edits) != {"schema", "basis", "ranges"} \
        or type(edits.get("schema")) is not int or edits["schema"] != 1 \
        or not isinstance(edits.get("ranges"), list):
        raise ValueError("Energy-threshold edits must be schema 1 or null.")
    if type(fragments) is not dict or type(fragments.get("echo")) is not dict:
        raise ValueError("The canonical fragment list is not ready.")
    if edits.get("basis") != fragments["echo"]:
        return None
    fragment_keys = {row["key"] for row in _fragment_rows(fragments)}
    protein_rows = _protein_rows(range_specs)
    current_ranges = {
        json.dumps(row, sort_keys=True, separators=(",", ":"))
        for row in protein_rows
    }
    seen_rows = set()
    for entry in edits["ranges"]:
        if type(entry) is not dict or set(entry) != {"row", "values"} \
            or type(entry.get("row")) is not dict or type(entry.get("values")) is not dict:
            raise ValueError("An energy-threshold range row is malformed.")
        row = entry["row"]
        if set(row) != {"chain", "start", "end"} or not isinstance(row.get("chain"), str) \
            or type(row.get("start")) is not int or type(row.get("end")) is not int \
            or row["end"] < row["start"]:
            raise ValueError("An energy-threshold range address is malformed.")
        signature = json.dumps(row, sort_keys=True, separators=(",", ":"))
        if signature in seen_rows:
            raise ValueError("Energy-threshold edits contain duplicate range rows.")
        seen_rows.add(signature)
        if signature not in current_ranges:
            continue
        for key, value in entry["values"].items():
            if key not in fragment_keys:
                raise ValueError(f"Energy threshold uses unknown fragment {key!r}.")
            if value is not None and (not isinstance(value, str) or _ENERGY_VALUE.fullmatch(value) is None
                                      or not math.isfinite(float(value))):
                raise ValueError(f"Energy threshold for {key} must be a finite number with at most two decimals.")
    return None


def make_interaction_matrix(
    fragments: object,
    clip: object,
    range_specs: object,
    selection: object,
    stackings: object,
) -> dict:
    """Build the empty S14 matrix with one ordered column per protein range."""
    fragment_rows = _fragment_rows(fragments)
    protein_rows = _protein_rows(range_specs)
    if type(clip) is not dict or type(clip.get("clip1")) is not int \
        or type(clip.get("clip2")) is not int or clip["clip1"] < 0 or clip["clip2"] < 0 \
        or clip["clip1"] + clip["clip2"] >= len(fragment_rows):
        raise ValueError("The current clip has an invalid extent.")
    if type(selection) is not dict or type(selection.get("selected")) is not dict \
        or any(type(flag) is not bool for flag in selection["selected"].values()):
        raise ValueError("The current stacking selection is malformed.")
    all_rows = [row["key"] for row in fragment_rows]
    kept_rows = all_rows[clip["clip1"]:len(all_rows) - clip["clip2"]]
    basis = sorted(key for key, selected in selection["selected"].items() if selected)
    anchors = candidate_anchors(fragments, selection, stackings, kept_rows)
    ranges = [f"range-{index}" for index in range(2, len(protein_rows) + 2)]
    if any(anchor["range"] not in ranges for anchor in anchors):
        raise ValueError("A selected anchor refers to a protein range outside the current range table.")
    return {
        "schema": 1,
        "basis": basis,
        "active": True,
        "rows": kept_rows,
        "ranges": ranges,
        "grid": [[False for _ in ranges] for _ in kept_rows],
    }


def reconcile_energy_threshold_edits(
    existing: object,
    fragments: object,
    range_specs: object,
) -> dict:
    """Reconcile S24 by addressed range rows, preserving unchanged entries."""
    fragment_rows = _fragment_rows(fragments)
    protein_rows = _protein_rows(range_specs)
    basis = fragments.get("echo")
    if type(basis) is not dict:
        raise ValueError("The canonical fragment basis is not ready.")
    previous: dict[str, dict] = {}
    if type(existing) is dict and existing.get("schema") == 1 and existing.get("basis") == basis \
        and isinstance(existing.get("ranges"), list):
        for entry in existing["ranges"]:
            if type(entry) is dict and type(entry.get("row")) is dict and type(entry.get("values")) is dict:
                signature = json.dumps(entry["row"], sort_keys=True, separators=(",", ":"))
                previous[signature] = entry["values"]
    keys = [row["key"] for row in fragment_rows]
    reconciled = []
    for row in protein_rows:
        signature = json.dumps(row, sort_keys=True, separators=(",", ":"))
        old_values = previous.get(signature, {})
        values = {key: old_values.get(key) if isinstance(old_values.get(key), (str, type(None))) else None
                  for key in keys}
        reconciled.append({"row": dict(row), "values": values})
    return {"schema": 1, "basis": dict(basis), "ranges": reconciled}


def candidate_anchors(fragments: object, selection: object, stackings: object, matrix_rows: list[str]) -> list[dict]:
    """Resolve selected stacking locks through T24 fragment residue identities."""
    fragment_rows = _fragment_rows(fragments)
    if type(selection) is not dict or selection.get("active") is not True \
        or type(selection.get("selected")) is not dict \
        or any(type(value) is not bool for value in selection["selected"].values()):
        raise ValueError("A current stacking selection is required.")
    stacking_rows = stackings.get("rows") if type(stackings) is dict else stackings
    if not isinstance(stacking_rows, list):
        raise ValueError("Fresh stacking rows are not ready.")
    selected = {key for key, flag in selection["selected"].items() if flag}
    by_resid: dict[int, list[tuple[str, int]]] = {}
    for fragment in fragment_rows:
        for nucleotide, residue in enumerate(fragment.get("residues", []), start=1):
            if type(residue) is dict and type(residue.get("resid")) is int:
                by_resid.setdefault(residue["resid"], []).append((fragment["key"], nucleotide))

    anchors: dict[tuple[int, int], dict] = {}
    for stacking in stacking_rows:
        if type(stacking) is not dict:
            raise ValueError("A stacking row is malformed.")
        try:
            rna_name = str(stacking["rna_resname"]).strip()
            rna_resid = stacking["rna_resid"]
            range_index = stacking["protein_range_pdb_id"]
            protein_resid = stacking["protein_resid"]
        except KeyError as exc:
            raise ValueError(f"A stacking row is missing {exc.args[0]!r}.") from exc
        raw_key = json.dumps([rna_name, rna_resid, range_index, protein_resid], separators=(",", ":"))
        if raw_key not in selected:
            continue
        if type(rna_resid) is not int or type(range_index) is not int or range_index < 2 \
            or type(protein_resid) is not int:
            raise ValueError("A selected stacking has an invalid residue address.")
        identity = (rna_resid, range_index)
        anchor = anchors.setdefault(identity, {
            "resid": rna_resid,
            "range": f"range-{range_index}",
            "stackings": [],
            "locks": [],
        })
        anchor["stackings"].append(raw_key)
        for fragment_key, nucleotide in by_resid.get(rna_resid, []):
            lock = {"fragment": fragment_key, "nucleotide": nucleotide}
            if lock not in anchor["locks"]:
                anchor["locks"].append(lock)

    row_set = set(matrix_rows)
    result = []
    for (resid, range_index), anchor in sorted(anchors.items(), key=lambda item: (item[0][0], item[0][1])):
        locks = [lock for lock in anchor["locks"] if lock["fragment"] in row_set]
        if not locks:
            raise ValueError(f"Selected anchor at RNA residue {resid} in range-{range_index} is clipped or unmapped.")
        result.append({**anchor, "locks": locks})
    return result


def recompute_interaction_matrix(matrix: object, anchors: list[dict], energetics: object) -> dict:
    """P5.4: rebuild one column per protein range while retaining every anchor lock."""
    if type(matrix) is not dict or matrix.get("schema") != 1 \
        or not isinstance(matrix.get("rows"), list) or not isinstance(matrix.get("ranges"), list):
        raise ValueError("The interaction matrix is malformed.")
    if type(energetics) is not dict or energetics.get("schema") != 1 \
        or not isinstance(energetics.get("ranges"), list):
        raise ValueError("Fresh energetics are not ready.")
    rows = matrix["rows"]
    ranges = matrix["ranges"]
    if not rows or not ranges or any(type(range_id) is not str or not range_id for range_id in ranges) \
        or len(set(ranges)) != len(ranges) or len(matrix.get("grid", [])) != len(rows):
        raise ValueError("The interaction matrix has an invalid extent.")
    if any(type(row) is not list or len(row) != len(ranges) or any(type(cell) is not bool for cell in row)
           for row in matrix["grid"]):
        raise ValueError("The interaction matrix has an invalid grid.")
    energy_by_range = _energy_index(energetics)
    row_positions = {key: index for index, key in enumerate(rows)}
    grid = [[False for _ in ranges] for _ in rows]
    if any(anchor.get("range") not in ranges for anchor in anchors):
        raise ValueError("A selected anchor refers to an unknown S14 protein range.")
    locks_by_range: dict[str, set[int]] = {range_id: set() for range_id in ranges}
    positions_by_anchor: list[set[int]] = []
    for anchor in anchors:
        range_id = anchor["range"]
        raw_locks = [lock.get("fragment") for lock in anchor.get("locks", [])]
        if not raw_locks or any(fragment not in row_positions for fragment in raw_locks):
            raise ValueError(f"Selected anchor {anchor.get('resid')} / {range_id} has no kept lock.")
        locks = {row_positions[fragment] for fragment in raw_locks}
        positions_by_anchor.append(locks)
        locks_by_range[range_id].update(locks)

    for anchor, locks in zip(anchors, positions_by_anchor):
        range_id = anchor["range"]
        for position in locks:
            if not _has_energy(energy_by_range.get(range_id, {}).get(rows[position])):
                raise ValueError(f"Selected anchor {anchor['resid']} / {range_id} has a scoring error at {rows[position]}.")
        left, right = min(locks), max(locks)
        other_locks = locks_by_range[range_id] - locks
        while left > 0:
            candidate = left - 1
            if candidate in other_locks or not _grows(energy_by_range.get(range_id, {}).get(rows[candidate])):
                break
            left = candidate
        while right + 1 < len(rows):
            candidate = right + 1
            if candidate in other_locks or not _grows(energy_by_range.get(range_id, {}).get(rows[candidate])):
                break
            right = candidate
        column = ranges.index(range_id)
        for position in range(left, right + 1):
            grid[position][column] = True
    result = dict(matrix)
    result["grid"] = grid
    return result


def _energy_index(energetics: object) -> dict[str, dict[str, tuple[dict, dict]]]:
    if type(energetics) is not dict or energetics.get("schema") != 1 \
        or not isinstance(energetics.get("ranges"), list):
        raise ValueError("Fresh energetics are not ready.")
    output: dict[str, dict[str, tuple[dict, dict]]] = {}
    for range_data in energetics["ranges"]:
        if type(range_data) is not dict or not isinstance(range_data.get("range"), str) \
            or not isinstance(range_data.get("rows"), list):
            raise ValueError("The energetics result contains a malformed range.")
        by_fragment = {}
        for row in range_data["rows"]:
            if type(row) is not dict or not isinstance(row.get("key"), str) \
                or type(row.get("reference")) is not dict or type(row.get("bestfit")) is not dict:
                raise ValueError("The energetics result contains a malformed fragment row.")
            by_fragment[row["key"]] = (row["reference"], row["bestfit"])
        output[range_data["range"]] = by_fragment
    return output


def _has_energy(pair: tuple[dict, dict] | None) -> bool:
    if pair is None:
        return False
    for item in pair:
        value = item.get("energy")
        if item.get("error") or type(value) not in (int, float) or not math.isfinite(value):
            return False
    return True


def _grows(pair: tuple[dict, dict] | None) -> bool:
    if not _has_energy(pair):
        return False
    return any(abs(float(item["energy"])) >= 3 for item in pair)


def validate_interaction_matrix(
    matrix: object,
    rows: list[str],
    ranges: list[str],
    anchors: list[dict],
    energetics: object,
    *,
    require_active: bool = True,
) -> dict:
    """Strict addressed-matrix validation for range columns, anchor locks, and runs."""
    if type(matrix) is not dict or set(matrix) != {"schema", "basis", "active", "rows", "ranges", "grid"}:
        raise ValueError("S14 schema check failed.")
    if type(matrix.get("schema")) is not int or matrix["schema"] != 1:
        raise ValueError("S14 schema check failed.")
    if type(matrix.get("active")) is not bool or (require_active and matrix["active"] is not True):
        raise ValueError("S14 active-state check failed.")
    if not isinstance(matrix.get("basis"), list) or any(not isinstance(key, str) or not key for key in matrix["basis"]) \
        or len(set(matrix["basis"])) != len(matrix["basis"]):
        raise ValueError("S14 selection basis check failed.")
    if matrix.get("rows") != rows:
        raise ValueError("S14 addressed-row order or clip extent check failed.")
    if any(type(range_id) is not str or not range_id for range_id in ranges) \
        or len(set(ranges)) != len(ranges) or matrix.get("ranges") != ranges:
        raise ValueError("S14 ordered protein-range check failed.")
    if any(anchor.get("range") not in ranges for anchor in anchors):
        raise ValueError("S14 selected anchor refers to an unknown protein range.")
    grid = matrix.get("grid")
    if not isinstance(grid, list) or len(grid) != len(rows):
        raise ValueError("S14 row-count check failed.")
    for row in grid:
        if not isinstance(row, list) or len(row) != len(ranges):
            raise ValueError("S14 rectangular-grid check failed.")
        if any(type(value) is not bool for value in row):
            raise ValueError("S14 boolean-cell check failed.")
    energy_by_range = _energy_index(energetics)
    row_position = {key: index for index, key in enumerate(rows)}
    anchors_by_range: dict[str, list[dict]] = {range_id: [] for range_id in ranges}
    for anchor in anchors:
        range_id = anchor.get("range")
        if not isinstance(range_id, str):
            raise ValueError(f"S14 lock check failed for anchor {anchor.get('resid')} in {range_id}.")
        if range_id not in anchors_by_range:
            raise ValueError(f"S14 anchor {anchor.get('resid')} refers to unknown range {range_id}.")
        column = ranges.index(range_id)
        locks_by_column: set[int] = set()
        for lock in anchor.get("locks", []):
            fragment = lock.get("fragment")
            if fragment not in row_position:
                raise ValueError(f"S14 clip removed a selected anchor lock {fragment}.")
            index = row_position[fragment]
            if grid[index][column] is not True:
                raise ValueError(f"S14 lock check failed for anchor {anchor.get('resid')} at {fragment} / {range_id}.")
            locks_by_column.add(index)
        if not locks_by_column:
            raise ValueError(f"S14 anchor {anchor.get('resid')} / {range_id} has no lock.")
        anchors_by_range[range_id].append(anchor)
    for column, range_id in enumerate(ranges):
        lock_positions = {row_position[lock["fragment"]]
                          for anchor in anchors_by_range[range_id]
                          for lock in anchor.get("locks", []) if lock.get("fragment") in row_position}
        index = 0
        while index < len(rows):
            if not grid[index][column]:
                index += 1
                continue
            first = index
            while index + 1 < len(rows) and grid[index + 1][column]:
                index += 1
            last = index
            if not any(first <= lock <= last for lock in lock_positions):
                raise ValueError(f"S14 true run {rows[first]}–{rows[last]} in {range_id} has no selected anchor.")
            for fragment in rows[first:last + 1]:
                if not _has_energy(energy_by_range.get(range_id, {}).get(fragment)):
                    owner = next((anchor for anchor in anchors_by_range[range_id]
                                  if any(lock.get("fragment") == fragment
                                         for lock in anchor.get("locks", []))), None)
                    if owner is not None:
                        raise ValueError(
                            f"Selected anchor {owner['resid']} / {range_id} has a scoring error at {fragment}."
                        )
                    raise ValueError(f"S14 selects {fragment} in {range_id}, which has a scoring error.")
            index += 1
    return matrix


def validate_interaction_matrix_structure(
    matrix: object,
    rows: list[str],
    ranges: list[str],
    *,
    require_active: bool = True,
) -> dict:
    """Validate the addressed S14 envelope without checking anchor admission.

    T15's backend contract ends at shape and address validation. Anchor locks,
    true-run ownership, and scoring-error admission remain frontend read checks
    and are repeated at export.
    """
    if type(matrix) is not dict or set(matrix) != {"schema", "basis", "active", "rows", "ranges", "grid"}:
        raise ValueError("S14 schema check failed.")
    if type(matrix.get("schema")) is not int or matrix["schema"] != 1:
        raise ValueError("S14 schema check failed.")
    if type(matrix.get("active")) is not bool or (require_active and matrix["active"] is not True):
        raise ValueError("S14 active-state check failed.")
    basis = matrix.get("basis")
    if not isinstance(basis, list) or any(not isinstance(key, str) or not key for key in basis) \
        or len(set(basis)) != len(basis):
        raise ValueError("S14 selection basis check failed.")
    if matrix.get("rows") != rows:
        raise ValueError("S14 addressed-row order or clip extent check failed.")
    if any(type(range_id) is not str or not range_id for range_id in ranges) \
        or len(set(ranges)) != len(ranges) or matrix.get("ranges") != ranges:
        raise ValueError("S14 ordered range check failed.")
    grid = matrix.get("grid")
    if not isinstance(grid, list) or len(grid) != len(rows):
        raise ValueError("S14 row-count check failed.")
    for row in grid:
        if type(row) is not list or len(row) != len(ranges):
            raise ValueError("S14 rectangular-grid check failed.")
        if any(type(value) is not bool for value in row):
            raise ValueError("S14 boolean-cell check failed.")
    return matrix


def matrix_to_columns(matrix: object, anchors: list[dict]) -> list[dict]:
    """Serialize the addressed runs once per protein range, retaining all run locks."""
    if type(matrix) is not dict or not isinstance(matrix.get("rows"), list) \
        or not isinstance(matrix.get("ranges"), list) or not isinstance(matrix.get("grid"), list):
        raise ValueError("The interaction matrix is malformed.")
    rows, ranges, grid = matrix["rows"], matrix["ranges"], matrix["grid"]
    if len(grid) != len(rows) or any(type(row) is not list or len(row) != len(ranges)
                                     or any(type(cell) is not bool for cell in row) for row in grid):
        raise ValueError("The interaction matrix is not a rectangular boolean grid.")
    if any(type(range_id) is not str or not range_id for range_id in ranges) \
        or len(set(ranges)) != len(ranges) or any(anchor.get("range") not in ranges for anchor in anchors):
        raise ValueError("S14 columns do not match the ordered protein ranges.")
    columns = []
    for range_id in ranges:
        column = ranges.index(range_id)
        runs = []
        position = 0
        while position < len(rows):
            if not grid[position][column]:
                position += 1
                continue
            first = position
            while position + 1 < len(rows) and grid[position + 1][column]:
                position += 1
            last = position
            run_anchors = []
            for anchor in anchors:
                if anchor.get("range") != range_id:
                    continue
                for lock in anchor.get("locks", []):
                    if first <= rows.index(lock["fragment"]) <= last:
                        item = {"fragment": lock["fragment"], "nucleotide": lock["nucleotide"]}
                        if item not in run_anchors:
                            run_anchors.append(item)
            runs.append({
                "first": rows[first],
                "last": rows[last],
                "rows": rows[first:last + 1],
                "anchors": run_anchors,
            })
            position += 1
        columns.append({"range": range_id, "runs": runs})
    return columns


def columns_to_matrix(
    columns: list[dict], rows: list[str], ranges: list[str], basis: list[str], *, active: bool = True,
) -> dict:
    """Restore the addressed boolean grid from matrix column runs."""
    if not isinstance(columns, list) or len(set(rows)) != len(rows):
        raise ValueError("Addressed matrix rows must be a unique list.")
    if any(type(range_id) is not str or not range_id for range_id in ranges) \
        or len(set(ranges)) != len(ranges) or len(columns) != len(ranges) or any(
        type(column) is not dict or column.get("range") != ranges[index]
        or not isinstance(column.get("runs"), list)
        for index, column in enumerate(columns)
    ):
        raise ValueError("Interaction-matrix columns do not match the ordered ranges.")
    row_position = {key: index for index, key in enumerate(rows)}
    grid = [[False for _ in ranges] for _ in rows]
    for column_index, (column, range_id) in enumerate(zip(columns, ranges)):
        occupied: set[int] = set()
        for run in column["runs"]:
            if type(run) is not dict or not isinstance(run.get("rows"), list) or not run["rows"]:
                raise ValueError(f"A run in {range_id} is malformed.")
            positions = [row_position.get(key, -1) for key in run["rows"]]
            if any(position < 0 for position in positions) or positions != list(range(positions[0], positions[-1] + 1)):
                raise ValueError(f"A run in {range_id} is not a contiguous addressed sequence.")
            if run.get("first") != run["rows"][0] or run.get("last") != run["rows"][-1]:
                raise ValueError(f"A run in {range_id} has inconsistent endpoints.")
            if any(position in occupied for position in positions):
                raise ValueError(f"Runs in {range_id} overlap.")
            occupied.update(positions)
            for position in positions:
                grid[position][column_index] = True
    return {"schema": 1, "basis": list(basis), "active": bool(active),
            "rows": list(rows), "ranges": list(ranges), "grid": grid}


def adapt_interaction_matrix(
    matrix: object,
    new_rows: list[str],
    anchors: list[dict],
    energetics: object,
) -> dict:
    """P5.5: preserve addressed cells across clip changes and grow touched ends."""
    if type(matrix) is not dict or not isinstance(matrix.get("rows"), list) \
        or not isinstance(matrix.get("ranges"), list) or not isinstance(matrix.get("grid"), list):
        raise ValueError("The existing interaction matrix is malformed; recompute it.")
    old_rows, ranges, old_grid = matrix["rows"], matrix["ranges"], matrix["grid"]
    if not isinstance(new_rows, list) or not new_rows or len(set(new_rows)) != len(new_rows):
        raise ValueError("The new clip has an invalid row extent; recompute the matrix.")
    if any(type(range_id) is not str or not range_id for range_id in ranges) \
        or len(set(ranges)) != len(ranges) or len(old_grid) != len(old_rows) or any(type(row) is not list or len(row) != len(ranges)
                                             or any(type(cell) is not bool for cell in row) for row in old_grid):
        raise ValueError("The existing interaction matrix has an invalid grid; recompute it.")
    new_position = {key: index for index, key in enumerate(new_rows)}
    if any(lock["fragment"] not in new_position for anchor in anchors for lock in anchor.get("locks", [])):
        raise ValueError("The new clip removes a selected anchor lock; return to stage 3.")
    if any(type(range_id) is not str or not range_id for range_id in ranges) \
        or len(set(ranges)) != len(ranges) or any(anchor.get("range") not in ranges for anchor in anchors):
        raise ValueError("The selected protein-range columns changed; recompute the matrix.")
    old_positions = [new_position[key] for key in old_rows if key in new_position]
    if not old_positions or old_positions != list(range(old_positions[0], old_positions[-1] + 1)):
        raise ValueError("The new clip changed fragment order; recompute the matrix.")
    energy_by_range = _energy_index(energetics)
    old_position = {key: index for index, key in enumerate(old_rows)}
    grid = [[False for _ in ranges] for _ in new_rows]
    for row_index, key in enumerate(new_rows):
        if key not in old_position:
            continue
        for column in range(len(ranges)):
            grid[row_index][column] = old_grid[old_position[key]][column]

    # Every kept candidate is a mandatory lock, including rows newly exposed
    # by a clip change and candidates whose energy is below the growth cutoff.
    for anchor in anchors:
        column = ranges.index(anchor["range"])
        for lock in anchor.get("locks", []):
            row = new_position[lock["fragment"]]
            # Keep the mandatory lock in the adapted state even when this
            # newly uncovered candidate has a scoring error. The subsequent
            # matrix read rejects it with the anchor-specific diagnostic.
            grid[row][column] = True

    if old_positions:
        first_new, last_new = old_positions[0], old_positions[-1]
        for column, range_id in enumerate(ranges):
            if old_grid[0][column]:
                position = first_new - 1
                while position >= 0 and _grows(energy_by_range.get(range_id, {}).get(new_rows[position])):
                    grid[position][column] = True
                    position -= 1
            if old_grid[-1][column]:
                position = last_new + 1
                while position < len(new_rows) and _grows(energy_by_range.get(range_id, {}).get(new_rows[position])):
                    grid[position][column] = True
                    position += 1
    result = dict(matrix)
    result["rows"] = list(new_rows)
    result["grid"] = grid
    return result


def interaction_handoff(
    matrix: object,
    fragments: object,
    selection: object,
    stackings: object,
    clip: object,
    range_specs: object,
) -> dict | None:
    """T15: validate the addressed matrix and build `interactions.json`."""
    # S14 is not in force during stages 3–4. ST §4.3–4.4 defines an inactive
    # marker as null for its consumers, so T15 must remain a sentinel then.
    if matrix is None or (type(matrix) is dict and matrix.get("active") is False):
        return None
    all_fragments = _fragment_rows(fragments)
    protein_rows = _protein_rows(range_specs)
    if type(clip) is not dict or type(clip.get("clip1")) is not int or type(clip.get("clip2")) is not int:
        raise ValueError("The current clip is invalid.")
    start, stop = clip["clip1"], len(all_fragments) - clip["clip2"]
    kept_rows = all_fragments[start:stop]
    expected_rows = [row["key"] for row in kept_rows]
    anchors = candidate_anchors(fragments, selection, stackings, expected_rows)
    unique_ranges = [f"range-{index}" for index in range(2, len(protein_rows) + 2)]
    expected_ranges = unique_ranges
    validate_interaction_matrix_structure(matrix, expected_rows, expected_ranges)
    grid = matrix["grid"]
    row_positions = {key: index for index, key in enumerate(expected_rows)}
    anchors_by_range: dict[str, list[dict]] = {range_id: [] for range_id in unique_ranges}
    columns_by_range: dict[str, list[int]] = {range_id: [column] for column, range_id in enumerate(unique_ranges)}
    for anchor in anchors:
        if anchor["range"] not in anchors_by_range:
            raise ValueError(f"Selected anchor refers to unknown {anchor['range']}.")
        anchors_by_range[anchor["range"]].append(anchor)

    runs = []
    interactions: dict[str, list[str]] = {}
    for range_id in unique_ranges:
        columns = columns_by_range[range_id]
        position = 0
        while position < len(expected_rows):
            if not any(grid[position][column] for column in columns):
                position += 1
                continue
            first = position
            while position + 1 < len(expected_rows) and any(
                grid[position + 1][column] for column in columns
            ):
                position += 1
            last = position
            run_anchors = []
            for anchor in anchors_by_range[range_id]:
                for lock in anchor["locks"]:
                    lock_position = row_positions[lock["fragment"]]
                    if first <= lock_position <= last:
                        candidate = {"fragment": lock["fragment"], "nucleotide": lock["nucleotide"]}
                        if candidate not in run_anchors:
                            run_anchors.append(candidate)
            for fragment_key in expected_rows[first:last + 1]:
                interactions.setdefault(fragment_key, [])
                if range_id not in interactions[fragment_key]:
                    interactions[fragment_key].append(range_id)
            runs.append({
                "range": range_id,
                "first_fragment": int(_FRAGMENT_KEY.fullmatch(expected_rows[first]).group(1)),
                "last_fragment": int(_FRAGMENT_KEY.fullmatch(expected_rows[last]).group(1)),
                "anchors": run_anchors,
            })
            position += 1

    kept_residues = []
    seen_residues = set()
    for fragment in kept_rows:
        for residue in fragment.get("residues", []):
            key = json.dumps(residue, sort_keys=True, separators=(",", ":"))
            if key not in seen_residues:
                seen_residues.add(key)
                kept_residues.append(dict(residue))
    return {
        "schema": 1,
        "clip1": start,
        "clip2": clip["clip2"],
        "kept_residues": kept_residues,
        "ranges": [
            {"index": index, **dict(row)} for index, row in enumerate(protein_rows, start=2)
        ],
        "runs": runs,
        "interactions": interactions,
    }
