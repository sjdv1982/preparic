"""Pure planner for the export's fragment-pool generator catalogue.

This module describes action graphs only.  It does not read Preparic mounts,
construct a Seamless Context, write files, or serialise YAML.
"""

from __future__ import annotations

import itertools
import json
import math
import re
from collections import defaultdict
from copy import deepcopy
from typing import Any


_POOL_REF_KEYS = ("input", "input1", "input2", "restrict_input", "mask_input")
_FRAG_KEY = re.compile(r"^frag-([1-9][0-9]*)-([ACGU]{2})$")


def _range_number(value: str) -> int:
    match = re.fullmatch(r"range-([1-9][0-9]*)", value)
    if not match:
        raise ValueError(f"Invalid protein range identifier: {value!r}.")
    return int(match.group(1))


def _fragment_rows(snapshot: dict) -> tuple[list[dict], dict[int, dict], dict[str, dict]]:
    rows = snapshot.get("kept_fragments")
    if not isinstance(rows, list):
        rows = snapshot.get("fragments")
    if not isinstance(rows, list) or not rows:
        raise ValueError("The export snapshot has no kept fragment rows.")
    ordered = sorted((dict(row) for row in rows), key=lambda row: row["number"])
    by_number: dict[int, dict] = {}
    by_key: dict[str, dict] = {}
    for row in ordered:
        if type(row.get("number")) is not int or not isinstance(row.get("key"), str) \
                or not isinstance(row.get("sequence"), str) or _FRAG_KEY.fullmatch(row["key"]) is None:
            raise ValueError("A fragment row is malformed.")
        if row["number"] in by_number or row["key"] in by_key:
            raise ValueError(f"The fragment map repeats {row['key']}.")
        if len(row.get("residues", [])) != 2:
            raise ValueError(f"{row['key']} must identify two RNA residues.")
        by_number[row["number"]] = row
        by_key[row["key"]] = row
    return ordered, by_number, by_key


def _key_number(key: str) -> int:
    match = _FRAG_KEY.fullmatch(key)
    if not match:
        raise ValueError(f"Invalid fragment key: {key!r}.")
    return int(match.group(1))


def _range_suffix(ranges_for_fragment: dict[int, set[str]], number: int,
                  mode: str, root_range: str | None) -> str:
    if mode == "unfixed" and root_range is not None \
            and len(ranges_for_fragment.get(number, set())) > 1:
        return f"-range{_range_number(root_range)}"
    return ""


def _final_pool(number: int, ranges_for_fragment: dict[int, set[str]],
                mode: str, root_range: str | None) -> str:
    return f"frag{number}{_range_suffix(ranges_for_fragment, number, mode, root_range)}"


def _intermediate(number: int, role: str, step: str,
                  ranges_for_fragment: dict[int, set[str]], mode: str,
                  root_range: str | None) -> str:
    return f"frag{number}-{role}{_range_suffix(ranges_for_fragment, number, mode, root_range)}-{step}"


def _score_stem(number: int, role: str, ranges_for_fragment: dict[int, set[str]],
                mode: str, root_range: str | None) -> str:
    return f"frag{number}-{role}{_range_suffix(ranges_for_fragment, number, mode, root_range)}"


def _merged_pool(number: int, ranges_for_fragment: dict[int, set[str]],
                 mode: str, root_range: str | None) -> str:
    return f"frag{number}-merged{_range_suffix(ranges_for_fragment, number, mode, root_range)}"


def _class_name(row: dict) -> str:
    rrm = row.get("rrm")
    if type(rrm) is not int:
        raise ValueError(f"Gaussian class row {row.get('key')!r} has no integer rrm.")
    result = f"rrm{rrm}"
    if rrm in (90, 95):
        if row.get("tilt30") is True:
            result += "-tilt30"
        if row.get("disp45") is True:
            result += "-disp45"
    return result


def _model(snapshot: dict) -> dict:
    value = snapshot.get("gaussian_model_json")
    if isinstance(value, bytes):
        value = value.decode("utf-8")
    if isinstance(value, str):
        value = json.loads(value)
    if type(value) is not dict or value.get("schema") != 1 \
            or type(value.get("thresholds")) is not dict:
        raise ValueError("The schema-1 Gaussian model is required by the generator planner.")
    return value


def _stacking_map(snapshot: dict) -> dict[str, dict]:
    result = {}
    for row in snapshot.get("stackings", []):
        if type(row) is not dict or not isinstance(row.get("key"), str):
            raise ValueError("A selected stacking row has no stable key.")
        if row["key"] in result:
            raise ValueError(f"The selected stacking table repeats {row['key']}.")
        result[row["key"]] = row
    return result


def normalise_candidate_anchors(snapshot: dict) -> list[dict]:
    """Return stable candidate anchors, with range and ordinal lock addresses."""
    _, _, by_key = _fragment_rows(snapshot)
    stackings = _stacking_map(snapshot)
    anchors = snapshot.get("candidate_anchors")
    if not isinstance(anchors, list):
        raise ValueError("The export snapshot has no candidate anchor list.")
    result = []
    seen = set()
    for source in anchors:
        if type(source) is not dict or type(source.get("resid")) is not int \
                or not isinstance(source.get("range"), str) \
                or not isinstance(source.get("stackings"), list) \
                or not isinstance(source.get("locks"), list):
            raise ValueError("A candidate anchor record is malformed.")
        _range_number(source["range"])
        if any(not isinstance(key, str) for key in source["stackings"]) \
                or len(source["stackings"]) != len(set(source["stackings"])):
            raise ValueError(f"Candidate anchor {source['resid']}/{source['range']} repeats a stacking key.")
        stack_keys = sorted(source["stackings"])
        if not stack_keys or any(key not in stackings for key in stack_keys):
            raise ValueError(f"Candidate anchor {source['resid']}/{source['range']} has unknown stackings.")
        for key in stack_keys:
            stacking = stackings[key]
            if stacking.get("range") != source["range"] or stacking.get("rna_resid") != source["resid"]:
                raise ValueError(f"Candidate anchor {source['resid']}/{source['range']} disagrees with {key}.")
        locks = []
        seen_locks = set()
        for lock in source["locks"]:
            if type(lock) is not dict or lock.get("fragment") not in by_key \
                    or type(lock.get("nucleotide")) is not int or lock["nucleotide"] not in (1, 2):
                raise ValueError(f"Candidate anchor {source['resid']}/{source['range']} has a malformed lock.")
            fragment = by_key[lock["fragment"]]
            lock_identity = (lock["fragment"], lock["nucleotide"])
            if lock_identity in seen_locks:
                raise ValueError(f"Candidate anchor {source['resid']}/{source['range']} repeats a lock.")
            seen_locks.add(lock_identity)
            ordinal = fragment["number"] + lock["nucleotide"] - 1
            residue = fragment["residues"][lock["nucleotide"] - 1]
            if residue.get("resid") != source["resid"]:
                raise ValueError(f"Candidate anchor {source['resid']} disagrees with {lock['fragment']}.")
            locks.append({"fragment": lock["fragment"], "number": fragment["number"],
                          "nucleotide": lock["nucleotide"], "ordinal": ordinal})
        locks.sort(key=lambda row: (row["number"], row["nucleotide"]))
        identity = (source["resid"], source["range"])
        if identity in seen:
            raise ValueError(f"Candidate anchor {identity!r} is duplicated.")
        seen.add(identity)
        result.append({"resid": source["resid"], "range": source["range"],
                       "stackings": stack_keys, "locks": locks})
    result.sort(key=lambda row: (min(lock["ordinal"] for lock in row["locks"]),
                                 _range_number(row["range"]), row["resid"]))
    return result


def _matrix(snapshot: dict) -> tuple[list[int], list[str], dict[int, dict[str, bool]]]:
    matrix = snapshot.get("interaction_matrix")
    if type(matrix) is not dict:
        raise ValueError("The export snapshot has no interaction matrix.")
    keys = matrix.get("rows")
    ranges = matrix.get("ranges")
    grid = matrix.get("grid")
    _, by_number, by_key = _fragment_rows(snapshot)
    if not isinstance(keys, list) or not isinstance(ranges, list) or not isinstance(grid, list) \
            or len(keys) != len(grid):
        raise ValueError("The interaction matrix is malformed.")
    numbers = []
    cells = {}
    for key, row in zip(keys, grid):
        if key not in by_key or type(row) is not list or len(row) != len(ranges) \
                or any(type(value) is not bool for value in row):
            raise ValueError("The interaction matrix has a malformed row.")
        number = by_key[key]["number"]
        if number in cells:
            raise ValueError(f"The interaction matrix repeats {key}.")
        numbers.append(number)
        cells[number] = dict(zip(ranges, row))
    if len(set(ranges)) != len(ranges):
        raise ValueError("The interaction matrix repeats a protein range.")
    for range_id in ranges:
        _range_number(range_id)
    return numbers, list(ranges), cells


def _reachable(snapshot: dict, anchors: list[dict], ranges: list[str], cells: dict[int, dict[str, bool]]) -> dict[str, set[int]]:
    row_numbers = sorted(cells)
    reachable = {range_id: set() for range_id in ranges}
    for anchor in anchors:
        for lock in anchor["locks"]:
            if anchor["range"] in reachable:
                reachable[anchor["range"]].add(lock["number"])
    # The interaction matrix has already been clipped.  A range walk stops at
    # its first false cell on either side of each anchored fragment.
    for range_id in ranges:
        seeds = list(reachable[range_id])
        for seed in seeds:
            reachable[range_id].add(seed)
            for direction in (-1, 1):
                number = seed + direction
                while number in cells and cells[number].get(range_id) is True:
                    reachable[range_id].add(number)
                    number += direction
    return reachable


def classify_modes(snapshot: dict) -> dict:
    """Compute fixed and range-rooted reachability before generating names."""
    anchors = normalise_candidate_anchors(snapshot)
    numbers, ranges, cells = _matrix(snapshot)
    reachable_by_range = _reachable(snapshot, anchors, ranges, cells)
    fixed = {lock["number"] for anchor in anchors for lock in anchor["locks"]}
    all_interacting = {number for number, row in cells.items() if any(row.values())}
    for seed in list(fixed):
        for direction in (-1, 1):
            number = seed + direction
            while number in all_interacting:
                fixed.add(number)
                number += direction
    ranges_for_fragment = defaultdict(set)
    for range_id, values in reachable_by_range.items():
        for number in values:
            ranges_for_fragment[number].add(range_id)
    return {
        "ranges": ranges,
        "row_numbers": numbers,
        "cells": cells,
        "reachable_fixed": fixed,
        "reachable_unfixed": reachable_by_range,
        "ranges_for_fragment": dict(ranges_for_fragment),
        "anchors": anchors,
    }


def assign_pool_names(snapshot: dict, classification: dict | None = None) -> dict:
    """Assign canonical final pool names for every reachable mode/range."""
    classification = classification or classify_modes(snapshot)
    names = {"fixed": {}, "unfixed": {}}
    ranges_for_fragment = classification["ranges_for_fragment"]
    for number in sorted(classification["reachable_fixed"]):
        names["fixed"][number] = _final_pool(number, ranges_for_fragment, "fixed", None)
    for range_id, values in classification["reachable_unfixed"].items():
        names["unfixed"][range_id] = {
            number: _final_pool(number, ranges_for_fragment, "unfixed", range_id)
            for number in sorted(values)
        }
    return names


def _class_index(snapshot: dict) -> tuple[dict[str, str], dict[str, float]]:
    model = _model(snapshot)
    class_by_key = {}
    for row in snapshot.get("class_table", []):
        if type(row) is not dict or not isinstance(row.get("key"), str):
            raise ValueError("A Gaussian class row is malformed.")
        if row["key"] in class_by_key:
            raise ValueError(f"The Gaussian class table repeats {row['key']}.")
        name = _class_name(row)
        spec = model["thresholds"].get(name)
        if type(spec) is not dict or type(spec.get("cutoff")) not in (int, float) \
                or not math.isfinite(float(spec["cutoff"])):
            raise ValueError(f"The Gaussian model has no finite cutoff for {name}.")
        class_by_key[row["key"]] = name
    return class_by_key, {name: float(spec["cutoff"]) for name, spec in model["thresholds"].items()
                          if type(spec) is dict and type(spec.get("cutoff")) in (int, float)}


def _base_thresholds(snapshot: dict) -> dict[tuple[int, str], tuple[Any, Any]]:
    result = {}
    for row in snapshot.get("threshold_bases", []):
        if type(row) is not dict:
            continue
        number = row.get("fragment")
        sequence = row.get("sequence")
        if type(number) is not int or not isinstance(sequence, str):
            continue
        for side, key in (("first", "nuc1"), ("second", "nuc2")):
            rounded = row.get(f"{key}_ovrmsd")
            raw = row.get(f"{key}_ovrmsd_raw")
            if rounded is not None:
                address = (number, side)
                candidate = (raw, rounded)
                if address in result and result[address] != candidate:
                    raise ValueError(f"Conflicting base threshold evidence for frag{number}/{side}.")
                result[address] = candidate
    return result


def _anchor_runs(snapshot: dict, classification: dict, mode: str,
                 root_range: str | None) -> list[dict]:
    anchors = classification["anchors"]
    by_ordinal: dict[int, list[dict]] = defaultdict(list)
    for anchor in anchors:
        if mode == "unfixed" and anchor["range"] != root_range:
            continue
        for lock in anchor["locks"]:
            by_ordinal[lock["ordinal"]].append({"anchor": anchor, "lock": lock})
    ordinals = sorted(by_ordinal)
    runs = []
    current = []
    for ordinal in ordinals:
        if current and ordinal != current[-1] + 1:
            runs.append(current)
            current = []
        current.append(ordinal)
    if current:
        runs.append(current)

    _, _, by_key = _fragment_rows(snapshot)
    result = []
    reachable = classification["reachable_fixed"] if mode == "fixed" \
        else classification["reachable_unfixed"].get(root_range, set())
    for run_ordinals in runs:
        by_frag: dict[int, list[dict]] = defaultdict(list)
        for ordinal in run_ordinals:
            for item in by_ordinal[ordinal]:
                by_frag[item["lock"]["number"]].append(item)
        frag_numbers = sorted(number for number in by_frag if number in reachable)
        if not frag_numbers:
            continue
        result.append({"ordinals": run_ordinals, "by_frag": by_frag,
                       "fragments": frag_numbers, "left": frag_numbers[0], "right": frag_numbers[-1],
                       "root_range": root_range})
    return result


def _action(name: str, parameters: dict, comments: list[str] | None = None) -> dict:
    return {"name": name, "parameters": parameters,
            "comments": list(comments or [])}


def _pool_fragment_number(pool: str) -> int | None:
    match = re.match(r"^frag([1-9][0-9]*)", pool)
    return int(match.group(1)) if match else None


def _pool_sequence(pool: str, fragment_rows: dict[int, dict]) -> str | None:
    number = _pool_fragment_number(pool)
    row = fragment_rows.get(number)
    return row.get("sequence") if row else None


def _recipe(path: str, mode: str, root_range: str | None, direction: str | None) -> dict:
    return {"path": path, "modes": [mode], "range": root_range,
            "direction": direction, "final_pools": [],
            "external_inputs": [], "actions": []}


def _append_action(recipe: dict, action: dict, output: str, local_outputs: set[str]) -> None:
    if output in local_outputs:
        raise ValueError(f"Generator {recipe['path']} produces local pool {output} twice.")
    local_outputs.add(output)
    recipe["actions"].append(action)


def _score_pool(recipe: dict, snapshot: dict, output_fragment: int, input_pool: str,
                final_output: str, mode: str, root_range: str | None,
                classif: dict, output_stem: str, local_outputs: set[str]) -> str:
    seq_row = next(row for row in _fragment_rows(snapshot)[0] if row["number"] == output_fragment)
    numeric_ranges = [range_id for range_id in classif["ranges"]
                      if classif["cells"].get(output_fragment, {}).get(range_id) is True]
    if mode == "unfixed":
        numeric_ranges = [root_range] if root_range in numeric_ranges else []
    threshold_by_range = {row["range"]: row["threshold"] for row in snapshot.get("energy_cells", [])
                          if row.get("fragment") == seq_row["key"]}
    numeric_ranges.sort(key=lambda range_id: (float(threshold_by_range[range_id]), _range_number(range_id)))
    if not numeric_ranges:
        raise ValueError(f"No scoring range is available for frag{output_fragment} in {mode} mode.")
    prior = input_pool
    for index, range_id in enumerate(numeric_ranges):
        is_last = index == len(numeric_ranges) - 1
        # Tandem-filter intermediates are named from the requested output
        # pool, then receive their scoring-step suffix (-e1, -e2, ...).
        out = final_output if is_last else f"{final_output}-e{index + 1}"
        params = {
            "action": "energy-filter", "input": prior, "fragment": output_fragment,
            "protein": range_id, "threshold": "auto", "epsilon": 15,
            "cascade": "auto", "min_cascade": "auto", "sequence": "auto",
            "exclude": "auto",
        }
        _append_action(recipe, _action(out, params), out, local_outputs)
        prior = out
    return final_output


def _stacking_order(stack_keys: list[str], stacking_map: dict[str, dict]) -> list[dict]:
    rows = [stacking_map[key] for key in stack_keys]
    return sorted(rows, key=lambda row: (_range_number(row["range"]), row["protein_resid"], row["key"]))


def _raw_anchor_pools(recipe: dict, snapshot: dict, lock_items: list[dict], fragment_number: int,
                      raw_step: str, mode: str, root_range: str | None, class_names: dict[str, str],
                      cutoffs: dict[str, float], stackings: dict[str, dict], local_outputs: set[str],
                      reference: bool = False) -> tuple[str, int, list[str]]:
    """Emit one raw pool per Gaussian stacking (or deduped reference candidate)."""
    ranges_for_fragment = classify_modes(snapshot)["ranges_for_fragment"]
    base = _intermediate(fragment_number, "anchor", raw_step, ranges_for_fragment, mode, root_range)
    entries = []
    for item in lock_items:
        anchor = item["anchor"]
        lock = item["lock"]
        for stacking_key in anchor["stackings"]:
            entries.append((anchor, lock, stacking_key))
    if reference:
        dedup = {}
        for anchor, lock, key in entries:
            address = (fragment_number, lock["nucleotide"])
            dedup.setdefault(address, (anchor, lock, key))
        entries = [dedup[key] for key in sorted(dedup)]
    else:
        entries.sort(key=lambda entry: (_range_number(entry[0]["range"]),
                                        stackings[entry[2]]["protein_resid"], entry[2]))

    # The pool is described by its candidate identity.  Repeated Gaussian
    # stackings are explicitly folded; reference candidates are already unique.
    raw_names = []
    if reference:
        for anchor, lock, key in entries:
            raw_name = base
            raw_names.append(raw_name)
            threshold_map = _base_thresholds(snapshot)
            threshold = threshold_map.get((fragment_number,
                                           "first" if lock["nucleotide"] == 1 else "second"))
            if threshold is None:
                raise ValueError(f"No base threshold is recorded for frag{fragment_number}/{lock['nucleotide']}.")
            raw_value, rounded = threshold
            if type(rounded) not in (int, float, str) or isinstance(rounded, bool):
                raise ValueError(f"Invalid reference threshold for frag{fragment_number}.")
            try:
                rounded_number = float(rounded)
            except (TypeError, ValueError, OverflowError) as exc:
                raise ValueError(f"Invalid reference threshold for frag{fragment_number}.") from exc
            if not math.isfinite(rounded_number):
                raise ValueError(f"Invalid reference threshold for frag{fragment_number}.")
            l_ov = math.sqrt(3.0)
            comments = []
            action_threshold: Any = rounded
            if rounded_number > l_ov:
                action_threshold = -1
                raw_text = str(raw_value)
                rounded_text = str(rounded)
                explanation = (f"Computed base ovRMSD is {raw_text} Å raw ({rounded_text} Å rounded); "
                               f"Alaric's supported ovRMSD limit is L_ov = {l_ov:.8f} Å. "
                               "The required rounded positive-control threshold exceeds the range Alaric "
                               "can search without misses. This action is kept on purpose as a disabled "
                               "positive control.")
                if type(raw_value) in (int, float) and float(raw_value) > l_ov:
                    explanation += " The best-fit pose itself lies beyond that limit."
                comments = [explanation]
            params = {
                "action": "anchor-refe", "fragment": fragment_number,
                "nucleotide": "first" if lock["nucleotide"] == 1 else "second",
                "ovrmsd": action_threshold, "reference": "reference.pdb",
                "sequence": "auto", "exclude": "auto",
            }
            _append_action(recipe, _action(raw_name, params, comments), raw_name, local_outputs)
    else:
        by_stacking_count = len(entries)
        if by_stacking_count == 0:
            raise ValueError(f"No anchor actions can produce {base}.")
        for index, (anchor, lock, key) in enumerate(entries, start=1):
            raw_name = base if by_stacking_count == 1 else f"{base}-{index}"
            raw_names.append(raw_name)
            class_name = class_names.get(key)
            if class_name is None or key not in stackings:
                raise ValueError(f"No Gaussian class exists for selected stacking {key}.")
            stack = stackings[key]
            params = {
                "action": "anchor", "fragment": fragment_number,
                "protein": stack["range"], "resid": stack["protein_resid"],
                "nucleotide": "first" if lock["nucleotide"] == 1 else "second",
                "precomputed_filter_name": "gaussian-rrm",
                "precomputed_filter_threshold": class_name,
                "sequence": "auto", "exclude": "auto",
            }
            _append_action(recipe, _action(raw_name, params), raw_name, local_outputs)
            recipe.setdefault("_raw_anchor_order", {})[raw_name] = (
                _range_number(stack["range"]), stack["protein_resid"], key,
                fragment_number, lock["nucleotide"],
            )

    if len(raw_names) == 1:
        return raw_names[0], 1, raw_names
    # Deterministic binary fold: earlier folds use numbered temporaries; the
    # unsuffixed raw/from pool is the last fold's output.
    current = raw_names[0]
    for index, next_name in enumerate(raw_names[1:], start=1):
        fold_step = f"ident-{index}" if raw_step == "raw" else f"{raw_step}-ident-{index}"
        output = _intermediate(fragment_number, "anchor", fold_step,
                               ranges_for_fragment, mode, root_range) \
            if index < len(raw_names) - 1 else base
        params = {"action": "identity", "input1": current, "input2": next_name}
        _append_action(recipe, _action(output, params), output, local_outputs)
        current = output
    return base, len(raw_names), raw_names


def _intersect(recipe: dict, first: str, second: str, output: str,
               local_outputs: set[str]) -> None:
    params = {"action": "identity", "input1": first, "input2": second}
    _append_action(recipe, _action(output, params), output, local_outputs)


def _grow_params(snapshot: dict, fragment: int, direction: str, input_pool: str,
                 restrict_pool: str | None = None) -> dict:
    params = {
        "action": "grow", "fragment": fragment, "direction": direction,
        "input": input_pool, "sequence": "auto", "exclude": "auto",
        "crmsd": "auto", "ovrmsd": "auto",
    }
    if restrict_pool is not None:
        params["restrict_input"] = restrict_pool
    return params


def _gen_path(kind: str, direction: str | None, mode: str, root_range: str | None,
              number_or_bounds: tuple[int, int] | int, *, has_single_range: bool,
              branch_ranges: set[str], source_tree: str | None = None) -> str:
    parts = [kind]
    if direction is not None:
        parts.append(direction)
    if not has_single_range:
        parts.append("fixed" if mode == "fixed" else "unfixed")
    if mode == "unfixed" and not has_single_range and len(branch_ranges) > 1:
        parts.append(root_range)
    if source_tree is not None:
        parts.append(source_tree)
    if isinstance(number_or_bounds, tuple):
        label = f"gen-frag{number_or_bounds[0]}-{number_or_bounds[1]}"
    else:
        label = f"gen-frag{number_or_bounds}"
    parts.append(label)
    return "/".join(parts)


def _anchor_recipes(snapshot: dict, classification: dict, mode: str, root_range: str | None,
                    tree: str, branch_ranges: dict | set[str]) -> list[dict]:
    runs = _anchor_runs(snapshot, classification, mode, root_range)
    class_names, cutoffs = _class_index(snapshot)
    stackings = _stacking_map(snapshot)
    range_for_frag = classification["ranges_for_fragment"]
    one_range = len(classification["ranges"]) == 1
    recipes = []
    for run in runs:
        run_ordinals = run["ordinals"]
        run_fragments = run["fragments"]
        if len(run_ordinals) == 1 and len(run_fragments) == 1:
            directions: list[str | None] = [None]
        else:
            directions = ["fwd", "bwd"]
        span = (run_fragments[0], run_fragments[-1])
        is_single_folder = directions == [None]
        for direction in directions:
            folder_label: tuple[int, int] | int = span if not is_single_folder else span[0]
            branch = branch_ranges.get((tree, direction), set()) if isinstance(branch_ranges, dict) \
                else branch_ranges
            path = _gen_path(tree, direction, mode, root_range, folder_label,
                             has_single_range=one_range, branch_ranges=branch)
            recipe = _recipe(path, mode, root_range, direction)
            local_outputs: set[str] = set()
            pool_by_fragment: dict[int, str] = {}
            raw_count_by_fragment: dict[int, int] = {}
            # Pool source description: inner fragments for K>=2 intersect
            # the two adjacent nucleotide candidates; outer fragments use
            # the one candidate that anchors them.
            k = len(run_ordinals)
            if k == 1:
                for fragment_number in run_fragments:
                    items = [item for ordinal in run_ordinals for item in run["by_frag"].get(fragment_number, [])
                             if item["lock"]["ordinal"] == ordinal]
                    target = _intermediate(fragment_number, "anchor", "raw", range_for_frag, mode, root_range)
                    raw, count, _ = _raw_anchor_pools(
                        recipe, snapshot, items, fragment_number, "raw", mode, root_range,
                        class_names, cutoffs, stackings, local_outputs, reference=(tree == "anchor-refe"),
                    )
                    raw_count_by_fragment[fragment_number] = 1 if tree == "anchor-refe" else count
                    if raw != target:
                        # Reference dedup can still return the canonical name;
                        # this guard keeps the output pool identity explicit.
                        raise ValueError(f"Anchor raw pool name {raw} does not match {target}.")
                    scored_name = _intermediate(fragment_number, "anchor", "scored", range_for_frag, mode, root_range)
                    final = _final_pool(fragment_number, range_for_frag, mode, root_range)
                    _score_pool(recipe, snapshot, fragment_number, raw, final if len(run_fragments) == 1 else scored_name,
                                mode, root_range, classification,
                                _score_stem(fragment_number, "anchor", range_for_frag, mode, root_range),
                                local_outputs)
                    pool_by_fragment[fragment_number] = final if len(run_fragments) == 1 else scored_name
                    if len(run_fragments) == 2:
                        # One anchored nucleotide interior: score both anchor
                        # pools, then construct the opposite direction below.
                        pool_by_fragment[fragment_number] = scored_name
                if len(run_fragments) == 2:
                    low, high = run_fragments
                    sequence = [high, low] if direction == "fwd" else [low, high]
                    current = pool_by_fragment[low if direction == "fwd" else high]
                    for target in sequence:
                        out = _final_pool(target, range_for_frag, mode, root_range)
                        target_scored = pool_by_fragment[target]
                        grow_direction = "fwd" if target > _pool_fragment_number(current) else "bwd"
                        grown = _intermediate(target, "anchor", "grown", range_for_frag, mode, root_range)
                        # R12 turns a single-stacking anchor restriction into
                        # grow + identity.  Combined anchor pools stay on the
                        # cheaper restrict_input path.
                        if raw_count_by_fragment[target] == 1:
                            _append_action(recipe, _action(grown, _grow_params(snapshot, target,
                                                                             grow_direction, current)),
                                           grown, local_outputs)
                            _intersect(recipe, grown, target_scored, out, local_outputs)
                        else:
                            params = _grow_params(snapshot, target, grow_direction, current, target_scored)
                            _append_action(recipe, _action(out, params), out, local_outputs)
                        pool_by_fragment[target] = out
                        current = out
            else:
                # Build raw candidate pools for every anchor nucleotide.  The
                # group-by-fragment records make the inner intersection and
                # outer anchoring explicit and deterministic.
                raw_for_ordinal_fragment: dict[tuple[int, int], tuple[str, int]] = {}
                inner_fragments = [number for number in run_fragments
                                   if number in set(run_ordinals[:-1])]
                for fragment_number in run_fragments:
                    items = run["by_frag"].get(fragment_number, [])
                    # One candidate endpoint for an outer fragment; two for
                    # an inner fragment.  Naming uses the anchored nucleotide
                    # ordinal (`from-N`), not its local first/second index.
                    ords = sorted({item["lock"]["ordinal"] for item in items})
                    for ordinal in ords:
                        candidate_items = [item for item in items if item["lock"]["ordinal"] == ordinal]
                        if fragment_number in inner_fragments:
                            step = f"from-{ordinal}"
                        else:
                            step = "raw"
                        raw, count, _ = _raw_anchor_pools(
                            recipe, snapshot, candidate_items, fragment_number, step, mode, root_range,
                            class_names, cutoffs, stackings, local_outputs, reference=(tree == "anchor-refe"),
                        )
                        raw_for_ordinal_fragment[(fragment_number, ordinal)] = (raw, 1 if tree == "anchor-refe" else count)
                # Score the conjunction of the two anchored nucleotide pools
                # for each inner fragment.
                scored_inner = {}
                for fragment_number in inner_fragments:
                    ordinal_sources = sorted(ordinal for (number, ordinal) in raw_for_ordinal_fragment
                                             if number == fragment_number)
                    if len(ordinal_sources) != 2:
                        raise ValueError(f"Inner fragment {fragment_number} must have two anchor candidates.")
                    raw_intersection = _intermediate(fragment_number, "anchor", "ident", range_for_frag,
                                                     mode, root_range)
                    first = raw_for_ordinal_fragment[(fragment_number, ordinal_sources[0])][0]
                    second = raw_for_ordinal_fragment[(fragment_number, ordinal_sources[1])][0]
                    _intersect(recipe, first, second, raw_intersection, local_outputs)
                    scored = _intermediate(fragment_number, "anchor", "scored", range_for_frag, mode, root_range)
                    if k == 2:
                        scored = _final_pool(fragment_number, range_for_frag, mode, root_range)
                    _score_pool(recipe, snapshot, fragment_number, raw_intersection, scored, mode, root_range,
                                classification,
                                _score_stem(fragment_number, "anchor", range_for_frag, mode, root_range),
                                local_outputs)
                    scored_inner[fragment_number] = scored
                # Construct the first and reverse pruning passes.  With K=2
                # the sole inner pool is already final and there are no passes.
                final_inner = dict(scored_inner)
                if k > 2:
                    ordered_inner = sorted(inner_fragments)
                    if direction == "fwd":
                        first_order = ordered_inner
                    else:
                        first_order = list(reversed(ordered_inner))
                    first_pass = {first_order[0]: scored_inner[first_order[0]]}
                    for index, number in enumerate(first_order[1:], start=1):
                        prior_number = first_order[index - 1]
                        target_final = number == first_order[-1]
                        output = _final_pool(number, range_for_frag, mode, root_range) if target_final \
                            else _intermediate(number, "anchor", "pass", range_for_frag, mode, root_range)
                        grow_direction = "fwd" if direction == "fwd" else "bwd"
                        grow_name = _intermediate(number, "anchor", "grow-pass", range_for_frag, mode, root_range)
                        grow_params = _grow_params(snapshot, number, grow_direction,
                                                   first_pass[prior_number], scored_inner[number])
                        _append_action(recipe, _action(output, grow_params), output, local_outputs)
                        first_pass[number] = output
                        final_inner[number] = output
                    # The end of the first pass is final.  Walk back through
                    # every other inner fragment, restricting against the
                    # previous pass pool (or the scored pool at the edge).
                    reverse_order = list(reversed(first_order[:-1]))
                    if direction == "bwd":
                        # First pass went downwards; reverse pass goes up.
                        reverse_order = list(reversed(first_order[:-1]))
                        back_direction = "fwd"
                    else:
                        back_direction = "bwd"
                    current_number = first_order[-1]
                    for number in reverse_order:
                        output = _final_pool(number, range_for_frag, mode, root_range)
                        restrict_pool = first_pass[number]
                        grow_params = _grow_params(snapshot, number, back_direction,
                                                   final_inner[current_number], restrict_pool)
                        _append_action(recipe, _action(output, grow_params), output, local_outputs)
                        final_inner[number] = output
                        current_number = number
                    final_inner[first_order[-1]] = first_pass[first_order[-1]]
                # Outer fragments are grown from final inner pools, joined to
                # their raw anchor candidates, then scored.
                low_anchor = run_ordinals[0]
                high_anchor = run_ordinals[-1]
                for fragment_number in run_fragments:
                    if fragment_number in inner_fragments:
                        continue
                    if fragment_number < min(inner_fragments):
                        neighbor = min(inner_fragments)
                        ordinal = low_anchor
                        grow_direction = "bwd"
                    else:
                        neighbor = max(inner_fragments)
                        ordinal = high_anchor
                        grow_direction = "fwd"
                    raw_pool, constraint_count = raw_for_ordinal_fragment[(fragment_number, ordinal)]
                    grown = _intermediate(fragment_number, "anchor", "grown", range_for_frag, mode, root_range)
                    final = _final_pool(fragment_number, range_for_frag, mode, root_range)
                    grow_input = final_inner[neighbor]
                    if constraint_count == 1:
                        _append_action(recipe, _action(grown, _grow_params(snapshot, fragment_number,
                                                                           grow_direction, grow_input)),
                                       grown, local_outputs)
                        ident_name = _intermediate(fragment_number, "anchor", "ident", range_for_frag,
                                                   mode, root_range)
                        _intersect(recipe, grown, raw_pool, ident_name, local_outputs)
                        _score_pool(recipe, snapshot, fragment_number, ident_name, final, mode, root_range,
                                    classification,
                                    _score_stem(fragment_number, "anchor", range_for_frag, mode, root_range),
                                    local_outputs)
                    else:
                        params = _grow_params(snapshot, fragment_number, grow_direction, grow_input, raw_pool)
                        scored = _intermediate(fragment_number, "anchor", "scored", range_for_frag, mode, root_range)
                        _append_action(recipe, _action(scored, params), scored, local_outputs)
                        _score_pool(recipe, snapshot, fragment_number, scored, final, mode, root_range,
                                    classification,
                                    _score_stem(fragment_number, "anchor", range_for_frag, mode, root_range),
                                    local_outputs)
                    pool_by_fragment[fragment_number] = final
                pool_by_fragment.update(final_inner)

            recipe["final_pools"] = [pool_by_fragment[number] for number in run_fragments]
            if tree == "anchor":
                _order_gaussian_anchor_actions(recipe)
            _finalize_recipe(recipe, snapshot, local_outputs)
            recipes.append(recipe)
    return recipes


def _finalize_recipe(recipe: dict, snapshot: dict, local_outputs: set[str]) -> None:
    """Attach external pool contracts and validate local dependency ordering."""
    _, by_number, by_key = _fragment_rows(snapshot)
    ranges_for_fragment = classify_modes(snapshot)["ranges_for_fragment"]
    local = set(local_outputs)
    external = {}
    produced_before = set()
    for action in recipe["actions"]:
        for key in _POOL_REF_KEYS:
            value = action["parameters"].get(key)
            if not isinstance(value, str):
                continue
            if value in local and value not in produced_before:
                # A local dependency must appear before its consumer.
                producer_index = next((i for i, candidate in enumerate(recipe["actions"])
                                       if candidate["name"] == value), None)
                current_index = recipe["actions"].index(action)
                if producer_index is None or producer_index >= current_index:
                    raise ValueError(f"Recipe {recipe['path']} has an unresolved local dependency {value}.")
            elif value not in local:
                number = _pool_fragment_number(value)
                row = by_number.get(number)
                if row is None:
                    raise ValueError(f"Recipe {recipe['path']} references unknown fragment pool {value}.")
                pool_range = re.search(r"-range([1-9][0-9]*)(?:-|$)", value)
                if pool_range:
                    frame = f"range-{pool_range.group(1)}"
                else:
                    possible_frames = ranges_for_fragment.get(number, set())
                    frame = next(iter(possible_frames)) if len(possible_frames) == 1 else "fixed"
                external[value] = {
                    "pool": value, "fragment": row["key"], "sequence": row["sequence"],
                    "frame": frame,
                    "producers": [],
                }
        produced_before.add(action["name"])
    recipe["external_inputs"] = sorted(external.values(), key=lambda row: row["pool"])
    validate_recipe(recipe)


def _order_gaussian_anchor_actions(recipe: dict) -> None:
    order = recipe.pop("_raw_anchor_order", {})
    if not order:
        return
    raw = [action for action in recipe["actions"] if action["name"] in order]
    other = [action for action in recipe["actions"] if action["name"] not in order]
    raw.sort(key=lambda action: order[action["name"]])
    recipe["actions"] = raw + other


def validate_recipe(recipe: dict) -> None:
    """Reject unsafe paths, duplicate actions, broken local DAGs and finals."""
    path = recipe.get("path")
    if not isinstance(path, str) or not path or path.startswith("/") or "\\" in path \
            or any(part in ("", ".", "..") for part in path.split("/")):
        raise ValueError(f"Generator path is unsafe: {path!r}.")
    actions = recipe.get("actions")
    if not isinstance(actions, list) or not actions:
        raise ValueError(f"Generator {path} has no actions.")
    names = set()
    if any(type(action) is not dict or not isinstance(action.get("name"), str)
           for action in actions):
        raise ValueError(f"Generator {path} has duplicate or unnamed actions.")
    local = {action["name"] for action in actions}
    if len(local) != len(actions):
        raise ValueError(f"Generator {path} has duplicate or unnamed actions.")
    if any(name in ("", ".", "..") or "/" in name or "\\" in name for name in local):
        raise ValueError(f"Generator {path} has an unsafe action directory name.")
    external_rows = recipe.get("external_inputs")
    if not isinstance(external_rows, list):
        raise ValueError(f"Generator {path} has no external-input list.")
    external_names = set()
    for entry in external_rows:
        if type(entry) is not dict or not isinstance(entry.get("pool"), str) \
                or not isinstance(entry.get("fragment"), str) \
                or not isinstance(entry.get("sequence"), str) \
                or not isinstance(entry.get("frame"), str) or not entry["frame"] \
                or not isinstance(entry.get("producers"), list):
            raise ValueError(f"Generator {path} has a malformed external-input contract.")
        if "/" in entry["pool"] or "\\" in entry["pool"]:
            raise ValueError(f"Generator {path} has an unsafe external pool name.")
        if entry["pool"] in external_names:
            raise ValueError(f"Generator {path} repeats external pool {entry['pool']}.")
        external_names.add(entry["pool"])
    used_external = set()
    available = set()
    for action in actions:
        if type(action) is not dict or type(action.get("parameters")) is not dict \
                or not isinstance(action.get("comments"), list):
            raise ValueError(f"Generator {path} has a malformed action record.")
        name = action["name"]
        if name in names:
            raise ValueError(f"Generator {path} repeats action {name}.")
        names.add(name)
        if not isinstance(action["parameters"].get("action"), str):
            raise ValueError(f"Action {name} has no supported action type.")
        for key in _POOL_REF_KEYS:
            value = action["parameters"].get(key)
            if not isinstance(value, str):
                continue
            if value in local and value not in available:
                raise ValueError(f"Action {name} has forward or cyclic dependency on {value}.")
            if value not in local and value not in external_names:
                raise ValueError(f"Action {name} uses undeclared external pool {value}.")
            if value not in local:
                used_external.add(value)
        available.add(name)
    finals = recipe.get("final_pools")
    if not isinstance(finals, list) or not finals or any(not isinstance(name, str) for name in finals) \
            or len(finals) != len(set(finals)) or any(name not in local for name in finals):
        raise ValueError(f"Generator {path} has an undeclared final pool.")
    if used_external != external_names:
        raise ValueError(f"Generator {path} declares unused external inputs.")


def _ordinary_grow_recipes(snapshot: dict, classification: dict, mode: str,
                           root_range: str | None, runs: list[dict], branch_ranges: set[str]) -> list[dict]:
    ranges_for_fragment = classification["ranges_for_fragment"]
    reachable = classification["reachable_fixed"] if mode == "fixed" \
        else classification["reachable_unfixed"].get(root_range, set())
    allowed_ranges = classification["ranges"] if mode == "fixed" else [root_range]
    rows, by_number, _ = _fragment_rows(snapshot)
    runs = sorted(runs, key=lambda row: row["left"])
    result = []
    one_range = len(classification["ranges"]) == 1
    for run_index, run in enumerate(runs):
        # Grow to the right until the next anchor run's left pool, or a matrix
        # gap/clip boundary.  The opposite direction is built symmetrically.
        right_stop = runs[run_index + 1]["left"] if run_index + 1 < len(runs) else None
        prior_pool = _final_pool(run["right"], ranges_for_fragment, mode, root_range)
        target = run["right"] + 1
        while target in reachable and (right_stop is None or target < right_stop):
            output = _final_pool(target, ranges_for_fragment, mode, root_range)
            grown = _intermediate(target, "fwd", "grown", ranges_for_fragment, mode, root_range)
            unique = _intermediate(target, "fwd", "unique", ranges_for_fragment, mode, root_range)
            score_site = _intermediate(target, "fwd", "score", ranges_for_fragment, mode, root_range)
            path = _gen_path("grow", "fwd", mode, root_range, target,
                             has_single_range=one_range,
                             branch_ranges=(branch_ranges.get(("grow", "fwd"), set())
                                            if isinstance(branch_ranges, dict) else branch_ranges))
            recipe = _recipe(path, mode, root_range, "fwd")
            local_outputs = set()
            params = _grow_params(snapshot, target, "fwd", prior_pool)
            _append_action(recipe, _action(grown, params), grown, local_outputs)
            _append_action(recipe, _action(unique, {"action": "identity", "input1": grown}), unique, local_outputs)
            _score_pool(recipe, snapshot, target, unique, output, mode, root_range,
                        classification, _score_stem(target, "fwd", ranges_for_fragment, mode, root_range),
                        local_outputs)
            recipe["final_pools"] = [output]
            _finalize_recipe(recipe, snapshot, local_outputs)
            result.append(recipe)
            prior_pool = output
            target += 1

        left_stop = runs[run_index - 1]["right"] if run_index > 0 else None
        prior_pool = _final_pool(run["left"], ranges_for_fragment, mode, root_range)
        target = run["left"] - 1
        while target in reachable and (left_stop is None or target > left_stop):
            output = _final_pool(target, ranges_for_fragment, mode, root_range)
            grown = _intermediate(target, "bwd", "grown", ranges_for_fragment, mode, root_range)
            unique = _intermediate(target, "bwd", "unique", ranges_for_fragment, mode, root_range)
            score_site = _intermediate(target, "bwd", "score", ranges_for_fragment, mode, root_range)
            path = _gen_path("grow", "bwd", mode, root_range, target,
                             has_single_range=one_range,
                             branch_ranges=(branch_ranges.get(("grow", "bwd"), set())
                                            if isinstance(branch_ranges, dict) else branch_ranges))
            recipe = _recipe(path, mode, root_range, "bwd")
            local_outputs = set()
            params = _grow_params(snapshot, target, "bwd", prior_pool)
            _append_action(recipe, _action(grown, params), grown, local_outputs)
            _append_action(recipe, _action(unique, {"action": "identity", "input1": grown}), unique, local_outputs)
            _score_pool(recipe, snapshot, target, unique, output, mode, root_range,
                        classification, _score_stem(target, "bwd", ranges_for_fragment, mode, root_range),
                        local_outputs)
            recipe["final_pools"] = [output]
            _finalize_recipe(recipe, snapshot, local_outputs)
            result.append(recipe)
            prior_pool = output
            target -= 1
    return result


def _anchor_constraint_count(snapshot: dict, run: dict, fragment: int,
                             mode: str, root_range: str | None, tree: str) -> int:
    candidates = []
    for item in run["by_frag"].get(fragment, []):
        for key in item["anchor"]["stackings"]:
            candidates.append(key)
    if tree == "anchor-refe":
        return len({(fragment, item["lock"]["nucleotide"]) for item in run["by_frag"].get(fragment, [])})
    return len(candidates)


def _merge_recipes(snapshot: dict, classification: dict, mode: str, root_range: str | None,
                   runs: list[dict], branch_ranges: set[str]) -> list[dict]:
    if len(runs) < 2:
        return []
    one_range = len(classification["ranges"]) == 1
    range_for_frag = classification["ranges_for_fragment"]
    cells = classification["cells"]
    result = []
    runs = sorted(runs, key=lambda row: row["left"])
    for left_run, right_run in zip(runs, runs[1:]):
        if left_run["right"] >= right_run["left"]:
            continue
        gap_numbers = range(left_run["right"], right_run["left"] + 1)
        if any(number not in cells or not (any(cells[number].values()) if mode == "fixed"
                                             else cells[number].get(root_range) is True)
               for number in gap_numbers):
            continue
        # In unfixed mode all runs already share root_range.  Fixed runs are
        # compatible by definition because all range PDBs share a frame.
        for direction in ("fwd", "bwd"):
            if direction == "fwd":
                targets = list(range(left_run["right"] + 1, right_run["left"] + 1))
                for target in targets:
                    # Every merge is an independent choice from the plain
                    # directional pools.  It never consumes an earlier
                    # `-merged` result or redirects the next edge.
                    input_pool = _final_pool(target - 1, range_for_frag, mode, root_range)
                    restrict_pool = _final_pool(target, range_for_frag, mode, root_range)
                    # For non-boundary restrictions, the opposite direction's
                    # ordinary grow pool has the same canonical name.
                    left_is_anchor = target == right_run["left"]
                    gaussian_count = _anchor_constraint_count(
                        snapshot, right_run, target, mode, root_range, "anchor") if left_is_anchor else 0
                    reference_count = _anchor_constraint_count(
                        snapshot, right_run, target, mode, root_range, "anchor-refe") if left_is_anchor else 0
                    if left_is_anchor and (gaussian_count == 1) != (reference_count == 1):
                        variants = [("anchor", gaussian_count), ("anchor-refe", reference_count)]
                    else:
                        variants = [(None, gaussian_count)]
                    for source_tree, count in variants:
                        path = _gen_path("merge", direction, mode, root_range, target,
                                         has_single_range=one_range,
                                         branch_ranges=(branch_ranges.get(("merge", direction), set())
                                                        if isinstance(branch_ranges, dict) else branch_ranges),
                                         source_tree=source_tree)
                        recipe = _recipe(path, mode, root_range, direction)
                        local_outputs = set()
                        final = _merged_pool(target, range_for_frag, mode, root_range)
                        grown = _intermediate(target, "merged", "grown", range_for_frag, mode, root_range)
                        if left_is_anchor and count == 1:
                            _append_action(recipe, _action(grown, _grow_params(snapshot, target, "fwd", input_pool)),
                                           grown, local_outputs)
                            _intersect(recipe, grown, restrict_pool, final, local_outputs)
                        else:
                            _append_action(recipe, _action(final, _grow_params(snapshot, target, "fwd",
                                                                                input_pool, restrict_pool)),
                                           final, local_outputs)
                        recipe["final_pools"] = [final]
                        _finalize_recipe(recipe, snapshot, local_outputs)
                        result.append(recipe)
            else:
                targets = list(range(right_run["left"] - 1, left_run["right"] - 1, -1))
                for target in targets:
                    input_pool = _final_pool(target + 1, range_for_frag, mode, root_range)
                    restrict_pool = _final_pool(target, range_for_frag, mode, root_range)
                    left_is_anchor = target == left_run["right"]
                    gaussian_count = _anchor_constraint_count(
                        snapshot, left_run, target, mode, root_range, "anchor") if left_is_anchor else 0
                    reference_count = _anchor_constraint_count(
                        snapshot, left_run, target, mode, root_range, "anchor-refe") if left_is_anchor else 0
                    if left_is_anchor and (gaussian_count == 1) != (reference_count == 1):
                        variants = [("anchor", gaussian_count), ("anchor-refe", reference_count)]
                    else:
                        variants = [(None, gaussian_count)]
                    for source_tree, count in variants:
                        path = _gen_path("merge", direction, mode, root_range, target,
                                         has_single_range=one_range,
                                         branch_ranges=(branch_ranges.get(("merge", direction), set())
                                                        if isinstance(branch_ranges, dict) else branch_ranges),
                                         source_tree=source_tree)
                        recipe = _recipe(path, mode, root_range, direction)
                        local_outputs = set()
                        final = _merged_pool(target, range_for_frag, mode, root_range)
                        grown = _intermediate(target, "merged", "grown", range_for_frag, mode, root_range)
                        if left_is_anchor and count == 1:
                            _append_action(recipe, _action(grown, _grow_params(snapshot, target, "bwd", input_pool)),
                                           grown, local_outputs)
                            _intersect(recipe, grown, restrict_pool, final, local_outputs)
                        else:
                            _append_action(recipe, _action(final, _grow_params(snapshot, target, "bwd",
                                                                                input_pool, restrict_pool)),
                                           final, local_outputs)
                        recipe["final_pools"] = [final]
                        _finalize_recipe(recipe, snapshot, local_outputs)
                        result.append(recipe)
    return result


def _common_conformer_recipes(snapshot: dict, classification: dict) -> list[dict]:
    result = []
    ranges_for_fragment = classification["ranges_for_fragment"]
    rows, by_number, _ = _fragment_rows(snapshot)
    for number in sorted(ranges_for_fragment):
        roots = sorted(ranges_for_fragment[number], key=_range_number)
        if len(roots) < 2:
            continue
        for first, second in itertools.combinations(roots, 2):
            path = f"common-conformer/gen-frag{number}-range-{_range_number(first)}-{_range_number(second)}"
            recipe = _recipe(path, "unfixed", None, None)
            recipe["range"] = None
            local_outputs = set()
            first_pool = _final_pool(number, ranges_for_fragment, "unfixed", first)
            second_pool = _final_pool(number, ranges_for_fragment, "unfixed", second)
            mask = f"frag{number}-range{_range_number(first)}-range{_range_number(second)}-mask"
            _append_action(recipe, _action(mask, {
                "action": "mask-common-conformer", "input1": first_pool, "input2": second_pool,
            }), mask, local_outputs)
            if len(roots) == 2:
                finals = [f"{first_pool}-common-conf", f"{second_pool}-common-conf"]
            else:
                finals = [f"frag{number}-range{_range_number(first)}-range{_range_number(second)}-common-conf",
                          f"frag{number}-range{_range_number(second)}-range{_range_number(first)}-common-conf"]
            for side, (source, final) in enumerate(zip((first_pool, second_pool), finals), start=1):
                _append_action(recipe, _action(final, {
                    "action": "apply-mask", "input": source, "mask_input": mask, "mask": side,
                }), final, local_outputs)
            recipe["final_pools"] = finals
            _finalize_recipe(recipe, snapshot, local_outputs)
            result.append(recipe)
    return result


def _recipe_signature(recipe: dict) -> str:
    """Signature used only for fixed/unfixed R7 hoisting."""
    actions = [{"name": row["name"], "parameters": row["parameters"], "comments": row["comments"]}
               for row in recipe["actions"]]
    external = [{key: row[key] for key in ("pool", "fragment", "sequence", "frame")}
                for row in recipe["external_inputs"]]
    value = {"direction": recipe["direction"], "final_pools": recipe["final_pools"],
             "external_inputs": external, "actions": actions}
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def _assign_paths_and_hoist(snapshot: dict, fixed: list[dict], unfixed: list[dict],
                            branch_counts: dict[tuple[str, str], set[str]]) -> list[dict]:
    all_ranges = len(classify_modes(snapshot)["ranges"])
    unfixed_by_sig: dict[str, list[dict]] = defaultdict(list)
    for recipe in unfixed:
        unfixed_by_sig[_recipe_signature(recipe)].append(recipe)
    result = []
    consumed = set()
    for recipe in fixed:
        signature = _recipe_signature(recipe)
        candidates = [candidate for candidate in unfixed_by_sig.get(signature, [])
                      if candidate["path"] not in consumed]
        if candidates:
            candidate = candidates[0]
            consumed.add(candidate["path"])
            merged = deepcopy(recipe)
            merged["modes"] = ["fixed", "unfixed"]
            merged["range"] = None
            # The fixed path is the shared R7 location.
            merged["path"] = re.sub(r"/(?:fixed|unfixed)(?:/range-[0-9]+)?/", "/", recipe["path"])
            result.append(merged)
        else:
            standalone = deepcopy(recipe)
            standalone["path"] = re.sub(r"/(?:fixed|unfixed)(?:/range-[0-9]+)?/", "/fixed/", recipe["path"])
            result.append(standalone)
    for recipe in unfixed:
        if recipe["path"] in consumed:
            continue
        standalone = deepcopy(recipe)
        direction = standalone["direction"] or "single"
        kind = standalone["path"].split("/")[0]
        branches = branch_counts.get((kind, direction), set())
        if all_ranges == 1:
            standalone["path"] = re.sub(r"/(?:fixed|unfixed)(?:/range-[0-9]+)?/", "/", standalone["path"])
            standalone["modes"] = ["fixed", "unfixed"]
            standalone["range"] = None
        else:
            mode_level = "/unfixed/"
            if mode_level not in standalone["path"]:
                standalone["path"] = standalone["path"].replace(f"/{direction}/", f"/{direction}/unfixed/", 1) \
                    if standalone["direction"] else standalone["path"].replace(f"{kind}/", f"{kind}/unfixed/", 1)
            if len(branches) > 1 and standalone["range"] is not None:
                current = f"/range-{_range_number(standalone['range'])}/"
                if current not in standalone["path"]:
                    standalone["path"] = standalone["path"].replace("/unfixed/", f"/unfixed{current}", 1)
        result.append(standalone)
    return result


def _populate_producers(recipes: list[dict]) -> None:
    producer_map: dict[str, set[str]] = defaultdict(set)
    for recipe in recipes:
        for pool in recipe["final_pools"]:
            producer_map[pool].add(recipe["path"])
    for recipe in recipes:
        path_parts = recipe["path"].split("/")
        source_tree = next((tree for tree in ("anchor", "anchor-refe")
                            if tree in path_parts), None) if path_parts[0] == "merge" else None
        target_number = _pool_fragment_number(recipe["final_pools"][0]) \
            if source_tree is not None and recipe["final_pools"] else None
        for external in recipe["external_inputs"]:
            producers = producer_map.get(external["pool"], set())
            # A tree-qualified boundary merge selects which anchor recipe may
            # provide its restriction pool. Its grow-side input remains a
            # normal external choice across every compatible producer.
            if source_tree is not None \
                    and _pool_fragment_number(external["pool"]) == target_number:
                producers = {path for path in producers if path.startswith(f"{source_tree}/")}
            external["producers"] = sorted(producers)
            if not external["producers"]:
                raise ValueError(f"External pool {external['pool']} has no possible producer.")


def _validate_catalogue(recipes: list[dict]) -> list[dict]:
    paths = set()
    member_keys = set()
    for recipe in recipes:
        validate_recipe(recipe)
        if recipe["path"] in paths:
            raise ValueError(f"Two generator recipes collide at {recipe['path']}.")
        paths.add(recipe["path"])
        for action in recipe["actions"]:
            member = f"{recipe['path']}/{action['name']}/alaric.yaml"
            if member in member_keys:
                raise ValueError(f"Two generators write action path {member}.")
            member_keys.add(member)
    ordered_members = sorted(member_keys)
    for previous, current in zip(ordered_members, ordered_members[1:]):
        if current.startswith(previous + "/"):
            raise ValueError(f"Generated archive key {previous} is a directory prefix of {current}.")
    ordered = sorted(recipes, key=lambda row: row["path"])
    return ordered


def plan_gens(snapshot: dict) -> list[dict]:
    """Plan the complete ordered R0–R12 generator catalogue for a snapshot."""
    if type(snapshot) is not dict or snapshot.get("schema") != 1:
        raise ValueError("A schema-1 export snapshot is required.")
    classification = classify_modes(snapshot)
    range_ids = classification["ranges"]
    one_range = len(range_ids) == 1
    ranges_for_fragment = classification["ranges_for_fragment"]

    # Branch roots decide whether R7 needs an explicit unfixed range folder for
    # every member of a (type, direction) branch.
    branch_ranges: dict[tuple[str, str | None], set[str]] = defaultdict(set)
    fixed_runs = _anchor_runs(snapshot, classification, "fixed", None)
    fixed_recipes = []
    for tree in ("anchor", "anchor-refe"):
        fixed_recipes.extend(_anchor_recipes(snapshot, classification, "fixed", None,
                                             tree, set()))
    fixed_recipes.extend(_ordinary_grow_recipes(snapshot, classification, "fixed", None,
                                                 fixed_runs, set()))
    fixed_recipes.extend(_merge_recipes(snapshot, classification, "fixed", None,
                                        fixed_runs, set()))

    unfixed_recipes = []
    for root_range in range_ids:
        runs = _anchor_runs(snapshot, classification, "unfixed", root_range)
        for tree in ("anchor", "anchor-refe"):
            unfixed_recipes.extend(_anchor_recipes(snapshot, classification, "unfixed", root_range,
                                                   tree, set()))
        unfixed_recipes.extend(_ordinary_grow_recipes(snapshot, classification, "unfixed", root_range,
                                                      runs, set()))
        unfixed_recipes.extend(_merge_recipes(snapshot, classification, "unfixed", root_range,
                                              runs, set()))

    # A range folder is used consistently only inside a branch populated by
    # more than one unfixed root.  Count actual emitted direction/type pairs,
    # then insert the level into each affected path before comparing recipes.
    for recipe in unfixed_recipes:
        kind = recipe["path"].split("/")[0]
        branch_ranges[(kind, recipe["direction"])].add(recipe["range"])
    if not one_range:
        for recipe in unfixed_recipes:
            branch = branch_ranges[(recipe["path"].split("/")[0], recipe["direction"])]
            if len(branch) > 1:
                range_part = f"range-{_range_number(recipe['range'])}"
                recipe["path"] = recipe["path"].replace("/unfixed/", f"/unfixed/{range_part}/", 1)

    # Give raw candidates structural paths first; R7 hoists only exact recipe
    # matches, preserving mutually exclusive direction and range variants.
    def seed_path(recipe: dict, mode: str) -> str:
        return recipe["path"]

    all_recipes = _assign_paths_and_hoist(snapshot, fixed_recipes, unfixed_recipes, branch_ranges)
    all_recipes.extend(_common_conformer_recipes(snapshot, classification))
    # One-range projects have identical fixed/unfixed recipes by contract.
    if one_range:
        for recipe in all_recipes:
            if recipe["path"].startswith("common-conformer/"):
                continue
            recipe["modes"] = ["fixed", "unfixed"]
            recipe["range"] = None
    _populate_producers(all_recipes)
    return _validate_catalogue(all_recipes)

