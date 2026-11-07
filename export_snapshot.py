"""Pure export-input normalisation, admission, warnings, and snapshot identity.

The functions in this module consume the named values pinned by the export
transformers.  They deliberately do not access a Seamless Context, deployment
mount, local installation, or browser state.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from decimal import Decimal, InvalidOperation
from typing import Any


EXPORTER_VERSION = 1
ANCHOR_MIN_CELL_FRACTION = 1.0 / 3.0
L_OV = math.sqrt(3.0)

_FRAGMENT_KEY = re.compile(r"^frag-([1-9][0-9]*)-([ACGU]{2})$")
_PAIR_DECIMAL = re.compile(r"^(?:0|[1-9][0-9]*)(?:\.[0-9]{1,2})?$")
_ENERGY_DECIMAL = re.compile(r"^-?(?:0|[1-9][0-9]*)(?:\.[0-9]{1,2})?$")
_STACKING_KEY_FIELDS = ("rna_resname", "rna_resid", "protein_range_pdb_id", "protein_resid")
_RING_PINS = (
    "ring_refe_a", "ring_refe_c", "ring_refe_phe", "ring_refe_his",
    "ring_refe_arg", "ring_refe_trp",
)
_VALIDATOR_PINS = (
    "validated_pairs", "validated_classes", "validated_energies",
)


class ExportNotReady(Exception):
    """Expected refusal to build a tree for a blocked or stale admission."""


class _InputProblem(ValueError):
    def __init__(self, code: str, address: str | None, message: str, remedy: str):
        super().__init__(message)
        self.code = code
        self.address = address
        self.message = message
        self.remedy = remedy


def _block(code: str, address: str | None, message: str, remedy: str) -> dict:
    return {"code": code, "address": address, "message": message, "remedy": remedy}


def _problem(code: str, address: str | None, message: str, remedy: str) -> None:
    raise _InputProblem(code, address, message, remedy)


def _is_int(value: Any) -> bool:
    return type(value) is int


def _is_finite_number(value: Any) -> bool:
    if type(value) not in (int, float):
        return False
    try:
        return math.isfinite(float(value))
    except (OverflowError, ValueError):
        return False


def _is_record(value: Any) -> bool:
    return type(value) is dict


def _canonical_json(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def _sha_bytes(value: bytes) -> dict:
    return {"sha256": hashlib.sha256(value).hexdigest(), "length": len(value)}


def _digest_value(value: Any) -> Any:
    if isinstance(value, bytes):
        return {"$bytes": _sha_bytes(value)}
    if isinstance(value, dict):
        return {str(key): _digest_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_digest_value(item) for item in value]
    if value is None or type(value) in (str, bool, int, float):
        if type(value) is float and not math.isfinite(value):
            raise ValueError("A snapshot cannot contain a non-finite number.")
        return value
    raise TypeError(f"Unsupported value in export snapshot: {type(value).__name__}.")


def snapshot_digest(snapshot: Any) -> str:
    """Return the canonical SHA-256 identity for a normalised snapshot."""
    payload = {"exporter_version": EXPORTER_VERSION, "snapshot": _digest_value(snapshot)}
    return hashlib.sha256(_canonical_json(payload)).hexdigest()


def _range_id(number: int) -> str:
    return f"range-{number}"


def _range_signature(row: dict) -> str:
    return json.dumps(row, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def _fragment_sort_key(row: dict) -> tuple[int, str]:
    return row["number"], row["sequence"]


def _parse_model(model_json: Any) -> dict:
    if isinstance(model_json, bytes):
        try:
            model_json = model_json.decode("utf-8")
        except UnicodeDecodeError:
            _problem("gaussian-model", None, "The Gaussian model is not UTF-8 JSON.",
                     "Repair the installed Gaussian model and recompute stage 6.")
    if not isinstance(model_json, str):
        _problem("gaussian-model", None, "The Gaussian model is unavailable.",
                 "Repair the installed Gaussian model and recompute stage 6.")
    try:
        value = json.loads(model_json)
    except json.JSONDecodeError:
        _problem("gaussian-model", None, "The Gaussian model is not valid JSON.",
                 "Repair the installed Gaussian model and recompute stage 6.")
    if not _is_record(value) or type(value.get("schema")) is not int or value["schema"] != 1 \
        or not isinstance(value.get("filter_set_version"), str) \
        or not _is_record(value.get("thresholds")):
        _problem("gaussian-model", None, "The Gaussian model must use schema 1 and declare its filters.",
                 "Install a schema-1 Gaussian model and recompute stage 6.")
    return value


def _selected_keys(selection: Any) -> list[str]:
    if not _is_record(selection) or type(selection.get("schema")) is not int \
        or selection.get("schema") != 1 or selection.get("active") is not True \
        or not _is_record(selection.get("selected")):
        _problem("stacking-selection", None, "The current stacking selection is not active schema 1.",
                 "Complete stage 3 and recompute the selection.")
    if any(type(flag) is not bool for flag in selection["selected"].values()):
        _problem("stacking-selection", None, "The stacking selection contains a non-boolean choice.",
                 "Repair the stage-3 selection and recompute stage 3.")
    return sorted(key for key, selected in selection["selected"].items() if selected)


def _rows_by_key(value: Any, field: str, address_field: str = "key") -> dict[str, dict]:
    rows = value.get(field) if _is_record(value) else None
    if not isinstance(rows, list):
        _problem("boundary-schema", None, f"The {field} boundary rows are missing.",
                 "Recompute the producer for this boundary.")
    result: dict[str, dict] = {}
    for row in rows:
        if not _is_record(row) or not isinstance(row.get(address_field), str):
            _problem("boundary-schema", None, f"A {field} row has no stable address.",
                     "Recompute the producer for this boundary.")
        address = row[address_field]
        if address in result:
            _problem("duplicate-address", address, f"The {field} boundary repeats {address}.",
                     "Recompute the producer and remove the duplicate row.")
        result[address] = row
    return result


def _validate_schema1(name: str, value: Any, address: str | None = None) -> None:
    if not _is_record(value) or type(value.get("schema")) is not int or value.get("schema") != 1:
        _problem("boundary-schema", address, f"{name} must use top-level schema 1.",
                 f"Recompute {name} with a schema-1 producer.")


def _validate_fragments(fragments: Any) -> tuple[list[dict], dict[str, dict], list[dict]]:
    _validate_schema1("the canonical fragment list", fragments)
    rows = fragments.get("fragments")
    residues = fragments.get("residues")
    adjacent = fragments.get("adjacent_pairs")
    if not isinstance(rows, list) or not rows or not isinstance(residues, list) \
        or not isinstance(adjacent, list) or not _is_record(fragments.get("echo")):
        _problem("fragment-map", None, "The canonical fragment map is incomplete.",
                 "Recompute the fragment list from the prepared RNA.")
    seen: set[str] = set()
    by_key: dict[str, dict] = {}
    previous = 0
    for row in rows:
        if not _is_record(row) or not isinstance(row.get("key"), str) \
            or not _is_int(row.get("number")) or not isinstance(row.get("sequence"), str) \
            or re.fullmatch(r"[ACGU]{2}", row.get("sequence", "")) is None:
            _problem("fragment-map", None, "A canonical fragment row is malformed.",
                     "Recompute the fragment list from the prepared RNA.")
        match = _FRAGMENT_KEY.fullmatch(row["key"])
        if not match or int(match.group(1)) != row["number"] or match.group(2) != row["sequence"]:
            _problem("fragment-map", row["key"], "A fragment key, ordinal, and sequence disagree.",
                     "Recompute the fragment list from the prepared RNA.")
        if row["key"] in seen or row["number"] <= previous:
            _problem("fragment-order", row["key"], "Fragment addresses must be unique and ordinally ordered.",
                     "Recompute the fragment list from the prepared RNA.")
        if not isinstance(row.get("residues"), list) or len(row["residues"]) != 2 \
            or any(not _is_record(residue) or not _is_int(residue.get("resid")) for residue in row["residues"]):
            _problem("fragment-map", row["key"], "A fragment must identify both RNA residues.",
                     "Recompute the fragment list from the prepared RNA.")
        seen.add(row["key"])
        previous = row["number"]
        by_key[row["key"]] = row

    expected_pairs = []
    by_number = {row["number"]: row for row in rows}
    for number, down in sorted(by_number.items()):
        up = by_number.get(number + 1)
        if up is None or down["residues"][1] != up["residues"][0]:
            continue
        expected_pairs.append({
            "down": number, "up": number + 1,
            "sequence": down["sequence"] + up["sequence"][1],
        })
    if adjacent != expected_pairs:
        _problem("fragment-pairs", None, "Adjacent fragment addresses or trinucleotide sequences disagree.",
                 "Recompute the canonical fragment list.")
    return rows, by_key, expected_pairs


def _validate_clip(clip: Any, fragments: list[dict], fragment_map: dict[str, dict],
                   expected_echo: dict) -> tuple[list[dict], set[int]]:
    if not _is_record(clip) or type(clip.get("schema")) is not int or clip.get("schema") != 1 \
        or not _is_int(clip.get("clip1")) or not _is_int(clip.get("clip2")) \
        or clip["clip1"] < 0 or clip["clip2"] < 0 \
        or clip["clip1"] + clip["clip2"] >= len(fragments):
        _problem("clip", None, "The current clip is not a valid schema-1 fragment interval.",
                 "Adjust the clip while keeping at least one fragment.")
    if clip.get("basis") != expected_echo:
        _problem("stale-clip", None, "The clip belongs to a different RNA range.",
                 "Reopen stage 4 and confirm the clip for the current RNA range.")
    kept = fragments[clip["clip1"]:len(fragments) - clip["clip2"]]
    kept_numbers = {row["number"] for row in kept}
    return kept, kept_numbers


def _validate_ranges(pdbcode: Any, biounit: Any, range_decision: Any,
                     range_specs: Any, structure_id: Any) -> tuple[list[dict], dict]:
    if not isinstance(pdbcode, str) or not pdbcode.strip():
        _problem("structure", None, "A non-empty structure code is required.",
                 "Select a structure and complete stage 1.")
    if not _is_int(biounit) or biounit < 1:
        _problem("structure", None, "The selected biological assembly is invalid.",
                 "Select a valid biological assembly in stage 1.")
    if not isinstance(structure_id, str) or structure_id.strip().upper() != pdbcode.strip().upper():
        _problem("structure-identity", None, "The stored structure does not match the selected code.",
                 "Download the selected structure again and recompute stage 1.")
    if not _is_record(range_specs) or set(range_specs) != {"schema", "pdbcode", "biounit", "rows"} \
        or type(range_specs.get("schema")) is not int \
        or range_specs.get("schema") != 1 or range_specs.get("pdbcode") != pdbcode \
        or range_specs.get("biounit") != biounit or not isinstance(range_specs.get("rows"), list):
        _problem("range-specs", None, "Normalised range specifications are not current schema 1.",
                 "Complete stage 2 for the current structure and assembly.")
    rows = range_specs["rows"]
    if len(rows) < 2 or any(
        not _is_record(row) or set(row) != {"chain", "start", "end"}
        or not isinstance(row.get("chain"), str) or not row["chain"]
        or not _is_int(row.get("start")) or not _is_int(row.get("end"))
        or row["end"] < row["start"] for row in rows
    ):
        _problem("range-specs", None, "At least one RNA range and one protein range are required.",
                 "Complete stage 2 with valid RNA and protein ranges.")
    if type(range_decision) is str:
        # The editable S4 cell is source text; range_specs is its validated,
        # normalized projection and supplies the stable decision identity.
        if not range_decision.strip():
            _problem("range-decision", None, "The range decision is empty.",
                     "Complete stage 2 with valid ranges.")
    elif _is_record(range_decision):
        if set(range_decision) != {"schema", "basis", "active", "rows"} \
            or type(range_decision.get("schema")) is not int or range_decision.get("schema") != 1 \
            or range_decision.get("active") is not True \
            or range_decision.get("basis") != {"pdbcode": pdbcode, "biounit": biounit} \
            or range_decision.get("rows") != rows:
            _problem("range-decision", None, "The range decision is not current schema 1 for this structure.",
                     "Reopen stage 2 and confirm the current structure ranges.")
    elif range_decision is None:
        _problem("range-decision", None, "The range decision is not active.",
                 "Complete stage 2 with valid ranges.")
    else:
        _problem("range-decision", None, "The range decision has an unsupported shape.",
                 "Repair the range decision and recompute stage 2.")
    return [dict(row) for row in rows], {"chain": rows[0]["chain"], "start": rows[0]["start"], "end": rows[0]["end"]}


def _validate_file_folders(inputs: dict, range_rows: list[dict], pdbcode: str,
                           rna_range: dict, fit_echo: Any) -> tuple[dict, dict, bytes]:
    b1, b3, rna = inputs.get("range_pdb_files"), inputs.get("range_pdb_files_aa"), inputs.get("rna_pdb")
    if not _is_record(b1) or not _is_record(b3):
        _problem("prepared-files", None, "Prepared range files are not available as folders.",
                 "Recompute structure preparation for the current ranges.")
    expected_b1 = [f"range-{index}.pdb" for index in range(1, len(range_rows) + 1)]
    expected_b3 = [f"range-{index}-aa.pdb" for index in range(1, len(range_rows) + 1)]
    if set(b1) != set(expected_b1) or any(not isinstance(b1.get(name), bytes) for name in expected_b1):
        _problem("prepared-files", None, "Prepared plain range files do not match the current ranges.",
                 "Recompute structure preparation for the current ranges.")
    if set(b3) != set(expected_b3) or any(not isinstance(b3.get(name), bytes) for name in expected_b3):
        _problem("prepared-files", None, "Prepared all-atom range files do not match the current ranges.",
                 "Recompute structure preparation for the current ranges.")
    exact_rna = b3["range-1-aa.pdb"]
    if not isinstance(rna, bytes) or rna != exact_rna:
        _problem("reference-bytes", "range-1", "The displayed RNA bytes differ from the prepared RNA range.",
                 "Recompute stage 3 from the current all-atom RNA file.")
    displayed = inputs.get("range_pdbs_aa")
    if not _is_record(displayed) or set(displayed) != {"schema", "pdb_code", "range_specs", "pdbs"} \
        or type(displayed.get("schema")) is not int \
        or displayed.get("schema") != 1 or displayed.get("pdb_code") != pdbcode.strip().lower() \
        or displayed.get("range_specs") != inputs.get("range_specs") \
        or not isinstance(displayed.get("pdbs"), list) or len(displayed["pdbs"]) != len(range_rows):
        _problem("displayed-files", None, "The displayed all-atom ranges are not current schema 1.",
                 "Recompute the displayed range files from the current preparation.")
    for index, item in enumerate(displayed["pdbs"], start=1):
        name = f"range-{index}-aa.pdb"
        try:
            decoded = b3[name].decode("utf-8")
        except UnicodeDecodeError:
            _problem("prepared-files", _range_id(index), "A prepared all-atom PDB is not UTF-8 text.",
                     "Recompute structure preparation for the current ranges.")
        if not _is_record(item) or item.get("range") != _range_id(index) \
            or item.get("name") != name or item.get("pdb") != decoded:
            _problem("displayed-files", _range_id(index),
                     "A displayed all-atom range does not match the exact prepared bytes.",
                     "Recompute the displayed range files from the current preparation.")
    exclude = inputs.get("exclude_pdb_code")
    if not isinstance(exclude, str) or not exclude or exclude != pdbcode.strip().lower():
        _problem("exclude-code", None, "The fit and Gaussian exclusion code must be the lower-cased structure code.",
                 "Recompute stage 4 and stage 6 with the current structure code.")
    expected_residues = None
    if _is_record(fit_echo):
        expected_residues = fit_echo.get("resids")
    expected_fit_echo = {
        "pdb_code": pdbcode.strip().lower(), "rna": rna_range,
        "resids": expected_residues,
    }
    if not _is_record(fit_echo) or fit_echo.get("pdb_code") != expected_fit_echo["pdb_code"] \
        or fit_echo.get("rna") != rna_range or not isinstance(expected_residues, list):
        _problem("fit-echo", "range-1", "The fit results do not echo the current RNA range and exclusion code.",
                 "Recompute best fits for the current prepared RNA.")
    return dict(b1), dict(b3), exact_rna


def _validate_pair_data(fragments: list[dict], kept_numbers: set[int], expected_pairs: list[dict],
                        computed: Any, edits: Any, clipped: Any) -> list[dict]:
    if not _is_record(computed) or type(computed.get("schema")) is not int or computed.get("schema") != 1 \
        or not isinstance(computed.get("pairs"), list) or not isinstance(computed.get("bases"), list):
        _problem("thresholds", None, "Computed pair and base thresholds are not schema 1.",
                 "Recompute stage 4 thresholds from the current best-fit tables.")
    if not _is_record(edits) or type(edits.get("schema")) is not int or edits.get("schema") != 1 \
        or edits.get("active") is not True or not isinstance(edits.get("pairs"), list):
        _problem("threshold-edits", None, "Stage 4 pair thresholds are not an active schema-1 decision.",
                 "Complete stage 4 and fill the pair thresholds for every kept pair.")
    computed_by_pair = {}
    for row in computed["pairs"]:
        if not _is_record(row) or not _is_int(row.get("down")) or not _is_int(row.get("up")):
            _problem("pair-address", None, "A computed pair threshold has an invalid address.",
                     "Recompute stage 4 thresholds.")
        address = f"frag-{row['down']}..frag-{row['up']}"
        key = (row["down"], row["up"])
        if key in computed_by_pair:
            _problem("pair-address", address, "A computed pair threshold is duplicated.",
                     "Recompute stage 4 thresholds.")
        if not _is_finite_number(row.get("ovrmsd_raw")) or not _is_finite_number(row.get("crmsd_raw")) \
            or row["ovrmsd_raw"] <= 0 or row["crmsd_raw"] <= 0:
            _problem("pair-value", address, "Raw pair thresholds must be positive finite numbers.",
                     "Recompute stage 4 thresholds from the current library values.")
        if any(not isinstance(row.get(name), str) or _PAIR_DECIMAL.fullmatch(row[name]) is None
               for name in ("ovrmsd", "crmsd")):
            _problem("pair-value", address, "Rounded pair thresholds must be positive finite values with at most two decimals.",
                     "Recompute stage 4 thresholds from the current library values.")
        computed_by_pair[key] = row
    expected_all = {(row["down"], row["up"]) for row in expected_pairs}
    if set(computed_by_pair) != expected_all:
        _problem("pair-address", None, "Computed pair thresholds do not match the current fragment map.",
                 "Recompute stage 4 thresholds.")
    bases = computed["bases"]
    if len(bases) != len(fragments):
        _problem("base-thresholds", None, "Computed base thresholds do not match the full fragment map.",
                 "Recompute stage 4 base thresholds.")
    for fragment, row in zip(fragments, bases):
        fields = ("nuc1_ovrmsd", "nuc2_ovrmsd")
        raw_fields = ("nuc1_ovrmsd_raw", "nuc2_ovrmsd_raw")
        if not _is_record(row) or row.get("fragment") != fragment["number"] \
            or row.get("sequence") != fragment["sequence"] \
            or any(not isinstance(row.get(name), str) or _PAIR_DECIMAL.fullmatch(row[name]) is None
                   for name in fields) \
            or any(not _is_finite_number(row.get(name)) or row[name] <= 0 for name in raw_fields):
            _problem("base-thresholds", fragment["key"], "Computed base thresholds are malformed.",
                     "Recompute stage 4 base thresholds for the current fragments.")
    edit_by_pair = {}
    for row in edits["pairs"]:
        if not _is_record(row) or not _is_int(row.get("down")) or not _is_int(row.get("up")):
            _problem("pair-address", None, "An edited pair threshold has an invalid address.",
                     "Repair the stage 4 threshold table.")
        key = (row["down"], row["up"])
        address = f"frag-{key[0]}..frag-{key[1]}"
        if key in edit_by_pair or key not in expected_all:
            _problem("pair-address", address, "An edited pair is duplicated or absent from the current fragment map.",
                     "Repair the stage 4 threshold table and recompute stage 4.")
        edit_by_pair[key] = row
    if set(edit_by_pair) != expected_all:
        _problem("pair-address", None, "Edited pair rows do not exactly match the current adjacent fragments.",
                 "Complete every pair row in stage 4.")

    result = []
    for pair in expected_pairs:
        key = (pair["down"], pair["up"])
        source, edited = computed_by_pair[key], edit_by_pair[key]
        address = f"frag-{key[0]}..frag-{key[1]}"
        if source.get("sequence") != pair["sequence"] or edited.get("sequence") != pair["sequence"]:
            _problem("pair-sequence", address, "The pair sequence does not match its adjacent fragments.",
                     "Recompute the fragment and threshold tables.")
        row = dict(pair)
        for name in ("ovrmsd", "crmsd"):
            value = edited.get(name)
            if key[0] in kept_numbers and key[1] in kept_numbers:
                if not isinstance(value, str) or _PAIR_DECIMAL.fullmatch(value) is None:
                    _problem("pair-threshold", address, f"Kept pair {name} must be a positive two-decimal value.",
                             "Fill the missing kept pair threshold in stage 4.")
                try:
                    number = Decimal(value)
                except InvalidOperation:
                    number = Decimal("NaN")
                if not number.is_finite() or number <= 0:
                    _problem("pair-threshold", address, f"Kept pair {name} must be positive and finite.",
                             "Correct the pair threshold in stage 4.")
                if name == "ovrmsd" and float(number) > L_OV:
                    _problem("pair-coverage", address, "The edited pair ovRMSD exceeds Alaric's guaranteed coverage limit.",
                             "Lower the pair ovRMSD threshold or revise the fit inputs; Proceed anyway cannot override this block.")
            row[name] = value
            raw_name = f"{name}_raw"
            row[raw_name] = source.get(raw_name)
        result.append(row)

    # B11 is a versioned addressed boundary.  Its exact byte PDB members are
    # retained below; these two metadata files are independently checked here.
    if not _is_record(clipped) or "thresholds.json" not in clipped or "clip.json" not in clipped:
        _problem("clipped-boundary", None, "The kept-fragment folder is missing its threshold or clip record.",
                 "Recompute the kept-fragment boundary from current stage-4 decisions.")
    if not isinstance(clipped["thresholds.json"], bytes) or not isinstance(clipped["clip.json"], bytes):
        _problem("clipped-boundary", None, "The kept-fragment metadata members must be bytes.",
                 "Recompute the kept-fragment boundary from current stage-4 decisions.")
    try:
        clipped_thresholds = json.loads(clipped["thresholds.json"])
        clipped_clip = json.loads(clipped["clip.json"])
    except (TypeError, json.JSONDecodeError, UnicodeDecodeError):
        _problem("clipped-boundary", None, "The kept-fragment metadata is not valid JSON.",
                 "Recompute the kept-fragment boundary from current stage-4 decisions.")
    if not _is_record(clipped_thresholds) or set(clipped_thresholds) != {"schema", "pairs"} \
        or type(clipped_thresholds.get("schema")) is not int \
        or clipped_thresholds["schema"] != 1 or not isinstance(clipped_thresholds.get("pairs"), list) \
        or not _is_record(clipped_clip) or set(clipped_clip) != {"schema", "clip1", "clip2", "kept_residues"} \
        or type(clipped_clip.get("schema")) is not int \
        or clipped_clip["schema"] != 1:
        _problem("clipped-boundary", None, "The kept-fragment metadata must use schema 1.",
                 "Recompute the kept-fragment boundary with a schema-1 producer.")
    kept_pairs = [row for row in result if row["down"] in kept_numbers and row["up"] in kept_numbers]
    if clipped_thresholds["pairs"] != [
        {key: row.get(key) for key in ("down", "up", "sequence", "ovrmsd", "crmsd")}
        for row in kept_pairs
    ]:
        _problem("clipped-boundary", None, "The kept-fragment pair record disagrees with edited thresholds.",
                 "Recompute the kept-fragment boundary from current stage-4 decisions.")
    return kept_pairs


def _stacking_key(row: dict) -> str:
    if any(key not in row for key in _STACKING_KEY_FIELDS):
        _problem("stacking-row", None, "A stacking row is missing identity fields.",
                 "Recompute stage 3 stacking detection.")
    return json.dumps([
        str(row["rna_resname"]).strip(), row["rna_resid"],
        row["protein_range_pdb_id"], row["protein_resid"],
    ], separators=(",", ":"), ensure_ascii=False)


def _canonical_rna_base(resname: Any) -> str | None:
    if not isinstance(resname, str):
        return None
    name = resname.strip().upper()
    if len(name) == 2 and name[0] in {"R", "D"}:
        name = name[1]
    return name if name in {"A", "C", "G", "U"} else None


def _candidate_anchors(fragments: Any, kept_fragments: list[dict],
                       selection: Any, stackings: Any, matrix_ranges: list[str],
                       expected_echo: dict) -> tuple[list[dict], list[dict]]:
    selected = _selected_keys(selection)
    _validate_schema1("the stacking table", stackings)
    if stackings.get("echo") != expected_echo:
        address = "range-1"
        actual_echo = stackings.get("echo")
        if _is_record(actual_echo) and actual_echo.get("rows") != expected_echo.get("rows"):
            actual_rows = actual_echo.get("rows")
            if isinstance(actual_rows, list):
                for index, expected_row in enumerate(expected_echo.get("rows", []), start=1):
                    if index > len(actual_rows) or actual_rows[index - 1] != expected_row:
                        address = _range_id(index)
                        break
                else:
                    if len(actual_rows) < len(expected_echo.get("rows", [])):
                        address = _range_id(len(actual_rows) + 1)
            else:
                address = "range-1"
        _problem("stackings-echo", address,
                 "Detected stackings do not echo the current structure identity, code, assembly, and range rows.",
                 "Recompute stage 3 stacking detection for the current structure and ranges.")
    stack_rows = stackings.get("rows")
    if not isinstance(stack_rows, list):
        _problem("stackings", None, "Fresh stacking rows are unavailable.",
                 "Recompute stage 3 stacking detection.")
    fragment_rows = fragments.get("fragments") if _is_record(fragments) else None
    if not isinstance(fragment_rows, list):
        _problem("fragment-map", None, "The candidate anchor helper needs the canonical fragment map.",
                 "Recompute the canonical fragment list for the prepared RNA.")
    by_key = {}
    for row in stack_rows:
        if not _is_record(row):
            _problem("stacking-row", None, "A stacking row is malformed.",
                     "Recompute stage 3 stacking detection.")
        key = _stacking_key(row)
        if key in by_key:
            _problem("stacking-row", key, "The stacking table repeats an identity.",
                     "Recompute stage 3 stacking detection.")
        by_key[key] = row
    missing = sorted(set(selected) - set(by_key))
    if missing:
        _problem("stacking-selection", missing[0], "A selected stacking is absent from the current stacking table.",
                 "Recompute stage 3 and restore a current selection.")
    selected_rows = []
    for key in selected:
        row = by_key[key]
        rna_resid = row.get("rna_resid")
        range_index = row.get("protein_range_pdb_id")
        protein_resid = row.get("protein_resid")
        if not _is_int(rna_resid) or not _is_int(range_index) or range_index < 2 \
            or not _is_int(protein_resid):
            _problem("stacking-address", key, "A selected stacking uses a residue identity unsupported by the action schema.",
                     "Correct the selected stacking or prepare a protein residue with an integer identifier.")
        range_id = _range_id(range_index)
        if range_id not in matrix_ranges:
            _problem("stacking-range", key, f"The selected stacking refers to absent {range_id}.",
                     "Reopen stage 3 and select a stacking from a current protein range.")
        observed_bases = {
            residue.get("base")
            for fragment in fragment_rows
            for residue in fragment.get("residues", [])
            if _is_record(residue) and residue.get("resid") == rna_resid
        }
        canonical_base = _canonical_rna_base(row.get("rna_resname"))
        if canonical_base is None or observed_bases != {canonical_base}:
            _problem("stacking-rna-identity", key,
                     "The selected stacking residue name disagrees with the canonical RNA residue identity.",
                     "Recompute stage 3 stacking detection from the current prepared RNA.")
        locks_all = []
        for fragment in fragment_rows:
            for nucleotide, residue in enumerate(fragment.get("residues", []), start=1):
                if _is_record(residue) and residue.get("resid") == rna_resid:
                    lock = {"fragment": fragment["key"], "nucleotide": nucleotide}
                    if lock not in locks_all:
                        locks_all.append(lock)
        selected_rows.append({
            "key": key, "rna_resname": str(row["rna_resname"]).strip(),
            "rna_resid": rna_resid, "range": range_id,
            "protein_resid": protein_resid,
            "locks": [dict(lock) for lock in locks_all if lock["fragment"] in {f["key"] for f in kept_fragments}],
            "candidate_locks": [dict(lock) for lock in locks_all],
        })
    try:
        import preparic_stage5

        anchors = preparic_stage5.candidate_anchors(
            fragments, selection, stackings, [row["key"] for row in kept_fragments],
        )
    except ValueError as exc:
        _problem("candidate-anchors", selected[0] if selected else None, str(exc),
                 "Recompute stage 3 stacking selection and stage 5 anchor locks.")
    for anchor in anchors:
        if anchor.get("range") not in matrix_ranges:
            _problem("stacking-range", anchor.get("stackings", [None])[0],
                     "A selected anchor has no current interaction-matrix range.",
                     "Recompute stage 3 selection for the current protein ranges.")
    return anchors, selected_rows


def _validate_interactions(matrix: Any, handoff: Any, fragments: list[dict],
                           kept_fragments: list[dict], anchors: list[dict],
                           range_specs: list[dict], clip: Any, selected_keys: list[str]) -> dict:
    _validate_schema1("the interaction matrix", matrix)
    ranges = [_range_id(index) for index in range(2, len(range_specs) + 1)]
    kept_keys = [row["key"] for row in kept_fragments]
    if matrix.get("active") is not True or matrix.get("rows") != kept_keys \
        or matrix.get("ranges") != ranges or matrix.get("basis") != selected_keys:
        _problem("interaction-matrix", None, "The stage-5 matrix is stale or has the wrong addressed rows.",
                 "Reopen stage 5 and recompute interactions for the current clip and selections.")
    grid = matrix.get("grid")
    if not isinstance(grid, list) or len(grid) != len(kept_keys) or any(
        not isinstance(row, list) or len(row) != len(ranges)
        or any(type(value) is not bool for value in row) for row in grid
    ):
        _problem("interaction-grid", None, "The interaction matrix grid has an invalid shape or cell type.",
                 "Recompute the stage-5 interaction matrix.")
    row_pos = {key: index for index, key in enumerate(kept_keys)}
    range_pos = {key: index for index, key in enumerate(ranges)}
    for anchor in anchors:
        col = range_pos[anchor["range"]]
        for lock in anchor["locks"]:
            row = row_pos[lock["fragment"]]
            if grid[row][col] is not True:
                _problem("anchor-lock", lock["fragment"],
                         f"The selected {anchor['range']} anchor lock is false in the current interaction matrix.",
                         "Correct the interaction matrix so every selected anchor lock is true.")
    if any(not any(row) for row in grid):
        row = kept_keys[next(index for index, cells in enumerate(grid) if not any(cells))]
        _problem("uninteracting-fragment", row, "A kept fragment does not interact with any protein range.",
                 "Revise the stage-5 matrix or clip away this fragment.")
    _validate_schema1("the interaction hand-off", handoff)
    expected_ranges = [
        {"index": index, **dict(row)} for index, row in enumerate(range_specs[1:], start=2)
    ]
    if handoff.get("clip1") != clip.get("clip1") or handoff.get("clip2") != clip.get("clip2") \
        or handoff.get("ranges") != expected_ranges:
        _problem("handoff", None, "The interaction hand-off does not match the current clip and protein ranges.",
                 "Recompute the stage-5 interaction hand-off.")
    if not isinstance(handoff.get("runs"), list) or not isinstance(handoff.get("kept_residues"), list) \
        or not _is_record(handoff.get("interactions")):
        _problem("handoff", None, "The interaction hand-off is missing runs, kept residues, or fragment membership.",
                 "Recompute the stage-5 interaction hand-off.")
    expected_interactions = {
        key: [range_id for col, range_id in enumerate(ranges) if grid[row][col]]
        for row, key in enumerate(kept_keys)
    }
    if handoff["interactions"] != expected_interactions:
        _problem("handoff-matrix", None, "The hand-off fragment/range memberships disagree with the current matrix.",
                 "Recompute the stage-5 matrix and hand-off.")
    for run in handoff["runs"]:
        if not _is_record(run) or run.get("range") not in ranges \
            or not _is_int(run.get("first_fragment")) or not _is_int(run.get("last_fragment")) \
            or run["last_fragment"] < run["first_fragment"] or not isinstance(run.get("anchors"), list):
            _problem("handoff-run", None, "A hand-off run has a malformed range or fragment interval.",
                     "Recompute the stage-5 interaction hand-off.")

    # Rebuild the exact ordered run partition from the current matrix and
    # selected locks.  The hand-off is a boundary echo, not an independent
    # source of interaction or anchor evidence.
    anchors_by_range: dict[str, list[dict]] = {range_id: [] for range_id in ranges}
    for anchor in anchors:
        if anchor["range"] not in anchors_by_range:
            _problem("handoff-anchor", None, "A selected anchor has no current matrix range.",
                     "Recompute stage 3 selection and the stage-5 interaction matrix.")
        anchors_by_range[anchor["range"]].append(anchor)
    expected_runs = []
    for column, range_id in enumerate(ranges):
        position = 0
        while position < len(kept_keys):
            if not grid[position][column]:
                position += 1
                continue
            first = position
            while position + 1 < len(kept_keys) and grid[position + 1][column]:
                position += 1
            last = position
            run_anchors = []
            for anchor in anchors_by_range[range_id]:
                for lock in anchor["locks"]:
                    lock_position = row_pos[lock["fragment"]]
                    if first <= lock_position <= last:
                        candidate = {"fragment": lock["fragment"], "nucleotide": lock["nucleotide"]}
                        if candidate not in run_anchors:
                            run_anchors.append(candidate)
            first_number = int(_FRAGMENT_KEY.fullmatch(kept_keys[first]).group(1))
            last_number = int(_FRAGMENT_KEY.fullmatch(kept_keys[last]).group(1))
            expected_runs.append({
                "range": range_id, "first_fragment": first_number,
                "last_fragment": last_number, "anchors": run_anchors,
            })
            position += 1

    fragment_by_number = {row["number"]: row["key"] for row in kept_fragments}

    def run_address(expected: Any, actual: Any) -> str:
        expected_anchors = expected.get("anchors", []) if _is_record(expected) else []
        actual_anchors = actual.get("anchors", []) if _is_record(actual) else []
        if isinstance(expected_anchors, list) and isinstance(actual_anchors, list):
            for index in range(max(len(expected_anchors), len(actual_anchors))):
                left = expected_anchors[index] if index < len(expected_anchors) else None
                right = actual_anchors[index] if index < len(actual_anchors) else None
                if left != right:
                    for candidate in (left, right):
                        if _is_record(candidate) and isinstance(candidate.get("fragment"), str):
                            return candidate["fragment"]
        for candidate in (expected, actual):
            if _is_record(candidate):
                number = candidate.get("first_fragment")
                if _is_int(number) and number in fragment_by_number:
                    return fragment_by_number[number]
        return kept_keys[0] if kept_keys else "range-1"

    actual_runs = handoff["runs"]
    if actual_runs != expected_runs:
        mismatch = next(
            (index for index in range(max(len(actual_runs), len(expected_runs)))
             if index >= len(actual_runs) or index >= len(expected_runs)
             or actual_runs[index] != expected_runs[index]),
            0,
        )
        expected_run = expected_runs[mismatch] if mismatch < len(expected_runs) else None
        actual_run = actual_runs[mismatch] if mismatch < len(actual_runs) else None
        expected_core = {key: expected_run.get(key) for key in ("range", "first_fragment", "last_fragment")} \
            if _is_record(expected_run) else None
        actual_core = {key: actual_run.get(key) for key in ("range", "first_fragment", "last_fragment")} \
            if _is_record(actual_run) else None
        code = "handoff-runs" if expected_core != actual_core else "handoff-anchors"
        _problem(code, run_address(expected_run, actual_run),
                 "The interaction hand-off runs or their anchor locks disagree with the current matrix and selected anchors.",
                 "Recompute the stage-5 interaction hand-off from the current matrix and selection.")

    expected_kept_residues = []
    seen_residues = set()
    for fragment in kept_fragments:
        for residue in fragment["residues"]:
            signature = json.dumps(residue, sort_keys=True, separators=(",", ":"))
            if signature not in seen_residues:
                seen_residues.add(signature)
                expected_kept_residues.append(dict(residue))
    if handoff["kept_residues"] != expected_kept_residues:
        differing = None
        actual_residues = handoff["kept_residues"]
        for index in range(max(len(expected_kept_residues), len(actual_residues))):
            expected_residue = expected_kept_residues[index] if index < len(expected_kept_residues) else None
            actual_residue = actual_residues[index] if index < len(actual_residues) else None
            if expected_residue != actual_residue:
                differing = expected_residue or actual_residue
                break
        address = "range-1"
        if _is_record(differing):
            for fragment in kept_fragments:
                if differing in fragment["residues"]:
                    address = fragment["key"]
                    break
        _problem("handoff-kept-residues", address,
                 "The interaction hand-off kept-residue list disagrees with the current clipped fragments.",
                 "Recompute the stage-5 interaction hand-off for the current clip.")
    return {
        "schema": 1, "rows": kept_keys, "ranges": ranges,
        "grid": [list(row) for row in grid], "handoff": handoff,
    }


def _validate_energy(energetics: Any, edits: Any, range_specs: list[dict],
                     range_specs_value: Any, fragments: Any, kept_fragments: list[dict],
                     interaction: dict, best_fits: Any) -> list[dict]:
    try:
        import preparic_stage5

        preparic_stage5.validate_energy_threshold_edits(edits, fragments, range_specs_value)
    except ValueError as exc:
        _problem("energy-edits", None, str(exc),
                 "Repair stage-5 energy thresholds for the current fragments and protein ranges.")
    if not _is_record(edits) or edits.get("basis") != fragments.get("echo"):
        _problem("energy-edit-echo", None, "Energy thresholds belong to a different RNA fragment basis.",
                 "Reopen stage 5 and edit thresholds for the current fragments.")
    _validate_schema1("the energetics result", energetics)
    if not _is_record(energetics) or not _is_record(energetics.get("echo")) \
        or not isinstance(energetics.get("ranges"), list):
        _problem("energetics", None, "Fresh best-fit energetics are unavailable.",
                 "Recompute stage 5 energetics for the current best fits.")
    if energetics["echo"].get("best_fits") != best_fits.get("echo") \
        or energetics["echo"].get("protein_ranges") != [dict(row) for row in range_specs[1:]]:
        _problem("energetics-echo", None, "Energetics do not echo the current fits and protein range rows.",
                 "Recompute stage 5 energetics for the current inputs.")
    by_range = {}
    for range_row in energetics["ranges"]:
        range_id = range_row.get("range") if _is_record(range_row) else None
        if not isinstance(range_id, str) or range_id in by_range or not isinstance(range_row.get("rows"), list):
            _problem("energetics-range", None, "An energetics range row is malformed or duplicated.",
                     "Recompute stage 5 energetics.")
        by_range[range_id] = range_row
    expected_ranges = interaction["ranges"]
    if set(by_range) != set(expected_ranges):
        _problem("energetics-range", None, "Energetics do not contain the current protein ranges.",
                 "Recompute stage 5 energetics.")
    edit_by_range = {}
    if not _is_record(edits) or type(edits.get("schema")) is not int or edits.get("schema") != 1 \
        or not isinstance(edits.get("ranges"), list):
        _problem("energy-edits", None, "Energy threshold edits must use schema 1.",
                 "Complete stage 5 energy thresholds for the current protein ranges.")
    protein_rows = range_specs[1:]
    for entry in edits["ranges"]:
        if not _is_record(entry) or set(entry) != {"row", "values"} \
            or not _is_record(entry.get("row")) or not _is_record(entry.get("values")):
            _problem("energy-edits", None, "An energy-threshold row is malformed.",
                     "Repair the stage-5 energy thresholds.")
        signature = _range_signature(entry["row"])
        if signature in edit_by_range:
            _problem("energy-edits", None, "Energy-threshold edits repeat a protein range row.",
                     "Repair the stage-5 energy thresholds.")
        edit_by_range[signature] = entry["values"]
    energy_out = []
    kept_keys = [row["key"] for row in kept_fragments]
    for row_index, fragment_key in enumerate(kept_keys):
        true_ranges = [range_id for col, range_id in enumerate(expected_ranges)
                       if interaction["grid"][row_index][col]]
        for range_id in true_ranges:
            numeric = int(range_id.split("-")[1])
            range_identity = protein_rows[numeric - 2]
            source = by_range[range_id]
            address = f"{fragment_key}/{range_id}"
            rows = {}
            for energy_row in source["rows"]:
                if not _is_record(energy_row) or not isinstance(energy_row.get("key"), str):
                    _problem("energy-result", address,
                             "An energetics range contains a malformed fragment address.",
                             "Recompute stage 5 energetics for the current fragment map.")
                if energy_row["key"] in rows:
                    _problem("energy-result", address,
                             f"Energetics for {range_id} repeat {energy_row['key']}.",
                             "Recompute stage 5 energetics without duplicate fragment rows.")
                rows[energy_row["key"]] = energy_row
            item = rows.get(fragment_key)
            bestfit = item.get("bestfit") if _is_record(item) else None
            energy = bestfit.get("energy") if _is_record(bestfit) else None
            error = bestfit.get("error") if _is_record(bestfit) else "missing best-fit energy"
            if error is not None or not _is_finite_number(energy):
                _problem("energy-result", address, f"The true interaction cell {range_id} has no finite best-fit energy.",
                         "Resolve the addressed energetics error and recompute stage 5.")
            edited = edit_by_range.get(_range_signature(range_identity), {}).get(fragment_key)
            if not isinstance(edited, str) or _ENERGY_DECIMAL.fullmatch(edited) is None:
                _problem("energy-threshold", address, f"The true interaction cell {range_id} has no valid edited threshold.",
                         "Fill the addressed energy threshold in stage 5.")
            try:
                threshold = Decimal(edited)
            except InvalidOperation:
                threshold = Decimal("NaN")
            if not threshold.is_finite():
                _problem("energy-threshold", address, f"The threshold for {range_id} is not finite.",
                         "Correct the addressed energy threshold in stage 5.")
            energy_out.append({
                "fragment": fragment_key, "range": range_id,
                "range_row": dict(range_identity), "energy": float(energy),
                "threshold": edited,
            })
    return energy_out


def _candidate_order(row: Any) -> tuple:
    if not _is_record(row):
        return ("", 0, 0)
    nucleotide = {"first": 1, "second": 2}.get(row.get("nucleotide"), 0)
    fragment = row.get("fragment") if _is_int(row.get("fragment")) else 0
    return (str(row.get("key", "")), fragment, nucleotide)


def _canonical_class_table(rows: list[dict]) -> list[dict]:
    return [dict(row) for row in sorted(rows, key=lambda row: row["key"])]


def _canonical_validation(value: dict) -> dict:
    result = dict(value)
    echo = dict(value["echo"])
    echo["classes"] = _canonical_class_table(echo["classes"])
    result["echo"] = echo
    rows = []
    for source in value["rows"]:
        row = dict(source)
        if isinstance(row.get("candidates"), list):
            row["candidates"] = sorted((dict(item) for item in row["candidates"] if _is_record(item)),
                                        key=_candidate_order)
        rows.append(row)
    result["rows"] = sorted(rows, key=lambda row: row.get("key", ""))
    return result


def _canonical_gaussian_scores(value: dict) -> dict:
    result = dict(value)
    for field in ("classes", "candidate_scores", "reference_scores"):
        rows = value.get(field)
        if isinstance(rows, list):
            result[field] = sorted((dict(row) for row in rows if _is_record(row)), key=_candidate_order)
    return result


def _canonical_bestfit_stackings(value: dict) -> dict:
    result = dict(value)
    result["rows"] = sorted((dict(row) for row in value["rows"] if _is_record(row)), key=_candidate_order)
    return result


def _validate_gaussian(inputs: dict, selected_keys: list[str], best_fits: Any,
                       clip: Any, range_specs: list[dict], model: dict) -> dict:
    class_marker = inputs.get("gaussian_classes")
    if not _is_record(class_marker) or type(class_marker.get("schema")) is not int \
        or class_marker.get("schema") != 1 or class_marker.get("active") is not True \
        or not isinstance(class_marker.get("rows"), list) or class_marker.get("basis") != selected_keys:
        _problem("gaussian-classes", None, "Gaussian classes are not active schema 1 for the selected stackings.",
                 "Complete stage 6 with one class for every selected stacking.")
    class_rows = class_marker["rows"]
    class_keys = [row.get("key") for row in class_rows if _is_record(row) and isinstance(row.get("key"), str)]
    if len(class_keys) != len(class_rows) or len(class_keys) != len(set(class_keys)) \
        or sorted(class_keys) != selected_keys:
        address = next((key for key in class_keys if class_keys.count(key) > 1), None)
        if address is None:
            address = next((key for key in selected_keys if key not in class_keys), None)
        _problem("gaussian-class-keys", address,
                 "Gaussian class keys must equal the selected stacking keys exactly once.",
                 "Complete exactly one class row for every selected stacking.")
    try:
        import preparic_stage6
        import alaric_gaussian_filter

        preparic_stage6.validate_gaussian_classes(
            class_marker, inputs.get("gaussian_model_json"),
            inputs.get("stacking_selection"), alaric_gaussian_filter,
        )
    except ValueError as exc:
        _problem("gaussian-class", None, str(exc),
                 "Repair stage-6 classes using a filter in the current Gaussian model.")
    class_table = [dict(row) for row in class_rows]
    model_names = set(model["thresholds"])
    validation = inputs.get("gaussian_validation")
    if not _is_record(validation) or type(validation.get("schema")) is not int \
        or validation.get("schema") != 1 or not isinstance(validation.get("rows"), list):
        _problem("gaussian-validation", None, "Fresh Gaussian validation results are unavailable.",
                 "Recompute the stage-6 Gaussian validation.")
    expected_validation_echo = {
        "classes": _canonical_class_table(class_rows), "clip": clip, "selected_keys": selected_keys,
        "best_fits": best_fits.get("echo"),
    }
    validation_echo = validation.get("echo")
    canonical_validation_echo = dict(validation_echo) if _is_record(validation_echo) else None
    if canonical_validation_echo is not None and isinstance(canonical_validation_echo.get("classes"), list):
        canonical_validation_echo["classes"] = sorted(
            (dict(row) for row in canonical_validation_echo["classes"] if _is_record(row)),
            key=lambda row: row.get("key", ""),
        )
    if canonical_validation_echo != expected_validation_echo:
        _problem("gaussian-validation-echo", None, "Gaussian validation does not echo the current classes, clip, selections, and fits.",
                 "Recompute stage 6 validation for the current decisions.")
    if validation.get("filter_set_version") != model.get("filter_set_version") \
        or not isinstance(validation.get("filter_names"), list) \
        or not all(type(name) is str for name in validation["filter_names"]) \
        or len(validation["filter_names"]) != len(set(validation["filter_names"])) \
        or set(validation["filter_names"]) != model_names:
        _problem("gaussian-validation-version", None, "Gaussian validation uses a different filter set from the current model.",
                 "Recompute stage 6 with the current Gaussian model.")
    check = inputs.get("gaussian_check")
    if not _is_record(check) or type(check.get("schema")) is not int or check.get("schema") != 1 \
        or not _is_finite_number(check.get("min_cell_fraction")) \
        or float(check.get("min_cell_fraction")) != ANCHOR_MIN_CELL_FRACTION:
        _problem("gaussian-cell-fraction", None,
                 "The Gaussian check minimum cell fraction differs from the canonical one-third default.",
                 "Run the Gaussian check with the canonical minimum cell fraction and recompute stage 6.")
    if check.get("filter_set_version") != model.get("filter_set_version"):
        _problem("gaussian-check-version", None, "The Gaussian check uses a different model version.",
                 "Recompute the Gaussian check with the current model.")
    scores = inputs.get("gaussian_scores")
    bestfit_stackings = inputs.get("bestfit_stackings")
    if not _is_record(bestfit_stackings) or type(bestfit_stackings.get("schema")) is not int \
        or bestfit_stackings.get("schema") != 1 or not isinstance(bestfit_stackings.get("rows"), list):
        _problem("bestfit-stackings", None, "Fresh best-fit stacking evidence is unavailable.",
                 "Recompute the stage-6 best-fit stacking detector.")
    if not _is_record(scores) or type(scores.get("schema")) is not int or scores.get("schema") != 1 \
        or not isinstance(scores.get("classes"), list) \
        or not isinstance(scores.get("candidate_scores"), list) \
        or not isinstance(scores.get("reference_scores"), list):
        _problem("gaussian-scores", None, "Fresh Gaussian score evidence is unavailable.",
                 "Recompute the stage-6 Gaussian scorer.")
    expected_scores_echo = {
        "bestfit_stackings": bestfit_stackings.get("echo"),
        "selected_keys": selected_keys, "clip": clip,
    }
    if scores.get("echo") != expected_scores_echo:
        _problem("gaussian-scores-echo", None, "Gaussian scores do not echo the current descriptors, selections, and clip.",
                 "Recompute stage 6 Gaussian scores.")
    bestfit_echo = bestfit_stackings.get("echo")
    if not _is_record(bestfit_echo) or bestfit_echo.get("best_fits") != best_fits.get("echo") \
        or bestfit_echo.get("protein_ranges") != [dict(row) for row in range_specs[1:]]:
        _problem("bestfit-stackings-echo", None, "Best-fit stacking evidence is stale for the current fits or protein ranges.",
                 "Recompute the stage-6 best-fit stacking detector.")
    if scores.get("filter_set_version") != model.get("filter_set_version") \
        or not isinstance(scores.get("filter_names"), list) \
        or not all(type(name) is str for name in scores["filter_names"]) \
        or len(scores["filter_names"]) != len(set(scores["filter_names"])) \
        or set(scores["filter_names"]) != model_names:
        _problem("gaussian-scores-version", None, "Gaussian score evidence uses a different model version.",
                 "Recompute stage 6 with the current Gaussian model.")
    score_classes = scores["classes"]
    score_class_keys = []
    for row in score_classes:
        if not _is_record(row) or not isinstance(row.get("key"), str) \
            or type(row.get("below_rrm95")) is not bool:
            _problem("gaussian-score-class", row.get("key") if _is_record(row) else None,
                     "A computed Gaussian class row is malformed.",
                     "Recompute stage 6 Gaussian scores.")
        score_class_keys.append(row["key"])
    if len(score_class_keys) != len(set(score_class_keys)) or set(score_class_keys) != set(selected_keys):
        missing = next((key for key in selected_keys if key not in score_class_keys), None)
        address = missing or next((key for key in score_class_keys if score_class_keys.count(key) > 1), None)
        _problem("gaussian-score-classes", address,
                 "Computed Gaussian classes must contain exactly one row per selected stacking.",
                 "Recompute stage 6 Gaussian scores for every selected stacking.")
    validation_keys = []
    for row in validation["rows"]:
        if not _is_record(row) or not isinstance(row.get("key"), str):
            _problem("gaussian-validation-row", None, "A Gaussian validation row has no class address.",
                     "Recompute stage 6 validation.")
        key = row["key"]
        if key not in selected_keys:
            _problem("gaussian-validation-row", key, "A Gaussian validation row has an unknown class address.",
                     "Recompute stage 6 validation for the current selected stackings.")
        if key in validation_keys:
            _problem("gaussian-validation-duplicate", key,
                     "Gaussian validation contains more than one row for this selected class.",
                     "Recompute stage 6 validation so each selected class has exactly one row.")
        validation_keys.append(key)
    if set(validation_keys) != set(selected_keys) or len(validation_keys) != len(selected_keys):
        missing = next((key for key in selected_keys if key not in validation_keys), None)
        _problem("gaussian-validation-missing", missing,
                 "Gaussian validation must contain exactly one row for every selected class.",
                 "Recompute stage 6 validation for all current selected stackings.")
    for row in validation["rows"]:
        if row.get("worst_state") not in {"green", "implied", "amber", "red", "error"}:
            _problem("gaussian-validation-state", row["key"], "A Gaussian validation row has an unknown result state.",
                     "Recompute stage 6 validation.")
        if row.get("worst_state") == "error":
            _problem("gaussian-job", row["key"], "The Gaussian check contains an addressed candidate error.",
                     "Resolve the Gaussian check error and recompute stage 6.")
    return {
        "class_table": _canonical_class_table(class_table),
        "validation": _canonical_validation(validation),
        "scores": _canonical_gaussian_scores(scores),
        "bestfit_stackings": _canonical_bestfit_stackings(bestfit_stackings),
        "check": check,
        "model": inputs.get("gaussian_model_json"),
        "library": inputs.get("fraglib_identities"),
        "rings": {name: inputs.get(name) for name in _RING_PINS},
        "min_cell_fraction": float(check["min_cell_fraction"]),
    }


def _validate_clipped_pdbs(clipped: Any, kept_fragments: list[dict], best_fits: Any) -> dict:
    best_fit_rows = best_fits.get("rows") if _is_record(best_fits) else None
    if not isinstance(best_fit_rows, list):
        _problem("best-fits", None, "Best-fit rows are unavailable.",
                 "Recompute best fits for the current prepared RNA.")
    fit_by_key = {row.get("key"): row for row in best_fit_rows if _is_record(row)}
    expected_paths = {"thresholds.json", "clip.json"}
    result = {}
    for fragment in kept_fragments:
        key = fragment["key"]
        number, sequence = fragment["number"], fragment["sequence"]
        ref_name = f"reference/frag-{number}-{sequence}-refe.pdb"
        fit_name = f"bestfit/frag-{number}-{sequence}-bestfit.pdb"
        for name in (ref_name, fit_name):
            value = clipped.get(name) if _is_record(clipped) else None
            if not isinstance(value, bytes):
                _problem("clipped-pdb", key, f"The kept-fragment folder is missing {name}.",
                         "Recompute the kept-fragment folder from current stage-4 results.")
            result[name] = value
        fitted = fit_by_key.get(key)
        if not _is_record(fitted) or not isinstance(fitted.get("pdb"), str) \
            or result[fit_name] != fitted["pdb"].encode("utf-8"):
            _problem("bestfit-pdb", key, "The clipped best-fit bytes do not match the current fitted pose.",
                     "Recompute the kept-fragment folder from current best-fit results.")
        expected_paths.update((ref_name, fit_name))
    if set(clipped) != expected_paths:
        _problem("clipped-pdb", None, "The kept-fragment folder contains missing, clipped, or unknown members.",
                 "Recompute the kept-fragment folder for the current clip.")
    result["thresholds.json"] = clipped["thresholds.json"]
    result["clip.json"] = clipped["clip.json"]
    return result


def _stage6_gate(inputs: dict) -> list[dict]:
    """Apply project-state's persisted-decision stage prefix before deep checks."""
    def not_ready(stage: int, address: str | None = None) -> list[dict]:
        return [_block("stage-gate", address,
                       f"Stage {stage} decisions are not ready for export.",
                       f"Complete stage {stage} with decisions for the current upstream inputs.")]

    pdbcode, biounit = inputs.get("pdbcode"), inputs.get("biounit")
    if not isinstance(pdbcode, str) or not pdbcode.strip() or not _is_int(biounit) or biounit < 1:
        return not_ready(1)
    code = pdbcode.strip().upper()
    specs = inputs.get("range_specs")
    if not _is_record(specs) or specs.get("schema") != 1 \
        or specs.get("pdbcode") != code or specs.get("biounit") != biounit \
        or not isinstance(specs.get("rows"), list) or not specs["rows"]:
        return not_ready(2, "range-1")
    range_rows = specs["rows"]
    if any(not _is_record(row) or set(row) != {"chain", "start", "end"}
           or not isinstance(row.get("chain"), str) or not row["chain"]
           or not _is_int(row.get("start")) or not _is_int(row.get("end"))
           or row["end"] < row["start"] for row in range_rows):
        return not_ready(2, "range-1")
    ranges_marker = inputs.get("ranges")
    if ranges_marker is None or (isinstance(ranges_marker, str) and not ranges_marker.strip()):
        return not_ready(2, "range-1")
    if _is_record(ranges_marker) and (
        set(ranges_marker) != {"schema", "basis", "active", "rows"}
        or ranges_marker.get("schema") != 1 or ranges_marker.get("active") is not True
        or ranges_marker.get("basis") != {"pdbcode": code, "biounit": biounit}
        or ranges_marker.get("rows") != range_rows
    ):
        return not_ready(2, "range-1")

    selection = inputs.get("stacking_selection")
    if not _is_record(selection) or set(selection) != {"schema", "basis", "active", "selected"} \
        or selection.get("schema") != 1 or selection.get("active") is not True \
        or not _is_record(selection.get("selected")) \
        or any(type(value) is not bool for value in selection["selected"].values()) \
        or any(not isinstance(key, str) or not key for key in selection["selected"]) \
        or selection.get("basis") != {"pdbcode": code, "biounit": biounit, "rows": range_rows}:
        return not_ready(3)
    selected = sorted(key for key, value in selection["selected"].items() if value)

    rna_range = range_rows[0]
    expected_rna_basis = {"pdbcode": code, "rna": rna_range}
    thresholds = inputs.get("threshold_edits")
    if not _is_record(thresholds) or set(thresholds) != {"schema", "basis", "active", "pairs"} \
        or thresholds.get("schema") != 1 or thresholds.get("active") is not True \
        or not isinstance(thresholds.get("pairs"), list) \
        or thresholds.get("basis") != expected_rna_basis:
        return not_ready(4, "range-1")
    clip = inputs.get("best_fit_clip")
    if not _is_record(clip) or set(clip) != {"schema", "basis", "clip1", "clip2"} \
        or clip.get("schema") != 1 or not _is_int(clip.get("clip1")) \
        or not _is_int(clip.get("clip2")) or clip["clip1"] < 0 or clip["clip2"] < 0 \
        or clip.get("basis") != expected_rna_basis:
        return not_ready(4, "range-1")
    pairs = thresholds["pairs"]
    if not pairs:
        fragments = inputs.get("fragments")
        rows = fragments.get("fragments") if _is_record(fragments) else None
        if not isinstance(rows, list) or len(rows) != 1 or not isinstance(rows[0].get("key"), str):
            return not_ready(4, "range-1")
        all_keys = [rows[0]["key"]]
    else:
        all_keys = []
        for index, pair in enumerate(pairs):
            if not _is_record(pair) or pair.get("down") != index + 1 \
                or pair.get("up") != index + 2 or not isinstance(pair.get("sequence"), str) \
                or re.fullmatch(r"[ACGU]{3}", pair["sequence"]) is None:
                return not_ready(4, "range-1")
            if index == 0:
                all_keys.append(f"frag-1-{pair['sequence'][:2]}")
            elif pair["sequence"][:2] != pairs[index - 1].get("sequence", "")[1:3]:
                return not_ready(4, f"frag-{index + 1}-{pair['sequence'][:2]}")
            else:
                all_keys.append(f"frag-{index + 1}-{pair['sequence'][:2]}")
        all_keys.append(f"frag-{len(pairs) + 1}-{pairs[-1]['sequence'][1:3]}")
    if clip["clip1"] + clip["clip2"] >= len(all_keys):
        return not_ready(4, "range-1")
    expected_matrix_rows = all_keys[clip["clip1"]:len(all_keys) - clip["clip2"]]
    expected_ranges = [_range_id(index) for index in range(2, len(range_rows) + 1)]
    matrix = inputs.get("interaction_matrix")
    if not _is_record(matrix) or set(matrix) != {"schema", "basis", "active", "rows", "ranges", "grid"} \
        or matrix.get("schema") != 1 or matrix.get("active") is not True \
        or matrix.get("basis") != selected or matrix.get("rows") != expected_matrix_rows \
        or matrix.get("ranges") != expected_ranges or not isinstance(matrix.get("grid"), list) \
        or len(matrix["grid"]) != len(expected_matrix_rows) \
        or any(not isinstance(row, list) or len(row) != len(expected_ranges)
               or any(type(cell) is not bool for cell in row) for row in matrix["grid"]):
        return not_ready(5, expected_matrix_rows[0] if expected_matrix_rows else "range-1")

    classes = inputs.get("gaussian_classes")
    if not _is_record(classes) or set(classes) != {"schema", "basis", "active", "rows"} \
        or classes.get("schema") != 1 or classes.get("active") is not True \
        or classes.get("basis") != selected or not isinstance(classes.get("rows"), list):
        return not_ready(6, selected[0] if selected else None)
    class_keys = [row.get("key") for row in classes["rows"] if _is_record(row)]
    if len(class_keys) != len(classes["rows"]) \
        or any(not isinstance(key, str) for key in class_keys) or sorted(class_keys) != selected \
        or len(class_keys) != len(set(class_keys)):
        return not_ready(6, next((key for key in selected if key not in class_keys), None))
    return []


def _build_snapshot(inputs: dict) -> tuple[dict | None, list[dict]]:
    blocks: list[dict] = []
    required = (
        "pdbcode", "biounit", "cif", "ranges", "range_specs", "structure_id",
        "exclude_pdb_code", "range_pdb_files", "range_pdb_files_aa", "range_pdbs_aa", "rna_pdb",
        "fragments", "best_fit_thresholds", "best_fits", "clipped_fragments",
        "stackings", "stacking_selection", "threshold_edits", "best_fit_clip",
        "interaction_matrix", "interactions", "energetics", "energy_threshold_edits",
        "bestfit_stackings", "gaussian_scores", "gaussian_validation", "gaussian_check",
        "gaussian_model_json", "fraglib_identities", *_RING_PINS, *_VALIDATOR_PINS,
    )
    absent = [name for name in required if name not in inputs]
    if absent:
        return None, [_block("missing-input", None, f"Required export input {absent[0]} is missing.",
                             "Wire every named export input and recompute admission.")]
    for name in _VALIDATOR_PINS:
        if inputs[name] is not None:
            blocks.append(_block("validator-result", None,
                                 f"The {name} validator did not return its successful null sentinel.",
                                 "Resolve the upstream validator failure before export."))

    pdbcode = inputs.get("pdbcode")
    biounit = inputs.get("biounit")
    range_rows, rna_range = _validate_ranges(pdbcode, biounit, inputs.get("ranges"),
                                              inputs.get("range_specs"), inputs.get("structure_id"))
    code_upper = pdbcode.strip().upper()
    code_lower = pdbcode.strip().lower()
    range_specs = {"schema": 1, "pdbcode": code_upper, "biounit": biounit,
                   "rows": [dict(row) for row in range_rows]}
    fragments_value = inputs.get("fragments")
    fragment_rows, fragment_map, expected_pairs = _validate_fragments(fragments_value)
    expected_fragment_echo = {"pdbcode": code_upper, "rna": rna_range}
    if fragments_value.get("echo") != expected_fragment_echo:
        _problem("fragment-echo", "range-1", "The fragment list is stale for the current structure and RNA range.",
                 "Recompute the canonical fragment list for the current RNA.")
    fit_echo = inputs.get("best_fits", {}).get("echo") if _is_record(inputs.get("best_fits")) else None
    b1, b3, reference_bytes = _validate_file_folders(inputs, range_rows, pdbcode, rna_range, fit_echo)
    expected_fit_echo = {
        "pdb_code": code_lower, "rna": rna_range,
        "resids": fragments_value.get("residues"),
    }
    if fit_echo != expected_fit_echo:
        _problem("fit-echo", "range-1", "Best-fit results do not match the current RNA residue identities.",
                 "Recompute best fits from the current prepared RNA.")
    _validate_schema1("computed best-fit thresholds", inputs["best_fit_thresholds"])
    if not isinstance(inputs["best_fit_thresholds"].get("pairs"), list) \
        or not isinstance(inputs["best_fit_thresholds"].get("bases"), list):
        _problem("thresholds", None, "Computed best-fit pair/base thresholds are incomplete.",
                 "Recompute stage 4 thresholds from the current best-fit tables.")
    if inputs["best_fit_thresholds"].get("echo") != fit_echo:
        _problem("threshold-echo", None, "Computed thresholds are stale for the current best fits.",
                 "Recompute stage 4 thresholds.")
    _validate_schema1("best-fit results", inputs["best_fits"])
    if not isinstance(inputs["best_fits"].get("rows"), list):
        _problem("best-fits", None, "Best-fit result rows are unavailable.",
                 "Recompute best fits for the current prepared RNA.")

    clip = inputs.get("best_fit_clip")
    expected_clip_basis = {"pdbcode": code_upper, "rna": rna_range}
    if not _is_record(clip) or clip.get("basis") != expected_clip_basis:
        _problem("clip-echo", None, "The clip does not echo the current structure and RNA range.",
                 "Reopen stage 4 and set a current clip.")
    kept_fragments, kept_numbers = _validate_clip(clip, fragment_rows, fragment_map, expected_clip_basis)
    kept_fragment_rows = [dict(row) for row in kept_fragments]
    selected = _selected_keys(inputs.get("stacking_selection"))
    selection_basis = {"pdbcode": code_upper, "biounit": biounit, "rows": range_rows}
    if inputs["stacking_selection"].get("basis") != selection_basis:
        _problem("selection-echo", None, "The stacking selection belongs to different structure or ranges.",
                 "Reopen stage 3 and restore a current stacking selection.")
    anchors, selected_stackings = _candidate_anchors(
        fragments_value, kept_fragments, inputs["stacking_selection"], inputs["stackings"],
        inputs.get("interaction_matrix", {}).get("ranges", []) if _is_record(inputs.get("interaction_matrix")) else [],
        {
            "structure_id": inputs["structure_id"], "pdbcode": code_upper,
            "biounit": biounit, "rows": range_rows,
        },
    )
    if inputs["threshold_edits"].get("basis") != expected_fragment_echo:
        _problem("threshold-echo", None, "Pair thresholds belong to a different RNA range.",
                 "Reopen stage 4 and recompute pair thresholds for the current RNA.")
    if inputs["threshold_edits"].get("active") is not True:
        _problem("threshold-edits", None, "Pair thresholds are not active for export.",
                 "Complete stage 4 pair thresholds.")
    try:
        import preparic_stage4

        preparic_stage4._validate_threshold_edits(inputs["threshold_edits"], fragments_value)
    except ValueError as exc:
        _problem("threshold-edits", None, str(exc),
                 "Repair stage-4 pair thresholds for the current fragment list.")
    kept_pairs = _validate_pair_data(fragment_rows, kept_numbers, expected_pairs,
                                     inputs["best_fit_thresholds"], inputs["threshold_edits"],
                                     inputs["clipped_fragments"])
    clipped_json = json.loads(inputs["clipped_fragments"]["clip.json"])
    if clipped_json.get("clip1") != clip.get("clip1") or clipped_json.get("clip2") != clip.get("clip2"):
        _problem("clip-boundary", None, "The kept-fragment boundary does not echo the current clip.",
                 "Recompute the kept-fragment boundary for the current clip.")
    expected_kept_residues = []
    seen_kept_residues = set()
    for fragment in kept_fragments:
        for residue in fragment["residues"]:
            identity = json.dumps(residue, sort_keys=True, separators=(",", ":"))
            if identity not in seen_kept_residues:
                seen_kept_residues.add(identity)
                expected_kept_residues.append(dict(residue))
    if clipped_json.get("kept_residues") != expected_kept_residues:
        _problem("clip-boundary-residues", "range-1",
                 "The kept-fragment clip record has stale RNA residue identities.",
                 "Recompute the kept-fragment boundary for the current fragment clip.")
    clipped_files = _validate_clipped_pdbs(inputs["clipped_fragments"], kept_fragments, inputs["best_fits"])

    model = _parse_model(inputs.get("gaussian_model_json"))
    class_rows = inputs.get("gaussian_classes", {}).get("rows", []) if _is_record(inputs.get("gaussian_classes")) else []
    class_keys = sorted(row.get("key") for row in class_rows if _is_record(row) and isinstance(row.get("key"), str))
    if class_keys != selected:
        _problem("gaussian-class-keys", None, "Gaussian class keys must equal the selected stacking keys.",
                 "Complete exactly one class row for every selected stacking.")
    if inputs["gaussian_classes"].get("basis") != selected:
        _problem("gaussian-class-basis", None, "Gaussian classes belong to a different stacking selection.",
                 "Reopen stage 6 and set classes for the current stacking selection.")

    interaction = _validate_interactions(
        inputs["interaction_matrix"], inputs["interactions"], fragment_rows,
        kept_fragments, anchors, range_rows, clip, selected,
    )
    energy_rows = _validate_energy(
        inputs["energetics"], inputs["energy_threshold_edits"], range_rows,
        inputs["range_specs"], fragments_value, kept_fragments, interaction, inputs["best_fits"],
    )
    gaussian = _validate_gaussian(inputs, selected, inputs["best_fits"], clip, range_rows, model)

    expected_selected_keys = sorted(row["key"] for row in selected_stackings)
    if expected_selected_keys != selected:
        _problem("selected-stackings", None, "Selected stacking identities could not be resolved exactly.",
                 "Recompute stage 3 and restore the current selected stackings.")

    snapshot = {
        "schema": 1,
        "structure": {"pdbcode": code_lower, "biounit": biounit, "structure_id": code_upper},
        "range_specs": range_specs,
        "exclude_pdb_code": code_lower,
        "range_pdb_files": b1,
        "range_pdb_files_aa": b3,
        "reference_pdb": reference_bytes,
        "reference_sha256": hashlib.sha256(reference_bytes).hexdigest(),
        "fragments": [dict(row) for row in fragment_rows],
        "kept_fragments": kept_fragment_rows,
        "clip": {"clip1": clip["clip1"], "clip2": clip["clip2"]},
        "kept_pairs": kept_pairs,
        "threshold_bases": [dict(row) for row in inputs["best_fit_thresholds"].get("bases", [])],
        "stacking_selection": selected,
        "stackings": selected_stackings,
        "candidate_anchors": anchors,
        "interaction_matrix": interaction,
        "handoff": inputs["interactions"],
        "energy_cells": energy_rows,
        "clipped_fragments": clipped_files,
        "class_table": gaussian["class_table"],
        "gaussian_validation": gaussian["validation"],
        "gaussian_scores": gaussian["scores"],
        "bestfit_stackings": gaussian["bestfit_stackings"],
        "gaussian_model_json": gaussian["model"],
        "fraglib_identities": gaussian["library"],
        "ring_references": gaussian["rings"],
        "minimum_cell_fraction": gaussian["min_cell_fraction"],
    }
    return snapshot, blocks


def normalise_export_inputs(inputs: Any) -> tuple[dict | None, list[dict]]:
    """Normalise named export pins and return addressed hard blocks."""
    if not _is_record(inputs):
        return None, [_block("missing-input", None, "Export inputs must be a named mapping.",
                             "Wire the export input pins and recompute admission.")]
    stage_blocks = _stage6_gate(inputs)
    if stage_blocks:
        return None, stage_blocks
    try:
        return _build_snapshot(inputs)
    except _InputProblem as exc:
        return None, [_block(exc.code, exc.address, exc.message, exc.remedy)]


def export_warnings(snapshot: Any) -> list[dict]:
    """Return ordered X5 warning reasons for a valid normalised snapshot."""
    if not _is_record(snapshot) or snapshot.get("schema") != 1:
        raise ValueError("export_warnings requires a valid schema-1 ExportSnapshot.")
    warnings = []
    validation = snapshot.get("gaussian_validation", {})
    gaussian_reasons: dict[str, list[str]] = {}
    for row in validation.get("rows", []) if _is_record(validation) else []:
        if _is_record(row) and row.get("worst_state") in {"implied", "amber", "red"} \
            and isinstance(row.get("key"), str):
            gaussian_reasons.setdefault(row["key"], []).append(
                f"the selected Gaussian filter has {row['worst_state']} validation results"
            )
    scores = snapshot.get("gaussian_scores", {})
    for row in scores.get("classes", []) if _is_record(scores) else []:
        if _is_record(row) and row.get("below_rrm95") is True and isinstance(row.get("key"), str):
            gaussian_reasons.setdefault(row["key"], []).append(
                "the computed Gaussian class is below RRM-95"
            )
    for key in sorted(gaussian_reasons):
        warnings.append({
            "kind": "gaussian", "address": key,
            "message": f"The selected Gaussian validation for {key} does not pass: "
                       + "; ".join(gaussian_reasons[key]) + ".",
        })
    for row in snapshot.get("kept_pairs", []):
        if not _is_record(row):
            continue
        address = f"frag-{row.get('down')}..frag-{row.get('up')}"
        for name in ("crmsd", "ovrmsd"):
            raw = row.get(f"{name}_raw")
            threshold = row.get(name)
            if isinstance(threshold, str) and _is_finite_number(raw):
                try:
                    if float(Decimal(threshold)) - float(raw) < 1e-5:
                        warnings.append({
                            "kind": "pair-boundary", "address": address,
                            "message": (f"The edited pair {name} threshold {threshold} is within 1e-5 "
                                        f"of raw value {raw}."),
                        })
                except (InvalidOperation, ValueError):
                    continue
    for row in snapshot.get("energy_cells", []):
        if not _is_record(row):
            continue
        try:
            if float(row["energy"]) >= float(Decimal(row["threshold"])):
                warnings.append({
                    "kind": "energy", "address": f"{row.get('fragment')}/{row.get('range')}",
                    "message": f"The best-fit energy for {row.get('fragment')} in {row.get('range')} does not pass the edited threshold.",
                })
        except (KeyError, InvalidOperation, TypeError, ValueError):
            continue
    return warnings


def build_admission(inputs: Any) -> dict:
    """Build the schema-1 admission report from named export pin values."""
    snapshot, blocks = normalise_export_inputs(inputs)
    if snapshot is None and not blocks:
        blocks = [_block("missing-snapshot", None,
                         "A normalised export snapshot could not be constructed.",
                         "Resolve the export input errors and recompute admission.")]
    if snapshot is None:
        return {"schema": 1, "ready": False, "blocks": blocks, "warnings": [], "digest": None}
    warnings = export_warnings(snapshot) if not blocks else []
    return {
        "schema": 1, "ready": not blocks, "blocks": blocks,
        "warnings": warnings, "digest": snapshot_digest(snapshot),
    }
