"""Run Preparic's PDB, biounit, CIF, and range shares during development.

Use ``python -i backend.py`` to leave the context available in an interactive
prompt; call ``close_context()`` when you are done.
"""

import atexit
import copy
import json
import math
import os
import re
import signal
import sys
import threading
from pathlib import Path

import seamless
import yaml
from seamless import Cell
from seamless_transformer import Transformer
from seamless_workflow import Context, shareserver
from alaric_dir import resolve_alaric_dir
from alaric_runtime import (
    best_fit_factory,
    best_fit_tables_transformer,
    energy_threshold_module_code,
    energetics_range_transformer,
    gaussian_check_transformer,
    fragment_library_identities,
    refe_fragments_transformer,
    threshold_module_code,
    thresholds_from_tables,
)
from range_pdbs import _pdbs_to_folder, _ranges_to_pdb
from detect_stacking import _detect_stacking


PROJECT_ROOT = Path(__file__).resolve().parent
MOUNT_DIRECTORY = PROJECT_ROOT / "deployment-mounts"
PDBCODE_PATH = MOUNT_DIRECTORY / "pdbcode.txt"
BIOUNIT_PATH = MOUNT_DIRECTORY / "biounit.json"
RANGES_PATH = MOUNT_DIRECTORY / "ranges.yaml"
CIF_PATH = MOUNT_DIRECTORY / "cif.cif"
STACKING_SELECTION_PATH = MOUNT_DIRECTORY / "stacking-selection.json"
THRESHOLD_EDITS_PATH = MOUNT_DIRECTORY / "best-fit-threshold-edits.json"
BEST_FIT_CLIP_PATH = MOUNT_DIRECTORY / "best-fit-clip.json"
INTERACTION_MATRIX_PATH = MOUNT_DIRECTORY / "interaction-matrix.json"
GAUSSIAN_CLASSES_PATH = MOUNT_DIRECTORY / "gaussian-classes.json"
ENERGY_THRESHOLD_EDITS_PATH = MOUNT_DIRECTORY / "energy-threshold-edits.json"
PDBS_PATH = MOUNT_DIRECTORY / "pdbs"
PDBS_AA_PATH = MOUNT_DIRECTORY / "pdbs-aa"
CLIPPED_FRAGMENTS_PATH = MOUNT_DIRECTORY / "clipped-fragments"

LOCAL_WORKER_COUNT = 8
_LOCAL_EXECUTION_CONFIGURED = False
_STAGE5_RANGE_SYNC_STOP: threading.Event | None = None
_STAGE5_RANGE_SYNC_THREAD: threading.Thread | None = None

_RANGE_KEYS = {"chain", "start", "end"}
_MARKER_KEYS = {
    "ranges": {"schema", "basis", "active", "rows"},
    "stacking_selection": {"schema", "basis", "active", "selected"},
    "threshold_edits": {"schema", "basis", "active", "pairs", "bases"},
    "best_fit_clip": {"schema", "basis", "clip1", "clip2"},
    "interaction_matrix": {"schema", "basis", "active", "rows", "ranges", "grid"},
    "gaussian_classes": {"schema", "basis", "active", "rows"},
    "energy_threshold_edits": {"schema", "basis", "ranges"},
}
_NO_FALLBACK = object()


def _exact_keys(value: object, keys: set[str]) -> bool:
    return type(value) is dict and set(value) == keys


def _valid_range_row(value: object) -> bool:
    return _exact_keys(value, _RANGE_KEYS) and isinstance(value["chain"], str) and bool(value["chain"]) \
        and type(value["start"]) is int and type(value["end"]) is int and value["end"] >= value["start"]


def _unique(values: list[object]) -> bool:
    try:
        return len(values) == len(set(values))
    except TypeError:
        return False


def _valid_project_basis(value: object, with_rows: bool = False) -> bool:
    keys = {"pdbcode", "biounit", "rows"} if with_rows else {"pdbcode", "biounit"}
    return _exact_keys(value, keys) and isinstance(value["pdbcode"], str) and bool(value["pdbcode"]) \
        and type(value["biounit"]) is int and value["biounit"] >= 1 \
        and (not with_rows or (isinstance(value["rows"], list) and len(value["rows"]) > 0
            and all(_valid_range_row(row) for row in value["rows"])
            and _unique([json.dumps(row, sort_keys=True) for row in value["rows"]])))


def _valid_rna_basis(value: object) -> bool:
    return _exact_keys(value, {"pdbcode", "rna"}) and isinstance(value["pdbcode"], str) \
        and bool(value["pdbcode"]) and _valid_range_row(value["rna"])


def _valid_decision(name: str, value: object) -> bool:
    if value is None:
        return name != "best_fit_clip"
    if not _exact_keys(value, _MARKER_KEYS[name]) or type(value.get("schema")) is not int or value["schema"] != 1:
        return False
    if name == "ranges":
        return _valid_project_basis(value["basis"]) and type(value["active"]) is bool \
            and isinstance(value["rows"], list) and len(value["rows"]) > 0 \
            and all(_valid_range_row(row) for row in value["rows"]) \
            and _unique([json.dumps(row, sort_keys=True) for row in value["rows"]])
    if name == "stacking_selection":
        selected = value["selected"]
        return _valid_project_basis(value["basis"], with_rows=True) and type(value["active"]) is bool \
            and type(selected) is dict and all(isinstance(key, str) and key for key in selected) \
            and all(type(flag) is bool for flag in selected.values())
    if name == "threshold_edits":
        if not _valid_rna_basis(value["basis"]) or type(value["active"]) is not bool \
            or type(value["pairs"]) is not list or type(value["bases"]) is not list:
            return False
        pair_ids = []
        for row in value["pairs"]:
            if not _exact_keys(row, {"down", "up", "sequence", "ovrmsd", "crmsd"}) \
                or type(row["down"]) is not int or row["down"] < 1 or type(row["up"]) is not int \
                or row["up"] != row["down"] + 1 or not isinstance(row["sequence"], str) \
                or re.fullmatch(r"[ACGU]{3}", row["sequence"]) is None:
                return False
            for key in ("ovrmsd", "crmsd"):
                if row[key] is not None and (not isinstance(row[key], str)
                    or re.fullmatch(r"(?:[1-9][0-9]*(?:\.[0-9]{1,2})?|0\.(?:[1-9][0-9]?|0[1-9]))", row[key]) is None):
                    return False
            pair_ids.append((row["down"], row["up"]))
        if not _unique(pair_ids):
            return False
        base_ids = []
        for row in value["bases"]:
            if not _exact_keys(row, {"fragment", "sequence", "nuc1_ovrmsd", "nuc2_ovrmsd"}) \
                or type(row["fragment"]) is not int or row["fragment"] < 1 \
                or not isinstance(row["sequence"], str) or re.fullmatch(r"[ACGU]{2}", row["sequence"]) is None:
                return False
            for key in ("nuc1_ovrmsd", "nuc2_ovrmsd"):
                if row[key] is not None and (not isinstance(row[key], str)
                    or re.fullmatch(r"(?:[1-9][0-9]*(?:\.[0-9]{1,2})?|0\.(?:[1-9][0-9]?|0[1-9]))", row[key]) is None):
                    return False
            base_ids.append(row["fragment"])
        return _unique(base_ids)
    if name == "best_fit_clip":
        if type(value["clip1"]) is not int or value["clip1"] < 0 or type(value["clip2"]) is not int or value["clip2"] < 0:
            return False
        return (value["basis"] is None and value["clip1"] == 0 and value["clip2"] == 0) or _valid_rna_basis(value["basis"])
    if name == "interaction_matrix":
        if type(value["active"]) is not bool or type(value["basis"]) is not list \
            or not all(isinstance(key, str) and key for key in value["basis"]) or not _unique(value["basis"]):
            return False
        if type(value["rows"]) is not list or not all(isinstance(key, str) and re.fullmatch(r"frag-[1-9][0-9]*-[ACGU]{2}", key) for key in value["rows"]) \
            or not _unique(value["rows"]):
            return False
        if type(value["ranges"]) is not list or not all(isinstance(key, str) and re.fullmatch(r"range-[1-9][0-9]*", key) for key in value["ranges"]) \
            or type(value["grid"]) is not list:
            return False
        return all(type(row) is list and len(row) == len(value["ranges"]) and all(type(cell) is bool for cell in row) for row in value["grid"])
    if name == "gaussian_classes":
        if type(value["active"]) is not bool or type(value["basis"]) is not list \
            or not all(isinstance(key, str) and key for key in value["basis"]) or not _unique(value["basis"]) \
            or type(value["rows"]) is not list:
            return False
        keys = []
        for row in value["rows"]:
            if type(row) is not dict or not isinstance(row.get("key"), str) or not row["key"] \
                or row.get("rrm") not in (50, 60, 70, 80, 90, 95) or type(row.get("rrm")) is not int:
                return False
            expected = {"key", "rrm"} if row["rrm"] < 90 else {"key", "rrm", "tilt30", "disp45"}
            if set(row) != expected or (row["rrm"] >= 90 and (type(row["tilt30"]) is not bool or type(row["disp45"]) is not bool)):
                return False
            keys.append(row["key"])
        return _unique(keys) and sorted(keys) == sorted(value["basis"])
    if name == "energy_threshold_edits":
        if not _valid_rna_basis(value["basis"]) or type(value["ranges"]) is not list:
            return False
        rows = []
        for entry in value["ranges"]:
            if not _exact_keys(entry, {"row", "values"}) or not _valid_range_row(entry["row"]) or type(entry["values"]) is not dict:
                return False
            rows.append(json.dumps(entry["row"], sort_keys=True))
            for fragment, threshold in entry["values"].items():
                if re.fullmatch(r"frag-[1-9][0-9]*-[ACGU]{2}", fragment) is None:
                    return False
                if threshold is not None:
                    if not isinstance(threshold, str) or re.fullmatch(
                        r"[+-]?(?:0|[1-9][0-9]*)(?:\.[0-9]{1,2})?", threshold
                    ) is None:
                        return False
                    try:
                        if not math.isfinite(float(threshold)):
                            return False
                    except ValueError:
                        return False
        return _unique(rows)
    return False


def _yaml_ranges_value(value: object) -> object:
    """Coerce BaseLoader scalars by field, preserving chain IDs like ``NO``."""
    if value in ("null", "Null", "NULL", "~"):
        return None
    if type(value) is not dict:
        return value
    result = copy.deepcopy(value)
    if "schema" in result and isinstance(result["schema"], str) and result["schema"].isdigit():
        result["schema"] = int(result["schema"])
    if isinstance(result.get("basis"), dict):
        biounit = result["basis"].get("biounit")
        if isinstance(biounit, str) and re.fullmatch(r"[0-9]+", biounit):
            result["basis"]["biounit"] = int(biounit)
    active = result.get("active")
    if isinstance(active, str) and active.lower() in ("true", "false"):
        result["active"] = active.lower() == "true"
    if isinstance(result.get("rows"), list):
        for row in result["rows"]:
            if not isinstance(row, dict):
                continue
            for key in ("start", "end"):
                if isinstance(row.get(key), str) and re.fullmatch(r"-?[0-9]+", row[key]):
                    row[key] = int(row[key])
    return result


def _default_marker(name: str) -> object:
    if name == "best_fit_clip":
        return {"schema": 1, "basis": None, "clip1": 0, "clip2": 0}
    return None


def load_marker_payload(path: Path, marker_name: str, fallback: object = _NO_FALLBACK) -> tuple[object, list[str]]:
    """Load a schema-1 mounted decision without migrating legacy content.

    Diagnostics are returned before startup normalization touches a malformed
    file.  The caller can log/present them, then replace that file with the
    undecided default before mounting it file-authoritatively.
    """
    if marker_name not in _MARKER_KEYS:
        raise ValueError(f"Unknown mounted decision {marker_name!r}.")
    if fallback is _NO_FALLBACK:
        fallback = _default_marker(marker_name)
    if not path.exists():
        return fallback, []
    raw = path.read_text(encoding="utf-8")
    try:
        if marker_name == "ranges":
            value = _yaml_ranges_value(yaml.load(raw, Loader=yaml.BaseLoader))
        else:
            value = json.loads(raw)
    except (UnicodeError, json.JSONDecodeError, yaml.YAMLError) as exc:
        return fallback, [f"{marker_name}: stored decision could not be parsed and was reset: {exc}"]
    if not _valid_decision(marker_name, value):
        schema = value.get("schema") if isinstance(value, dict) else None
        if schema is not None and schema != 1:
            reason = f"unsupported schema {schema!r}"
        elif value is not None and isinstance(value, dict) and "basis" not in value:
            reason = "stored decision has no basis, so its upstream project cannot be established"
        else:
            reason = "stored decision does not match the schema-1 contract"
        return fallback, [f"{marker_name}: {reason}; reset to undecided"]
    return value, []


def _write_marker_default(path: Path, marker_name: str, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if marker_name == "ranges":
        raw = yaml.safe_dump(value, sort_keys=False)
    else:
        raw = json.dumps(value, indent=2, sort_keys=True) + "\n"
    path.write_text(raw, encoding="utf-8")


def _mounted_decision(context: Context, name: str, path: Path, celltype: str = "plain") -> Cell:
    fallback = _default_marker(name)
    payload, diagnostics = load_marker_payload(path, name, fallback)
    for diagnostic in diagnostics:
        print(diagnostic, file=sys.stderr, flush=True)
        _write_marker_default(path, name, fallback)
    cell = Cell(celltype=celltype)
    if celltype == "yaml":
        cell.set(yaml.safe_dump(payload, sort_keys=False))
    else:
        cell.set(payload)
    setattr(context, name, cell)
    cell.mount(path, mode="rw", authority="file")
    return cell


def configure_local_execution(spawn_fn=None) -> int:
    """Select the approved eight-worker local spawn runtime exactly once.

    No configuration loader is called here: this backend has no persistent
    Seamless result cache, hashserver, database, or remote execution service.
    """
    global _LOCAL_EXECUTION_CONFIGURED
    if _LOCAL_EXECUTION_CONFIGURED:
        return LOCAL_WORKER_COUNT
    if spawn_fn is None:
        from seamless_transformer.worker import has_spawned, spawn
        if not has_spawned():
            spawn_fn = spawn
        else:
            _LOCAL_EXECUTION_CONFIGURED = True
            return LOCAL_WORKER_COUNT
    spawn_fn(LOCAL_WORKER_COUNT)
    _LOCAL_EXECUTION_CONFIGURED = True
    return LOCAL_WORKER_COUNT


def _initial_pdb_code(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8").strip().upper()


def _initial_biounit_index(path: Path) -> int:
    if not path.exists():
        return 1
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return 1
    return value if type(value) is int and value >= 1 else 1


def _initial_range_specs(marker: object) -> dict | None:
    """Small, deterministic T0 projection used before the live transformer."""
    return _parse_ranges(marker, _initial_pdb_code(PDBCODE_PATH), _initial_biounit_index(BIOUNIT_PATH))


def _parse_ranges(ranges: object, pdbcode: str, biounit: int) -> dict | None:
    """T0: normalize the active range decision for backend consumers."""
    import re
    import yaml

    empty_ranges = {"schema": 1, "pdbcode": pdbcode, "biounit": biounit, "rows": []}
    if ranges is None:
        return empty_ranges
    if isinstance(ranges, str):
        try:
            ranges = yaml.load(ranges, Loader=yaml.BaseLoader)
        except yaml.YAMLError:
            return None
        # BaseLoader preserves scalar spellings as strings, so an undecided
        # YAML cell mounted as ``null\n...`` arrives here as the literal
        # string ``null`` rather than Python None. Treat YAML null documents
        # as the valid empty S5 envelope before validating the marker mapping.
        if ranges is None or (isinstance(ranges, str) and ranges.strip() in {"null", "Null", "NULL", "~"}):
            return empty_ranges
        if type(ranges) is dict:
            schema = ranges.get("schema")
            if isinstance(schema, str) and re.fullmatch(r"[0-9]+", schema):
                ranges["schema"] = int(schema)
            basis = ranges.get("basis")
            if type(basis) is dict:
                raw_biounit = basis.get("biounit")
                if isinstance(raw_biounit, str) and re.fullmatch(r"[0-9]+", raw_biounit):
                    basis["biounit"] = int(raw_biounit)
            active = ranges.get("active")
            if isinstance(active, str) and active.lower() in ("true", "false"):
                ranges["active"] = active.lower() == "true"
            rows = ranges.get("rows")
            if isinstance(rows, list):
                for row in rows:
                    if type(row) is not dict:
                        continue
                    for key in ("start", "end"):
                        value = row.get(key)
                        if isinstance(value, str) and re.fullmatch(r"-?[0-9]+", value):
                            row[key] = int(value)
    if type(ranges) is not dict or set(ranges) != {"schema", "basis", "active", "rows"} \
        or type(ranges.get("schema")) is not int or ranges.get("schema") != 1 \
        or type(ranges.get("active")) is not bool:
        return None
    basis = ranges.get("basis")
    rows = ranges.get("rows")
    # Keep this transformer self-contained. Seamless may serialize the
    # function into a worker process, where module globals such as
    # ``_valid_range_row`` are not attached automatically.
    def valid_row(row: object) -> bool:
        return type(row) is dict and set(row) == {"chain", "start", "end"} \
            and isinstance(row["chain"], str) and bool(row["chain"]) \
            and type(row["start"]) is int and type(row["end"]) is int \
            and row["end"] >= row["start"]

    if type(basis) is not dict or set(basis) != {"pdbcode", "biounit"} \
        or not isinstance(basis.get("pdbcode"), str) or basis.get("pdbcode") != pdbcode \
        or type(basis.get("biounit")) is not int or basis.get("biounit") != biounit \
        or not isinstance(rows, list) or not all(valid_row(row) for row in rows):
        return None
    if ranges["active"] is False:
        return empty_ranges
    if not rows:
        return empty_ranges
    addresses = [(row["chain"], row["start"], row["end"]) for row in rows]
    if len(addresses) != len(set(addresses)):
        return None
    return {"schema": 1, "pdbcode": pdbcode, "biounit": biounit, "rows": [dict(row) for row in rows]}


def _structure_id(cif: str) -> str:
    """T26: read the entry id from the persisted CIF used by T1."""
    if not isinstance(cif, str) or not cif.strip():
        raise ValueError("The stored CIF is empty; download the structure again.")
    from io import StringIO
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    values = MMCIF2Dict(StringIO(cif)).get("_entry.id")
    if isinstance(values, list) and len(values) == 1:
        value = values[0]
    elif isinstance(values, str):
        value = values
    else:
        raise ValueError("The stored CIF has no single _entry.id value.")
    if not isinstance(value, str) or not value.strip() or value in {".", "?"}:
        raise ValueError("The stored CIF has no usable _entry.id value.")
    return value.strip().upper()


def _lower_code(pdbcode: str) -> str:
    return pdbcode.strip().lower() if isinstance(pdbcode, str) and pdbcode.strip() else ""


def _protein_pdbs_from_aareduce(pdb_files: dict) -> dict:
    if not isinstance(pdb_files, dict):
        return {}
    return {
        name: value for name, value in pdb_files.items()
        if isinstance(name, str) and name.startswith("range-") and name.endswith("-aa.pdb")
        and name != "range-1-aa.pdb" and isinstance(value, bytes)
    }


def _first_rna_range(range_specs: dict) -> dict | None:
    if not isinstance(range_specs, dict) or range_specs.get("schema") != 1:
        return None
    rows = range_specs.get("rows")
    return dict(rows[0]) if isinstance(rows, list) and rows and isinstance(rows[0], dict) else None


def _fit_input_gate(rna_pdb: bytes | None, selection: object) -> bytes | None:
    """Leave T5 unwired at the content level until S8 is active."""
    rna_pdb = getattr(rna_pdb, "content", rna_pdb)
    if not isinstance(rna_pdb, bytes) or not isinstance(selection, dict) \
        or selection.get("active") is not True:
        return None
    return rna_pdb


def _stackings_envelope(
    range_pdbs: list[str],
    ring_ref_a: bytes,
    ring_ref_c: bytes,
    ring_ref_phe: bytes,
    ring_ref_his: bytes,
    ring_ref_arg: bytes,
    ring_ref_trp: bytes,
    range_specs: dict,
    structure_id: str,
) -> dict:
    import preparic_detect_stacking

    if not isinstance(range_specs, dict) or range_specs.get("schema") != 1:
        return None
    # Seamless passes bytes pins as Buffer instances to transformers.
    # Unwrap their exact payloads before calling the detector.
    ring_references = {
        "refe-A.pdb": getattr(ring_ref_a, "content", ring_ref_a),
        "refe-C.pdb": getattr(ring_ref_c, "content", ring_ref_c),
        "refe-PHE.pdb": getattr(ring_ref_phe, "content", ring_ref_phe),
        "refe-HIS.pdb": getattr(ring_ref_his, "content", ring_ref_his),
        "refe-ARG.pdb": getattr(ring_ref_arg, "content", ring_ref_arg),
        "refe-TRP.pdb": getattr(ring_ref_trp, "content", ring_ref_trp),
    }
    rows = preparic_detect_stacking._detect_stacking(range_pdbs, ring_references=ring_references)
    return {
        "schema": 1,
        "echo": {
            "structure_id": structure_id,
            "pdbcode": range_specs.get("pdbcode"),
            "biounit": range_specs.get("biounit"),
            "rows": range_specs.get("rows", []),
        },
        "rows": rows,
    }


def _share_aareduce_pdbs(pdb_files: dict, range_specs: dict, exclude_pdb_code: str) -> dict:
    if not isinstance(range_specs, dict) or not range_specs.get("rows") \
        or not isinstance(exclude_pdb_code, str) or not exclude_pdb_code.strip():
        return None
    import preparic_stage4
    return preparic_stage4._share_range_pdbs_aa(pdb_files, range_specs, exclude_pdb_code)


def _stage4_fit_echo(rna_pdb: bytes, pdbcode: str, rna_range: dict) -> dict:
    import preparic_stage4
    rna_pdb = getattr(rna_pdb, "content", rna_pdb)
    return preparic_stage4._fit_echo(rna_pdb, pdbcode, rna_range)


def _stage4_list_fragments(
    refe_fragments: dict, rna_pdb: bytes, pdbcode: str,
    rna_range: dict,
) -> dict:
    import preparic_stage4
    rna_pdb = getattr(rna_pdb, "content", rna_pdb)
    return preparic_stage4._list_fragments(refe_fragments, rna_pdb, pdbcode, rna_range)


def _stage4_reference_fragment_pdbs(refe_fragments: dict, fragments: dict) -> dict:
    import preparic_stage4
    return preparic_stage4._reference_fragment_pdbs(refe_fragments, fragments)


def _stage4_parse_best_fits(best_fit: dict, fit_echo: dict) -> dict:
    import preparic_stage4
    return preparic_stage4._parse_best_fits(best_fit, fit_echo)


def _stage4_validate_threshold_edits(edits: object, fragments: object) -> None:
    import preparic_stage4
    return preparic_stage4._validate_threshold_edits(edits, fragments)


def _stage4_clip_fragments(
    best_fits: dict, refe_fragments: dict, thresholds: object,
    clip: dict, fragments: dict, selection: object, stackings: object,
) -> dict:
    import preparic_stage4
    if selection is None or (type(selection) is dict and selection.get("active") is False):
        return {}
    candidates = preparic_stage4._selected_candidate_numbers(fragments, selection, stackings)
    return preparic_stage4._clip_fragments(
        best_fits, refe_fragments, thresholds, clip, fragments, candidates
    )


def _stage5_energetics_ligands(refe_fragments: dict, best_fit: dict, fragments: dict) -> dict:
    import preparic_stage5
    return preparic_stage5.energetics_ligands(refe_fragments, best_fit, fragments)


def _stage5_merge_energetics(
    range_results: dict, range_specs: dict, best_fits: dict, fragments: dict,
) -> dict:
    import json
    import preparic_stage5
    decoded = {}
    for name, payload in range_results.items():
        if not isinstance(name, str) or not name.endswith(".json") or not isinstance(payload, bytes):
            continue
        try:
            decoded[name[:-5]] = json.loads(payload.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise ValueError(f"Energetics result {name} is not valid UTF-8 JSON: {exc}") from exc
    return preparic_stage5.merge_energetics(decoded, range_specs, best_fits, fragments)


def _reconcile_stage5_range_scorers(context: Context, alaric_dir: Path, stop: threading.Event) -> None:
    """Keep one T13 graph node and one result path for each current S5 row."""
    installed: dict[int, str] = {}
    while not stop.wait(0.25):
        try:
            specs = context.range_specs.value
        except Exception:
            continue
        rows = specs.get("rows") if type(specs) is dict and specs.get("schema") == 1 else None
        wanted = set(range(2, len(rows) + 1)) if isinstance(rows, list) and len(rows) >= 2 else set()

        for index in sorted(wanted - installed.keys()):
            name = f"energetics_range_{index}"
            try:
                setattr(context, name, energetics_range_transformer(alaric_dir))
                scorer = getattr(context, name)
                scorer.pins.ligands = context.energetics_ligands.result
                scorer.pins.receptor = context.protein_pdbs.result[f"range-{index}-aa.pdb"]
                scorer.celltypes.ligands = "folder"
                scorer.celltypes.receptor = "bytes"
                scorer.celltypes.result = "folder"
                context.energetics_range_results[f"range-{index}.json"] = scorer.result["energetics.json"]
                installed[index] = name
            except Exception:
                # Keep the desired index uninstalled so a transient Context
                # ingress failure is retried on the next reconciliation pass.
                continue

        for index in sorted(installed.keys() - wanted):
            name = installed[index]
            try:
                del context.energetics_range_results[f"range-{index}.json"]
                delattr(context, name)
                del installed[index]
            except Exception:
                continue


def _stage5_energy_thresholds(energetics: dict, best_fits: dict) -> dict:
    import energy_threshold
    import preparic_stage5
    return preparic_stage5.energy_thresholds(energetics, best_fits, energy_threshold)


def _stage5_validate_energy_threshold_edits(edits: object, fragments: object, range_specs: object) -> None:
    import preparic_stage5
    return preparic_stage5.validate_energy_threshold_edits(edits, fragments, range_specs)


def _stage5_interaction_handoff(
    matrix: object,
    fragments: object,
    selection: object,
    stackings: object,
    clip: object,
    range_specs: object,
) -> dict | None:
    import preparic_stage5
    return preparic_stage5.interaction_handoff(
        matrix, fragments, selection, stackings, clip, range_specs
    )


def _stage6_detect_stacking_bestfit(
    best_fits: object, proteins: object, fragments: object, range_specs: object,
    ring_ref_a: object, ring_ref_c: object, ring_ref_phe: object,
    ring_ref_his: object, ring_ref_arg: object, ring_ref_trp: object,
) -> dict | None:
    import preparic_stage6
    return preparic_stage6.detect_stacking_bestfit(
        best_fits, proteins, fragments, range_specs,
        ring_ref_a, ring_ref_c, ring_ref_phe, ring_ref_his, ring_ref_arg, ring_ref_trp,
    )


def _stage6_gaussian_check_jobs(
    best_fits: object, fragments: object, stackings: object,
    selection: object, clip: object,
) -> dict | None:
    import preparic_stage6
    return preparic_stage6.gaussian_check_jobs(best_fits, fragments, stackings, selection, clip)


def _stage6_validate_gaussian_classes(
    marker: object, model_json: object, selection: object,
) -> None:
    import preparic_stage6
    import alaric_gaussian_filter
    return preparic_stage6.validate_gaussian_classes(
        marker, model_json, selection, alaric_gaussian_filter,
    )


def _stage6_semantic_gaussian_identity(model_json: object) -> str:
    import preparic_stage6
    return preparic_stage6.semantic_gaussian_identity(model_json)


def _stage6_jobs_json(jobs_result: object) -> str | None:
    import json

    if jobs_result is None:
        return None
    if type(jobs_result) is not dict or not isinstance(jobs_result.get("jobs"), list):
        raise ValueError("T20 Gaussian job list is malformed.")
    return json.dumps(jobs_result["jobs"], sort_keys=True, separators=(",", ":"))


def _stage6_read_gaussian_check(result_folder: object) -> dict:
    import preparic_stage6
    return preparic_stage6.read_gaussian_check_result(result_folder)


def _stage6_merge_gaussian_validation(
    check_result: object, jobs_result: object, class_marker: object,
    selection: object, clip: object, best_fits: object, model_json: object,
    class_validation: object,
) -> dict | None:
    import preparic_stage6
    import alaric_gaussian_filter
    if class_validation is not None:
        raise ValueError("Gaussian class validation did not return its success sentinel.")
    return preparic_stage6.merge_gaussian_validation(
        check_result, jobs_result, class_marker, selection, clip, best_fits,
        model_json, alaric_gaussian_filter,
    )


def _stage6_validate_refe_anchors(
    threshold_edits: object, best_fit_thresholds: object, fragments: object,
    stackings: object, selection: object, clip: object, best_fits: object,
) -> dict | None:
    import preparic_stage6
    return preparic_stage6.validate_refe_anchors(
        threshold_edits, best_fit_thresholds, fragments, stackings, selection, clip, best_fits,
    )


def _stage6_gaussian_scores(
    bestfit_stackings: object, reference_stackings: object, selection: object,
    clip: object, fragments: object, model_json: object,
) -> dict | None:
    import preparic_stage6
    import alaric_gaussian_filter
    return preparic_stage6.gaussian_scores(
        bestfit_stackings, reference_stackings, selection, clip, fragments,
        model_json, alaric_gaussian_filter,
    )


def _add_stacking_pipeline(
    context: Context, ring_references: dict[str, Cell], mount_directory: Path | None = None
) -> None:
    """Wire range PDBs and live mounted ring references into T4."""
    if mount_directory is None:
        mount_directory = MOUNT_DIRECTORY
    context.detect_stacking = _stackings_envelope
    context.detect_stacking.modules.preparic_detect_stacking = {
        "type": "interpreted",
        "language": "python",
        "code": (PROJECT_ROOT / "detect_stacking.py").read_text(encoding="utf-8"),
    }
    context.detect_stacking.celltypes.range_pdbs = "plain"
    context.detect_stacking.celltypes.ring_ref_a = "bytes"
    context.detect_stacking.celltypes.ring_ref_c = "bytes"
    context.detect_stacking.celltypes.ring_ref_phe = "bytes"
    context.detect_stacking.celltypes.ring_ref_his = "bytes"
    context.detect_stacking.celltypes.ring_ref_arg = "bytes"
    context.detect_stacking.celltypes.ring_ref_trp = "bytes"
    context.detect_stacking.celltypes.range_specs = "plain"
    context.detect_stacking.celltypes.structure_id = "plain"
    context.detect_stacking.celltypes.result = "plain"
    context.detect_stacking.pins.range_pdbs = context.range_pdbs
    context.detect_stacking.pins.ring_ref_a = ring_references["A"]
    context.detect_stacking.pins.ring_ref_c = ring_references["C"]
    context.detect_stacking.pins.ring_ref_phe = ring_references["PHE"]
    context.detect_stacking.pins.ring_ref_his = ring_references["HIS"]
    context.detect_stacking.pins.ring_ref_arg = ring_references["ARG"]
    context.detect_stacking.pins.ring_ref_trp = ring_references["TRP"]
    context.detect_stacking.pins.range_specs = context.range_specs
    context.detect_stacking.pins.structure_id = context.structure_id_result

    context.stackings = Cell(celltype="plain")
    context.stackings = context.detect_stacking.result
    context.stackings.mount(
        mount_directory / "stackings.json", mode="w", authority="cell"
    )


def build_context() -> Context:
    """Build the durable project Context backing the Vite development proxy."""
    global ctx
    configure_local_execution()
    alaric_dir = resolve_alaric_dir()
    alaric_root = alaric_dir.parent
    if shareserver.url is None:
        configured_host = shareserver.host
        host = os.environ.get("SEAMLESS_SHARE_HOST")
        if host is None:
            host = "127.0.0.1" if configured_host == "0.0.0.0" else configured_host
        port = int(os.environ.get("SEAMLESS_SHARE_PORT", shareserver.port))
        shareserver.configure(host=host, port=port)

    MOUNT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    context = Context()
    try:
        context.shares.namespace = "preparic"

        context.pdbcode = Cell(celltype="text")
        context.pdbcode.set(_initial_pdb_code(PDBCODE_PATH))
        context.pdbcode.mount(PDBCODE_PATH, mode="rw", authority="file")

        context.biounit = Cell(celltype="plain")
        context.biounit.set(_initial_biounit_index(BIOUNIT_PATH))
        context.biounit.mount(BIOUNIT_PATH, mode="rw", authority="file")

        context.cif = Cell(celltype="text")
        context.cif.set("")
        context.cif.mount(CIF_PATH, mode="rw", authority="file")

        _mounted_decision(context, "ranges", RANGES_PATH, celltype="yaml")

        context.parse_ranges = _parse_ranges
        context.parse_ranges.pins.ranges = context.ranges
        context.parse_ranges.pins.pdbcode = context.pdbcode
        context.parse_ranges.pins.biounit = context.biounit
        context.parse_ranges.celltypes.result = "plain"
        context.range_specs = context.parse_ranges.result

        context.structure_id = _structure_id
        context.structure_id.pins.cif = context.cif
        context.structure_id.celltypes.result = "plain"
        context.structure_id_result = context.structure_id.result

        for name in ("A", "C", "PHE", "HIS", "ARG", "TRP"):
            attribute = f"ring_refe_{name.lower()}"
            setattr(context, attribute, Cell(celltype="bytes"))
            reference = getattr(context, attribute)
            reference.mount(alaric_root / "ring-refe" / f"refe-{name}.pdb", mode="r")

        _mounted_decision(context, "stacking_selection", STACKING_SELECTION_PATH)
        _mounted_decision(context, "threshold_edits", THRESHOLD_EDITS_PATH)
        _mounted_decision(context, "best_fit_clip", BEST_FIT_CLIP_PATH)
        _mounted_decision(context, "interaction_matrix", INTERACTION_MATRIX_PATH)
        _mounted_decision(context, "gaussian_classes", GAUSSIAN_CLASSES_PATH)
        _mounted_decision(context, "energy_threshold_edits", ENERGY_THRESHOLD_EDITS_PATH)

        context.ranges_to_pdb = _ranges_to_pdb
        context.ranges_to_pdb.pins.cif = context.cif
        context.ranges_to_pdb.pins.ranges = context.range_specs
        context.range_pdbs = Cell(celltype="plain")
        context.range_pdbs = context.ranges_to_pdb.result

        context.pdbs_to_folder = _pdbs_to_folder
        context.pdbs_to_folder.celltypes.pdb_texts = "plain"
        context.pdbs_to_folder.celltypes.result = "folder"
        context.pdbs_to_folder.pins.pdb_texts = context.range_pdbs
        context.range_pdb_files = context.pdbs_to_folder.result
        context.range_pdb_files.mount(PDBS_PATH, mode="w")

        context.aareduce_files = Cell(celltype="folder")
        context.aareduce_files.mount(PROJECT_ROOT / "aareduce", mode="r")

        aareduce_ranges = Transformer("bash")
        aareduce_ranges.code = (
            PROJECT_ROOT / "range-pdbs-aareduce.sh"
        ).read_text()
        aareduce_ranges.environment.set_conda_env("preparic")
        context.aareduce_ranges = aareduce_ranges
        context.aareduce_ranges.pins.pdbs = context.range_pdb_files
        context.aareduce_ranges.pins.aareduce = context.aareduce_files
        context.aareduce_ranges.celltypes.pdbs = "folder"
        context.aareduce_ranges.celltypes.aareduce = "folder"
        context.aareduce_ranges.celltypes.result = "folder"
        context.range_pdb_files_aa = context.aareduce_ranges.result
        context.range_pdb_files_aa.mount(PDBS_AA_PATH, mode="w")

        _add_stacking_pipeline(context, {
            "A": context.ring_refe_a,
            "C": context.ring_refe_c,
            "PHE": context.ring_refe_phe,
            "HIS": context.ring_refe_his,
            "ARG": context.ring_refe_arg,
            "TRP": context.ring_refe_trp,
        })

        # B4/B5/B6/B7 are plain backend-only inputs shared by the Alaric steps.
        context.rna_pdb = context.range_pdb_files_aa["range-1-aa.pdb"]

        context.exclude_pdb_code = _lower_code
        context.exclude_pdb_code.pins.pdbcode = context.pdbcode
        context.exclude_pdb_code.celltypes.result = "text"

        context.protein_pdbs = _protein_pdbs_from_aareduce
        context.protein_pdbs.pins.pdb_files = context.range_pdb_files_aa
        context.protein_pdbs.celltypes.pdb_files = "folder"
        context.protein_pdbs.celltypes.result = "folder"

        identities = fragment_library_identities(alaric_dir)
        context.fraglib_identities = Cell(celltype="text")
        context.fraglib_identities.set(json.dumps(identities, sort_keys=True, separators=(",", ":")))

        context.first_rna_range = _first_rna_range
        context.first_rna_range.pins.range_specs = context.range_specs
        context.first_rna_range.celltypes.result = "plain"

        # T6/T24 are independent of the stage-3 selection and of all fit work.
        context.refe_fragments = refe_fragments_transformer(alaric_dir)
        context.refe_fragments.pins.rna_pdb = context.rna_pdb

        stage4_module = {
            "type": "interpreted",
            "language": "python",
            "code": (PROJECT_ROOT / "stage4.py").read_text(encoding="utf-8"),
        }
        context.list_fragments = _stage4_list_fragments
        context.list_fragments.modules.preparic_stage4 = stage4_module
        context.list_fragments.pins.refe_fragments = context.refe_fragments.result
        context.list_fragments.pins.rna_pdb = context.rna_pdb
        context.list_fragments.pins.pdbcode = context.pdbcode
        context.list_fragments.pins.rna_range = context.first_rna_range.result
        context.list_fragments.celltypes.refe_fragments = "folder"
        context.list_fragments.celltypes.rna_pdb = "bytes"
        context.list_fragments.celltypes.pdbcode = "text"
        context.list_fragments.celltypes.rna_range = "plain"
        context.list_fragments.celltypes.result = "plain"
        context.fragments = context.list_fragments.result

        context.reference_fragment_pdbs_transformer = _stage4_reference_fragment_pdbs
        context.reference_fragment_pdbs_transformer.modules.preparic_stage4 = stage4_module
        context.reference_fragment_pdbs_transformer.pins.refe_fragments = context.refe_fragments.result
        context.reference_fragment_pdbs_transformer.pins.fragments = context.fragments
        context.reference_fragment_pdbs_transformer.celltypes.refe_fragments = "folder"
        context.reference_fragment_pdbs_transformer.celltypes.fragments = "plain"
        context.reference_fragment_pdbs_transformer.celltypes.result = "plain"
        context.reference_fragment_pdbs = context.reference_fragment_pdbs_transformer.result

        context.share_range_pdbs_aa = _share_aareduce_pdbs
        context.share_range_pdbs_aa.modules.preparic_stage4 = stage4_module
        context.share_range_pdbs_aa.pins.pdb_files = context.range_pdb_files_aa
        context.share_range_pdbs_aa.pins.range_specs = context.range_specs
        context.share_range_pdbs_aa.pins.exclude_pdb_code = context.exclude_pdb_code
        context.share_range_pdbs_aa.celltypes.pdb_files = "folder"
        context.share_range_pdbs_aa.celltypes.range_specs = "plain"
        context.share_range_pdbs_aa.celltypes.exclude_pdb_code = "text"
        context.share_range_pdbs_aa.celltypes.result = "plain"
        context.range_pdbs_aa = context.share_range_pdbs_aa.result

        # T5 is wired to a None-or-bytes gate, so it does no library work while
        # S8 is null/inactive. Its Alaric factory keeps only project inputs as
        # unbound pins, and its code/template inputs remain checksum-pinned.
        context.fit_rna_pdb = _fit_input_gate
        context.fit_rna_pdb.pins.rna_pdb = context.rna_pdb
        context.fit_rna_pdb.pins.selection = context.stacking_selection
        context.fit_rna_pdb.celltypes.rna_pdb = "bytes"
        context.fit_rna_pdb.celltypes.selection = "plain"
        context.fit_rna_pdb.celltypes.result = "bytes"
        fraglib_config = Path.home() / ".alaric" / "fraglib.yaml"
        fit_factory = best_fit_factory(alaric_dir)
        context.best_fit = fit_factory.seamlessify_refe_best_fit(
            fraglib_config=fraglib_config,
            conda_environment="alaric",
        )
        context.best_fit.pins.rna_pdb = context.fit_rna_pdb.result
        context.best_fit.pins.exclude_pdb_code = context.exclude_pdb_code
        context.best_fit.pins.fraglib_identities = context.fraglib_identities

        # T7 uses the exact T5 inputs and still runs in the Alaric environment.
        context.best_fit_tables = best_fit_tables_transformer(alaric_dir, fraglib_config)
        context.best_fit_tables.pins.best_fit = context.best_fit.result
        context.best_fit_tables.pins.rna_pdb = context.rna_pdb
        context.best_fit_tables.pins.exclude_pdb_code = context.exclude_pdb_code
        context.best_fit_tables.pins.fraglib_identities = context.fraglib_identities

        context.fit_echo = _stage4_fit_echo
        context.fit_echo.modules.preparic_stage4 = stage4_module
        context.fit_echo.pins.rna_pdb = context.rna_pdb
        context.fit_echo.pins.pdbcode = context.pdbcode
        context.fit_echo.pins.rna_range = context.first_rna_range.result
        context.fit_echo.celltypes.rna_pdb = "bytes"
        context.fit_echo.celltypes.pdbcode = "text"
        context.fit_echo.celltypes.rna_range = "plain"
        context.fit_echo.celltypes.result = "plain"

        context.best_fit_thresholds_transformer = thresholds_from_tables
        context.best_fit_thresholds_transformer.pins.tables = context.best_fit_tables.result
        context.best_fit_thresholds_transformer.pins.fit_echo = context.fit_echo.result
        context.best_fit_thresholds_transformer.modules.refe_thresholds, context.best_fit_thresholds_transformer.modules.threshold_rounding = threshold_module_code(alaric_dir)
        context.best_fit_thresholds_transformer.celltypes.tables = "folder"
        context.best_fit_thresholds_transformer.celltypes.fit_echo = "plain"
        context.best_fit_thresholds_transformer.celltypes.result = "plain"
        context.best_fit_thresholds = context.best_fit_thresholds_transformer.result

        context.parse_best_fits = _stage4_parse_best_fits
        context.parse_best_fits.modules.preparic_stage4 = stage4_module
        context.parse_best_fits.pins.best_fit = context.best_fit.result
        context.parse_best_fits.pins.fit_echo = context.fit_echo.result
        context.parse_best_fits.celltypes.best_fit = "folder"
        context.parse_best_fits.celltypes.fit_echo = "plain"
        context.parse_best_fits.celltypes.result = "plain"
        context.best_fits = context.parse_best_fits.result

        context.validate_threshold_edits = _stage4_validate_threshold_edits
        context.validate_threshold_edits.modules.preparic_stage4 = stage4_module
        context.validate_threshold_edits.pins.edits = context.threshold_edits
        context.validate_threshold_edits.pins.fragments = context.fragments
        context.validate_threshold_edits.celltypes.edits = "plain"
        context.validate_threshold_edits.celltypes.fragments = "plain"

        context.clip_fragments = _stage4_clip_fragments
        context.clip_fragments.modules.preparic_stage4 = stage4_module
        context.clip_fragments.pins.best_fits = context.best_fits
        context.clip_fragments.pins.refe_fragments = context.refe_fragments.result
        context.clip_fragments.pins.thresholds = context.threshold_edits
        context.clip_fragments.pins.clip = context.best_fit_clip
        context.clip_fragments.pins.fragments = context.fragments
        context.clip_fragments.pins.selection = context.stacking_selection
        context.clip_fragments.pins.stackings = context.stackings
        context.clip_fragments.celltypes.best_fits = "plain"
        context.clip_fragments.celltypes.refe_fragments = "folder"
        context.clip_fragments.celltypes.thresholds = "plain"
        context.clip_fragments.celltypes.clip = "plain"
        context.clip_fragments.celltypes.fragments = "plain"
        context.clip_fragments.celltypes.selection = "plain"
        context.clip_fragments.celltypes.stackings = "plain"
        context.clip_fragments.celltypes.result = "folder"
        context.clipped_fragments = context.clip_fragments.result
        context.clipped_fragments.mount(CLIPPED_FRAGMENTS_PATH, mode="w")

        # W4/W5 energetics and interaction hand-off. T12 deliberately reads
        # the full T6/T5 folders and has no clip input, so clip edits do not
        # redo fitting or ligand preparation. The live S5-row reconciler adds
        # one receptor-pinned T13 scorer per protein range; T14 restores order.
        stage5_module = {
            "type": "interpreted",
            "language": "python",
            "code": (PROJECT_ROOT / "stage5.py").read_text(encoding="utf-8"),
        }
        context.energetics_ligands = _stage5_energetics_ligands
        context.energetics_ligands.modules.preparic_stage5 = stage5_module
        context.energetics_ligands.pins.refe_fragments = context.refe_fragments.result
        context.energetics_ligands.pins.best_fit = context.best_fit.result
        context.energetics_ligands.pins.fragments = context.fragments
        context.energetics_ligands.celltypes.refe_fragments = "folder"
        context.energetics_ligands.celltypes.best_fit = "folder"
        context.energetics_ligands.celltypes.fragments = "plain"
        context.energetics_ligands.celltypes.result = "folder"

        # T13 scorers are installed and removed by the S5-row reconciler below.
        # Each node is pinned to one receptor file, so a range edit invalidates
        # only that node while a ligand-folder change invalidates every scorer.
        context.energetics_range_results = Cell(celltype="folder")
        context.merge_energetics = _stage5_merge_energetics
        context.merge_energetics.modules.preparic_stage5 = stage5_module
        context.merge_energetics.pins.range_results = context.energetics_range_results
        context.merge_energetics.pins.range_specs = context.range_specs
        context.merge_energetics.pins.best_fits = context.best_fits
        context.merge_energetics.pins.fragments = context.fragments
        context.merge_energetics.celltypes.range_results = "folder"
        context.merge_energetics.celltypes.range_specs = "plain"
        context.merge_energetics.celltypes.best_fits = "plain"
        context.merge_energetics.celltypes.fragments = "plain"
        context.merge_energetics.celltypes.result = "plain"
        context.energetics = context.merge_energetics.result

        energy_threshold_code = energy_threshold_module_code(alaric_dir)
        context.energy_thresholds_transformer = _stage5_energy_thresholds
        context.energy_thresholds_transformer.modules.preparic_stage5 = stage5_module
        context.energy_thresholds_transformer.modules.energy_threshold = energy_threshold_code
        context.energy_thresholds_transformer.pins.energetics = context.energetics
        context.energy_thresholds_transformer.pins.best_fits = context.best_fits
        context.energy_thresholds_transformer.celltypes.energetics = "plain"
        context.energy_thresholds_transformer.celltypes.best_fits = "plain"
        context.energy_thresholds_transformer.celltypes.result = "plain"
        context.energy_thresholds = context.energy_thresholds_transformer.result

        context.validate_energy_threshold_edits = _stage5_validate_energy_threshold_edits
        context.validate_energy_threshold_edits.modules.preparic_stage5 = stage5_module
        context.validate_energy_threshold_edits.pins.edits = context.energy_threshold_edits
        context.validate_energy_threshold_edits.pins.fragments = context.fragments
        context.validate_energy_threshold_edits.pins.range_specs = context.range_specs
        context.validate_energy_threshold_edits.celltypes.edits = "plain"
        context.validate_energy_threshold_edits.celltypes.fragments = "plain"
        context.validate_energy_threshold_edits.celltypes.range_specs = "plain"
        context.validate_energy_threshold_edits.celltypes.result = "plain"

        context.interaction_handoff = _stage5_interaction_handoff
        context.interaction_handoff.modules.preparic_stage5 = stage5_module
        context.interaction_handoff.pins.matrix = context.interaction_matrix
        context.interaction_handoff.pins.fragments = context.fragments
        context.interaction_handoff.pins.selection = context.stacking_selection
        context.interaction_handoff.pins.stackings = context.stackings
        context.interaction_handoff.pins.clip = context.best_fit_clip
        context.interaction_handoff.pins.range_specs = context.range_specs
        for pin in ("matrix", "fragments", "selection", "stackings", "clip", "range_specs"):
            context.interaction_handoff.celltypes[pin] = "plain"
        context.interaction_handoff.celltypes.result = "plain"
        context.interactions = context.interaction_handoff.result

        # W7 stage 6: keep A7 private to the backend and publish only the
        # scores, computed classes and validation payloads consumed by the UI.
        context.gaussian_model_json = Cell(celltype="bytes")
        context.gaussian_model_json.mount(
            alaric_root / "alaric" / "gaussian_stacking_filter.json", mode="r"
        )
        gaussian_model_json = context.gaussian_model_json
        stage6_module = {
            "type": "interpreted",
            "language": "python",
            "code": (PROJECT_ROOT / "stage6.py").read_text(encoding="utf-8"),
        }
        gaussian_filter_module = {
            "type": "interpreted",
            "language": "python",
            "code": (alaric_root / "alaric" / "gaussian_filter.py").read_text(encoding="utf-8"),
        }
        detector_module = {
            "type": "interpreted",
            "language": "python",
            "code": (PROJECT_ROOT / "detect_stacking.py").read_text(encoding="utf-8"),
        }
        ring_cells = {
            "ring_ref_a": context.ring_refe_a,
            "ring_ref_c": context.ring_refe_c,
            "ring_ref_phe": context.ring_refe_phe,
            "ring_ref_his": context.ring_refe_his,
            "ring_ref_arg": context.ring_refe_arg,
            "ring_ref_trp": context.ring_refe_trp,
        }

        context.detect_stacking_bestfit = _stage6_detect_stacking_bestfit
        context.detect_stacking_bestfit.modules.preparic_stage6 = stage6_module
        context.detect_stacking_bestfit.modules.preparic_detect_stacking = detector_module
        for pin, cell in {
            "best_fits": context.best_fits,
            "proteins": context.protein_pdbs,
            "fragments": context.fragments,
            "range_specs": context.range_specs,
            **ring_cells,
        }.items():
            context.detect_stacking_bestfit.pins[pin] = cell
        context.detect_stacking_bestfit.celltypes.best_fits = "plain"
        context.detect_stacking_bestfit.celltypes.proteins = "folder"
        context.detect_stacking_bestfit.celltypes.fragments = "plain"
        context.detect_stacking_bestfit.celltypes.range_specs = "plain"
        for pin in ring_cells:
            context.detect_stacking_bestfit.celltypes[pin] = "bytes"
        context.detect_stacking_bestfit.celltypes.result = "plain"
        context.bestfit_stackings = context.detect_stacking_bestfit.result

        context.validate_gaussian_classes = _stage6_validate_gaussian_classes
        context.validate_gaussian_classes.modules.preparic_stage6 = stage6_module
        context.validate_gaussian_classes.modules.alaric_gaussian_filter = gaussian_filter_module
        context.validate_gaussian_classes.pins.marker = context.gaussian_classes
        context.validate_gaussian_classes.pins.model_json = gaussian_model_json
        context.validate_gaussian_classes.pins.selection = context.stacking_selection
        context.validate_gaussian_classes.celltypes.marker = "plain"
        context.validate_gaussian_classes.celltypes.model_json = "bytes"
        context.validate_gaussian_classes.celltypes.selection = "plain"
        context.validate_gaussian_classes.celltypes.result = "plain"

        context.gaussian_check_jobs = _stage6_gaussian_check_jobs
        context.gaussian_check_jobs.modules.preparic_stage6 = stage6_module
        context.gaussian_check_jobs.pins.best_fits = context.best_fits
        context.gaussian_check_jobs.pins.fragments = context.fragments
        context.gaussian_check_jobs.pins.stackings = context.stackings
        context.gaussian_check_jobs.pins.selection = context.stacking_selection
        context.gaussian_check_jobs.pins.clip = context.best_fit_clip
        for pin in ("best_fits", "fragments", "stackings", "selection", "clip"):
            context.gaussian_check_jobs.celltypes[pin] = "plain"
        context.gaussian_check_jobs.celltypes.result = "plain"

        context.gaussian_semantic_identity = _stage6_semantic_gaussian_identity
        context.gaussian_semantic_identity.modules.preparic_stage6 = stage6_module
        context.gaussian_semantic_identity.pins.model_json = gaussian_model_json
        context.gaussian_semantic_identity.celltypes.model_json = "bytes"
        context.gaussian_semantic_identity.celltypes.result = "text"

        context.gaussian_check_job_json = _stage6_jobs_json
        context.gaussian_check_job_json.pins.jobs_result = context.gaussian_check_jobs.result
        context.gaussian_check_job_json.celltypes.jobs_result = "plain"
        context.gaussian_check_job_json.celltypes.result = "text"

        context.gaussian_check_transformer = gaussian_check_transformer(
            alaric_dir, Path.home() / ".alaric" / "fraglib.yaml"
        )
        context.gaussian_check_transformer.pins.jobs_json = context.gaussian_check_job_json.result
        context.gaussian_check_transformer.pins.best_fit = context.best_fit.result
        context.gaussian_check_transformer.pins.proteins = context.protein_pdbs
        context.gaussian_check_transformer.pins.exclude_code = context.exclude_pdb_code
        context.gaussian_check_transformer.pins.fraglib_identities = context.fraglib_identities
        context.gaussian_check_transformer.pins.semantic_identity = context.gaussian_semantic_identity.result

        context.gaussian_check = _stage6_read_gaussian_check
        context.gaussian_check.modules.preparic_stage6 = stage6_module
        context.gaussian_check.pins.result_folder = context.gaussian_check_transformer.result
        context.gaussian_check.celltypes.result_folder = "folder"
        context.gaussian_check.celltypes.result = "plain"

        context.gaussian_scores_transformer = _stage6_gaussian_scores
        context.gaussian_scores_transformer.modules.preparic_stage6 = stage6_module
        context.gaussian_scores_transformer.modules.alaric_gaussian_filter = gaussian_filter_module
        context.gaussian_scores_transformer.pins.bestfit_stackings = context.bestfit_stackings
        context.gaussian_scores_transformer.pins.reference_stackings = context.stackings
        context.gaussian_scores_transformer.pins.selection = context.stacking_selection
        context.gaussian_scores_transformer.pins.clip = context.best_fit_clip
        context.gaussian_scores_transformer.pins.fragments = context.fragments
        context.gaussian_scores_transformer.pins.model_json = gaussian_model_json
        for pin in ("bestfit_stackings", "reference_stackings", "selection", "clip", "fragments"):
            context.gaussian_scores_transformer.celltypes[pin] = "plain"
        context.gaussian_scores_transformer.celltypes.model_json = "bytes"
        context.gaussian_scores_transformer.celltypes.result = "plain"
        context.gaussian_scores = context.gaussian_scores_transformer.result

        context.merge_gaussian_validation = _stage6_merge_gaussian_validation
        context.merge_gaussian_validation.modules.preparic_stage6 = stage6_module
        context.merge_gaussian_validation.modules.alaric_gaussian_filter = gaussian_filter_module
        context.merge_gaussian_validation.pins.check_result = context.gaussian_check
        context.merge_gaussian_validation.pins.jobs_result = context.gaussian_check_jobs.result
        context.merge_gaussian_validation.pins.class_marker = context.gaussian_classes
        context.merge_gaussian_validation.pins.selection = context.stacking_selection
        context.merge_gaussian_validation.pins.clip = context.best_fit_clip
        context.merge_gaussian_validation.pins.best_fits = context.best_fits
        context.merge_gaussian_validation.pins.model_json = gaussian_model_json
        context.merge_gaussian_validation.pins.class_validation = context.validate_gaussian_classes
        for pin in ("check_result", "jobs_result", "class_marker", "selection", "clip", "best_fits"):
            context.merge_gaussian_validation.celltypes[pin] = "plain"
        context.merge_gaussian_validation.celltypes.model_json = "bytes"
        context.merge_gaussian_validation.celltypes.class_validation = "plain"
        context.merge_gaussian_validation.celltypes.result = "plain"
        context.gaussian_validation = context.merge_gaussian_validation.result

        context.validate_refe_anchors = _stage6_validate_refe_anchors
        context.validate_refe_anchors.modules.preparic_stage6 = stage6_module
        context.validate_refe_anchors.pins.threshold_edits = context.threshold_edits
        context.validate_refe_anchors.pins.best_fit_thresholds = context.best_fit_thresholds
        context.validate_refe_anchors.pins.fragments = context.fragments
        context.validate_refe_anchors.pins.stackings = context.stackings
        context.validate_refe_anchors.pins.selection = context.stacking_selection
        context.validate_refe_anchors.pins.clip = context.best_fit_clip
        context.validate_refe_anchors.pins.best_fits = context.best_fits
        for pin in ("threshold_edits", "best_fit_thresholds", "fragments", "stackings", "selection", "clip", "best_fits"):
            context.validate_refe_anchors.celltypes[pin] = "plain"
        context.validate_refe_anchors.celltypes.result = "plain"
        context.refe_validation = context.validate_refe_anchors.result

        context.pdbcode.share(readonly=False)
        context.biounit.share(readonly=False)
        context.cif.share(path="cif.cif", readonly=False)
        context.ranges.share(readonly=False)
        context.range_specs.share(readonly=True)
        context.stackings.share()
        context.range_pdbs_aa.share(path="range_pdbs_aa")
        context.fragments.share(path="fragments")
        context.reference_fragment_pdbs.share()
        context.best_fit_thresholds.mount(
            MOUNT_DIRECTORY / "best-fit-thresholds.json", mode="w", authority="cell"
        )
        context.best_fit_thresholds.share()
        context.best_fits.share()
        context.stacking_selection.share(readonly=False)
        context.threshold_edits.share(readonly=False)
        context.best_fit_clip.share(readonly=False)
        context.interaction_matrix.share(readonly=False)
        context.gaussian_classes.share(readonly=False)
        context.energy_threshold_edits.share(readonly=False)
        context.energetics.share()
        context.energy_thresholds.mount(
            MOUNT_DIRECTORY / "computed-energy-thresholds.json", mode="w", authority="cell"
        )
        context.energy_thresholds.share()
        context.interactions.mount(
            MOUNT_DIRECTORY / "interactions.json", mode="w", authority="cell"
        )
        context.interactions.share()
        context.bestfit_stackings.share()
        context.gaussian_validation.mount(
            MOUNT_DIRECTORY / "gaussian-validation.json", mode="w", authority="cell"
        )
        context.gaussian_validation.share()
        context.refe_validation.mount(
            MOUNT_DIRECTORY / "refe-validation.json", mode="w", authority="cell"
        )
        context.refe_validation.share()
        context.gaussian_scores.share()
        context.structure_id_result.share(path="structure_id")
        context.mounts.sync(timeout=60)
        global _STAGE5_RANGE_SYNC_STOP, _STAGE5_RANGE_SYNC_THREAD
        if _STAGE5_RANGE_SYNC_STOP is not None:
            _STAGE5_RANGE_SYNC_STOP.set()
        _STAGE5_RANGE_SYNC_STOP = threading.Event()
        _STAGE5_RANGE_SYNC_THREAD = threading.Thread(
            target=_reconcile_stage5_range_scorers,
            args=(context, alaric_dir, _STAGE5_RANGE_SYNC_STOP),
            name="Preparic-stage5-range-reconciler",
            daemon=True,
        )
        _STAGE5_RANGE_SYNC_THREAD.start()
        ctx = context
        return context
    except BaseException:
        try:
            context.close()
        finally:
            seamless.close()
        raise


def close_context() -> None:
    """Close the module-level context and its share server."""
    global ctx, _STAGE5_RANGE_SYNC_STOP, _STAGE5_RANGE_SYNC_THREAD
    if _STAGE5_RANGE_SYNC_STOP is not None:
        _STAGE5_RANGE_SYNC_STOP.set()
    if _STAGE5_RANGE_SYNC_THREAD is not None and _STAGE5_RANGE_SYNC_THREAD is not threading.current_thread():
        _STAGE5_RANGE_SYNC_THREAD.join(timeout=2)
    _STAGE5_RANGE_SYNC_STOP = None
    _STAGE5_RANGE_SYNC_THREAD = None
    if ctx is not None:
        ctx.close()
        ctx = None
    seamless.close()


def _interactive() -> bool:
    """True under ``python -i`` and under ``ipython -i`` (which leaves sys.flags alone)."""
    if sys.flags.interactive:
        return True
    ipython = sys.modules.get("IPython")
    return ipython is not None and ipython.get_ipython() is not None


def _close_on_sigterm(signum, frame) -> None:
    """Close the context, then die from SIGTERM as if it had not been handled.

    Raising SystemExit instead makes prompt_toolkit (IPython) print a traceback.
    """
    close_context()
    signal.signal(signum, signal.SIG_DFL)
    os.kill(os.getpid(), signum)


ctx: Context | None = None


if __name__ == "__main__":
    ctx = build_context()
    # Close when the interpreter exits, e.g. through exit() in a shell. These
    # run before seamless's own atexit hook, which would warn. They are
    # registered as objects, since IPython clears this module's globals before
    # atexit hooks run; both are safe to call again after close_context().
    atexit.register(seamless.close)
    atexit.register(ctx.close)
    signal.signal(signal.SIGTERM, _close_on_sigterm)
    print(f"PDB code share: {ctx.pdbcode.share.url}")
    print(f"Biounit share: {ctx.biounit.share.url}")
    print(f"CIF share: {ctx.cif.share.url}")
    print(f"Ranges share: {ctx.ranges.share.url}")
    if not _interactive():
        try:
            input("Press Enter to stop the Seamless context... ")
        except EOFError:
            pass
        finally:
            close_context()
