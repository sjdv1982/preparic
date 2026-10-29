"""Run Preparic's PDB, biounit, CIF, and range shares during development.

Use ``python -i backend.py`` to leave the context available in an interactive
prompt; call ``close_context()`` when you are done.
"""

import atexit
import copy
import json
import os
import re
import signal
import sys
from pathlib import Path

import seamless
import yaml
from seamless import Cell
from seamless_transformer import Transformer
from seamless_workflow import Context, shareserver
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

LOCAL_WORKER_COUNT = 8
_LOCAL_EXECUTION_CONFIGURED = False

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
            or not _unique(value["ranges"]) or type(value["grid"]) is not list:
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
                if threshold is not None and (not isinstance(threshold, str)
                    or re.fullmatch(r"-?(?:0|[1-9][0-9]*)\.[0-9]{2}", threshold) is None):
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


def _add_stacking_pipeline(
    context: Context, mount_directory: Path | None = None
) -> None:
    """Wire the pre-aareduce range PDBs into a content-bound stacking step."""
    if mount_directory is None:
        mount_directory = MOUNT_DIRECTORY
    context.detect_stacking = _detect_stacking
    context.detect_stacking.celltypes.range_pdbs = "plain"
    context.detect_stacking.celltypes.result = "plain"
    context.detect_stacking.pins.range_pdbs = context.range_pdbs

    context.stackings = Cell(celltype="plain")
    context.stackings = context.detect_stacking.result
    context.stackings.mount(
        mount_directory / "stackings.json", mode="w", authority="cell"
    )


def build_context() -> Context:
    """Build the durable project Context backing the Vite development proxy."""
    global ctx
    configure_local_execution()
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
        _add_stacking_pipeline(context)

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

        context.structure_id = _structure_id
        context.structure_id.pins.cif = context.cif
        context.structure_id.celltypes.result = "plain"
        context.structure_id_result = context.structure_id.result

        context.pdbcode.share(readonly=False)
        context.biounit.share(readonly=False)
        context.cif.share(path="cif.cif", readonly=False)
        context.ranges.share(readonly=False)
        context.range_specs.share(readonly=True)
        context.stackings.share()
        context.stacking_selection.share(readonly=False)
        context.threshold_edits.share(readonly=False)
        context.best_fit_clip.share(readonly=False)
        context.interaction_matrix.share(readonly=False)
        context.gaussian_classes.share(readonly=False)
        context.energy_threshold_edits.share(readonly=False)
        context.structure_id_result.share(path="structure_id")
        context.mounts.sync(timeout=10)
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
    global ctx
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
