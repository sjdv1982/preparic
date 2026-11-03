"""Build sandboxed Alaric transformations from the single ALARIC_DIR setting.

This module imports no installed Alaric package. Factory source is loaded by
file path, while scientific CLI code and helper modules are checksum-pinned or
attached to Seamless transformations for execution on the worker.
"""

from __future__ import annotations

import importlib.util
import shlex
from pathlib import Path
from types import ModuleType

from alaric_dir import resolve_alaric_repo_root


def _load_source(path: Path, module_name: str) -> ModuleType:
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load Alaric factory source {path}.")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def best_fit_factory(alaric_dir: str | Path):
    """Load the stdlib/Seamless-only T5 factory without importing `alaric`."""
    root = resolve_alaric_repo_root(alaric_dir)
    return _load_source(root / "util" / "seamlessify_refe_best_fit.py", "preparic_refe_best_fit_factory")


def fragment_library_identities(alaric_dir: str | Path) -> dict:
    """Use Alaric's delivered checksum-sidecar reader through its source file."""
    root = resolve_alaric_repo_root(alaric_dir)
    support = _load_source(
        root / "util" / "refe_best_fit_seamless_support.py",
        "preparic_refe_best_fit_support",
    )
    return support.identities(root / "fraglib")


def pinned_cli_transformer(
    alaric_dir: str | Path,
    deps_filename: str,
    command: str,
    dynamic_pins: dict[str, str],
    *,
    conda_environment: str = "alaric",
    meta_file: Path | None = None,
):
    """Build a CLI transformer with flat checksum pins and a declared layout.

    The DEPS entries are pinned by content, then reconstructed under the
    relative `util/`, `alaric/`, `fraglib/`, and `attract-jax/` layout expected
    by Alaric's command-line tools. Runtime project data remains on unbound pins.
    """
    from seamless import Checksum
    from seamless_transformer import Transformer
    from seamless_transformer.cmd.file_load import files_to_checksums

    root = resolve_alaric_repo_root(alaric_dir)
    util_dir = root / "util"
    deps_path = util_dir / deps_filename
    entries = deps_path.read_text(encoding="utf-8").splitlines()
    if not entries or any(not line or line != line.strip() or line.startswith("#") for line in entries):
        raise ValueError(f"{deps_filename} must contain one relative dependency per line.")

    # The manifest is source-side factory input only. Stage the eight named
    # files themselves so a worker never needs an unstaged .DEPS.txt file.
    declared = [(util_dir / entry).resolve() for entry in entries]
    destinations = {}
    sources = {}
    for path in declared:
        try:
            relative = path.relative_to(root)
        except ValueError as exc:
            raise ValueError(f"Dependency {path} is outside Alaric's repository root.") from exc
        if relative.parts[0] not in {"util", "alaric", "fraglib", "attract-jax"}:
            raise ValueError(f"Unexpected Alaric dependency location {relative}.")
        if relative.parts[0] == "fraglib" and len(relative.parts) > 1:
            destination = Path(*relative.parts)
        else:
            destination = relative
        destinations[path.name] = destination
        sources[path.name] = path
    if len(destinations) != len(declared):
        raise ValueError(f"{deps_filename} contains duplicate dependency basenames.")
    pin_names = set(destinations)
    if pin_names & set(dynamic_pins):
        raise ValueError("Static and project input pins have colliding names.")

    checksums, _ = files_to_checksums(
        [str(path) for path in declared],
        directories={},
        max_upload_files=None,
        max_upload_size=None,
        auto_confirm="yes",
    )
    stage = [
        "set -euo pipefail",
        "mkdir -p work/util work/alaric work/fraglib/templates",
    ]
    for pin_name, destination in sorted(destinations.items()):
        target = Path("work") / destination
        stage.append(f"mkdir -p {shlex.quote(str(target.parent))}")
        stage.append(f"cp {shlex.quote(pin_name)} {shlex.quote(str(target))}")
    stage.extend(command.splitlines())
    transformer = Transformer("bash")
    transformer.code = "\n".join(stage)
    for pin_name, path in destinations.items():
        transformer.celltypes[pin_name] = "bytes"
        transformer.args[pin_name].checksum = Checksum(checksums[str(sources[pin_name])])
    for pin_name, celltype in dynamic_pins.items():
        transformer.celltypes[pin_name] = celltype
    if meta_file is not None:
        meta_file = Path(meta_file).expanduser().resolve()
        meta_checksums, _ = files_to_checksums(
            [str(meta_file)], directories={}, max_upload_files=None,
            max_upload_size=None, auto_confirm="yes",
        )
        meta_pin = "META__FILE__fraglib.yaml"
        transformer.celltypes[meta_pin] = "bytes"
        transformer.args[meta_pin].checksum = Checksum(meta_checksums[str(meta_file)])
    transformer.celltypes.result = "folder"
    transformer.environment.set_conda_env(conda_environment)
    return transformer


def refe_fragments_transformer(alaric_dir: str | Path):
    command = """cp rna_pdb work/rna.pdb
cd work
python -B util/refe-fragments.py rna.pdb --templates fraglib/templates --output ../RESULT"""
    return pinned_cli_transformer(
        alaric_dir,
        "refe-fragments.DEPS",
        command,
        {"rna_pdb": "bytes"},
    )


def best_fit_tables_transformer(alaric_dir: str | Path, fraglib_config: Path):
    command = """cp rna_pdb work/rna.pdb
cp fraglib_identities work/fraglib_identities.json
cp fraglib.yaml work/fraglib.yaml
cp -a best_fit work/best_fit
cd work
python -B util/refe-best-fit-tables.py best_fit/best-fits.tsv rna.pdb \
  --exclude-pdb-code "$(cat ../exclude_pdb_code)" \
  --fraglib-identities fraglib_identities.json \
  --fraglib-config fraglib.yaml \
  --templates-dir fraglib/templates \
  --output ../RESULT"""
    return pinned_cli_transformer(
        alaric_dir,
        "refe-best-fit-tables.DEPS",
        command,
        {
            "rna_pdb": "bytes",
            "best_fit": "folder",
            "exclude_pdb_code": "str",
            "fraglib_identities": "str",
        },
        meta_file=fraglib_config,
    )


def energetics_range_transformer(alaric_dir: str | Path):
    """Build one independent T13 scorer for a single protein range."""
    command = """mkdir -p RESULT
cp -a ligands work/ligands
cp receptor work/receptor-aa.pdb
cd work
python -B util/energetics_seamless_support.py score ligands receptor-aa.pdb
cp energies.json ../RESULT/energetics.json
"""
    return pinned_cli_transformer(
        alaric_dir,
        "energetics_seamless_support.py.DEPS.txt",
        command,
        {"ligands": "folder", "receptor": "bytes"},
        conda_environment="alaric",
    )


def energy_threshold_module_code(alaric_dir: str | Path) -> tuple[dict, dict]:
    """Return Alaric's stdlib rounding module and its decimal helper."""
    root = resolve_alaric_repo_root(alaric_dir)
    return {
        "type": "interpreted",
        "language": "python",
        "code": (root / "alaric" / "energy_threshold.py").read_text(encoding="utf-8"),
    }, {
        "type": "interpreted",
        "language": "python",
        "code": (root / "alaric" / "threshold_rounding.py").read_text(encoding="utf-8"),
    }


def threshold_module_code(alaric_dir: str | Path) -> tuple[dict, dict]:
    root = resolve_alaric_repo_root(alaric_dir)
    thresholds = (root / "util" / "refe-thresholds.py").read_text(encoding="utf-8")
    rounding = (root / "alaric" / "threshold_rounding.py").read_text(encoding="utf-8")
    return (
        {"type": "interpreted", "language": "python", "code": thresholds},
        {"type": "interpreted", "language": "python", "code": rounding},
    )


def thresholds_from_tables(tables: dict[str, bytes], fit_echo: dict) -> dict:
    """Invoke the attached Alaric T8 module over T7's staged TSV bytes."""
    import tempfile
    from pathlib import Path

    import refe_thresholds

    if not isinstance(tables, dict) or not isinstance(tables.get("pairs.tsv"), bytes) \
        or not isinstance(tables.get("bases.tsv"), bytes):
        raise ValueError("The best-fit table folder is incomplete.")
    with tempfile.TemporaryDirectory(prefix="preparic-t8-") as temporary:
        directory = Path(temporary)
        pairs_path = directory / "pairs.tsv"
        bases_path = directory / "bases.tsv"
        pairs_path.write_bytes(tables["pairs.tsv"])
        bases_path.write_bytes(tables["bases.tsv"])
        result = refe_thresholds.build_thresholds(pairs_path, bases_path)
        result["echo"] = fit_echo
        return result
