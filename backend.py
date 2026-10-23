"""Run Preparic's PDB, biounit, CIF, and range shares during development.

Use ``python -i backend.py`` to leave the context available in an interactive
prompt; call ``close_context()`` when you are done.
"""

import atexit
import json
import os
import signal
import sys
from pathlib import Path

import seamless
import yaml
from seamless import Cell
from seamless_transformer import Transformer
from seamless_workflow import Context, shareserver
from range_pdbs import _pdbs_to_folder, _ranges_to_pdb


PROJECT_ROOT = Path(__file__).resolve().parent
MOUNT_DIRECTORY = PROJECT_ROOT / "deployment-mounts"
PDBCODE_PATH = MOUNT_DIRECTORY / "pdbcode.txt"
BIOUNIT_PATH = MOUNT_DIRECTORY / "biounit.json"
RANGES_PATH = MOUNT_DIRECTORY / "ranges.yaml"
PDBS_PATH = MOUNT_DIRECTORY / "pdbs"
PDBS_AA_PATH = MOUNT_DIRECTORY / "pdbs-aa"


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


def _initial_ranges(path: Path) -> dict:
    if not path.exists():
        return {"ranges": []}
    try:
        # BaseLoader avoids YAML 1.1 coercions such as chain id "NO" to False.
        value = yaml.load(path.read_text(encoding="utf-8"), Loader=yaml.BaseLoader)
    except (OSError, yaml.YAMLError):
        return {"ranges": []}
    if isinstance(value, dict) and isinstance(value.get("ranges"), list):
        return value
    return {"ranges": []}


def build_context() -> Context:
    """Build the live context backing the Vite development proxy."""
    global ctx
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

        initial_ranges = _initial_ranges(RANGES_PATH)
        context.ranges = Cell(celltype="yaml")
        context.ranges.set(yaml.safe_dump(initial_ranges, sort_keys=False))
        context.ranges.mount(RANGES_PATH, mode="rw", authority="file")

        # The UI keeps its mounted YAML file, while the transformer receives
        # the corresponding structured value through this plain cell.
        context.range_specs = Cell(celltype="plain")
        context.range_specs.set(initial_ranges)

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

        context.pdbcode.share(readonly=False)
        context.biounit.share(readonly=False)
        # The dotted path disables the browser client's automatic CIF read.
        context.cif.share(path="cif.cif", readonly=False)
        context.ranges.share(readonly=False)
        context.range_specs.share(readonly=False)
        context.range_pdbs.share()
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
