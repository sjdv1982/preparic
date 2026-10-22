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
from seamless_workflow import Context, shareserver


PROJECT_ROOT = Path(__file__).resolve().parent
MOUNT_DIRECTORY = PROJECT_ROOT / "deployment-mounts"
PDBCODE_PATH = MOUNT_DIRECTORY / "pdbcode.txt"
BIOUNIT_PATH = MOUNT_DIRECTORY / "biounit.json"
RANGES_PATH = MOUNT_DIRECTORY / "ranges.yaml"


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
        context.pdbcode.mount(PDBCODE_PATH, mode="rw", authority="cell")

        context.biounit = Cell(celltype="plain")
        context.biounit.set(_initial_biounit_index(BIOUNIT_PATH))
        context.biounit.mount(BIOUNIT_PATH, mode="rw", authority="cell")

        context.cif = Cell(celltype="text")
        context.cif.set("")

        initial_ranges = _initial_ranges(RANGES_PATH)
        context.ranges = Cell(celltype="yaml")
        context.ranges.set(yaml.safe_dump(initial_ranges, sort_keys=False))
        context.ranges.mount(RANGES_PATH, mode="rw", authority="cell")

        # The UI keeps its mounted YAML file, while the transformer receives
        # the corresponding structured value through this plain cell.
        context.range_specs = Cell(celltype="plain")
        context.range_specs.set(initial_ranges)

        context.ranges_to_pdb = _ranges_to_pdb
        context.ranges_to_pdb.pins.cif = context.cif
        context.ranges_to_pdb.pins.ranges = context.range_specs
        context.range_pdbs = Cell(celltype="plain")
        context.range_pdbs = context.ranges_to_pdb.result

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
