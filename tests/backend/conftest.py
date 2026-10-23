"""Self-contained inputs shared by direct and standalone Seamless tests."""

from pathlib import Path

import pytest
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA = Path(__file__).parent / "data" / "1b7f"


@pytest.fixture(scope="session")
def cif_text():
    return (DATA / "structure.cif").read_text()


@pytest.fixture
def range_specs():
    return yaml.safe_load((DATA / "ranges.yaml").read_text())


@pytest.fixture(scope="session")
def expected_pdbs():
    return [(DATA / "pdbs" / f"range-{i}.pdb").read_text() for i in (1, 2, 3)]


@pytest.fixture(scope="session")
def aareduce_files():
    root = PROJECT_ROOT / "aareduce"
    return {
        path.relative_to(root).as_posix(): path.read_bytes()
        for path in root.rglob("*")
        if path.is_file() and "__pycache__" not in path.parts
    }


@pytest.fixture(scope="session")
def seamless_runtime():
    # Select local execution explicitly; never initialize a user's cluster config.
    import seamless
    from seamless_config.select import (
        get_execution, get_persistent, select_execution, select_persistent,
    )

    previous = get_execution(), get_persistent()
    select_execution("process")
    select_persistent(False)
    yield
    seamless.close()
    select_execution(previous[0])
    select_persistent(previous[1])


@pytest.fixture
def aareduce_folder(seamless_runtime, aareduce_files):
    from seamless import Buffer, Cell

    # A folder cell holds an index. Keep explicit owners of the file buffers
    # alive for the test; no filesystem mounts or workflow are involved.
    buffers = {name: Buffer(data) for name, data in aareduce_files.items()}
    cell = Cell(celltype="folder")
    cell.set({name: buffer.get_checksum() for name, buffer in buffers.items()})
    yield cell
