"""Generated PDB backups survive leaf eviction without restarting the backend."""

import gc
from pathlib import Path
import subprocess
import sys
import time

from seamless import Buffer, Cell
from seamless.caching.buffer_cache import get_buffer_cache
from seamless_transformer import Transformer
from seamless_workflow import Context
from seamless_workflow.attachments.fs.service import FileSystemService

from alaric_dir import resolve_alaric_dir, resolve_alaric_repo_root
from result_backups import mount_result_backup
from stage4 import _list_fragments, _parse_best_fits


def _clip_from_backed_folder(best_fits, refe_fragments, thresholds, clip, fragments):
    import preparic_stage4
    return preparic_stage4._clip_fragments(
        best_fits, refe_fragments, thresholds, clip, fragments,
    )


def _index(contents):
    index = {}
    for name, content in contents.items():
        buffer = Buffer(content)
        buffer.tempref()
        index[name] = buffer.get_checksum()
    return index


def _evict_unowned_leaves(checksums):
    cache = get_buffer_cache()
    for checksum in checksums:
        assert cache.refholder_counts.get(checksum, 0) == 0
        entry = cache.strong_cache.get(checksum)
        if entry is not None and entry.tempref is not None:
            entry.tempref.last_refreshed -= 60
    del entry
    gc.collect()
    cache.run_eviction_once()
    gc.collect()
    assert all(cache.get(checksum) is None for checksum in checksums)


def _wait_for_backup(cell):
    deadline = time.monotonic() + 10
    while time.monotonic() < deadline:
        if cell.mount.status["in_sync"]:
            return
        time.sleep(.01)
    raise AssertionError(cell.mount.status)


def test_t5_t6_backups_recover_evicted_1b7f_leaves_during_repeated_clips(tmp_path, monkeypatch):
    # Readback would add leaf leases and hide the missing-file fallback. Wait
    # only for the acknowledged writes, then prove the leaves leave memory.
    monkeypatch.setattr(FileSystemService, "poll_interval", 100)
    root = resolve_alaric_repo_root(resolve_alaric_dir())
    fixture = root / "tests" / "1b7f"
    reference_dir = tmp_path / "reference-fragments"
    subprocess.run([
        sys.executable, str(root / "util" / "refe-fragments.py"),
        str(fixture / "pdbs" / "rna.pdb"), "--templates", str(root / "fraglib" / "templates"),
        "--output", str(reference_dir),
    ], check=True, capture_output=True)
    refe_files = {path.name: path.read_bytes() for path in reference_dir.glob("*.pdb")}
    fit_files = {"best-fits.tsv": (fixture / "data" / "refe-best-fit.tsv").read_bytes()}
    fit_files.update({f"pdbs/{path.name}": path.read_bytes()
                      for path in (fixture / "best-fits").glob("*-bestfit.pdb")})
    refe_index, fit_index = _index(refe_files), _index(fit_files)
    assert refe_index["frag-10-UU-refe.pdb"].hex() == "c489a9be57b0047cb542295748b21c998d0eb022d2e02e15c0084eb5f5583fe4"
    assert refe_index["frag-11-UU-refe.pdb"].hex() == "906b4b421fa7a373efe658c87171a51fca0cad73bbdd419961e5b48b773d14e2"
    fit_leaf = fit_index["pdbs/frag-11-UU-bestfit.pdb"]
    assert fit_leaf.hex() == "78b6ee3bb009151bb4f0e0f725d85a41337e11a5ac9ab0325bd5da097c5e0b7c"
    fragments = _list_fragments(refe_files, (fixture / "pdbs" / "rna.pdb").read_bytes(),
                                "1B7F", {"chain": "B", "start": 1, "end": 12})
    basis = fragments["echo"]
    thresholds = {
        "schema": 1, "basis": basis, "active": True,
        "pairs": [{**row, "ovrmsd": "0.25", "crmsd": "0.25"}
                  for row in fragments["adjacent_pairs"]],
    }
    context = Context()
    try:
        context.refe_source = Cell(celltype="deepfolder")
        context.refe_source.set(refe_index)
        context.fit_source = Cell(celltype="deepfolder")
        context.fit_source.set(fit_index)
        refe_backup = mount_result_backup(context, "refe_fragments_disk_backup", context.refe_source)
        fit_backup = mount_result_backup(context, "best_fit_disk_backup", context.fit_source)
        for name in ("refe_fragments_disk_backup", "best_fit_disk_backup"):
            cell = getattr(context, name)
            assert cell.celltype == "folder"
            assert cell.mount.spec.mode == "w"
            assert cell.mount.spec.authority == "cell"
            assert cell.mount.spec.persistent is False
            _wait_for_backup(cell)

        context.best_fits = _parse_best_fits(fit_files, basis)
        context.fragments = fragments
        context.thresholds = thresholds
        context.clip = {"schema": 1, "basis": basis, "clip1": 1, "clip2": 1}
        context.clip_fragments = _clip_from_backed_folder
        context.clip_fragments.modules.preparic_stage4 = {
            "type": "interpreted", "language": "python",
            "code": (Path(__file__).resolve().parents[2] / "stage4.py").read_text(),
        }
        for pin in ("best_fits", "fragments", "thresholds", "clip"):
            context.clip_fragments.pins[pin] = getattr(context, pin)
            context.clip_fragments.celltypes[pin] = "plain"
        context.clip_fragments.celltypes.refe_fragments = "folder"
        context.clip_fragments.pins.refe_fragments = context.refe_source
        context.clip_fragments.celltypes.result = "folder"

        boundary = Transformer("bash")
        boundary.code = "sha256sum best_fit/pdbs/frag-11-UU-bestfit.pdb > RESULT"
        boundary.celltypes.best_fit = "folder"
        boundary.celltypes.revision = "plain"
        boundary.celltypes.result = "text"
        context.gaussian_boundary = boundary
        context.gaussian_revision = 0
        context.gaussian_boundary.pins.best_fit = context.fit_source
        context.gaussian_boundary.pins.revision = context.gaussian_revision
        assert context.clip_fragments.compute(timeout=10) is not None
        assert context.gaussian_boundary.compute(timeout=10) is not None
        assert context.gaussian_boundary.result.value.startswith(fit_leaf.hex())
        source_checksums = (context.refe_source.checksum, context.fit_source.checksum)

        for revision, clip2 in enumerate((0, 1), 1):
            _evict_unowned_leaves([*refe_index.values(), *fit_index.values()])
            # Resolving the exact missing hashes reads the mounted backups.
            assert refe_index["frag-11-UU-refe.pdb"].resolve().content == refe_files["frag-11-UU-refe.pdb"]
            assert fit_leaf.resolve().content == fit_files["pdbs/frag-11-UU-bestfit.pdb"]
            context.clip.set({"schema": 1, "basis": basis, "clip1": 1, "clip2": clip2})
            context.gaussian_revision.set(revision)
            assert context.clip_fragments.compute(timeout=10) is not None
            assert context.clip_fragments.exception is None
            assert context.gaussian_boundary.compute(timeout=10) is not None
            assert context.gaussian_boundary.result.value.startswith(fit_leaf.hex())
            assert (context.refe_source.checksum, context.fit_source.checksum) == source_checksums
    finally:
        context.close()
    assert not refe_backup.exists()
    assert not fit_backup.exists()
