"""Backend-lifetime disk backups for long-lived generated folder inputs."""

from pathlib import Path
from tempfile import TemporaryDirectory
import weakref


def mount_result_backup(context, name: str, source) -> Path:
    """Mount a folder view without changing its producer or consumer edges.

    The write-only mount cannot adopt an old result on startup. Its files are
    removed by Context.close(); the finalizer also cleans up if the Context is
    discarded before an explicit close.
    """
    directory = TemporaryDirectory(prefix=f"preparic-{name}-")
    try:
        setattr(context, name, source.as_celltype("folder"))
        getattr(context, name).mount(
            directory.name, mode="w", authority="cell", persistent=False,
        )
    except BaseException:
        directory.cleanup()
        raise
    weakref.finalize(context, directory.cleanup)
    return Path(directory.name)
