"""Resolve the one configured Alaric package location used by Preparic."""

from __future__ import annotations

import os
from pathlib import Path


def resolve_alaric_dir(override: str | os.PathLike[str] | None = None) -> Path:
    """Return the configured ``alaric/`` package directory.

    The precedence is an explicit caller override, ``ALARIC_DIR``, then the
    documented ``~/alaric/alaric`` default. The resolved path also makes
    symlinked installations produce one stable source for all downstream
    Alaric paths.
    """
    configured = override if override is not None else os.environ.get("ALARIC_DIR")
    candidate = Path(configured).expanduser() if configured is not None else Path.home() / "alaric" / "alaric"
    resolved = candidate.resolve()
    if not resolved.is_dir():
        source = "--alaric-dir/ALARIC_DIR" if configured is not None else "the default location"
        raise FileNotFoundError(
            f"Alaric package directory {resolved} does not exist ({source}). "
            "Set ALARIC_DIR to the repository's alaric/ directory."
        )
    return resolved


def resolve_alaric_repo_root(override: str | os.PathLike[str] | None = None) -> Path:
    """Return the repository root containing the configured package folder."""
    return resolve_alaric_dir(override).parent
