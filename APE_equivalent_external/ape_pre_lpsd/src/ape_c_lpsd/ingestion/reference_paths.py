"""Resolve the read-only reference project's vendored Dewesoft interface.

This adapter deliberately does not alter the vendor implementation. It only
supplies its documented paths from the configured ``reference_project`` root.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


class ReferenceDewesoftPathError(RuntimeError):
    """The configured reference project cannot supply its Dewesoft wrapper."""


@dataclass(frozen=True)
class ReferenceDewesoftPaths:
    reference_project_root: Path
    wrapper_path: Path
    library_path: Path

    def source_kwargs(self) -> dict[str, str]:
        """Keyword arguments accepted by unchanged reference source classes."""
        return {"wrapper_path": str(self.wrapper_path), "library_path": str(self.library_path)}


def resolve_reference_dewesoft_paths(reference_project_root: Path) -> ReferenceDewesoftPaths:
    """Resolve and preflight the vendor paths relative to one configured root."""
    root = Path(reference_project_root).expanduser().resolve()
    wrapper = root / "vendor" / "dewesoft" / "5.0.8" / "DWDataReaderHeader.py"
    library = root / "vendor" / "dewesoft" / "5.0.8" / "DWDataReaderLib64.dylib"
    if not root.is_dir() or not wrapper.is_file() or not library.is_file():
        raise ReferenceDewesoftPathError(
            "Dewesoft reference-path preflight failed. "
            f"Configured reference_project_root: {root}. "
            f"Expected Dewesoft wrapper path: {wrapper}. "
            f"Expected Dewesoft library path: {library}."
        )
    return ReferenceDewesoftPaths(root, wrapper, library)
