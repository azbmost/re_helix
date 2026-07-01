"""Shared Tk GUI helpers for bundled re_helix tools."""

from __future__ import annotations

from pathlib import Path
from typing import Optional


def _repo_root_for(source_file: str) -> Path:
    source_dir = Path(source_file).resolve().parent
    if source_dir.name == "re_helix_lib":
        return source_dir.parent
    return source_dir


def apply_optional_icon(root, source_file: str) -> Optional[Path]:
    """Apply the optional re_helix icon to a Tk root, if an icon file exists."""
    try:
        import tkinter as tk
    except Exception:
        return None

    source_dir = Path(source_file).resolve().parent
    repo_root = _repo_root_for(source_file)
    candidates = [
        repo_root / "assets" / "icon.png",
        source_dir / "assets" / "icon.png",
        source_dir / "icon.png",
        repo_root / "re_helix_lib" / "icon.png",
    ]

    seen = set()
    for icon_path in candidates:
        if icon_path in seen:
            continue
        seen.add(icon_path)
        if not icon_path.exists():
            continue
        try:
            icon_image = tk.PhotoImage(master=root, file=str(icon_path))
            root.iconphoto(True, icon_image)
            root._re_helix_icon = icon_image
            root._re_helix_icon_path = str(icon_path)
            return icon_path
        except Exception:
            continue
    return None
