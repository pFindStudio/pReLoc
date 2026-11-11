"""Locate bundled assets independently of the caller's working directory."""

from pathlib import Path
import sys


def resource_path(relative_path):
    path = Path(relative_path)
    if path.is_absolute():
        return str(path)
    root = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent))
    return str(root / path)
