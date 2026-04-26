"""Shared pytest fixtures."""

from __future__ import annotations

import sys
from pathlib import Path

# Ensure the project root (parent of `src/`) is on sys.path so tests can
# import the `src` package directly.
_ROOT = Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
