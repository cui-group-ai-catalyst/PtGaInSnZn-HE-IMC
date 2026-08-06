"""Path resolution for the Figure 3 analysis package.

All paths resolve relative to this file's parent directory (the Figure3/
root) so the package is self-contained on any machine. Set the environment
variable FIG3_ROOT to relocate the whole package, or FIG3_OUTPUT to redirect
generated files (useful for keeping canonical bundled files untouched).
"""

from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(os.environ.get("FIG3_ROOT", Path(__file__).resolve().parent))
ANALYSIS = ROOT / "analysis"
BUILD = ROOT / "build"
DATA = ROOT / "data"
SOURCE_DATA = ROOT / "source_data"
STATIC = ROOT / "static"
INTERACTIVE_3D = ROOT / "interactive_3d"
OUTPUT = Path(os.environ.get("FIG3_OUTPUT", ROOT))
