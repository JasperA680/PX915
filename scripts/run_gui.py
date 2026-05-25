"""Launch the PyQt5 GUI for the PX915 traffic network simulator.

Run from the repo root:
    .venv/bin/python scripts/run_gui.py
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from python.gui.main_window import main

if __name__ == "__main__":
    main()
