"""Offscreen smoke test for the PyQt5 GUI.

Cycles through every preset, runs a short simulation via the RunnerThread,
asserts the plot panel populates without exceptions.

Run from repo root:
    QT_QPA_PLATFORM=offscreen .venv/bin/python tests/test_gui_smoke.py
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

# Default to offscreen unless the caller overrode it.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

from PyQt5.QtWidgets import QApplication
from PyQt5.QtCore import QEventLoop

from python.gui.main_window import MainWindow


BINARY = ROOT / "build" / "run_network"


def main():
    if not BINARY.exists():
        sys.exit(f"missing binary {BINARY} — run `make run_network` first")

    app = QApplication.instance() or QApplication(["gui_smoke"])

    win = MainWindow(binary=BINARY, output_dir=ROOT / "data" / "output" / "gui_smoke")

    for preset in ("single_lane", "crossroads", "roundabout", "town"):
        print(f"--- {preset} ---")
        win.preset_combo.setCurrentText(preset)
        assert win._current_spec is not None, "preset did not build a spec"
        assert win._current_layout is not None, "preset did not build a layout"

        # Shrink the run for the test.
        win.param_form.n_steps.setValue(150)
        win.param_form.lane_length.setValue(12)

        loop = QEventLoop()
        state = {"ok": False, "err": None, "path": None}

        def on_ok(p):
            state["ok"] = True
            state["path"] = p
            loop.quit()

        def on_err(msg):
            state["err"] = msg
            loop.quit()

        win._on_run()
        assert win._runner is not None
        win._runner.finished_ok.connect(on_ok)
        win._runner.finished_error.connect(on_err)
        loop.exec_()

        assert state["err"] is None, f"{preset}: runner error: {state['err']}"
        assert state["ok"], f"{preset}: did not finish OK"
        assert len(win.plot_panel.tab_density.ax.lines) > 0, \
            f"{preset}: density tab has no lines"

        print(f"  PASS: {state['path']}")

    print("ALL OK")


if __name__ == "__main__":
    main()
