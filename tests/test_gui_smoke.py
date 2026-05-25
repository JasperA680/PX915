"""Offscreen smoke test for the PyQt5 GUI.

Cycles through every preset on the CA tab, runs a short simulation via the
RunnerThread, asserts the plot panel populates without exceptions.  Also
asserts the PDE placeholder tab is reachable.

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

from PyQt5.QtWidgets import QApplication, QLabel, QPushButton
from PyQt5.QtCore import QEventLoop

from python.gui.main_window import MainWindow
from python.gui.pde_plot_panel import PDEPlotPanel
from python.gui.pde_param_form import PDEParamForm


BINARY = ROOT / "build" / "run_network"


def run_one_preset(win: MainWindow, preset: str):
    print(f"--- {preset} ---")
    tab = win.ca_tab
    tab.preset_combo.setCurrentText(preset)
    assert tab._current_spec is not None, "preset did not build a spec"
    assert tab._current_layout is not None, "preset did not build a layout"
    # The network widget should have a drawn collection immediately,
    # before any sim has run.
    assert len(tab.network_widget._ax.collections) > 0, \
        f"{preset}: network preview is empty before run"

    tab.param_form.n_steps.setValue(150)

    loop = QEventLoop()
    state = {"ok": False, "err": None, "path": None}

    def on_ok(p):
        state["ok"] = True
        state["path"] = p
        loop.quit()

    def on_err(msg):
        state["err"] = msg
        loop.quit()

    tab._on_run()
    assert tab._runner is not None
    tab._runner.finished_ok.connect(on_ok)
    tab._runner.finished_error.connect(on_err)
    loop.exec_()

    assert state["err"] is None, f"{preset}: runner error: {state['err']}"
    assert state["ok"], f"{preset}: did not finish OK"
    # Density tab defaults to nothing selected (user must pick which roads /
    # lanes to display) — tick every road so the plot draws something we
    # can validate.
    from PyQt5.QtCore import Qt
    rs = tab.plot_panel.tab_density.road_selector
    for i in range(rs.count()):
        rs.model().item(i).setCheckState(Qt.Checked)
    rs.selectionChanged.emit(rs.checked_indices())
    assert len(tab.plot_panel.tab_density.ax.lines) > 0, \
        f"{preset}: density tab has no lines after ticking all roads"
    print(f"  PASS: {state['path']}")


def main():
    if not BINARY.exists():
        sys.exit(f"missing binary {BINARY} — run `make run_network` first")

    app = QApplication.instance() or QApplication(["gui_smoke"])

    win = MainWindow(binary=BINARY, output_dir=ROOT / "data" / "output" / "gui_smoke")

    # Tabs exist in the expected order.
    assert win.tabs.count() == 2, f"expected 2 tabs, got {win.tabs.count()}"
    assert win.tabs.tabText(0) == "Cellular Automaton (CA)"
    assert win.tabs.tabText(1) == "PDE Continuum Model"

    # PDE tab is wired up: form, plot panel, and Run button present.
    win.tabs.setCurrentIndex(1)
    assert isinstance(win.pde_tab.param_form, PDEParamForm), \
        "PDE tab missing PDEParamForm"
    assert isinstance(win.pde_tab.plot_panel, PDEPlotPanel), \
        "PDE tab missing PDEPlotPanel"
    assert win.pde_tab.run_button.text() == "Run", \
        "PDE tab missing Run button"
    assert win.pde_tab.plot_panel.count() == 4, \
        f"PDE plot panel should have 4 tabs, got {win.pde_tab.plot_panel.count()}"
    # Multi-lane-only tabs (Lane densities, Mass conservation) start disabled.
    for i in (2, 3):
        assert not win.pde_tab.plot_panel.isTabEnabled(i), \
            f"PDE plot tab {i} should be disabled before a run"
    print("--- PDE tab: PASS ---")
    win.tabs.setCurrentIndex(0)

    # Default model is NS.
    assert win.ca_tab.model_combo.currentText() == "NS"

    for preset in ("single_lane", "two_lane", "t_junction",
                   "crossroads", "roundabout", "town"):
        run_one_preset(win, preset)

    # TASEP should be greyed out for non-single-lane presets.
    win.ca_tab.preset_combo.setCurrentText("two_lane")
    tasep_idx = win.ca_tab.model_combo.findText("TASEP")
    tasep_item = win.ca_tab.model_combo.model().item(tasep_idx)
    assert tasep_item is not None and not (tasep_item.flags() & 0x20), \
        "TASEP should be disabled in the model combo when two_lane is selected"
    assert win.ca_tab.model_combo.currentText() == "NS", \
        "model should auto-switch to NS when TASEP is disabled"
    print("--- TASEP grey-out: PASS ---")

    # TASEP is now implemented Fortran-side; check that single_lane + TASEP
    # actually completes (this used to assert an error).
    print("--- TASEP single_lane ---")
    win.ca_tab.preset_combo.setCurrentText("single_lane")
    win.ca_tab.model_combo.setCurrentText("TASEP")
    win.ca_tab.param_form.n_steps.setValue(50)

    loop = QEventLoop()
    state = {"err": None, "ok": False}

    def on_ok(p):
        state["ok"] = True
        loop.quit()

    def on_err(msg):
        state["err"] = msg
        loop.quit()

    win.ca_tab._on_run()
    win.ca_tab._runner.finished_ok.connect(on_ok)
    win.ca_tab._runner.finished_error.connect(on_err)
    loop.exec_()

    assert state["err"] is None, f"TASEP single_lane run failed: {state['err']}"
    assert state["ok"], "TASEP single_lane did not finish OK"
    print(f"  PASS: TASEP single_lane completed")

    print("ALL OK")


if __name__ == "__main__":
    main()
