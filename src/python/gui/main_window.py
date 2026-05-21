"""Main window: top-level CA / PDE tabs, status bar, run log."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QComboBox,
    QPushButton, QProgressBar, QLabel, QStatusBar, QPlainTextEdit, QTabWidget,
    QDockWidget, QScrollArea,
)

from python.road_network import PRESETS, NetworkSpec, LayoutSpec
from python.io import load_network_netcdf
from python.pde_runner import load_pde_netcdf

from python.gui.param_form import ParamForm
from python.gui.network_widget import NetworkWidget
from python.gui.runner_thread import RunnerThread
from python.gui.plot_panel import PlotPanel
from python.gui.pde_param_form import PDEParamForm, PRESETS_ORDER, CUSTOM_LABEL
from python.gui.pde_runner_thread import PDERunnerThread
from python.gui.pde_plot_panel import PDEPlotPanel
from python.gui.road_edit_dialog import RoadEditDialog


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_BINARY = REPO_ROOT / "build" / "run_network"
DEFAULT_OUTDIR = REPO_ROOT / "data" / "output" / "gui"
DEFAULT_PDE_BINARY = REPO_ROOT / "build" / "pde_solver"
DEFAULT_PDE_OUTDIR = REPO_ROOT / "data" / "output" / "gui_pde"

# Available CA update rules.  Only "NS" is implemented physics-side; "TASEP"
# is plumbed through the JSON config but the Fortran driver errors out on it.
CA_MODELS = ("NS", "TASEP")


class CATab(QWidget):
    """Cellular-automaton tab: preset selector, model toggle, network preview, run + plots."""

    log = pyqtSignal(str)
    status = pyqtSignal(str)
    progress = pyqtSignal(int)

    def __init__(self, binary: Path, output_dir: Path, parent=None):
        super().__init__(parent)
        self._binary = binary
        self._output_dir = output_dir
        self._runner: Optional[RunnerThread] = None
        self._current_spec: Optional[NetworkSpec] = None
        self._current_layout: Optional[LayoutSpec] = None
        # Per-road parameter overrides set via the click-to-edit dialog.
        # Keyed by road_id; reset when the preset combo changes.
        self._road_overrides: dict = {}

        # --- Toolbar row: preset + model + run ---
        toolbar = QWidget()
        toolbar_row = QHBoxLayout(toolbar)
        toolbar_row.addWidget(QLabel("Preset:"))
        self.preset_combo = QComboBox()
        for name in PRESETS:
            self.preset_combo.addItem(name)
        toolbar_row.addWidget(self.preset_combo)
        toolbar_row.addSpacing(12)
        toolbar_row.addWidget(QLabel("Model:"))
        self.model_combo = QComboBox()
        for m in CA_MODELS:
            self.model_combo.addItem(m)
        toolbar_row.addWidget(self.model_combo)
        toolbar_row.addSpacing(12)
        self.run_button = QPushButton("Run")
        self.run_button.setMinimumWidth(120)
        toolbar_row.addWidget(self.run_button)
        toolbar_row.addStretch(1)

        # --- Left: params (scrollable, capped height) + network preview (fills rest) ---
        # The form has ~10 controls (~450 px tall natural).  Wrapping it in a
        # QScrollArea with maximumHeight=340 lets the network widget claim the
        # rest of the vertical space without the form ever pushing it out.
        self.param_form = ParamForm()
        form_scroll = QScrollArea()
        form_scroll.setWidget(self.param_form)
        form_scroll.setWidgetResizable(True)
        form_scroll.setMaximumHeight(340)
        form_scroll.setFrameShape(QScrollArea.NoFrame)

        self.network_widget = NetworkWidget()
        self.network_widget.setMinimumSize(360, 360)

        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.addWidget(form_scroll)
        left_layout.addWidget(self.network_widget, stretch=1)

        # --- Right: plot panel ---
        self.plot_panel = PlotPanel()

        # --- Central layout ---
        body = QSplitter(Qt.Horizontal)
        body.addWidget(left)
        body.addWidget(self.plot_panel)
        body.setStretchFactor(0, 0)
        body.setStretchFactor(1, 1)
        body.setSizes([520, 800])

        outer = QVBoxLayout(self)
        outer.addWidget(toolbar)
        outer.addWidget(body, stretch=1)

        # --- Signals ---
        self.preset_combo.currentTextChanged.connect(self._on_preset_changed_ca)
        self.model_combo.currentTextChanged.connect(lambda _m: None)  # no rebuild needed
        self.param_form.params_changed.connect(self._rebuild_spec)
        self.run_button.clicked.connect(self._on_run)
        self.network_widget.road_clicked.connect(self._on_road_clicked)

        # Bootstrap with the first preset.
        self._rebuild_spec()

    def _on_preset_changed_ca(self, name: str):
        """Reset per-road overrides whenever the preset changes, then rebuild."""
        self._road_overrides = {}
        self._rebuild_spec()

    def _rebuild_spec(self, *_args):
        name = self.preset_combo.currentText()
        if name not in PRESETS:
            return
        kwargs = self.param_form.preset_kwargs(name)
        try:
            spec, layout = PRESETS[name](**kwargs)
        except TypeError as exc:
            self.log.emit(f"preset error: {exc}")
            return
        # Apply any per-road overrides set via the click-to-edit dialog.
        for rid, ov in self._road_overrides.items():
            road = next((r for r in spec.roads if r.id == rid), None)
            if road is None:
                continue
            for ln in road.lanes:
                if ov.get("alpha") is not None and getattr(ln, "open_in", False):
                    ln.alpha = ov["alpha"]
                if ov.get("beta") is not None and getattr(ln, "open_out", False):
                    ln.beta = ov["beta"]
            if ov.get("length") is not None:
                road.length = ov["length"]
        self._current_spec = spec
        self._current_layout = layout
        self.network_widget.set_network(spec, layout)
        self.status.emit(
            f"preset {name}: {len(spec.roads)} roads, {len(spec.junctions)} junctions"
        )

    def _on_road_clicked(self, road_id: int):
        """Open the per-road editor dialog and apply any changes."""
        if self._current_spec is None:
            return
        road = next((r for r in self._current_spec.roads if r.id == road_id), None)
        if road is None:
            return

        # Read current α/β from the road's lanes.
        alpha = next((ln.alpha for ln in road.lanes if getattr(ln, "open_in", False)), 0.5)
        beta  = next((ln.beta  for ln in road.lanes if getattr(ln, "open_out", False)), 0.5)
        has_alpha = any(getattr(ln, "open_in",  False) for ln in road.lanes)
        has_beta  = any(getattr(ln, "open_out", False) for ln in road.lanes)
        length = getattr(road, "length", 20)

        dlg = RoadEditDialog(road_id, alpha, beta, length,
                             has_alpha, has_beta, parent=self)
        if dlg.exec_() == RoadEditDialog.Accepted:
            vals = dlg.values()
            self._road_overrides[road_id] = vals
            self._rebuild_spec()
            self.log.emit(
                f"R{road_id} overridden: α={vals['alpha']}  β={vals['beta']}  L={vals['length']}"
            )

    def _on_run(self):
        if self._current_spec is None or self._current_layout is None:
            self.log.emit("no spec to run")
            return
        if not self._binary.exists():
            self.log.emit(f"binary not found: {self._binary} — try `make run_network`")
            return
        if self._runner is not None and self._runner.isRunning():
            return

        params = self.param_form.sim_params()
        params.model = self.model_combo.currentText()
        preset_name = self.preset_combo.currentText()
        out_dir = self._output_dir / preset_name
        self.run_button.setEnabled(False)
        self.progress.emit(0)
        self.status.emit(f"running {preset_name} model={params.model} ({params.n_steps} steps)...")
        self.log.emit(f"--- run {preset_name} model={params.model} ---")

        self._runner = RunnerThread(
            self._current_spec, params, self._current_layout,
            output_dir=out_dir, binary=self._binary, parent=self,
        )
        self._runner.progress.connect(self._on_runner_progress)
        self._runner.finished_ok.connect(self._on_finished_ok)
        self._runner.finished_error.connect(self._on_finished_error)
        self._runner.start()

    def _on_runner_progress(self, frac: float, msg: str):
        if frac >= 0:
            self.progress.emit(int(round(100 * frac)))
        if msg:
            self.log.emit(msg)

    def _on_finished_ok(self, nc_path: str):
        self.progress.emit(100)
        self.status.emit(f"done — {nc_path}")
        self.run_button.setEnabled(True)
        self.log.emit(f"loaded {nc_path}")
        try:
            result = load_network_netcdf(nc_path)
            self.plot_panel.show_result(result)
        except Exception as exc:
            self.log.emit(f"plot error: {exc}")

    def _on_finished_error(self, err: str):
        self.status.emit("error")
        self.run_button.setEnabled(True)
        self.log.emit(f"ERROR: {err}")

    def stop_running(self):
        if self._runner is not None and self._runner.isRunning():
            self._runner.requestInterruption()
            self._runner.wait(2000)


class PDETab(QWidget):
    """LWR continuum PDE tab: presets + parameter form + plot panel."""

    log = pyqtSignal(str)
    status = pyqtSignal(str)
    progress = pyqtSignal(int)

    def __init__(self, binary: Path, output_dir: Path, parent=None):
        super().__init__(parent)
        self._binary = binary
        self._output_dir = output_dir
        self._runner: Optional[PDERunnerThread] = None
        self._suppress_custom_revert = False

        # --- Toolbar: preset + run ---
        toolbar = QWidget()
        toolbar_row = QHBoxLayout(toolbar)
        toolbar_row.addWidget(QLabel("Preset:"))
        self.preset_combo = QComboBox()
        for name in PRESETS_ORDER:
            self.preset_combo.addItem(name)
        toolbar_row.addWidget(self.preset_combo)
        toolbar_row.addSpacing(12)
        self.run_button = QPushButton("Run")
        self.run_button.setMinimumWidth(120)
        toolbar_row.addWidget(self.run_button)
        toolbar_row.addStretch(1)

        # --- Left: parameter form in a capped scroll area ---
        self.param_form = PDEParamForm()
        form_scroll = QScrollArea()
        form_scroll.setWidget(self.param_form)
        form_scroll.setWidgetResizable(True)
        form_scroll.setMaximumHeight(440)
        form_scroll.setFrameShape(QScrollArea.NoFrame)

        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.addWidget(form_scroll)
        left_layout.addStretch(1)

        # --- Right: plot panel ---
        self.plot_panel = PDEPlotPanel()

        body = QSplitter(Qt.Horizontal)
        body.addWidget(left)
        body.addWidget(self.plot_panel)
        body.setStretchFactor(0, 0)
        body.setStretchFactor(1, 1)
        body.setSizes([520, 800])

        outer = QVBoxLayout(self)
        outer.addWidget(toolbar)
        outer.addWidget(body, stretch=1)

        # --- Signals ---
        self.preset_combo.currentTextChanged.connect(self._on_preset_changed)
        self.param_form.params_changed.connect(self._on_form_edited)
        self.param_form.params_changed.connect(self._sync_multilane_tabs)
        self.run_button.clicked.connect(self._on_run)

        # Match the initial form state.
        self._sync_multilane_tabs()

    def _sync_multilane_tabs(self):
        self.plot_panel.set_multilane_enabled(self.param_form.n_lanes.value() > 1)

    def _on_preset_changed(self, name: str):
        if not name:
            return
        self._suppress_custom_revert = True
        try:
            self.param_form.apply_preset(name)
        finally:
            self._suppress_custom_revert = False
        if name != CUSTOM_LABEL:
            self.status.emit(f"preset: {name}")

    def _on_form_edited(self):
        if self._suppress_custom_revert:
            return
        if self.preset_combo.currentText() != CUSTOM_LABEL:
            self.preset_combo.blockSignals(True)
            self.preset_combo.setCurrentText(CUSTOM_LABEL)
            self.preset_combo.blockSignals(False)

    def _on_run(self):
        if not self._binary.exists():
            self.log.emit(f"binary not found: {self._binary} — try `make pde`")
            return
        if self._runner is not None and self._runner.isRunning():
            return
        try:
            params = self.param_form.pde_params()
        except ValueError as exc:
            self.log.emit(f"param error: {exc}")
            return

        preset = self.preset_combo.currentText() or "custom"
        safe = preset.replace(" ", "_").replace(":", "").replace("/", "_").strip("()") or "custom"
        out_path = self._output_dir / f"{safe}.nc"

        self.run_button.setEnabled(False)
        self.progress.emit(0)
        self.status.emit(
            f"running PDE preset={preset} M={params['M']} n_steps={params['n_steps']}..."
        )
        self.log.emit(f"--- run PDE preset={preset} ---")

        self._runner = PDERunnerThread(params, out_path, self._binary, parent=self)
        self._runner.progress.connect(self._on_runner_progress)
        self._runner.finished_ok.connect(self._on_finished_ok)
        self._runner.finished_error.connect(self._on_finished_error)
        self._runner.start()

    def _on_runner_progress(self, frac: float, msg: str):
        if frac >= 0:
            self.progress.emit(int(round(100 * frac)))
        if msg:
            self.log.emit(msg)

    def _on_finished_ok(self, nc_path: str):
        self.progress.emit(100)
        self.status.emit(f"done — {nc_path}")
        self.run_button.setEnabled(True)
        self.log.emit(f"loaded {nc_path}")
        try:
            data = load_pde_netcdf(nc_path)
            self.plot_panel.show_result(data)
        except Exception as exc:
            self.log.emit(f"plot error: {exc}")

    def _on_finished_error(self, err: str):
        self.status.emit("error")
        self.run_button.setEnabled(True)
        self.log.emit(f"ERROR: {err}")

    def stop_running(self):
        if self._runner is not None and self._runner.isRunning():
            self._runner.requestInterruption()
            self._runner.wait(2000)


class MainWindow(QMainWindow):
    def __init__(self, binary: os.PathLike = DEFAULT_BINARY, output_dir: os.PathLike = DEFAULT_OUTDIR,
                 pde_binary: os.PathLike = DEFAULT_PDE_BINARY,
                 pde_output_dir: os.PathLike = DEFAULT_PDE_OUTDIR):
        super().__init__()
        self.setWindowTitle("PX915 Traffic Network Simulator")
        self.resize(1400, 850)

        # --- Top tabs: CA / PDE ---
        self.tabs = QTabWidget()
        self.tabs.setStyleSheet("""
            QTabBar::tab {
                min-width: 220px; min-height: 38px;
                font-size: 13px; font-weight: bold;
                padding: 4px 16px;
            }
            QTabBar::tab:selected { background: #cce4ff; }
            QTabBar::tab:!selected { background: #e8e8e8; }
        """)
        self.ca_tab = CATab(Path(binary), Path(output_dir))
        self.pde_tab = PDETab(Path(pde_binary), Path(pde_output_dir))
        self.tabs.addTab(self.ca_tab, "Cellular Automaton (CA)")
        self.tabs.addTab(self.pde_tab, "PDE Continuum Model")
        self.setCentralWidget(self.tabs)

        # --- Status bar (shared across tabs) ---
        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setMinimumWidth(220)
        self.status_label = QLabel("ready")
        status = QStatusBar()
        status.addPermanentWidget(self.status_label, 1)
        status.addPermanentWidget(self.progress)
        self.setStatusBar(status)

        # --- Log dock (shared across tabs) ---
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(500)
        self.log.setMaximumHeight(140)
        log_holder = QWidget()
        ll = QVBoxLayout(log_holder)
        ll.setContentsMargins(2, 2, 2, 2)
        ll.addWidget(self.log)
        dock = QDockWidget("Run log", self)
        dock.setWidget(log_holder)
        dock.setFeatures(QDockWidget.DockWidgetMovable | QDockWidget.DockWidgetFloatable)
        self.addDockWidget(Qt.BottomDockWidgetArea, dock)

        # Wire tab signals to the shared widgets.
        self.ca_tab.log.connect(self.log.appendPlainText)
        self.ca_tab.status.connect(self.status_label.setText)
        self.ca_tab.progress.connect(self.progress.setValue)
        self.pde_tab.log.connect(self.log.appendPlainText)
        self.pde_tab.status.connect(self.status_label.setText)
        self.pde_tab.progress.connect(self.progress.setValue)

    def closeEvent(self, event):
        # Ask each running tab to wind down before the window closes.
        self.ca_tab.stop_running()
        self.pde_tab.stop_running()
        super().closeEvent(event)
