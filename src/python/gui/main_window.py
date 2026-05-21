"""Main window: top-level CA / PDE tabs, status bar, run log."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

from PyQt5.QtCore import Qt, QThread, pyqtSignal
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
from python.gui.road_edit_dialog import RoadEditDialog, JunctionEditDialog


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_BINARY = REPO_ROOT / "build" / "run_network"
DEFAULT_OUTDIR = REPO_ROOT / "data" / "output" / "gui"
DEFAULT_PDE_BINARY = REPO_ROOT / "build" / "pde_solver"
DEFAULT_PDE_OUTDIR = REPO_ROOT / "data" / "output" / "gui_pde"

# Available CA update rules.  Only "NS" is implemented physics-side; "TASEP"
# is plumbed through the JSON config but the Fortran driver errors out on it.
CA_MODELS = ("NS", "TASEP")

# Presets where TASEP physics is appropriate (single 1D chain only).
TASEP_OK_PRESETS = {"single_lane"}


class FDSweepThread(QThread):
    """Run the pure-Python TASEP fundamental_diagram() sweep off the UI thread."""

    finished_ok = pyqtSignal(object, object, int, int)   # rho, J, L, n_points
    finished_error = pyqtSignal(str)

    def __init__(self, L: int = 50, n_points: int = 30, parent=None):
        super().__init__(parent)
        self.L = L
        self.n_points = n_points

    def run(self):
        try:
            from python.analysis import fundamental_diagram
            rho, J = fundamental_diagram(L=self.L, n_points=self.n_points)
            self.finished_ok.emit(rho, J, self.L, self.n_points)
        except Exception as exc:
            self.finished_error.emit(str(exc))


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
        self._fd_runner: Optional[FDSweepThread] = None
        self._current_spec: Optional[NetworkSpec] = None
        self._current_layout: Optional[LayoutSpec] = None
        # Per-road parameter overrides set via the click-to-edit dialog.
        # Keyed by road_id; reset when the preset combo changes.
        self._road_overrides: dict = {}
        # Per-junction routing overrides (dict[junction_id -> list[list[float]]]).
        self._junction_overrides: dict = {}

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
        toolbar_row.addSpacing(12)
        self.fd_button = QPushButton("Run FD sweep")
        self.fd_button.setMinimumWidth(140)
        self.fd_button.setToolTip(
            "Run a TASEP fundamental-diagram sweep (α with β=1, β with α=1) "
            "using pure-Python analysis.fundamental_diagram(). Populates the "
            "Fundamental Diagram plot tab when finished."
        )
        toolbar_row.addWidget(self.fd_button)
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
        self.fd_button.clicked.connect(self._on_run_fd)
        self.network_widget.road_clicked.connect(self._on_road_clicked)
        self.network_widget.junction_clicked.connect(self._on_junction_clicked)

        # Bootstrap with the first preset (also greys out TASEP if needed).
        self._sync_model_availability()
        self._rebuild_spec()

    def _on_preset_changed_ca(self, name: str):
        """Reset per-road and per-junction overrides on preset change, then rebuild."""
        self._road_overrides = {}
        self._junction_overrides = {}
        self._sync_model_availability()
        self._rebuild_spec()

    def _sync_model_availability(self):
        """Grey out TASEP for any preset where it is not physically meaningful.

        TASEP is a strict 1D single-lane chain — anything with multiple lanes
        or junctions falls outside its definition.
        """
        name = self.preset_combo.currentText()
        tasep_ok = name in TASEP_OK_PRESETS
        # Find the TASEP entry in the combo and toggle its enabled flag.
        idx = self.model_combo.findText("TASEP")
        if idx >= 0:
            model = self.model_combo.model()
            item = model.item(idx)
            if item is not None:
                # Qt.ItemIsEnabled = 32; remove or add it to control greying.
                if tasep_ok:
                    item.setFlags(item.flags() | Qt.ItemIsEnabled)
                else:
                    item.setFlags(item.flags() & ~Qt.ItemIsEnabled)
        # If TASEP is now disabled but currently selected, switch to NS.
        if not tasep_ok and self.model_combo.currentText() == "TASEP":
            ns_idx = self.model_combo.findText("NS")
            if ns_idx >= 0:
                self.model_combo.setCurrentIndex(ns_idx)

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
        # Apply any per-junction routing overrides.
        for jid, routes in self._junction_overrides.items():
            j = next((j for j in spec.junctions if j.id == jid), None)
            if j is None:
                continue
            j.routes = [list(row) for row in routes]
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

    def _on_junction_clicked(self, junction_id: int):
        """Open the junction routing editor and apply any changes."""
        if self._current_spec is None:
            return
        j = next((j for j in self._current_spec.junctions if j.id == junction_id), None)
        if j is None or not j.routes:
            self.log.emit(f"J{junction_id}: no routing matrix to edit")
            return

        dlg = JunctionEditDialog(junction_id, j.routes, parent=self)
        if dlg.exec_() == JunctionEditDialog.Accepted:
            new_routes = dlg.routes()
            self._junction_overrides[junction_id] = new_routes
            self._rebuild_spec()
            summary = "  ".join(
                "[" + ", ".join(f"{v:.2f}" for v in row) + "]"
                for row in new_routes
            )
            self.log.emit(f"J{junction_id} routes overridden: {summary}")

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

    # ------------------------------------------------------------------
    # Fundamental-diagram sweep (pure-Python TASEP)
    # ------------------------------------------------------------------

    def _on_run_fd(self):
        if self._fd_runner is not None and self._fd_runner.isRunning():
            return
        L = 50
        n_points = 30
        self.fd_button.setEnabled(False)
        self.status.emit(f"running FD sweep L={L} ({2 * n_points} points)…")
        self.log.emit(f"--- FD sweep L={L} points={n_points} ---")
        self._fd_runner = FDSweepThread(L=L, n_points=n_points, parent=self)
        self._fd_runner.finished_ok.connect(self._on_fd_ok)
        self._fd_runner.finished_error.connect(self._on_fd_error)
        self._fd_runner.start()

    def _on_fd_ok(self, rho, J, L, n_points):
        self.fd_button.setEnabled(True)
        self.status.emit(f"FD sweep done — {len(rho)} points")
        self.log.emit(
            f"FD sweep complete: L={L}, {len(rho)} points  "
            f"(ρ range {float(rho.min()):.2f}–{float(rho.max()):.2f}, "
            f"J range {float(J.min()):.3f}–{float(J.max()):.3f})"
        )
        self.plot_panel.set_fd_sweep(rho, J, title=f"Fundamental diagram  (L={L})")

    def _on_fd_error(self, err: str):
        self.fd_button.setEnabled(True)
        self.status.emit("FD sweep error")
        self.log.emit(f"FD ERROR: {err}")

    # ------------------------------------------------------------------

    def _on_finished_ok(self, nc_path: str):
        self.progress.emit(100)
        self.status.emit(f"done — {nc_path}")
        self.run_button.setEnabled(True)
        self.log.emit(f"loaded {nc_path}")
        try:
            result = load_network_netcdf(nc_path)
            self.plot_panel.show_result(result)
            self._log_diagnostics(result)
        except Exception as exc:
            self.log.emit(f"plot error: {exc}")

    def _log_diagnostics(self, result):
        """Print per-road steady-state ρ̄ and J̄ to the run log."""
        import numpy as np
        n_steps, n_roads = result.road_density.shape
        # Average over the second half of the run to skip the transient.
        half = max(1, n_steps // 2)
        lane_road = np.asarray(result.lane_road_id)

        self.log.emit("--- diagnostics (steady-state averages, last half of run) ---")
        for r in range(n_roads):
            rho_bar = float(result.road_density[half:, r].mean())
            n_lanes = max(1, int(np.sum(lane_road == r + 1)))
            J_bar = float(result.road_exits[half:, r].mean()) / n_lanes
            total_in  = int(result.road_entries[:, r].sum())
            total_out = int(result.road_exits[:, r].sum())
            self.log.emit(
                f"  R{r + 1}: <ρ>={rho_bar:.3f}  <J>={J_bar:.3f}  "
                f"(in={total_in}, out={total_out}, lanes={n_lanes})"
            )

    def _on_finished_error(self, err: str):
        self.status.emit("error")
        self.run_button.setEnabled(True)
        self.log.emit(f"ERROR: {err}")

    def stop_running(self):
        if self._runner is not None and self._runner.isRunning():
            self._runner.requestInterruption()
            self._runner.wait(2000)
        if self._fd_runner is not None and self._fd_runner.isRunning():
            self._fd_runner.requestInterruption()
            self._fd_runner.wait(2000)


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
        # Scope styling to *this* widget only (use objectName selector) so it
        # doesn't cascade into the inner plot-panel tab widgets.
        self.tabs.setObjectName("mainTabs")
        self.tabs.setStyleSheet("""
            QTabWidget#mainTabs > QTabBar::tab {
                min-width: 220px; min-height: 30px;
                font-size: 13px; font-weight: bold;
                padding: 4px 16px;
            }
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
