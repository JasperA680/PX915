"""LWR continuum PDE tab: presets + parameter form + plot panel."""

from __future__ import annotations

from pathlib import Path

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QComboBox,
    QPushButton, QLabel, QScrollArea,
)

from python.pde_runner import load_pde_netcdf

from python.gui.pde_param_form import PDEParamForm, PRESETS_ORDER, CUSTOM_LABEL
from python.gui.pde_plot_panel import PDEPlotPanel
from python.gui.runners import PDERunnerThread
from python.gui._run_tab import RunTab
from python.gui._common import ALWAYS_VISIBLE_SCROLLBAR_QSS


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_PDE_BINARY = REPO_ROOT / "build" / "pde_solver"
DEFAULT_PDE_OUTDIR = REPO_ROOT / "data" / "output" / "gui_pde"


class PDETab(RunTab):
    """LWR continuum PDE tab: presets + parameter form + plot panel."""

    def __init__(self, binary: Path, output_dir: Path, parent=None):
        super().__init__(parent)
        self._binary = binary
        self._output_dir = output_dir
        self._suppress_custom_revert = False

        # --- Toolbar: preset + run ---
        toolbar = QWidget()
        toolbar_row = QHBoxLayout(toolbar)
        toolbar_row.addWidget(QLabel("Preset:"))
        self.preset_combo = QComboBox()
        for name in PRESETS_ORDER:
            self.preset_combo.addItem(name)
        self.preset_combo.setToolTip(
            "Pre-configured PDE scenario.  Picking one populates the\n"
            "parameter form below; editing any field afterwards flips the\n"
            "selector to '(custom)' so the run is clearly user-defined.\n"
            "Single-lane presets cover Riemann shocks and rarefactions;\n"
            "multilane presets demonstrate independent vs coupled lanes,\n"
            "mass conservation, and fast/slow overtaking lanes."
        )
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
        # Keep the vertical scrollbar visible at all times so users get a
        # static visual cue that the form has more content below the fold.
        # macOS auto-hides native scrollbars regardless of ScrollBarAlwaysOn;
        # applying a QScrollBar stylesheet kicks Qt off the native style so
        # AlwaysOn is actually honoured.
        form_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        form_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        form_scroll.setStyleSheet(ALWAYS_VISIBLE_SCROLLBAR_QSS)

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
        self.param_form.params_changed.connect(self._mark_stale)
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
        self._mark_stale()

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
        self._mark_running(f"{preset} / M={params['M']} / n_steps={params['n_steps']}")
        self._runner.start()

    def _on_finished_ok(self, nc_path: str):
        self.progress.emit(100)
        self.status.emit(f"done — {nc_path}")
        self.run_button.setEnabled(True)
        self.log.emit(f"loaded {nc_path}")
        try:
            data = load_pde_netcdf(nc_path)
            self.plot_panel.show_result(data)
            preset = self.preset_combo.currentText() or "custom"
            attrs = data.get("attrs", {})
            n_lanes = int(data.get("n_lanes", 1))
            summary = (
                f"{preset} / M={attrs.get('M','?')} / n_steps={attrs.get('n_steps','?')}"
                f" / lanes={n_lanes}"
            )
            self._mark_fresh(summary)
        except Exception as exc:
            self.log.emit(f"plot error: {exc}")
            self._mark_stale()
