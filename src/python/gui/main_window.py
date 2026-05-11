"""Main window: preset selector, parameter form, network preview, run button, plot tabs."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QComboBox,
    QPushButton, QProgressBar, QLabel, QStatusBar, QPlainTextEdit,
)

from python.road_network import PRESETS, NetworkSpec, LayoutSpec
from python.io import load_network_netcdf

from python.gui.param_form import ParamForm
from python.gui.network_widget import NetworkWidget
from python.gui.runner_thread import RunnerThread
from python.gui.plot_panel import PlotPanel


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_BINARY = REPO_ROOT / "build" / "run_network"
DEFAULT_OUTDIR = REPO_ROOT / "data" / "output" / "gui"


class MainWindow(QMainWindow):
    def __init__(self, binary: os.PathLike = DEFAULT_BINARY, output_dir: os.PathLike = DEFAULT_OUTDIR):
        super().__init__()
        self.setWindowTitle("PX915 Traffic Network Simulator")
        self.resize(1300, 800)
        self._binary = Path(binary)
        self._output_dir = Path(output_dir)
        self._runner: Optional[RunnerThread] = None
        self._current_spec: Optional[NetworkSpec] = None
        self._current_layout: Optional[LayoutSpec] = None

        # --- Toolbar row: preset + run ---
        toolbar = QWidget()
        toolbar_row = QHBoxLayout(toolbar)
        toolbar_row.addWidget(QLabel("Preset:"))
        self.preset_combo = QComboBox()
        for name in PRESETS:
            self.preset_combo.addItem(name)
        toolbar_row.addWidget(self.preset_combo)
        self.run_button = QPushButton("Run")
        self.run_button.setMinimumWidth(120)
        toolbar_row.addWidget(self.run_button)
        toolbar_row.addStretch(1)

        # --- Left: params + network preview ---
        self.param_form = ParamForm()
        self.network_widget = NetworkWidget()
        left = QSplitter(Qt.Vertical)
        left.addWidget(self.param_form)
        left.addWidget(self.network_widget)
        left.setStretchFactor(0, 0)
        left.setStretchFactor(1, 1)

        # --- Right: plot panel ---
        self.plot_panel = PlotPanel()

        # --- Central layout ---
        body = QSplitter(Qt.Horizontal)
        body.addWidget(left)
        body.addWidget(self.plot_panel)
        body.setStretchFactor(0, 0)
        body.setStretchFactor(1, 1)

        central = QWidget()
        cl = QVBoxLayout(central)
        cl.addWidget(toolbar)
        cl.addWidget(body, stretch=1)
        self.setCentralWidget(central)

        # --- Status bar + progress + log ---
        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setMinimumWidth(220)
        self.status_label = QLabel("ready")
        status = QStatusBar()
        status.addPermanentWidget(self.status_label, 1)
        status.addPermanentWidget(self.progress)
        self.setStatusBar(status)

        # --- Log dock at bottom ---
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(500)
        log_holder = QWidget()
        ll = QVBoxLayout(log_holder)
        ll.setContentsMargins(2, 2, 2, 2)
        ll.addWidget(QLabel("Run log"))
        ll.addWidget(self.log)
        from PyQt5.QtWidgets import QDockWidget
        dock = QDockWidget("Log", self)
        dock.setWidget(log_holder)
        self.addDockWidget(Qt.BottomDockWidgetArea, dock)

        # --- Connect signals ---
        self.preset_combo.currentTextChanged.connect(self._on_preset_changed)
        self.param_form.params_changed.connect(self._rebuild_spec)
        self.run_button.clicked.connect(self._on_run)

        # Bootstrap with the first preset.
        self._on_preset_changed(self.preset_combo.currentText())

    def _on_preset_changed(self, name: str):
        self._rebuild_spec()

    def _rebuild_spec(self):
        name = self.preset_combo.currentText()
        if name not in PRESETS:
            return
        kwargs = self.param_form.preset_kwargs(name)
        try:
            spec, layout = PRESETS[name](**kwargs)
        except TypeError as exc:
            self.log.appendPlainText(f"preset error: {exc}")
            return
        self._current_spec = spec
        self._current_layout = layout
        self.network_widget.set_network(spec, layout)
        self.status_label.setText(
            f"preset {name}: {len(spec.roads)} roads, {len(spec.junctions)} junctions"
        )

    def _on_run(self):
        if self._current_spec is None or self._current_layout is None:
            self.log.appendPlainText("no spec to run")
            return
        if not self._binary.exists():
            self.log.appendPlainText(
                f"binary not found: {self._binary} — try `make run_network`"
            )
            return
        if self._runner is not None and self._runner.isRunning():
            return

        params = self.param_form.sim_params()
        preset_name = self.preset_combo.currentText()
        out_dir = self._output_dir / preset_name
        self.run_button.setEnabled(False)
        self.progress.setValue(0)
        self.status_label.setText(f"running {preset_name} ({params.n_steps} steps)...")
        self.log.appendPlainText(f"--- run {preset_name} ---")

        self._runner = RunnerThread(
            self._current_spec, params, self._current_layout,
            output_dir=out_dir, binary=self._binary, parent=self,
        )
        self._runner.progress.connect(self._on_progress)
        self._runner.finished_ok.connect(self._on_finished_ok)
        self._runner.finished_error.connect(self._on_finished_error)
        self._runner.start()

    def _on_progress(self, frac: float, msg: str):
        if frac >= 0:
            self.progress.setValue(int(round(100 * frac)))
        if msg:
            self.log.appendPlainText(msg)

    def _on_finished_ok(self, nc_path: str):
        self.progress.setValue(100)
        self.status_label.setText(f"done — {nc_path}")
        self.run_button.setEnabled(True)
        self.log.appendPlainText(f"loaded {nc_path}")
        try:
            result = load_network_netcdf(nc_path)
            self.plot_panel.show_result(result)
        except Exception as exc:
            self.log.appendPlainText(f"plot error: {exc}")

    def _on_finished_error(self, err: str):
        self.status_label.setText("error")
        self.run_button.setEnabled(True)
        self.log.appendPlainText(f"ERROR: {err}")

    def closeEvent(self, event):
        if self._runner is not None and self._runner.isRunning():
            self._runner.requestInterruption()
            self._runner.wait(2000)
        super().closeEvent(event)
