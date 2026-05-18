"""Parameter form for the LWR PDE solver: grid, model, and lanes."""

from __future__ import annotations

from typing import Dict, List, Optional

from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import (
    QWidget, QFormLayout, QSpinBox, QDoubleSpinBox, QGroupBox, QVBoxLayout,
    QComboBox, QLineEdit,
)


# Single-lane Riemann presets (mirrors scripts/run_pde_model.py).
SINGLE_LANE_PRESETS: Dict[str, dict] = {
    "shock":             dict(rho_left_bc=0.1, rho_right_bc=0.7),
    "rarefaction":       dict(rho_left_bc=0.7, rho_right_bc=0.3),
    "shock-critical":    dict(rho_left_bc=0.2, rho_right_bc=0.8),
    "sonic-rarefaction": dict(rho_left_bc=0.8, rho_right_bc=0.2),
}

# Multi-lane scenarios (mirrors scripts/run_multilane_pde.py).
MULTILANE_PRESETS: Dict[str, dict] = {
    "multilane: independent": dict(
        n_lanes=2, lane_change_rate=0.0,
        v_max=1.0, rho_max=1.0,
        rho_left_bc=0.7, rho_right_bc=0.2,
        ic_type="staggered", flux_type="godunov", bc_type="periodic",
        M=200, n_steps=800,
    ),
    "multilane: coupled": dict(
        n_lanes=2, lane_change_rate=0.5,
        v_max=1.0, rho_max=1.0,
        rho_left_bc=0.7, rho_right_bc=0.2,
        ic_type="staggered", flux_type="godunov", bc_type="periodic",
        M=200, n_steps=800,
    ),
    "multilane: conservation": dict(
        n_lanes=2, lane_change_rate=0.5,
        v_max=1.0, rho_max=1.0,
        rho_left_bc=0.3, rho_right_bc=0.7,
        ic_type="sine", flux_type="godunov", bc_type="periodic",
        M=200, n_steps=800,
    ),
    "multilane: fast/slow": dict(
        n_lanes=2, lane_change_rate=0.5,
        v_max=1.0, rho_max=1.0,
        rho_left_bc=0.45, rho_right_bc=0.45,
        v_max_lanes=[1.0, 1.5],
        ic_type="constant", flux_type="godunov", bc_type="periodic",
        M=200, n_steps=1200,
    ),
}

CUSTOM_LABEL = "(custom)"
PRESETS_ORDER = [CUSTOM_LABEL] + list(SINGLE_LANE_PRESETS) + list(MULTILANE_PRESETS)


def _parse_csv_floats(text: str) -> Optional[List[float]]:
    s = text.strip()
    if not s:
        return None
    parts = [p.strip() for p in s.split(",") if p.strip()]
    try:
        return [float(p) for p in parts]
    except ValueError:
        raise ValueError(f"could not parse '{text}' as comma-separated floats")


def _csv_str(val) -> str:
    if val is None:
        return ""
    if isinstance(val, (list, tuple)):
        return ",".join(str(v) for v in val)
    return str(val)


class PDEParamForm(QWidget):
    """Editable PDE parameters that map directly to ``pde_runner.run_pde`` kwargs."""

    params_changed = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        outer = QVBoxLayout(self)

        # Grid / time
        grid_group = QGroupBox("Grid and time")
        gform = QFormLayout(grid_group)
        self.M = QSpinBox()
        self.M.setRange(20, 2000); self.M.setValue(200)
        self.n_steps = QSpinBox()
        self.n_steps.setRange(10, 10000); self.n_steps.setValue(500)
        gform.addRow("M (cells)", self.M)
        gform.addRow("n_steps", self.n_steps)
        outer.addWidget(grid_group)

        # Model
        model_group = QGroupBox("Model")
        mform = QFormLayout(model_group)
        self.v_max = QDoubleSpinBox()
        self.v_max.setRange(0.01, 10.0); self.v_max.setSingleStep(0.05); self.v_max.setValue(1.0); self.v_max.setDecimals(3)
        self.rho_max = QDoubleSpinBox()
        self.rho_max.setRange(0.01, 10.0); self.rho_max.setSingleStep(0.05); self.rho_max.setValue(1.0); self.rho_max.setDecimals(3)
        self.ic_type = QComboBox()
        self.ic_type.addItems(["riemann", "constant", "gaussian", "sine", "staggered"])
        self.flux_type = QComboBox()
        self.flux_type.addItems(["lf", "godunov", "newell"])
        self.bc_type = QComboBox()
        self.bc_type.addItems(["open", "periodic", "sponge"])
        self.rho_left_bc = QDoubleSpinBox()
        self.rho_left_bc.setRange(0.0, 10.0); self.rho_left_bc.setSingleStep(0.05); self.rho_left_bc.setValue(0.1); self.rho_left_bc.setDecimals(3)
        self.rho_right_bc = QDoubleSpinBox()
        self.rho_right_bc.setRange(0.0, 10.0); self.rho_right_bc.setSingleStep(0.05); self.rho_right_bc.setValue(0.9); self.rho_right_bc.setDecimals(3)
        mform.addRow("v_max", self.v_max)
        mform.addRow("rho_max", self.rho_max)
        mform.addRow("ic_type", self.ic_type)
        mform.addRow("flux_type", self.flux_type)
        mform.addRow("bc_type", self.bc_type)
        mform.addRow("rho_left_bc", self.rho_left_bc)
        mform.addRow("rho_right_bc", self.rho_right_bc)
        outer.addWidget(model_group)

        # Lanes
        lanes_group = QGroupBox("Lanes")
        lform = QFormLayout(lanes_group)
        self.n_lanes = QSpinBox()
        self.n_lanes.setRange(1, 4); self.n_lanes.setValue(1)
        self.lane_change_rate = QDoubleSpinBox()
        self.lane_change_rate.setRange(0.0, 5.0); self.lane_change_rate.setSingleStep(0.1)
        self.lane_change_rate.setValue(0.0); self.lane_change_rate.setDecimals(3)
        self.v_max_lanes = QLineEdit()
        self.v_max_lanes.setPlaceholderText("blank = broadcast v_max; e.g. 1.0,1.5")
        self.rho_max_lanes = QLineEdit()
        self.rho_max_lanes.setPlaceholderText("blank = broadcast rho_max")
        lform.addRow("n_lanes", self.n_lanes)
        lform.addRow("lane_change_rate", self.lane_change_rate)
        lform.addRow("v_max_lanes", self.v_max_lanes)
        lform.addRow("rho_max_lanes", self.rho_max_lanes)
        outer.addWidget(lanes_group)
        outer.addStretch(1)

        self._update_lane_enabled()

        for w in (self.M, self.n_steps, self.v_max, self.rho_max,
                  self.rho_left_bc, self.rho_right_bc,
                  self.n_lanes, self.lane_change_rate):
            w.valueChanged.connect(self.params_changed.emit)
        for w in (self.ic_type, self.flux_type, self.bc_type):
            w.currentTextChanged.connect(lambda _t: self.params_changed.emit())
        for w in (self.v_max_lanes, self.rho_max_lanes):
            w.textChanged.connect(lambda _t: self.params_changed.emit())
        self.n_lanes.valueChanged.connect(lambda _v: self._update_lane_enabled())

    def _update_lane_enabled(self):
        enabled = self.n_lanes.value() > 1
        self.lane_change_rate.setEnabled(enabled)
        self.v_max_lanes.setEnabled(enabled)
        self.rho_max_lanes.setEnabled(enabled)

    def apply_preset(self, name: str) -> None:
        """Programmatically set every widget for a named preset.

        ``(custom)`` is a no-op. Signals are blocked during the bulk set and
        ``params_changed`` is emitted exactly once at the end.
        """
        if name == CUSTOM_LABEL:
            return
        preset = SINGLE_LANE_PRESETS.get(name) or MULTILANE_PRESETS.get(name)
        if preset is None:
            return

        defaults = dict(
            M=200, n_steps=500,
            v_max=1.0, rho_max=1.0,
            rho_left_bc=0.1, rho_right_bc=0.9,
            ic_type="riemann", flux_type="lf", bc_type="open",
            n_lanes=1, lane_change_rate=0.0,
            v_max_lanes=None, rho_max_lanes=None,
        )
        merged = {**defaults, **preset}

        widgets = [self.M, self.n_steps, self.v_max, self.rho_max,
                   self.rho_left_bc, self.rho_right_bc,
                   self.n_lanes, self.lane_change_rate,
                   self.ic_type, self.flux_type, self.bc_type,
                   self.v_max_lanes, self.rho_max_lanes]
        for w in widgets:
            w.blockSignals(True)
        try:
            self.M.setValue(int(merged["M"]))
            self.n_steps.setValue(int(merged["n_steps"]))
            self.v_max.setValue(float(merged["v_max"]))
            self.rho_max.setValue(float(merged["rho_max"]))
            self.rho_left_bc.setValue(float(merged["rho_left_bc"]))
            self.rho_right_bc.setValue(float(merged["rho_right_bc"]))
            self.ic_type.setCurrentText(str(merged["ic_type"]))
            self.flux_type.setCurrentText(str(merged["flux_type"]))
            self.bc_type.setCurrentText(str(merged["bc_type"]))
            self.n_lanes.setValue(int(merged["n_lanes"]))
            self.lane_change_rate.setValue(float(merged["lane_change_rate"]))
            self.v_max_lanes.setText(_csv_str(merged["v_max_lanes"]))
            self.rho_max_lanes.setText(_csv_str(merged["rho_max_lanes"]))
        finally:
            for w in widgets:
                w.blockSignals(False)
        self._update_lane_enabled()
        self.params_changed.emit()

    def pde_params(self) -> dict:
        n_lanes = int(self.n_lanes.value())
        v_lanes = _parse_csv_floats(self.v_max_lanes.text())
        rho_lanes = _parse_csv_floats(self.rho_max_lanes.text())
        if v_lanes is not None and len(v_lanes) != n_lanes:
            raise ValueError(
                f"v_max_lanes has {len(v_lanes)} value(s) but n_lanes={n_lanes}"
            )
        if rho_lanes is not None and len(rho_lanes) != n_lanes:
            raise ValueError(
                f"rho_max_lanes has {len(rho_lanes)} value(s) but n_lanes={n_lanes}"
            )
        return dict(
            M=int(self.M.value()),
            n_steps=int(self.n_steps.value()),
            v_max=float(self.v_max.value()),
            rho_max=float(self.rho_max.value()),
            rho_left_bc=float(self.rho_left_bc.value()),
            rho_right_bc=float(self.rho_right_bc.value()),
            ic_type=self.ic_type.currentText(),
            flux_type=self.flux_type.currentText(),
            bc_type=self.bc_type.currentText(),
            n_lanes=n_lanes,
            lane_change_rate=float(self.lane_change_rate.value()),
            v_max_lanes=v_lanes,
            rho_max_lanes=rho_lanes,
        )
