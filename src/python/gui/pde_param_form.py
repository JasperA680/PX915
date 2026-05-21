"""Parameter form for the LWR PDE solver: grid, model, and lanes."""

from __future__ import annotations

import sys
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

_BIG_INT = 2_000_000_000
_BIG_FLOAT = 1e12


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
        self.M.setRange(2, _BIG_INT)
        self.M.setValue(200)
        self.M.setToolTip(
            "M: Number of spatial grid cells.\n"
            "More cells = finer resolution but slower simulation.\n"
            "Also determines the grid spacing dx = 1/M and the CFL timestep dt."
        )

        self.n_steps = QSpinBox()
        self.n_steps.setRange(1, _BIG_INT)
        self.n_steps.setValue(500)
        self.n_steps.setToolTip(
            "n_steps: Number of time steps to simulate.\n"
            "Changing this updates the Total time T display automatically."
        )

        # t_max is derived from n_steps and the CFL dt estimate; both are editable
        self.t_max = QDoubleSpinBox()
        self.t_max.setRange(0.001, _BIG_FLOAT)
        self.t_max.setSingleStep(0.5)
        self.t_max.setDecimals(4)
        self.t_max.setToolTip(
            "T: Total physical simulation time (seconds).\n"
            "Changing this updates n_steps automatically using dt ≈ 0.9 / (M × v_max).\n"
            "Increase T if the boundary flow has not converged."
        )

        gform.addRow("Spatial cells M", self.M)
        gform.addRow("Steps n_steps", self.n_steps)
        gform.addRow("Total time T", self.t_max)
        outer.addWidget(grid_group)

        # Model
        model_group = QGroupBox("Model")
        mform = QFormLayout(model_group)

        self.v_max = QDoubleSpinBox()
        self.v_max.setRange(0.001, _BIG_FLOAT)
        self.v_max.setSingleStep(0.05)
        self.v_max.setValue(1.0)
        self.v_max.setDecimals(3)
        self.v_max.setToolTip(
            "v_max: Free-flow speed — maximum vehicle speed at zero density\n"
            "(Greenshields model: v(ρ) = v_max × (1 − ρ/ρ_max)).\n"
            "Also affects the CFL timestep and total time T."
        )

        self.rho_max = QDoubleSpinBox()
        self.rho_max.setRange(0.001, _BIG_FLOAT)
        self.rho_max.setSingleStep(0.05)
        self.rho_max.setValue(1.0)
        self.rho_max.setDecimals(3)
        self.rho_max.setToolTip(
            "rho_max (ρ_max): Jam density — maximum number of vehicles per unit length.\n"
            "Traffic is fully stopped (v=0) when ρ = ρ_max."
        )

        self.ic_type = QComboBox()
        self.ic_type.addItems(["riemann", "constant", "gaussian", "sine", "staggered"])
        self.ic_type.setToolTip(
            "ic_type: Shape of the initial density profile ρ(x, t=0).\n"
            "  riemann   — step function (rho_left on left, rho_right on right)\n"
            "  constant  — uniform density everywhere\n"
            "  gaussian  — Gaussian bump in the centre\n"
            "  sine      — sinusoidal wave\n"
            "  staggered — alternating values per lane (multilane only)"
        )

        self.flux_type = QComboBox()
        self.flux_type.addItems(["lf", "godunov", "newell"])
        self.flux_type.setToolTip(
            "flux_type: Numerical flux scheme used to solve the LWR PDE.\n"
            "  lf      — Lax–Friedrichs (more diffusive, more stable)\n"
            "  godunov — Godunov (exact Riemann solver, sharper shocks)\n"
            "  newell  — Newell–Daganzo triangular fundamental diagram"
        )

        self.bc_type = QComboBox()
        self.bc_type.addItems(["open", "periodic", "sponge"])
        self.bc_type.setToolTip(
            "bc_type: How the domain boundaries are handled.\n"
            "  open     — fixed inflow/outflow densities (rho_left_bc, rho_right_bc)\n"
            "  periodic — right boundary wraps to left (no net flow)\n"
            "  sponge   — absorbing boundary (density relaxes to bc value)"
        )

        self.rho_left_bc = QDoubleSpinBox()
        self.rho_left_bc.setRange(0.0, _BIG_FLOAT)
        self.rho_left_bc.setSingleStep(0.05)
        self.rho_left_bc.setValue(0.1)
        self.rho_left_bc.setDecimals(3)
        self.rho_left_bc.setToolTip(
            "rho_left_bc: Density imposed at the left boundary (open/sponge BC).\n"
            "Low value = sparse inflow (free-flow); high value = congested inflow."
        )

        self.rho_right_bc = QDoubleSpinBox()
        self.rho_right_bc.setRange(0.0, _BIG_FLOAT)
        self.rho_right_bc.setSingleStep(0.05)
        self.rho_right_bc.setValue(0.9)
        self.rho_right_bc.setDecimals(3)
        self.rho_right_bc.setToolTip(
            "rho_right_bc: Density imposed at the right boundary (open/sponge BC).\n"
            "High value = downstream congestion (creates a shock wave travelling left)."
        )

        mform.addRow("Free-flow speed v_max", self.v_max)
        mform.addRow("Max density ρ_max", self.rho_max)
        mform.addRow("Initial condition", self.ic_type)
        mform.addRow("Flux scheme", self.flux_type)
        mform.addRow("Boundary condition", self.bc_type)
        mform.addRow("Left BC density", self.rho_left_bc)
        mform.addRow("Right BC density", self.rho_right_bc)
        outer.addWidget(model_group)

        # Lanes
        lanes_group = QGroupBox("Lanes")
        lform = QFormLayout(lanes_group)

        self.n_lanes = QSpinBox()
        self.n_lanes.setRange(1, _BIG_INT)
        self.n_lanes.setValue(1)
        self.n_lanes.setToolTip(
            "n_lanes: Number of parallel traffic lanes.\n"
            "Set > 1 to enable the multilane LWR model with lateral density exchange."
        )

        self.lane_change_rate = QDoubleSpinBox()
        self.lane_change_rate.setRange(0.0, _BIG_FLOAT)
        self.lane_change_rate.setSingleStep(0.1)
        self.lane_change_rate.setValue(0.0)
        self.lane_change_rate.setDecimals(3)
        self.lane_change_rate.setToolTip(
            "lane_change_rate: Rate coefficient for lateral density exchange between adjacent lanes.\n"
            "0 = lanes evolve independently; higher values couple the lanes more strongly."
        )

        self.v_max_lanes = QLineEdit()
        self.v_max_lanes.setPlaceholderText("blank = broadcast v_max; e.g. 1.0,1.5")
        self.v_max_lanes.setToolTip(
            "v_max_lanes: Comma-separated free-flow speed for each lane.\n"
            "Must contain exactly n_lanes values, or leave blank to use v_max for all lanes.\n"
            "Example (2 lanes, fast overtaking lane): 1.0,1.5"
        )

        self.rho_max_lanes = QLineEdit()
        self.rho_max_lanes.setPlaceholderText("blank = broadcast rho_max")
        self.rho_max_lanes.setToolTip(
            "rho_max_lanes: Comma-separated jam density for each lane.\n"
            "Must contain exactly n_lanes values, or leave blank to use rho_max for all lanes."
        )

        lform.addRow("Number of lanes", self.n_lanes)
        lform.addRow("Lane-change rate", self.lane_change_rate)
        lform.addRow("Per-lane v_max (CSV)", self.v_max_lanes)
        lform.addRow("Per-lane ρ_max (CSV)", self.rho_max_lanes)
        outer.addWidget(lanes_group)
        outer.addStretch(1)

        self._update_lane_enabled()

        # Sync t_max display from initial n_steps
        self._sync_t_from_steps()

        # Connect signals — n_steps/t_max are cross-linked; M and v_max also
        # update t_max (since dt depends on them).
        self.n_steps.valueChanged.connect(self._on_n_steps_changed)
        self.t_max.valueChanged.connect(self._on_t_max_changed)
        self.M.valueChanged.connect(self._on_grid_changed)
        self.v_max.valueChanged.connect(self._on_grid_changed)

        for w in (self.rho_max, self.rho_left_bc, self.rho_right_bc,
                  self.n_lanes, self.lane_change_rate):
            w.valueChanged.connect(self.params_changed.emit)
        for w in (self.ic_type, self.flux_type, self.bc_type):
            w.currentTextChanged.connect(lambda _t: self.params_changed.emit())
        for w in (self.v_max_lanes, self.rho_max_lanes):
            w.textChanged.connect(lambda _t: self.params_changed.emit())
        self.n_lanes.valueChanged.connect(lambda _v: self._update_lane_enabled())

    # ------------------------------------------------------------------
    # t_max / n_steps cross-link helpers
    # ------------------------------------------------------------------

    def _dt_estimate(self) -> float:
        """CFL-based dt estimate: 0.9 / (M * v_max)."""
        return 0.9 / max(1, self.M.value()) / max(1e-9, self.v_max.value())

    def _sync_t_from_steps(self):
        t = self.n_steps.value() * self._dt_estimate()
        self.t_max.blockSignals(True)
        self.t_max.setValue(round(t, 4))
        self.t_max.blockSignals(False)

    def _sync_steps_from_t(self):
        n = max(1, round(self.t_max.value() / max(1e-12, self._dt_estimate())))
        self.n_steps.blockSignals(True)
        self.n_steps.setValue(n)
        self.n_steps.blockSignals(False)

    def _on_n_steps_changed(self):
        self._sync_t_from_steps()
        self.params_changed.emit()

    def _on_t_max_changed(self):
        self._sync_steps_from_t()
        self.params_changed.emit()

    def _on_grid_changed(self):
        # M or v_max changed: dt shifts, so update the t_max display while
        # keeping n_steps fixed (the user set an explicit step count).
        self._sync_t_from_steps()
        self.params_changed.emit()

    # ------------------------------------------------------------------

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

        widgets = [self.M, self.n_steps, self.t_max,
                   self.v_max, self.rho_max,
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
        # Update t_max display from the newly-set n_steps (signals were blocked).
        self._sync_t_from_steps()
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
