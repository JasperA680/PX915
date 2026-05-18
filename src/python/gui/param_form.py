"""Parameter form: global knobs for SimParams + preset constructor kwargs."""

from __future__ import annotations

from typing import Dict

from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import (
    QWidget, QFormLayout, QSpinBox, QDoubleSpinBox, QGroupBox, QVBoxLayout,
)

from python.road_network import SimParams


class ParamForm(QWidget):
    """Editable parameters: SimParams + preset-specific kwargs (alpha, beta, p_left ...)."""

    params_changed = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        outer = QVBoxLayout(self)

        # SimParams group.
        sim_group = QGroupBox("Simulation")
        form = QFormLayout(sim_group)
        self.n_steps = QSpinBox()
        self.n_steps.setRange(10, 100000)
        self.n_steps.setValue(2000)
        self.lane_length = QSpinBox()
        self.lane_length.setRange(5, 1000)
        self.lane_length.setValue(20)
        self.rng_seed = QSpinBox()
        self.rng_seed.setRange(0, 2_000_000_000)
        self.rng_seed.setValue(42)
        self.v_max = QSpinBox()
        self.v_max.setRange(1, 20)
        self.v_max.setValue(5)
        self.p_slow = QDoubleSpinBox()
        self.p_slow.setRange(0.0, 1.0)
        self.p_slow.setSingleStep(0.05)
        self.p_slow.setValue(0.2)
        form.addRow("n_steps", self.n_steps)
        form.addRow("lane_length", self.lane_length)
        form.addRow("rng_seed", self.rng_seed)
        form.addRow("v_max", self.v_max)
        form.addRow("p_slow", self.p_slow)
        outer.addWidget(sim_group)

        # Boundary / routing group.
        boundary = QGroupBox("Boundaries and routing")
        bform = QFormLayout(boundary)
        self.alpha = QDoubleSpinBox()
        self.alpha.setRange(0.0, 1.0); self.alpha.setSingleStep(0.05); self.alpha.setValue(0.4)
        self.beta = QDoubleSpinBox()
        self.beta.setRange(0.0, 1.0); self.beta.setSingleStep(0.05); self.beta.setValue(0.5)
        self.p_left = QDoubleSpinBox()
        self.p_left.setRange(0.0, 1.0); self.p_left.setSingleStep(0.05); self.p_left.setValue(0.25)
        self.p_right = QDoubleSpinBox()
        self.p_right.setRange(0.0, 1.0); self.p_right.setSingleStep(0.05); self.p_right.setValue(0.25)
        self.p_exit = QDoubleSpinBox()
        self.p_exit.setRange(0.0, 1.0); self.p_exit.setSingleStep(0.05); self.p_exit.setValue(0.25)
        bform.addRow("alpha (entry rate)", self.alpha)
        bform.addRow("beta (exit rate)", self.beta)
        bform.addRow("p_left (crossroad)", self.p_left)
        bform.addRow("p_right (crossroad)", self.p_right)
        bform.addRow("p_exit (roundabout)", self.p_exit)
        outer.addWidget(boundary)

        outer.addStretch(1)

        for w in (self.n_steps, self.lane_length, self.rng_seed, self.v_max, self.p_slow,
                  self.alpha, self.beta, self.p_left, self.p_right, self.p_exit):
            w.valueChanged.connect(self.params_changed.emit)

    def sim_params(self) -> SimParams:
        return SimParams(
            n_steps=self.n_steps.value(),
            v_max=self.v_max.value(),
            p_slow=float(self.p_slow.value()),
            rng_seed=self.rng_seed.value(),
        )

    def preset_kwargs(self, preset: str) -> Dict[str, float]:
        """Kwargs to pass to the preset constructor (excluding lane_length L)."""
        L = self.lane_length.value()
        a = float(self.alpha.value())
        b = float(self.beta.value())
        if preset == "single_lane":
            return {"L": L, "alpha": a, "beta": b}
        if preset == "crossroads":
            return {"L": L, "alpha": a, "beta": b,
                    "p_left":  float(self.p_left.value()),
                    "p_right": float(self.p_right.value())}
        if preset == "roundabout":
            return {"L_arm": L, "L_ring": max(5, L - 3), "alpha": a, "beta": b,
                    "p_exit": float(self.p_exit.value())}
        if preset == "town":
            return {"L": L, "alpha": a, "beta": b,
                    "p_left":  float(self.p_left.value()),
                    "p_right": float(self.p_right.value())}
        return {}
