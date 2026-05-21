"""Parameter form: global knobs for SimParams.

Per-road α, β, length are set via the click-to-edit dialog on the network
preview, and junction routing probabilities use sensible defaults — so the
boundary/routing group is no longer part of this form.
"""

from __future__ import annotations

from typing import Dict

from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import (
    QWidget, QFormLayout, QSpinBox, QDoubleSpinBox, QGroupBox, QVBoxLayout,
)

from python.road_network import SimParams

_BIG_INT = 2_000_000_000


class ParamForm(QWidget):
    """Editable SimParams (model-wide knobs)."""

    params_changed = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(parent)
        outer = QVBoxLayout(self)

        # SimParams group.
        sim_group = QGroupBox("Simulation")
        form = QFormLayout(sim_group)

        self.n_steps = QSpinBox()
        self.n_steps.setRange(10, _BIG_INT)
        self.n_steps.setValue(2000)
        self.n_steps.setToolTip(
            "n_steps: Number of CA time steps to simulate.\n"
            "More steps = longer run time but more time to reach steady state."
        )

        # Per-road length is set via the click-to-edit dialog on the
        # network preview; initial spec construction uses this default.
        self._default_lane_length = 20

        # Default α/β/routing probabilities used during initial spec
        # construction.  Override per-road via the click-to-edit dialog.
        self._default_alpha   = 0.4
        self._default_beta    = 0.5
        self._default_p_left  = 0.25
        self._default_p_right = 0.25
        self._default_p_exit  = 0.25

        self.rng_seed = QSpinBox()
        self.rng_seed.setRange(0, _BIG_INT)
        self.rng_seed.setValue(42)
        self.rng_seed.setToolTip(
            "rng_seed: Random number seed for reproducibility.\n"
            "Set to 0 for a random seed each run."
        )

        self.v_max = QSpinBox()
        self.v_max.setRange(1, _BIG_INT)
        self.v_max.setValue(5)
        self.v_max.setToolTip(
            "v_max: Maximum vehicle speed in cells per timestep (Nagel–Schreckenberg model).\n"
            "Vehicles accelerate toward this limit each step."
        )

        self.p_slow = QDoubleSpinBox()
        self.p_slow.setRange(0.0, 1.0)
        self.p_slow.setSingleStep(0.05)
        self.p_slow.setValue(0.2)
        self.p_slow.setToolTip(
            "p_slow: Dawdle (randomisation) probability.\n"
            "At each step a vehicle decelerates by 1 cell/step with this probability\n"
            "(Nagel–Schreckenberg noise term). 0 = no noise, 1 = maximum noise."
        )

        form.addRow("Steps", self.n_steps)
        form.addRow("RNG seed", self.rng_seed)
        form.addRow("Max speed (cells/step)", self.v_max)
        form.addRow("Dawdle probability", self.p_slow)
        outer.addWidget(sim_group)

        outer.addStretch(1)

        for w in (self.n_steps, self.rng_seed, self.v_max, self.p_slow):
            w.valueChanged.connect(self.params_changed.emit)

    def sim_params(self) -> SimParams:
        return SimParams(
            n_steps=self.n_steps.value(),
            v_max=self.v_max.value(),
            p_slow=float(self.p_slow.value()),
            rng_seed=self.rng_seed.value(),
        )

    def preset_kwargs(self, preset: str) -> Dict[str, float]:
        """Default kwargs for a preset constructor.

        Per-road α/β/length are overridden after construction via the
        network preview's click-to-edit dialog.
        """
        L = self._default_lane_length
        a = self._default_alpha
        b = self._default_beta
        if preset == "single_lane":
            return {"L": L, "alpha": a, "beta": b}
        if preset == "two_lane":
            return {"L": L, "alpha": a, "beta": b}
        if preset == "crossroads":
            return {"L": L, "alpha": a, "beta": b,
                    "p_left":  self._default_p_left,
                    "p_right": self._default_p_right}
        if preset == "roundabout":
            return {"L_arm": L, "L_ring": max(5, L - 3), "alpha": a, "beta": b,
                    "p_exit": self._default_p_exit}
        if preset == "town":
            return {"L": L, "alpha": a, "beta": b,
                    "p_left":  self._default_p_left,
                    "p_right": self._default_p_right}
        return {}
