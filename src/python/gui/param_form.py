"""Parameter form: global knobs for SimParams + preset constructor kwargs."""

from __future__ import annotations

from typing import Dict

from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import (
    QWidget, QFormLayout, QSpinBox, QDoubleSpinBox, QGroupBox, QVBoxLayout,
)

from python.road_network import SimParams

_BIG_INT = 2_000_000_000


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
        self.n_steps.setRange(10, _BIG_INT)
        self.n_steps.setValue(2000)
        self.n_steps.setToolTip(
            "n_steps: Number of CA time steps to simulate.\n"
            "More steps = longer run time but more time to reach steady state."
        )

        self.lane_length = QSpinBox()
        self.lane_length.setRange(1, _BIG_INT)
        self.lane_length.setValue(20)
        self.lane_length.setToolTip(
            "L: Number of cells per road segment.\n"
            "Each cell can hold at most one vehicle."
        )

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
        form.addRow("Road length (cells)", self.lane_length)
        form.addRow("RNG seed", self.rng_seed)
        form.addRow("Max speed (cells/step)", self.v_max)
        form.addRow("Dawdle probability", self.p_slow)
        outer.addWidget(sim_group)

        # Boundary / routing group.
        boundary = QGroupBox("Boundaries and routing")
        bform = QFormLayout(boundary)

        self.alpha = QDoubleSpinBox()
        self.alpha.setRange(0.0, 1.0)
        self.alpha.setSingleStep(0.05)
        self.alpha.setValue(0.4)
        self.alpha.setToolTip(
            "alpha (α): Entry rate.\n"
            "Probability per timestep that a new vehicle enters at an open road boundary.\n"
            "0 = no entry, 1 = always enter when space is available."
        )

        self.beta = QDoubleSpinBox()
        self.beta.setRange(0.0, 1.0)
        self.beta.setSingleStep(0.05)
        self.beta.setValue(0.5)
        self.beta.setToolTip(
            "beta (β): Exit rate.\n"
            "Probability per timestep that a vehicle at an open road exit leaves.\n"
            "0 = never leave, 1 = always leave when at the exit cell."
        )

        self.p_left = QDoubleSpinBox()
        self.p_left.setRange(0.0, 1.0)
        self.p_left.setSingleStep(0.05)
        self.p_left.setValue(0.25)
        self.p_left.setToolTip(
            "p_left: Left-turn probability at a crossroads junction.\n"
            "Fraction of vehicles that turn left when passing through the junction."
        )

        self.p_right = QDoubleSpinBox()
        self.p_right.setRange(0.0, 1.0)
        self.p_right.setSingleStep(0.05)
        self.p_right.setValue(0.25)
        self.p_right.setToolTip(
            "p_right: Right-turn probability at a crossroads junction.\n"
            "Fraction of vehicles that turn right. "
            "Straight-ahead probability = 1 − p_left − p_right."
        )

        self.p_exit = QDoubleSpinBox()
        self.p_exit.setRange(0.0, 1.0)
        self.p_exit.setSingleStep(0.05)
        self.p_exit.setValue(0.25)
        self.p_exit.setToolTip(
            "p_exit: Roundabout exit probability.\n"
            "Probability that a vehicle leaves the roundabout ring road at each arm."
        )

        bform.addRow("Entry rate α", self.alpha)
        bform.addRow("Exit rate β", self.beta)
        bform.addRow("Turn left prob.", self.p_left)
        bform.addRow("Turn right prob.", self.p_right)
        bform.addRow("Roundabout exit prob.", self.p_exit)
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
