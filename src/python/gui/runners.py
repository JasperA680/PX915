"""QThread wrappers for the GUI's background work.

Each class owns a single off-UI job and surfaces its lifecycle through
``progress`` / ``finished_ok`` / ``finished_error`` signals so the main
window can drive its status bar and stale indicators without ever
touching the worker directly.

- ``RunnerThread``    — Fortran network/CA simulation
- ``PDERunnerThread`` — Fortran LWR PDE solver
- ``FDSweepThread``   — pure-Python fundamental-diagram sweep
"""

from __future__ import annotations

import os
from pathlib import Path

from PyQt5.QtCore import QThread, pyqtSignal

from python.road_network import NetworkSpec, SimParams, LayoutSpec
from python.run_simulation import run_simulation
from python.pde_runner import run_pde


class RunnerThread(QThread):
    """Background driver for the CA / network Fortran binary."""

    progress = pyqtSignal(float, str)        # frac in [0,1], text
    finished_ok = pyqtSignal(str)            # nc_path
    finished_error = pyqtSignal(str)

    def __init__(self, spec: NetworkSpec, params: SimParams, layout: LayoutSpec,
                 output_dir: os.PathLike, binary: os.PathLike, parent=None):
        super().__init__(parent)
        self._spec = spec
        self._params = params
        self._layout = layout
        self._output_dir = Path(output_dir)
        self._binary = Path(binary)

    def run(self) -> None:
        try:
            nc_path = run_simulation(
                self._spec, self._params, self._layout,
                output_dir=self._output_dir, binary=self._binary,
                progress=lambda frac, msg: self.progress.emit(frac, msg),
            )
        except Exception as exc:
            self.finished_error.emit(str(exc))
            return
        self.finished_ok.emit(str(nc_path))


class PDERunnerThread(QThread):
    """Background driver for the LWR PDE Fortran solver."""

    progress = pyqtSignal(float, str)
    finished_ok = pyqtSignal(str)
    finished_error = pyqtSignal(str)

    def __init__(self, params: dict, output_path: os.PathLike, exe: os.PathLike, parent=None):
        super().__init__(parent)
        self._params = params
        self._output_path = Path(output_path)
        self._exe = Path(exe)

    def run(self) -> None:
        try:
            self.progress.emit(0.0, "starting PDE solver...")
            self._output_path.parent.mkdir(parents=True, exist_ok=True)
            run_pde(self._params, output_path=self._output_path, exe=self._exe)
        except Exception as exc:
            self.finished_error.emit(str(exc))
            return
        self.progress.emit(1.0, "")
        self.finished_ok.emit(str(self._output_path))


class FDSweepThread(QThread):
    """Run a pure-Python fundamental-diagram sweep off the UI thread.

    Supports both ``TASEP`` (α/β phase-diagram sweep on an open chain) and
    ``NS`` (density sweep on a periodic ring; takes v_max and p_slow).
    """

    # rho, J, model_name, info_dict
    finished_ok = pyqtSignal(object, object, str, dict)
    finished_error = pyqtSignal(str)

    def __init__(self, model: str = "TASEP", L: int = 50, n_points: int = 30,
                 v_max: int = 5, p_slow: float = 0.2, parent=None):
        super().__init__(parent)
        self.model = model
        self.L = L
        self.n_points = n_points
        self.v_max = v_max
        self.p_slow = p_slow

    def run(self) -> None:
        try:
            if self.model == "NS":
                from python.analysis import fundamental_diagram_ns
                rho, J = fundamental_diagram_ns(
                    L=self.L, n_points=self.n_points,
                    v_max=self.v_max, p_slow=self.p_slow,
                )
            else:
                from python.analysis import fundamental_diagram
                rho, J = fundamental_diagram(L=self.L, n_points=self.n_points)
            info = dict(L=self.L, n_points=self.n_points,
                        v_max=self.v_max, p_slow=self.p_slow)
            self.finished_ok.emit(rho, J, self.model, info)
        except Exception as exc:
            self.finished_error.emit(str(exc))
