"""QThread wrappers for the GUI's background work.

Each class owns a single off-UI job and surfaces its lifecycle through
``progress`` / ``finished_ok`` / ``finished_error`` signals so the main
window can drive its status bar and stale indicators without ever
touching the worker directly.

- ``RunnerThread``    — Fortran network/CA simulation
- ``PDERunnerThread`` — Fortran LWR PDE solver
- ``FDSweepThread``   — Fortran fundamental-diagram sweep
"""

from __future__ import annotations

import os
from pathlib import Path

from PyQt5.QtCore import QThread, pyqtSignal

from python.road_network import NetworkSpec, SimParams, LayoutSpec
from python.run_simulation import run_simulation
from python.pde_runner import run_pde
from python.fd_runner import run_fd_sweep, load_fd_netcdf


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
    """Run the Fortran fundamental-diagram sweep off the UI thread.

    Supports both ``TASEP`` (α/β phase-diagram sweep on an open chain) and
    ``NS`` (density sweep on a periodic ring; uses v_max and p_slow).
    Shells out to ``build/fd_sweep`` and reads the resulting NetCDF.
    """

    # rho, J, model_name, info_dict
    finished_ok = pyqtSignal(object, object, str, dict)
    finished_error = pyqtSignal(str)

    def __init__(self, model: str = "TASEP", L: int = 50, n_points: int = 30,
                 v_max: int = 5, p_slow: float = 0.2,
                 output_path: os.PathLike = None,
                 exe: os.PathLike = None,
                 parent=None):
        super().__init__(parent)
        self.model = model
        self.L = L
        self.n_points = n_points
        self.v_max = v_max
        self.p_slow = p_slow
        self._output_path = Path(output_path) if output_path else None
        self._exe = Path(exe) if exe else None

    def run(self) -> None:
        try:
            kwargs = dict(
                model=self.model, L=self.L, n_points=self.n_points,
                v_max=self.v_max, p_slow=self.p_slow,
            )
            if self._output_path is not None:
                kwargs["output_path"] = self._output_path
            if self._exe is not None:
                kwargs["exe"] = self._exe
            out = run_fd_sweep(**kwargs)
            data = load_fd_netcdf(out)
            info = dict(
                L=self.L, n_points=self.n_points,
                v_max=self.v_max, p_slow=self.p_slow,
                n_burnin=int(data["n_burnin"]),
                n_measure=int(data["n_measure"]),
            )
            self.finished_ok.emit(data["rho"], data["J"], self.model, info)
        except Exception as exc:
            self.finished_error.emit(str(exc))
