"""QThread wrapper around run_simulation so the GUI stays responsive."""

from __future__ import annotations

import os
from pathlib import Path

from PyQt5.QtCore import QThread, pyqtSignal

from python.road_network import NetworkSpec, SimParams, LayoutSpec
from python.run_simulation import run_simulation


class RunnerThread(QThread):
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
        def cb(frac: float, msg: str) -> None:
            self.progress.emit(frac, msg)

        try:
            nc_path = run_simulation(
                self._spec, self._params, self._layout,
                output_dir=self._output_dir, binary=self._binary,
                progress=cb,
            )
        except Exception as exc:
            self.finished_error.emit(str(exc))
            return
        self.finished_ok.emit(str(nc_path))
