"""QThread wrapper around pde_runner.run_pde so the GUI stays responsive."""

from __future__ import annotations

import os
from pathlib import Path

from PyQt5.QtCore import QThread, pyqtSignal

from python.pde_runner import run_pde


class PDERunnerThread(QThread):
    progress = pyqtSignal(float, str)        # frac in [0,1], text
    finished_ok = pyqtSignal(str)            # nc_path
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
