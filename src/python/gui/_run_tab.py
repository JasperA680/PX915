"""Common lifecycle scaffolding for the CA and PDE tabs.

Both tabs follow the same skeleton: log/status/progress signals, a single
background ``QThread`` runner whose result toggles ``has_results``, and a
trio of fresh/stale/running visual markers driven by a
``StatusCorner`` on the plot panel + a stale dot on the Run button.
"""

from __future__ import annotations

from typing import Optional

from PyQt5.QtCore import QThread, pyqtSignal
from PyQt5.QtWidgets import QPushButton, QWidget

from python.gui._common import set_button_stale


class RunTab(QWidget):
    """Base class for the CA and PDE tabs.

    Subclasses are expected to assign ``self.run_button`` and
    ``self.plot_panel`` during ``__init__`` (the helpers below read both).
    """

    log = pyqtSignal(str)
    status = pyqtSignal(str)
    progress = pyqtSignal(int)

    _RUN_LABEL = "Run"   # subclasses may override

    def __init__(self, parent=None):
        super().__init__(parent)
        self._runner: Optional[QThread] = None
        self._has_results: bool = False
        self._last_run_summary: str = ""
        # Subclasses set these before any _mark_* is called.
        self.run_button: Optional[QPushButton] = None
        self.plot_panel = None

    # ------------------------------------------------------------------
    # fresh / stale / running state
    # ------------------------------------------------------------------

    def _mark_stale(self) -> None:
        if self._has_results:
            self.plot_panel.set_status("stale", self._last_run_summary)
        else:
            self.plot_panel.set_status(None)
        set_button_stale(self.run_button, self._RUN_LABEL, self._has_results)

    def _mark_running(self, summary: str) -> None:
        self.plot_panel.set_status("running", summary)
        set_button_stale(self.run_button, self._RUN_LABEL, False)

    def _mark_fresh(self, summary: str) -> None:
        self._last_run_summary = summary
        self._has_results = True
        self.plot_panel.set_status("fresh", summary)
        set_button_stale(self.run_button, self._RUN_LABEL, False)

    # ------------------------------------------------------------------
    # runner-thread plumbing
    # ------------------------------------------------------------------

    def _on_runner_progress(self, frac: float, msg: str) -> None:
        if frac >= 0:
            self.progress.emit(int(round(100 * frac)))
        if msg:
            self.log.emit(msg)

    def _on_finished_error(self, err: str) -> None:
        self.status.emit("error")
        self.run_button.setEnabled(True)
        self.log.emit(f"ERROR: {err}")
        self._mark_stale()

    def stop_running(self) -> None:
        if self._runner is not None and self._runner.isRunning():
            self._runner.requestInterruption()
            self._runner.wait(2000)
