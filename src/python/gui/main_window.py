"""Main window: top-level CA / PDE tabs, shared status bar, run log, and entry point.

The CA and PDE tabs each own their own runner thread, plot panel, and stale
indicator; this module just stitches them under a common shell, wires the
status bar / log dock to both, and provides the ``main()`` entry point.
"""

from __future__ import annotations

import os
import sys

from pathlib import Path

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout,
    QProgressBar, QLabel, QStatusBar, QPlainTextEdit, QTabWidget, QDockWidget,
)

from python.gui.ca_tab import CATab, DEFAULT_BINARY, DEFAULT_OUTDIR
from python.gui.pde_tab import PDETab, DEFAULT_PDE_BINARY, DEFAULT_PDE_OUTDIR


class MainWindow(QMainWindow):
    def __init__(self,
                 binary: os.PathLike = DEFAULT_BINARY,
                 output_dir: os.PathLike = DEFAULT_OUTDIR,
                 pde_binary: os.PathLike = DEFAULT_PDE_BINARY,
                 pde_output_dir: os.PathLike = DEFAULT_PDE_OUTDIR):
        super().__init__()
        self.setWindowTitle("PX915 Traffic Network Simulator")
        self.resize(1400, 850)

        # --- Top tabs: CA / PDE ---
        self.tabs = QTabWidget()
        # Scope styling to *this* widget only (use objectName selector) so it
        # doesn't cascade into the inner plot-panel tab widgets.
        self.tabs.setObjectName("mainTabs")
        self.tabs.setStyleSheet("""
            QTabWidget#mainTabs > QTabBar::tab {
                min-width: 220px; min-height: 30px;
                font-size: 13px; font-weight: bold;
                padding: 4px 16px;
            }
        """)
        self.ca_tab = CATab(Path(binary), Path(output_dir))
        self.pde_tab = PDETab(Path(pde_binary), Path(pde_output_dir))
        self.tabs.addTab(self.ca_tab, "Cellular Automaton (CA)")
        self.tabs.addTab(self.pde_tab, "PDE Continuum Model")
        self.setCentralWidget(self.tabs)

        # --- Status bar (shared across tabs) ---
        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setMinimumWidth(220)
        self.status_label = QLabel("ready")
        status = QStatusBar()
        status.addPermanentWidget(self.status_label, 1)
        status.addPermanentWidget(self.progress)
        self.setStatusBar(status)

        # --- Log dock (shared across tabs) ---
        self.log = QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(500)
        self.log.setMaximumHeight(140)
        log_holder = QWidget()
        ll = QVBoxLayout(log_holder)
        ll.setContentsMargins(2, 2, 2, 2)
        ll.addWidget(self.log)
        dock = QDockWidget("Run log", self)
        dock.setWidget(log_holder)
        dock.setFeatures(QDockWidget.DockWidgetMovable | QDockWidget.DockWidgetFloatable)
        self.addDockWidget(Qt.BottomDockWidgetArea, dock)

        # Wire tab signals to the shared widgets.
        for tab in (self.ca_tab, self.pde_tab):
            tab.log.connect(self.log.appendPlainText)
            tab.status.connect(self.status_label.setText)
            tab.progress.connect(self.progress.setValue)

    def closeEvent(self, event):
        # Ask each running tab to wind down before the window closes.
        self.ca_tab.stop_running()
        self.pde_tab.stop_running()
        super().closeEvent(event)


def main():
    """PyQt5 entry point — used by scripts/run_gui.py."""
    app = QApplication.instance() or QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
