"""PyQt5 application entry point."""

from __future__ import annotations

import sys

from PyQt5.QtWidgets import QApplication

from python.gui.main_window import MainWindow


def main():
    app = QApplication.instance() or QApplication(sys.argv)
    win = MainWindow()
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
