"""Shared GUI primitives used by both the CA and PDE panels.

Lives alongside the rest of ``python.gui`` and is private to it; nothing
outside the GUI package should import from here.
"""

from __future__ import annotations

from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import (
    QWidget, QHBoxLayout, QLabel, QComboBox, QPushButton,
)


# Generous numeric ranges for QSpinBox / QDoubleSpinBox; matches the values
# the param forms have used since the first GUI revision.
BIG_INT = 2_000_000_000
BIG_FLOAT = 1e12


# Stylesheet applied to QScrollArea instances that wrap the parameter forms.
# macOS auto-hides native scrollbars regardless of ScrollBarAlwaysOn, so the
# scrollbar can't act as a "there is more below" cue. Applying any
# QScrollBar stylesheet kicks Qt off the native macOS style and onto its own
# renderer, which honours AlwaysOn properly. The visual is a slim neutral-grey
# bar that fades darker on hover.
ALWAYS_VISIBLE_SCROLLBAR_QSS = """
QScrollBar:vertical {
    background: rgba(0, 0, 0, 0.04);
    width: 10px;
    margin: 0;
    border-radius: 5px;
}
QScrollBar::handle:vertical {
    background: rgba(0, 0, 0, 0.30);
    min-height: 24px;
    border-radius: 5px;
}
QScrollBar::handle:vertical:hover {
    background: rgba(0, 0, 0, 0.50);
}
QScrollBar::add-line:vertical,
QScrollBar::sub-line:vertical {
    height: 0;
    background: none;
}
QScrollBar::add-page:vertical,
QScrollBar::sub-page:vertical {
    background: none;
}
"""


def set_button_stale(button: QPushButton, base_label: str, stale: bool) -> None:
    """Style ``button`` to signal whether its last-produced output is stale.

    Stale state appends a ●, sets bold + amber text; clearing restores the
    plain label and default stylesheet.
    """
    if stale:
        button.setText(f"{base_label}  ●")
        button.setStyleSheet("color: #b08040; font-weight: bold;")
    else:
        button.setText(base_label)
        button.setStyleSheet("")


class StatusCorner(QWidget):
    """Top-right corner indicator for a QTabWidget.

    Shows one of four states: ``fresh`` (green dot), ``stale`` (amber),
    ``running`` (blue), or empty (grey).  Install via
    ``QTabWidget.setCornerWidget(self, Qt.TopRightCorner)``.
    """

    _PALETTE = {
        "fresh":   ("#4a8a6a", "● up to date"),
        "stale":   ("#b08040", "● stale"),
        "running": ("#5a7fad", "● running…"),
    }
    _EMPTY = ("#9a9a9a", "● no results")

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(8, 0, 14, 4)
        self._label = QLabel()
        layout.addWidget(self._label)
        self.set_status(None)

    def set_status(self, kind, summary: str = "") -> None:
        """Update the indicator.  ``summary`` is accepted for API parity."""
        color, text = self._PALETTE.get(kind, self._EMPTY)
        self._label.setText(text)
        self._label.setStyleSheet(
            f"color: {color}; font-size: 10px; font-weight: 500;"
        )


class LaneSelector(QWidget):
    """Horizontal 'Lane:' label + combo, hidden when only one lane is present.

    Combo items are labelled ``Lane 1`` … ``Lane N`` and carry their 0-based
    index as user data.  Emits ``lane_changed(int)`` (0-based) on user
    interaction; ``set_lanes`` is silent.
    """

    lane_changed = pyqtSignal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(QLabel("Lane:"))
        self._combo = QComboBox()
        self._combo.setToolTip(
            "Which lane of the selected road to display.  Only shown when "
            "the data has more than one lane."
        )
        layout.addWidget(self._combo)
        self._combo.currentIndexChanged.connect(self._emit)
        self.setVisible(False)

    def set_lanes(self, n_lanes: int) -> None:
        self._combo.blockSignals(True)
        self._combo.clear()
        for i in range(n_lanes):
            self._combo.addItem(f"Lane {i + 1}", i)
        self._combo.blockSignals(False)
        self.setVisible(n_lanes > 1)

    def current_lane(self) -> int:
        data = self._combo.currentData()
        return 0 if data is None else int(data)

    def _emit(self, _idx: int) -> None:
        self.lane_changed.emit(self.current_lane())
