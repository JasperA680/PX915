"""Dialog for editing per-road parameters (α, β, length).

Opens when the user clicks a road in the NetworkWidget preview.
"""

from __future__ import annotations

from PyQt5.QtWidgets import (
    QDialog, QDialogButtonBox, QFormLayout, QDoubleSpinBox, QSpinBox,
    QVBoxLayout, QLabel,
)

_BIG_INT = 2_000_000_000
_BIG_FLOAT = 1e12


class RoadEditDialog(QDialog):
    """Pop-up editor for a single road's open-boundary rates and length."""

    def __init__(self, road_id: int, alpha: float, beta: float, length: int,
                 has_alpha: bool, has_beta: bool, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Edit Road R{road_id}")
        self.setMinimumWidth(300)

        form = QFormLayout()

        if has_alpha:
            self.alpha_box = QDoubleSpinBox()
            self.alpha_box.setRange(0.0, 1.0)
            self.alpha_box.setSingleStep(0.05)
            self.alpha_box.setDecimals(3)
            self.alpha_box.setValue(alpha)
            form.addRow("Entry rate α (0–1):", self.alpha_box)
        else:
            self.alpha_box = None
            form.addRow(QLabel("Entry rate α:"), QLabel("N/A (no open entry boundary)"))

        if has_beta:
            self.beta_box = QDoubleSpinBox()
            self.beta_box.setRange(0.0, 1.0)
            self.beta_box.setSingleStep(0.05)
            self.beta_box.setDecimals(3)
            self.beta_box.setValue(beta)
            form.addRow("Exit rate β (0–1):", self.beta_box)
        else:
            self.beta_box = None
            form.addRow(QLabel("Exit rate β:"), QLabel("N/A (no open exit boundary)"))

        self.length_box = QSpinBox()
        self.length_box.setRange(1, _BIG_INT)
        self.length_box.setValue(max(1, int(length)))
        form.addRow("Length (cells):", self.length_box)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        outer = QVBoxLayout(self)
        outer.addLayout(form)
        outer.addWidget(buttons)

    def values(self) -> dict:
        return dict(
            alpha=float(self.alpha_box.value()) if self.alpha_box else None,
            beta=float(self.beta_box.value()) if self.beta_box else None,
            length=int(self.length_box.value()),
        )
