"""Dialogs for editing per-road and per-junction parameters.

Opens when the user clicks a road or junction in the NetworkWidget preview.
"""

from __future__ import annotations

from PyQt5.QtWidgets import (
    QDialog, QDialogButtonBox, QFormLayout, QDoubleSpinBox, QSpinBox,
    QVBoxLayout, QHBoxLayout, QLabel, QGridLayout, QGroupBox,
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


class JunctionEditDialog(QDialog):
    """Editor for a junction's routing matrix.

    The matrix has ``n_in`` rows (one per inbound leg) and ``n_out`` columns
    (one per outbound leg).  Each row should sum to 1.  We don't auto-enforce
    that — the user can normalise via the "Normalise rows" button or just type
    values that sum to 1.
    """

    def __init__(self, junction_id: int, routes, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Edit Junction J{junction_id} routing")
        self.setMinimumWidth(360)
        self._n_in  = len(routes)
        self._n_out = len(routes[0]) if routes else 0
        self._boxes: list = []   # [[spinbox, …], …] in row-major order

        info = QLabel(
            f"Routing matrix: {self._n_in} inbound × {self._n_out} outbound legs.\n"
            "Each row is the probability of an incoming vehicle taking each\n"
            "outbound leg.  Rows should sum to 1.0."
        )
        info.setWordWrap(True)

        grid_group = QGroupBox("Routes  (rows = incoming leg, cols = outgoing leg)")
        grid = QGridLayout(grid_group)
        # Column headers
        for c in range(self._n_out):
            grid.addWidget(QLabel(f"→ out {c + 1}"), 0, c + 1)
        for r in range(self._n_in):
            grid.addWidget(QLabel(f"in {r + 1}"), r + 1, 0)
            row_boxes = []
            for c in range(self._n_out):
                sb = QDoubleSpinBox()
                sb.setRange(0.0, 1.0)
                sb.setSingleStep(0.05)
                sb.setDecimals(3)
                sb.setValue(float(routes[r][c]))
                grid.addWidget(sb, r + 1, c + 1)
                row_boxes.append(sb)
            self._boxes.append(row_boxes)

        normalise_btn = QDialogButtonBox(QDialogButtonBox.NoButton)
        from PyQt5.QtWidgets import QPushButton
        norm_pb = QPushButton("Normalise rows to sum=1")
        norm_pb.clicked.connect(self._normalise)
        normalise_btn.addButton(norm_pb, QDialogButtonBox.ActionRole)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        outer = QVBoxLayout(self)
        outer.addWidget(info)
        outer.addWidget(grid_group)
        outer.addWidget(normalise_btn)
        outer.addWidget(buttons)

    def _normalise(self):
        for row in self._boxes:
            total = sum(b.value() for b in row)
            if total <= 0:
                # Set uniform 1/n distribution.
                v = 1.0 / max(1, len(row))
                for b in row:
                    b.blockSignals(True)
                    b.setValue(v)
                    b.blockSignals(False)
            else:
                for b in row:
                    b.blockSignals(True)
                    b.setValue(b.value() / total)
                    b.blockSignals(False)

    def routes(self) -> list:
        """Return the edited routing matrix as a list of lists of floats."""
        return [[float(b.value()) for b in row] for row in self._boxes]
