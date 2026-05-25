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

    Rows = inbound legs (one per inbound lane); columns = outbound legs.
    Each header is labelled with the road number the leg attaches to (e.g.
    ``from R3`` / ``→ R2``) so the user can match the matrix back to the
    network preview.  Cells where the inbound and outbound legs share the
    same road are U-turns — forced to 0 and read-only.  Each row should sum
    to 1; click "Normalise rows" to make it so (locked cells are excluded).
    """

    def __init__(self, junction_id: int, routes,
                 in_road_ids=None, out_road_ids=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Edit Junction J{junction_id} routing")
        self.setMinimumWidth(360)
        self._n_in  = len(routes)
        self._n_out = len(routes[0]) if routes else 0
        if in_road_ids is None:
            in_road_ids = list(range(1, self._n_in + 1))
        if out_road_ids is None:
            out_road_ids = list(range(1, self._n_out + 1))
        self._in_road_ids = list(in_road_ids)
        self._out_road_ids = list(out_road_ids)
        # locked_cells[(r, c)] = True when cell is a U-turn (same road).
        self._locked = {
            (r, c): (self._in_road_ids[r] == self._out_road_ids[c])
            for r in range(self._n_in) for c in range(self._n_out)
        }
        self._boxes: list = []   # [[spinbox, …], …] in row-major order

        info = QLabel(
            f"Routing matrix: {self._n_in} inbound × {self._n_out} outbound legs.\n"
            "Each row is the probability of an incoming vehicle exiting via\n"
            "each outbound leg.  U-turns (same-road cells) are locked at 0.\n"
            "Rows should sum to 1.0 across the editable cells."
        )
        info.setWordWrap(True)

        grid_group = QGroupBox("Routes  (rows = arriving from, cols = leaving to)")
        grid = QGridLayout(grid_group)
        # Top-left corner placeholder.
        grid.addWidget(QLabel(""), 0, 0)
        # Column headers — labelled by the road each outbound leg leads to.
        for c in range(self._n_out):
            grid.addWidget(QLabel(f"→ R{self._out_road_ids[c]}"), 0, c + 1)
        for r in range(self._n_in):
            grid.addWidget(QLabel(f"from R{self._in_road_ids[r]}"), r + 1, 0)
            row_boxes = []
            for c in range(self._n_out):
                sb = QDoubleSpinBox()
                sb.setRange(0.0, 1.0)
                sb.setSingleStep(0.05)
                sb.setDecimals(3)
                if self._locked[(r, c)]:
                    sb.setValue(0.0)
                    sb.setEnabled(False)
                    sb.setToolTip("U-turn (same road) — locked at 0")
                else:
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
        # OK auto-normalises so submitted rows are guaranteed to sum to 1.
        buttons.accepted.connect(self._on_accept)
        buttons.rejected.connect(self.reject)

        outer = QVBoxLayout(self)
        outer.addWidget(info)
        outer.addWidget(grid_group)
        outer.addWidget(normalise_btn)
        outer.addWidget(buttons)

    def _on_accept(self):
        self._normalise()
        self.accept()

    def _normalise(self):
        """Normalise each row so the editable cells sum to 1; locked cells stay 0."""
        for r, row in enumerate(self._boxes):
            editable = [(c, b) for c, b in enumerate(row) if not self._locked[(r, c)]]
            if not editable:
                continue
            total = sum(b.value() for _, b in editable)
            if total <= 0:
                v = 1.0 / len(editable)
                for _, b in editable:
                    b.blockSignals(True); b.setValue(v); b.blockSignals(False)
            else:
                for _, b in editable:
                    b.blockSignals(True)
                    b.setValue(b.value() / total)
                    b.blockSignals(False)

    def routes(self) -> list:
        """Return the edited routing matrix as a list of lists of floats.

        Locked (U-turn) cells are always returned as 0.0.
        """
        out = []
        for r, row in enumerate(self._boxes):
            out.append([0.0 if self._locked[(r, c)] else float(b.value())
                        for c, b in enumerate(row)])
        return out
