"""Tabbed panel that shows density, currents, space-time, and network heatmap."""

from __future__ import annotations

from typing import Optional

from PyQt5.QtCore import Qt, pyqtSignal, QEvent
from PyQt5.QtGui import QStandardItem, QStandardItemModel
from PyQt5.QtWidgets import (
    QTabWidget, QWidget, QVBoxLayout, QHBoxLayout,
    QComboBox, QSlider, QLabel, QListView,
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from python.io import NetworkResult
from python.visualisation import (
    plot_network_density,
    plot_network_currents,
    plot_network_spacetime,
    plot_network_layout,
)


class _CanvasTab(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._layout = QVBoxLayout(self)
        self._layout.setContentsMargins(2, 2, 2, 2)
        self.figure = Figure(figsize=(6, 4))
        self.canvas = FigureCanvasQTAgg(self.figure)
        self._layout.addWidget(self.canvas)
        self.ax = self.figure.add_subplot(111)

    def clear(self):
        self.ax.clear()
        self.canvas.draw_idle()


# ---------------------------------------------------------------------------
# Checkable multi-select combo box
# ---------------------------------------------------------------------------

class _CheckableComboBox(QComboBox):
    """QComboBox whose popup items are individually checkable.

    Keeps the popup open while the user ticks boxes; closes on click-away.
    Emits ``selectionChanged(list[int])`` with 0-based checked indices on
    every change.
    """

    selectionChanged = pyqtSignal(list)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setEditable(True)
        self.lineEdit().setReadOnly(True)
        self.setModel(QStandardItemModel(0, 1, self))
        self.setView(QListView(self))
        # Eat MouseButtonRelease on the viewport so the popup stays open
        # while the user is ticking boxes.
        self.view().viewport().installEventFilter(self)
        self.view().pressed.connect(self._on_item_pressed)

    def eventFilter(self, source, event):
        if (source is self.view().viewport()
                and event.type() == QEvent.MouseButtonRelease):
            return True
        return super().eventFilter(source, event)

    def _on_item_pressed(self, index):
        item = self.model().itemFromIndex(index)
        new_state = Qt.Unchecked if item.checkState() == Qt.Checked else Qt.Checked
        item.setCheckState(new_state)
        self._update_display()
        self.selectionChanged.emit(self.checked_indices())

    def populate(self, labels: list, all_checked: bool = True):
        """Replace all items with the given labels (all checked by default)."""
        self.model().clear()
        state = Qt.Checked if all_checked else Qt.Unchecked
        for label in labels:
            item = QStandardItem(label)
            item.setFlags(Qt.ItemIsEnabled)
            item.setCheckState(state)
            self.model().appendRow(item)
        self._update_display()

    def checked_indices(self) -> list:
        return [i for i in range(self.model().rowCount())
                if self.model().item(i).checkState() == Qt.Checked]

    def _update_display(self):
        n = self.model().rowCount()
        checked = self.checked_indices()
        if not n:
            self.lineEdit().setText("")
        elif len(checked) == n:
            self.lineEdit().setText("All roads")
        elif not checked:
            self.lineEdit().setText("No roads selected")
        else:
            self.lineEdit().setText(f"{len(checked)} / {n} roads")


# ---------------------------------------------------------------------------
# Density and Currents tabs with per-road checkable selector
# ---------------------------------------------------------------------------

class _DensityTab(_CanvasTab):
    def __init__(self, parent=None):
        super().__init__(parent)
        row = QHBoxLayout()
        row.addWidget(QLabel("Roads:"))
        self.road_selector = _CheckableComboBox()
        self.road_selector.setMinimumWidth(140)
        row.addWidget(self.road_selector)
        row.addStretch(1)
        self._layout.insertLayout(0, row)
        self._result: Optional[NetworkResult] = None
        self.road_selector.selectionChanged.connect(self._refresh)

    def set_result(self, result: NetworkResult):
        self._result = result
        n_roads = result.road_density.shape[1]
        self.road_selector.selectionChanged.disconnect(self._refresh)
        self.road_selector.populate([f"R{r + 1}" for r in range(n_roads)])
        self.road_selector.selectionChanged.connect(self._refresh)
        self._refresh()

    def _refresh(self):
        if self._result is None:
            return
        selected = self.road_selector.checked_indices()
        self.ax.clear()
        plot_network_density(self._result, ax=self.ax,
                             road_ids=selected if selected else None)
        self.canvas.draw_idle()


class _CurrentsTab(_CanvasTab):
    def __init__(self, parent=None):
        super().__init__(parent)
        row = QHBoxLayout()
        row.addWidget(QLabel("Roads:"))
        self.road_selector = _CheckableComboBox()
        self.road_selector.setMinimumWidth(140)
        row.addWidget(self.road_selector)
        row.addStretch(1)
        self._layout.insertLayout(0, row)
        self._result: Optional[NetworkResult] = None
        self.road_selector.selectionChanged.connect(self._refresh)

    def set_result(self, result: NetworkResult):
        self._result = result
        n_roads = result.road_exits.shape[1]
        self.road_selector.selectionChanged.disconnect(self._refresh)
        self.road_selector.populate([f"R{r + 1}" for r in range(n_roads)])
        self.road_selector.selectionChanged.connect(self._refresh)
        self._refresh()

    def _refresh(self):
        if self._result is None:
            return
        selected = self.road_selector.checked_indices()
        self.ax.clear()
        plot_network_currents(self._result, ax=self.ax,
                              road_ids=selected if selected else None)
        self.canvas.draw_idle()


# ---------------------------------------------------------------------------
# Space-time tab (Road + Lane selectors)
# ---------------------------------------------------------------------------

class _SpacetimeTab(_CanvasTab):
    def __init__(self, parent=None):
        super().__init__(parent)
        import numpy as np
        self._np = np

        row = QHBoxLayout()
        row.addWidget(QLabel("Road:"))
        self.road_combo = QComboBox()
        row.addWidget(self.road_combo)
        row.addSpacing(12)
        self._lane_label = QLabel("Lane:")
        self.lane_combo = QComboBox()
        row.addWidget(self._lane_label)
        row.addWidget(self.lane_combo)
        row.addStretch(1)
        self._layout.insertLayout(0, row)
        self._result: Optional[NetworkResult] = None
        self.road_combo.currentIndexChanged.connect(self._on_road_changed)
        self.lane_combo.currentIndexChanged.connect(self._refresh)
        self._lane_label.setVisible(False)
        self.lane_combo.setVisible(False)

    def set_result(self, result: NetworkResult):
        self._result = result
        self.road_combo.blockSignals(True)
        self.road_combo.clear()
        n_roads = result.road_density.shape[1]
        for r in range(1, n_roads + 1):
            self.road_combo.addItem(f"R{r}", r)
        self.road_combo.blockSignals(False)
        self._on_road_changed()

    def _on_road_changed(self):
        """Repopulate lane combo for the newly selected road, then refresh."""
        if self._result is None:
            return
        rid = self.road_combo.currentData() or 1
        lane_idxs = self._np.where(self._result.lane_road_id == int(rid))[0]
        n_lanes = len(lane_idxs)

        self.lane_combo.blockSignals(True)
        self.lane_combo.clear()
        for i in range(n_lanes):
            self.lane_combo.addItem(f"Lane {i + 1}", i)
        self.lane_combo.blockSignals(False)

        show_lane = n_lanes > 1
        self._lane_label.setVisible(show_lane)
        self.lane_combo.setVisible(show_lane)

        self._refresh()

    def _refresh(self):
        if self._result is None:
            return
        self.ax.clear()
        rid = self.road_combo.currentData() or 1
        lane_idx = self.lane_combo.currentData()
        if lane_idx is None:
            lane_idx = 0
        plot_network_spacetime(self._result, road_id=int(rid),
                               lane_index=int(lane_idx), ax=self.ax)
        self.canvas.draw_idle()


# ---------------------------------------------------------------------------
# Network heatmap tab
# ---------------------------------------------------------------------------

class _HeatmapTab(_CanvasTab):
    def __init__(self, parent=None):
        super().__init__(parent)
        row = QHBoxLayout()
        row.addWidget(QLabel("Timestep:"))
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setMinimum(0)
        self.slider.setMaximum(0)
        self.label = QLabel("0")
        self.slider.valueChanged.connect(self._on_slide)
        row.addWidget(self.slider)
        row.addWidget(self.label)
        self._layout.insertLayout(0, row)
        self._result: Optional[NetworkResult] = None

    def set_result(self, result: NetworkResult):
        self._result = result
        self.slider.blockSignals(True)
        self.slider.setMinimum(0)
        self.slider.setMaximum(max(0, result.occupancy.shape[0] - 1))
        self.slider.setValue(result.occupancy.shape[0] - 1)
        self.slider.blockSignals(False)
        self._on_slide(self.slider.value())

    def _on_slide(self, value: int):
        if self._result is None:
            return
        self.label.setText(str(value + 1))
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        plot_network_layout(self._result, ax=self.ax, occupancy_t=value)
        self.canvas.draw_idle()


# ---------------------------------------------------------------------------
# Container panel
# ---------------------------------------------------------------------------

class PlotPanel(QTabWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.tab_density = _DensityTab()
        self.tab_currents = _CurrentsTab()
        self.tab_spacetime = _SpacetimeTab()
        self.tab_heatmap = _HeatmapTab()
        self.addTab(self.tab_density,   "Density")
        self.addTab(self.tab_currents,  "Currents")
        self.addTab(self.tab_spacetime, "Space-time")
        self.addTab(self.tab_heatmap,   "Network heatmap")

    def clear(self):
        self.tab_density.clear()
        self.tab_currents.clear()
        self.tab_spacetime.clear()
        self.tab_heatmap.clear()

    def show_result(self, result: NetworkResult):
        self.tab_density.set_result(result)
        self.tab_currents.set_result(result)
        self.tab_spacetime.set_result(result)
        self.tab_heatmap.set_result(result)
