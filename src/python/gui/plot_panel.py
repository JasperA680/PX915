"""Tabbed panel that shows density, currents, space-time, and network heatmap."""

from __future__ import annotations

from typing import Optional

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QTabWidget, QWidget, QVBoxLayout, QHBoxLayout, QComboBox, QSlider, QLabel
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


class _SpacetimeTab(_CanvasTab):
    def __init__(self, parent=None):
        super().__init__(parent)
        row = QHBoxLayout()
        row.addWidget(QLabel("Road:"))
        self.road_combo = QComboBox()
        row.addWidget(self.road_combo)
        row.addStretch(1)
        self._layout.insertLayout(0, row)
        self._result: Optional[NetworkResult] = None
        self.road_combo.currentIndexChanged.connect(self._refresh)

    def set_result(self, result: NetworkResult):
        self._result = result
        self.road_combo.blockSignals(True)
        self.road_combo.clear()
        n_roads = result.road_density.shape[1]
        for r in range(1, n_roads + 1):
            self.road_combo.addItem(f"R{r}", r)
        self.road_combo.blockSignals(False)
        self._refresh()

    def _refresh(self):
        if self._result is None:
            return
        self.ax.clear()
        rid = self.road_combo.currentData() or 1
        plot_network_spacetime(self._result, road_id=int(rid), ax=self.ax)
        self.canvas.draw_idle()


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
        self.label.setText(str(value))
        self.ax.clear()
        plot_network_layout(self._result, ax=self.ax, occupancy_t=value)
        self.canvas.draw_idle()


class PlotPanel(QTabWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.tab_density = _CanvasTab()
        self.tab_currents = _CanvasTab()
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
        self.tab_density.ax.clear()
        plot_network_density(result, ax=self.tab_density.ax)
        self.tab_density.canvas.draw_idle()

        self.tab_currents.ax.clear()
        plot_network_currents(result, ax=self.tab_currents.ax)
        self.tab_currents.canvas.draw_idle()

        self.tab_spacetime.set_result(result)
        self.tab_heatmap.set_result(result)
