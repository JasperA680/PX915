"""Tabbed panel for LWR PDE plots.

Six tabs:
  Space-time, Snapshots, Boundary flow — populated for any result.
  Per-lane space-time, Lane densities, Mass conservation — enabled only
  when n_lanes > 1.

Each tab clears the whole Figure before re-plotting so colorbars and
legends don't pile up across runs.
"""

from __future__ import annotations

from typing import Optional

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QTabWidget, QHBoxLayout, QLabel, QComboBox, QSlider

from python.gui.plot_panel import _CanvasTab
from python.pde_visualisation import (
    plot_pde_spacetime,
    plot_pde_snapshots,
    plot_pde_flow,
    plot_lane_densities,
)


def _reset_axes(tab: _CanvasTab):
    """Drop every artist (axes, colorbars) and start with a fresh single axes."""
    tab.figure.clear()
    tab.ax = tab.figure.add_subplot(111)


class _SpacetimeTab(_CanvasTab):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._data: Optional[dict] = None

    def set_data(self, data: dict):
        self._data = data
        _reset_axes(self)
        plot_pde_spacetime(data, ax=self.ax)
        self.canvas.draw_idle()


class _FlowTab(_CanvasTab):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._data: Optional[dict] = None

    def set_data(self, data: dict):
        self._data = data
        _reset_axes(self)
        plot_pde_flow(data, ax=self.ax)
        self.canvas.draw_idle()


class _SnapshotsTab(_CanvasTab):
    def __init__(self, parent=None):
        super().__init__(parent)
        row = QHBoxLayout()
        row.addWidget(QLabel("Snapshots:"))
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setRange(2, 12)
        self.slider.setValue(6)
        self.label = QLabel("6")
        self.slider.valueChanged.connect(self._on_slide)
        row.addWidget(self.slider)
        row.addWidget(self.label)
        row.addStretch(1)
        self._layout.insertLayout(0, row)
        self._data: Optional[dict] = None

    def set_data(self, data: dict):
        self._data = data
        self._refresh()

    def _on_slide(self, value: int):
        self.label.setText(str(value))
        self._refresh()

    def _refresh(self):
        if self._data is None:
            return
        _reset_axes(self)
        plot_pde_snapshots(self._data, n_snapshots=self.slider.value(), ax=self.ax)
        self.canvas.draw_idle()


class _PerLaneTab(_CanvasTab):
    def __init__(self, parent=None):
        super().__init__(parent)
        row = QHBoxLayout()
        row.addWidget(QLabel("Lane:"))
        self.lane_combo = QComboBox()
        row.addWidget(self.lane_combo)
        row.addStretch(1)
        self._layout.insertLayout(0, row)
        self._data: Optional[dict] = None
        self.lane_combo.currentIndexChanged.connect(self._refresh)

    def set_data(self, data: dict):
        self._data = data
        self.lane_combo.blockSignals(True)
        self.lane_combo.clear()
        for i in range(int(data["n_lanes"])):
            self.lane_combo.addItem(f"Lane {i + 1}", i)
        self.lane_combo.blockSignals(False)
        self._refresh()

    def _refresh(self):
        if self._data is None:
            return
        lane_idx = self.lane_combo.currentData()
        if lane_idx is None:
            lane_idx = 0
        lane_idx = int(lane_idx)
        density = self._data["density_per_lane"][:, lane_idx, :]
        x = self._data["x"]
        time = self._data["time"]
        rho_max = float(self._data["attrs"].get("rho_max", 1.0))

        _reset_axes(self)
        im = self.ax.imshow(
            density,
            aspect="auto",
            origin="lower",
            cmap="viridis",
            vmin=0,
            vmax=rho_max,
            extent=[float(x[0]), float(x[-1]), float(time[0]), float(time[-1])],
        )
        self.figure.colorbar(im, ax=self.ax, label="ρ")
        self.ax.set_xlabel("Position x")
        self.ax.set_ylabel("Time t")
        self.ax.set_title(f"Per-lane density — Lane {lane_idx + 1}")
        self.canvas.draw_idle()


class _LaneDensitiesTab(_CanvasTab):
    def __init__(self, parent=None):
        super().__init__(parent)
        row = QHBoxLayout()
        row.addWidget(QLabel("x position:"))
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setRange(0, 0)
        self.label = QLabel("0.000")
        self.slider.valueChanged.connect(self._on_slide)
        row.addWidget(self.slider)
        row.addWidget(self.label)
        row.addStretch(1)
        self._layout.insertLayout(0, row)
        self._data: Optional[dict] = None

    def set_data(self, data: dict):
        self._data = data
        x = data["x"]
        self.slider.blockSignals(True)
        self.slider.setRange(0, len(x) - 1)
        self.slider.setValue(len(x) // 2)
        self.slider.blockSignals(False)
        self._on_slide(self.slider.value())

    def _on_slide(self, value: int):
        if self._data is None:
            return
        x = self._data["x"]
        idx = max(0, min(int(value), len(x) - 1))
        x_pos = float(x[idx])
        self.label.setText(f"{x_pos:.3f}")
        _reset_axes(self)
        plot_lane_densities(self._data, x_pos=x_pos, ax=self.ax)
        self.canvas.draw_idle()


class _MassTab(_CanvasTab):
    """Inline mass conservation plot.

    pde_visualisation.plot_total_mass imports `analysis` from sys.path with
    ``src/python`` on it; the GUI prepends ``src/`` instead, so we compute the
    one-line mass deviation here rather than relying on that import.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self._data: Optional[dict] = None

    def set_data(self, data: dict):
        self._data = data
        x = data["x"]
        time = data["time"]
        dx = float(data["attrs"].get("dx", x[1] - x[0]))
        mass = data["density_per_lane"].sum(axis=(1, 2)) * dx
        deviation = mass - mass[0]
        variation = float(deviation.max() - deviation.min())

        _reset_axes(self)
        self.ax.plot(time, deviation, color="steelblue", linewidth=1.2)
        self.ax.axhline(
            0.0, color="tomato", linestyle="--", linewidth=1,
            label=f"zero (initial mass = {float(mass[0]):.4f}), range = {variation:.2e}",
        )
        self.ax.set_xlabel("Time t")
        self.ax.set_ylabel("Δ Total mass")
        self.ax.set_title("Mass conservation check  (deviation from initial mass)")
        self.ax.legend(fontsize=9)
        self.canvas.draw_idle()


# Indices for multi-lane-only tabs (set in addTab order below).
_MULTILANE_TAB_INDICES = (3, 4, 5)


class PDEPlotPanel(QTabWidget):
    """Six-tab plot panel for PDE simulation results.

    Call ``show_result(data)`` with the dict returned by ``load_pde_netcdf``.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.tab_spacetime = _SpacetimeTab()
        self.tab_snapshots = _SnapshotsTab()
        self.tab_flow = _FlowTab()
        self.tab_per_lane = _PerLaneTab()
        self.tab_lane_dens = _LaneDensitiesTab()
        self.tab_mass = _MassTab()
        self.addTab(self.tab_spacetime, "Space-time")
        self.addTab(self.tab_snapshots, "Snapshots")
        self.addTab(self.tab_flow, "Boundary flow")
        self.addTab(self.tab_per_lane, "Per-lane space-time")
        self.addTab(self.tab_lane_dens, "Lane densities")
        self.addTab(self.tab_mass, "Mass conservation")
        for i in _MULTILANE_TAB_INDICES:
            self.setTabEnabled(i, False)

    def set_multilane_enabled(self, enabled: bool):
        """Enable/disable the three multi-lane tabs without touching plot data."""
        for i in _MULTILANE_TAB_INDICES:
            self.setTabEnabled(i, enabled)

    def show_result(self, data: dict):
        self.tab_spacetime.set_data(data)
        self.tab_snapshots.set_data(data)
        self.tab_flow.set_data(data)

        multi = int(data.get("n_lanes", 1)) > 1
        self.set_multilane_enabled(multi)
        if multi:
            self.tab_per_lane.set_data(data)
            self.tab_lane_dens.set_data(data)
            self.tab_mass.set_data(data)
