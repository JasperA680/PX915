"""Tabbed panel for LWR PDE plots.

Four tabs:
  Space-time, Snapshots — populated for any result.
  Lane densities, Mass conservation — enabled only when n_lanes > 1.

On multilane results the Space-time and Snapshots tabs each show a lane
selector so individual lane density fields can be explored.

Each tab clears the whole Figure before re-plotting so colorbars and
legends don't pile up across runs.
"""

from __future__ import annotations

from typing import Optional

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QTabWidget, QHBoxLayout, QLabel, QSlider

from python.gui.plot_panel import _CanvasTab
from python.gui._common import StatusCorner, LaneSelector
from python.pde_visualisation import (
    plot_pde_spacetime,
    plot_pde_snapshots,
    plot_lane_densities,
)


class _SpacetimeTab(_CanvasTab):
    """Space-time heatmap.

    Single-lane: shows the overall density field.
    Multi-lane:  shows a per-lane selector so each lane's density can be
                 inspected individually (matches what the old _PerLaneTab did).
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        row = QHBoxLayout()
        self.lane_selector = LaneSelector()
        row.addWidget(self.lane_selector)
        row.addStretch(1)
        self._layout.insertLayout(0, row)
        self._data: Optional[dict] = None
        self.lane_selector.lane_changed.connect(self._refresh)

    def set_data(self, data: dict):
        self._data = data
        n_lanes = int(data.get("n_lanes", 1))
        self.lane_selector.set_lanes(n_lanes)
        self._refresh()

    def _refresh(self, _=None):
        if self._data is None:
            return
        self._reset_axes()
        n_lanes = int(self._data.get("n_lanes", 1))
        if n_lanes > 1:
            lane_idx = self.lane_selector.current_lane()
            density = self._data["density_per_lane"][:, lane_idx, :]
            x = self._data["x"]
            time = self._data["time"]
            rho_max = float(self._data["attrs"].get("rho_max", 1.0))
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
            attrs = self._data["attrs"]
            ic = attrs.get("ic_type", "")
            flux = attrs.get("flux_type", "")
            self.ax.set_title(f"Space-time density — Lane {lane_idx + 1}  (IC: {ic},  flux: {flux})")
        else:
            plot_pde_spacetime(self._data, ax=self.ax)
        self.canvas.draw_idle()


class _SnapshotsTab(_CanvasTab):
    """Density snapshots at evenly-spaced times.

    On multilane results a lane selector appears so each lane can be explored.
    """

    def __init__(self, parent=None):
        super().__init__(parent)

        lane_row = QHBoxLayout()
        self.lane_selector = LaneSelector()
        lane_row.addWidget(self.lane_selector)
        lane_row.addStretch(1)
        self._layout.insertLayout(0, lane_row)

        snap_row = QHBoxLayout()
        snap_row.addWidget(QLabel("Snapshots:"))
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setRange(2, 12)
        self.slider.setValue(6)
        self.label = QLabel("6")
        self.slider.valueChanged.connect(self._on_slide)
        snap_row.addWidget(self.slider)
        snap_row.addWidget(self.label)
        snap_row.addStretch(1)
        self._layout.insertLayout(1, snap_row)

        self._data: Optional[dict] = None
        self.lane_selector.lane_changed.connect(self._refresh)

    def set_data(self, data: dict):
        self._data = data
        n_lanes = int(data.get("n_lanes", 1))
        self.lane_selector.set_lanes(n_lanes)
        self._refresh()

    def _on_slide(self, value: int):
        self.label.setText(str(value))
        self._refresh()

    def _refresh(self, _=None):
        if self._data is None:
            return
        self._reset_axes()
        n_lanes = int(self._data.get("n_lanes", 1))
        if n_lanes > 1:
            lane_idx = self.lane_selector.current_lane()
            # Build a temporary single-lane data view for the selected lane
            lane_data = dict(self._data)
            lane_data["density"] = self._data["density_per_lane"][:, lane_idx, :]
            attrs = dict(self._data["attrs"])
            attrs["ic_type"] = f'{attrs.get("ic_type", "")} (Lane {lane_idx + 1})'
            lane_data["attrs"] = attrs
            plot_pde_snapshots(lane_data, n_snapshots=self.slider.value(), ax=self.ax)
        else:
            plot_pde_snapshots(self._data, n_snapshots=self.slider.value(), ax=self.ax)
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
        self._reset_axes()
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

        self._reset_axes()
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


class PDEPlotPanel(QTabWidget):
    """Four-tab plot panel for PDE simulation results.

    Tab order: Space-time (0), Snapshots (1),
               Lane densities (2, multilane only), Mass conservation (3, multilane only).

    Call ``show_result(data)`` with the dict returned by ``load_pde_netcdf``.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.tab_spacetime = _SpacetimeTab()
        self.tab_snapshots = _SnapshotsTab()
        self.tab_lane_dens = _LaneDensitiesTab()
        self.tab_mass = _MassTab()
        self.addTab(self.tab_spacetime, "Space-time")
        self.addTab(self.tab_snapshots, "Snapshots")
        self.addTab(self.tab_lane_dens, "Lane densities")
        self.addTab(self.tab_mass, "Mass conservation")
        # Lookups are done lazily so they survive tab reordering.
        self._multilane_tabs = (self.tab_lane_dens, self.tab_mass)
        self.set_multilane_enabled(False)

        self._status_corner = StatusCorner()
        self.setCornerWidget(self._status_corner, Qt.TopRightCorner)

    def set_status(self, kind, summary: str = ""):
        self._status_corner.set_status(kind, summary)

    def set_multilane_enabled(self, enabled: bool):
        """Enable/disable the two multi-lane tabs without touching plot data."""
        for tab in self._multilane_tabs:
            self.setTabEnabled(self.indexOf(tab), bool(enabled))

    def show_result(self, data: dict):
        self.tab_spacetime.set_data(data)
        self.tab_snapshots.set_data(data)

        multi = int(data.get("n_lanes", 1)) > 1
        self.set_multilane_enabled(multi)
        if multi:
            self.tab_lane_dens.set_data(data)
            self.tab_mass.set_data(data)
