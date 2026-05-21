"""Tabbed panel for LWR PDE plots.

Five tabs:
  Space-time, Snapshots, Boundary flow — populated for any result.
  Lane densities, Mass conservation — enabled only when n_lanes > 1.

On multilane results the Space-time and Snapshots tabs each show a lane
selector so individual lane density fields can be explored.

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
    """Space-time heatmap.

    Single-lane: shows the overall density field.
    Multi-lane:  shows a per-lane selector so each lane's density can be
                 inspected individually (matches what the old _PerLaneTab did).
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        row = QHBoxLayout()
        row.addWidget(QLabel("Lane:"))
        self.lane_combo = QComboBox()
        row.addWidget(self.lane_combo)
        row.addStretch(1)
        self._lane_row_widget = row   # kept so we can show/hide it
        self._layout.insertLayout(0, row)
        self._data: Optional[dict] = None
        self.lane_combo.currentIndexChanged.connect(self._refresh)
        # Hide lane selector until multilane data arrives
        self._set_lane_row_visible(False)

    def _set_lane_row_visible(self, visible: bool):
        for i in range(self._lane_row_widget.count()):
            item = self._lane_row_widget.itemAt(i)
            if item and item.widget():
                item.widget().setVisible(visible)

    def set_data(self, data: dict):
        self._data = data
        n_lanes = int(data.get("n_lanes", 1))
        if n_lanes > 1:
            self.lane_combo.blockSignals(True)
            self.lane_combo.clear()
            for i in range(n_lanes):
                self.lane_combo.addItem(f"Lane {i + 1}", i)
            self.lane_combo.blockSignals(False)
            self._set_lane_row_visible(True)
        else:
            self._set_lane_row_visible(False)
        self._refresh()

    def _refresh(self):
        if self._data is None:
            return
        _reset_axes(self)
        n_lanes = int(self._data.get("n_lanes", 1))
        if n_lanes > 1:
            lane_idx = self.lane_combo.currentData()
            if lane_idx is None:
                lane_idx = 0
            lane_idx = int(lane_idx)
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
    """Density snapshots at evenly-spaced times.

    On multilane results a lane selector appears so each lane can be explored.
    """

    def __init__(self, parent=None):
        super().__init__(parent)

        # Lane selector row (shown only for multilane)
        lane_row = QHBoxLayout()
        lane_row.addWidget(QLabel("Lane:"))
        self.lane_combo = QComboBox()
        lane_row.addWidget(self.lane_combo)
        lane_row.addStretch(1)
        self._lane_row = lane_row
        self._layout.insertLayout(0, lane_row)
        self._set_lane_row_visible(False)

        # Snapshot count slider
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
        self.lane_combo.currentIndexChanged.connect(self._refresh)

    def _set_lane_row_visible(self, visible: bool):
        for i in range(self._lane_row.count()):
            item = self._lane_row.itemAt(i)
            if item and item.widget():
                item.widget().setVisible(visible)

    def set_data(self, data: dict):
        self._data = data
        n_lanes = int(data.get("n_lanes", 1))
        if n_lanes > 1:
            self.lane_combo.blockSignals(True)
            self.lane_combo.clear()
            for i in range(n_lanes):
                self.lane_combo.addItem(f"Lane {i + 1}", i)
            self.lane_combo.blockSignals(False)
            self._set_lane_row_visible(True)
        else:
            self._set_lane_row_visible(False)
        self._refresh()

    def _on_slide(self, value: int):
        self.label.setText(str(value))
        self._refresh()

    def _refresh(self):
        if self._data is None:
            return
        _reset_axes(self)
        n_lanes = int(self._data.get("n_lanes", 1))
        if n_lanes > 1:
            lane_idx = self.lane_combo.currentData()
            if lane_idx is None:
                lane_idx = 0
            lane_idx = int(lane_idx)
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
_MULTILANE_TAB_INDICES = (3, 4)


class PDEPlotPanel(QTabWidget):
    """Five-tab plot panel for PDE simulation results.

    Tab order: Space-time (0), Snapshots (1), Boundary flow (2),
               Lane densities (3, multilane only), Mass conservation (4, multilane only).

    Call ``show_result(data)`` with the dict returned by ``load_pde_netcdf``.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self.tab_spacetime = _SpacetimeTab()
        self.tab_snapshots = _SnapshotsTab()
        self.tab_flow = _FlowTab()
        self.tab_lane_dens = _LaneDensitiesTab()
        self.tab_mass = _MassTab()
        self.addTab(self.tab_spacetime, "Space-time")
        self.addTab(self.tab_snapshots, "Snapshots")
        self.addTab(self.tab_flow, "Boundary flow")
        self.addTab(self.tab_lane_dens, "Lane densities")
        self.addTab(self.tab_mass, "Mass conservation")
        for i in _MULTILANE_TAB_INDICES:
            self.setTabEnabled(i, False)

        # Stale-result status indicator (top-right corner of the tab bar).
        from PyQt5.QtCore import Qt
        from PyQt5.QtWidgets import QLabel, QWidget, QHBoxLayout
        holder = QWidget()
        hbox = QHBoxLayout(holder)
        hbox.setContentsMargins(8, 0, 14, 4)
        self._status = QLabel()
        hbox.addWidget(self._status)
        self.setCornerWidget(holder, Qt.TopRightCorner)
        self.set_status(None)

    def set_status(self, kind, summary: str = ""):
        """Update the corner indicator with a minimal label.  ``summary`` accepted for API parity."""
        if kind == "fresh":
            color, text = "#4a8a6a", "● up to date"
        elif kind == "stale":
            color, text = "#b08040", "● stale"
        elif kind == "running":
            color, text = "#5a7fad", "● running…"
        else:
            color, text = "#9a9a9a", "● no results"
        self._status.setText(text)
        self._status.setStyleSheet(
            f"color: {color}; font-size: 10px; font-weight: 500;"
        )

    def set_multilane_enabled(self, enabled: bool):
        """Enable/disable the two multi-lane tabs without touching plot data."""
        for i in _MULTILANE_TAB_INDICES:
            self.setTabEnabled(i, enabled)

    def show_result(self, data: dict):
        self.tab_spacetime.set_data(data)
        self.tab_snapshots.set_data(data)
        self.tab_flow.set_data(data)

        multi = int(data.get("n_lanes", 1)) > 1
        self.set_multilane_enabled(multi)
        if multi:
            self.tab_lane_dens.set_data(data)
            self.tab_mass.set_data(data)
