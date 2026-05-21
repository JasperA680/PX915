"""Tabbed panel that shows density, currents, space-time, and network heatmap."""

from __future__ import annotations

from typing import Optional

from PyQt5.QtCore import Qt, pyqtSignal, QEvent
from PyQt5.QtGui import QStandardItem, QStandardItemModel
from PyQt5.QtWidgets import (
    QTabWidget, QWidget, QVBoxLayout, QHBoxLayout,
    QComboBox, QSlider, QLabel, QListView, QFrame,
)
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from python.io import NetworkResult
from python.visualisation import (
    plot_network_density,
    plot_network_spacetime,
    plot_network_layout,
    plot_fundamental_diagram,
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


class _FundamentalDiagramTab(_CanvasTab):
    """Plot of J vs ρ from a parameter sweep (TASEP or NS).

    A single run is one point in J–ρ space; the full diagram requires a sweep:
    TASEP uses an α/β sweep on an open chain, NS uses a density sweep on a
    periodic ring.  Use the "Run FD sweep" button — the Model combo decides
    which sweep runs.  Plots use ``plot_fundamental_diagram`` (also used by
    scripts/run_fundamental_diagram.py).
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self._show_placeholder()

    def _show_placeholder(self):
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.text(
            0.5, 0.5,
            "Fundamental diagram\n\n"
            "A J–ρ diagram requires a parameter sweep, not a single run.\n\n"
            "Click 'Run FD sweep' on the toolbar — the active Model decides\n"
            "which sweep is performed:\n"
            "  • NS:    periodic-ring density sweep (uses v_max + p_slow)\n"
            "  • TASEP: open-chain α/β phase-diagram sweep",
            ha="center", va="center",
            transform=ax.transAxes,
            fontsize=11, color="gray",
        )
        ax.set_axis_off()
        self.ax = ax
        self.canvas.draw_idle()

    def set_sweep(self, rho, J, title: str = None,
                  model: str = 'TASEP', v_max: int = 1):
        """Populate the FD scatter from a parameter sweep."""
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        plot_fundamental_diagram(rho, J, ax=ax, title=title,
                                 model=model, v_max=v_max)
        self.ax = ax
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
    # Index of the Fundamental Diagram tab; disabled until a sweep populates it.
    FD_TAB_INDEX = 3

    def __init__(self, parent=None):
        super().__init__(parent)
        self.tab_density = _DensityTab()
        self.tab_spacetime = _SpacetimeTab()
        self.tab_heatmap = _HeatmapTab()
        self.tab_fd = _FundamentalDiagramTab()
        self.addTab(self.tab_density,   "Density")
        self.addTab(self.tab_spacetime, "Space-time")
        self.addTab(self.tab_heatmap,   "Network heatmap")
        self.addTab(self.tab_fd,        "Fundamental diagram")
        # FD tab disabled until a sweep populates it.
        self.setTabEnabled(self.FD_TAB_INDEX, False)

        # Stale-result status indicator (top-right corner of the tab bar).
        # Wrap in a holder so we can give it real padding without the QTabBar
        # eating it.  Padding pushes the label up + away from the tabs so it
        # doesn't overlap the tab text.
        holder = QWidget()
        hbox = QHBoxLayout(holder)
        hbox.setContentsMargins(8, 0, 14, 4)
        self._status = QLabel()
        hbox.addWidget(self._status)
        self.setCornerWidget(holder, Qt.TopRightCorner)
        self.set_status(None)   # initial "no results yet"

    # ------------------------------------------------------------------
    # Stale / fresh status indicator
    # ------------------------------------------------------------------

    def set_status(self, kind: str, summary: str = ""):
        """Update the corner indicator with a minimal label.

        kind: 'fresh', 'stale', 'running', or None / 'empty'.
        ``summary`` is accepted for API parity but no longer displayed.
        """
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

    def clear(self):
        self.tab_density.clear()
        self.tab_spacetime.clear()
        self.tab_heatmap.clear()

    def show_result(self, result: NetworkResult):
        self.tab_density.set_result(result)
        self.tab_spacetime.set_result(result)
        self.tab_heatmap.set_result(result)

    def set_fd_sweep(self, rho, J, title: str = None,
                     model: str = 'TASEP', v_max: int = 1):
        """Populate the FD tab with a sweep result and enable it."""
        self.tab_fd.set_sweep(rho, J, title=title, model=model, v_max=v_max)
        self.setTabEnabled(self.FD_TAB_INDEX, True)
        self.setTabText(self.FD_TAB_INDEX, "Fundamental diagram")
        self.setCurrentIndex(self.FD_TAB_INDEX)

    def set_fd_stale(self, stale: bool):
        """Annotate the FD tab title to flag its data as out-of-date."""
        if not self.isTabEnabled(self.FD_TAB_INDEX):
            return  # no data yet → no point flagging staleness
        label = "Fundamental diagram"
        if stale:
            label = "● Fundamental diagram (stale)"
        self.setTabText(self.FD_TAB_INDEX, label)

    def reset_fd(self):
        """Restore the FD tab placeholder and disable it (e.g. on preset switch)."""
        self.tab_fd._show_placeholder()
        self.setTabEnabled(self.FD_TAB_INDEX, False)
        self.setTabText(self.FD_TAB_INDEX, "Fundamental diagram")
