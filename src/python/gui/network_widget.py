"""Matplotlib-in-Qt canvas that draws a network spec preview.

Pure delegation to visualisation.plot_network_spec — the spec preview shares
its geometry-drawing code with the post-run heatmap (plot_network_layout).

The widget wraps both a NavigationToolbar (zoom/pan) and the canvas in a
QWidget so callers can embed it like any other widget.  Clicking near a road
emits ``road_clicked(road_id: int)`` so CATab can open an editor dialog.
"""

from __future__ import annotations

import math
from typing import Optional

from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import QVBoxLayout, QWidget
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg,
    NavigationToolbar2QT,
)
from matplotlib.figure import Figure

from python.road_network import NetworkSpec, LayoutSpec
from python.visualisation import plot_network_spec


def _point_to_segment_dist(px, py, x1, y1, x2, y2) -> float:
    """Euclidean distance from point (px,py) to segment (x1,y1)–(x2,y2)."""
    dx, dy = x2 - x1, y2 - y1
    seg_len_sq = dx * dx + dy * dy
    if seg_len_sq < 1e-12:
        return math.hypot(px - x1, py - y1)
    t = max(0.0, min(1.0, ((px - x1) * dx + (py - y1) * dy) / seg_len_sq))
    return math.hypot(px - (x1 + t * dx), py - (y1 + t * dy))


class NetworkWidget(QWidget):
    """Network preview canvas with zoom/pan toolbar, road/junction click support."""

    road_clicked     = pyqtSignal(int)   # road_id when user clicks near a road
    junction_clicked = pyqtSignal(int)   # junction_id when user clicks near a junction

    def __init__(self, parent=None):
        super().__init__(parent)
        fig = Figure(figsize=(5, 5))
        self.canvas = FigureCanvasQTAgg(fig)
        self._ax = fig.add_subplot(111)
        self._toolbar = NavigationToolbar2QT(self.canvas, self)
        self._spec: Optional[NetworkSpec] = None
        self._layout_spec: Optional[LayoutSpec] = None
        self._road_endpoints: dict = {}   # rid -> ((x1,y1), (x2,y2))
        self._junction_xy: dict = {}      # jid -> (x, y)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._toolbar)
        layout.addWidget(self.canvas, stretch=1)

        self.canvas.mpl_connect("button_press_event", self._on_canvas_click)

    def set_network(self, spec: NetworkSpec, layout: LayoutSpec):
        self._spec = spec
        self._layout_spec = layout
        self._road_endpoints = dict(layout.road_endpoints)
        self._junction_xy = dict(layout.junctions)
        self._ax.clear()
        plot_network_spec(spec, layout, ax=self._ax, alpha_beta_labels=True)
        self.canvas.draw_idle()

    def _on_canvas_click(self, event):
        """Emit road_clicked or junction_clicked, whichever is closest."""
        if event.inaxes is not self._ax:
            return
        # Only act on plain left-click (not when zoom/pan tool is active).
        if self._toolbar.mode:
            return
        if not (self._road_endpoints or self._junction_xy):
            return

        px, py = event.xdata, event.ydata
        xlim = self._ax.get_xlim()
        ylim = self._ax.get_ylim()
        diag = math.hypot(xlim[1] - xlim[0], ylim[1] - ylim[0])
        road_threshold     = diag * 0.08    # 8 % of axes diagonal
        junction_threshold = diag * 0.05    # 5 % — tighter so it doesn't
                                            # steal clicks from nearby roads

        best_rid, best_road_dist = None, float("inf")
        for rid, (p1, p2) in self._road_endpoints.items():
            d = _point_to_segment_dist(px, py, p1[0], p1[1], p2[0], p2[1])
            if d < best_road_dist:
                best_road_dist = d
                best_rid = rid

        best_jid, best_jdist = None, float("inf")
        for jid, (jx, jy) in self._junction_xy.items():
            d = math.hypot(px - jx, py - jy)
            if d < best_jdist:
                best_jdist = d
                best_jid = jid

        # If a junction is the closest hit AND within its threshold, emit that.
        if (best_jid is not None and best_jdist <= junction_threshold
                and best_jdist <= best_road_dist):
            self.junction_clicked.emit(int(best_jid))
        elif best_rid is not None and best_road_dist <= road_threshold:
            self.road_clicked.emit(int(best_rid))
