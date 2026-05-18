"""Matplotlib-in-Qt canvas that draws a network spec preview.

Pure delegation to visualisation.plot_network_spec — the spec preview shares
its geometry-drawing code with the post-run heatmap (plot_network_layout).
"""

from __future__ import annotations

from typing import Optional

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from python.road_network import NetworkSpec, LayoutSpec
from python.visualisation import plot_network_spec


class NetworkWidget(FigureCanvasQTAgg):
    def __init__(self, parent=None):
        fig = Figure(figsize=(5, 5))
        super().__init__(fig)
        self.setParent(parent)
        self._ax = fig.add_subplot(111)
        self._spec: Optional[NetworkSpec] = None
        self._layout: Optional[LayoutSpec] = None

    def set_network(self, spec: NetworkSpec, layout: LayoutSpec):
        self._spec = spec
        self._layout = layout
        self._ax.clear()
        plot_network_spec(spec, layout, ax=self._ax, alpha_beta_labels=True)
        self.draw_idle()
