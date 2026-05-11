"""Matplotlib-in-Qt canvas that draws a network spec preview."""

from __future__ import annotations

from typing import Optional, Tuple

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.collections import LineCollection
from matplotlib.figure import Figure

from python.road_network import NetworkSpec, LayoutSpec


class NetworkWidget(FigureCanvasQTAgg):
    """Renders a NetworkSpec + LayoutSpec.  Same visuals as
    visualisation.plot_network_layout but driven by the spec directly
    (not from a NetworkResult)."""

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
        self._redraw()

    def _redraw(self):
        ax = self._ax
        ax.clear()
        if self._spec is None or self._layout is None:
            self.draw_idle()
            return

        spec, layout = self._spec, self._layout

        # Roads
        rid_list = sorted(layout.road_endpoints.keys())
        segs = [list(layout.road_endpoints[rid]) for rid in rid_list]
        lc = LineCollection(segs, colors='steelblue', linewidths=3)
        ax.add_collection(lc)

        # Road labels at midpoint, with a small annotation for open boundaries.
        road_by_id = {r.id: r for r in spec.roads}
        for rid in rid_list:
            e1, e2 = layout.road_endpoints[rid]
            mx, my = 0.5 * (e1[0] + e2[0]), 0.5 * (e1[1] + e2[1])
            label = f"R{rid}"
            r = road_by_id[rid]
            # Annotate alpha/beta if the road has open boundaries.
            alpha_set = any(ln.open_in for ln in r.lanes)
            beta_set  = any(ln.open_out for ln in r.lanes)
            if alpha_set or beta_set:
                a = next((ln.alpha for ln in r.lanes if ln.open_in), 0.0)
                b = next((ln.beta  for ln in r.lanes if ln.open_out), 0.0)
                bits = []
                if alpha_set: bits.append(f"α={a:.2f}")
                if beta_set:  bits.append(f"β={b:.2f}")
                label = f"R{rid}\n" + " ".join(bits)
            ax.text(mx, my, label, fontsize=7, color='dimgray',
                    ha='center', va='center', backgroundcolor='white')

        # Junctions
        if layout.junctions:
            jx = [p[0] for p in layout.junctions.values()]
            jy = [p[1] for p in layout.junctions.values()]
            ax.scatter(jx, jy, s=200, c='black', zorder=3)
            for jid, (x, y) in layout.junctions.items():
                ax.text(x, y, f"J{jid}", color='white', fontsize=8,
                        ha='center', va='center', zorder=4, fontweight='bold')
            # If this is a 4-way crossroad, show turning probabilities as a
            # small label slightly above each junction.
            for j in spec.junctions:
                if j.n_in == 4 and j.routes and j.id in layout.junctions:
                    x, y = layout.junctions[j.id]
                    row = j.routes[0]
                    p_left  = row[1 % 4]
                    p_right = row[3 % 4]
                    p_straight = row[2 % 4]
                    ax.text(x, y + 0.1, f"L={p_left:.2f} S={p_straight:.2f} R={p_right:.2f}",
                            fontsize=6, color='darkblue', ha='center', va='bottom')

        # Open-boundary squares
        for r in spec.roads:
            e1, e2 = layout.road_endpoints[r.id]
            if r.end_junction[0] == 0 and any(ln.open_in or ln.open_out for ln in r.lanes):
                ax.plot(*e1, marker='s', color='tomato', markersize=8, zorder=5)
            if r.end_junction[1] == 0 and any(ln.open_in or ln.open_out for ln in r.lanes):
                ax.plot(*e2, marker='s', color='tomato', markersize=8, zorder=5)

        # Sizing.
        all_pts = [pt for ends in layout.road_endpoints.values() for pt in ends] + list(layout.junctions.values())
        if all_pts:
            xs = [p[0] for p in all_pts]; ys = [p[1] for p in all_pts]
            span = max(max(xs) - min(xs), max(ys) - min(ys), 1.0)
            pad = 0.3 * span
            ax.set_xlim(min(xs) - pad, max(xs) + pad)
            ax.set_ylim(min(ys) - pad, max(ys) + pad)
        ax.set_aspect('equal')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("Network preview")
        self.draw_idle()
