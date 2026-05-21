from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.collections import LineCollection


def plot_spacetime(data, ax=None, title=None):
    """Space-time diagram: sites on y-axis, time on x-axis, occupied=black."""
    history = data['history']   # shape (L, n_steps)
    L, n_steps = history.shape

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4))
    else:
        fig = ax.figure

    # history is (site, time); imshow expects (row, col) = (y, x)
    ax.imshow(
        history,
        aspect='auto',
        origin='lower',
        cmap='binary',
        interpolation='nearest',
        extent=[1, n_steps, 1, L],
    )
    ax.set_xlabel('Time step')
    ax.set_ylabel('Site')
    ax.set_title(title or f'Space-time diagram  (α={data["alpha"]}, β={data["beta"]})')
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    return fig, ax


def plot_density(data, ax=None, title=None):
    """Mean density ρ = N/L vs time step."""
    density = data['density']
    n_steps = len(density)
    steps = np.arange(1, n_steps + 1)

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 3))
    else:
        fig = ax.figure

    ax.plot(steps, density, color='steelblue', linewidth=1.2)
    # Theoretical steady-state density for min-current phase: 0.5 when alpha=beta=0.5
    rho_mean = float(np.mean(density))
    ax.axhline(rho_mean, color='tomato', linestyle='--', linewidth=1,
               label=f'time-mean ρ = {rho_mean:.3f}')
    ax.set_xlabel('Time step')
    ax.set_ylabel('Density ρ')
    ax.set_ylim(0, 1)
    ax.set_title(title or f'Density vs time  (α={data["alpha"]}, β={data["beta"]})')
    ax.legend(fontsize=9)
    return fig, ax


def plot_current(data, ax=None, title=None):
    """Particle exits (current) vs time step."""
    current = data['current']
    n_steps = len(current)
    steps = np.arange(1, n_steps + 1)

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 3))
    else:
        fig = ax.figure

    ax.bar(steps, current, color='steelblue', width=0.8)
    j_mean = float(np.mean(current))
    ax.axhline(j_mean, color='tomato', linestyle='--', linewidth=1,
               label=f'mean J = {j_mean:.3f}')
    ax.set_xlabel('Time step')
    ax.set_ylabel('Exits per step')
    ax.set_title(title or f'Current vs time  (α={data["alpha"]}, β={data["beta"]})')
    ax.legend(fontsize=9)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    return fig, ax


def plot_fundamental_diagram(rho_vals, J_vals, ax=None, title=None):
    """Flow J vs density rho with the theoretical J = rho(1-rho) parabola."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))
    else:
        fig = ax.figure

    rho_theory = np.linspace(0, 1, 300)
    J_theory = np.minimum(rho_theory, 1 - rho_theory)

    ax.plot(rho_theory, J_theory, 'k--',
        linewidth=1.5, label='J = min(ρ, 1−ρ) theory')
    ax.scatter(rho_vals, J_vals, s=18, color='steelblue',
               alpha=0.8, zorder=3, label='simulation')
    ax.set_xlabel('Density ρ')
    ax.set_ylabel('Current J')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 0.8)
    ax.set_title(title or 'Fundamental diagram')
    ax.legend(fontsize=9)
    return fig, ax


def plot_summary(data, save_path=None):
    """Three-panel summary figure: space-time, density, current."""
    fig = plt.figure(figsize=(12, 8))
    gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.35)

    ax_st = fig.add_subplot(gs[0, :])
    ax_rho = fig.add_subplot(gs[1, 0])
    ax_j = fig.add_subplot(gs[1, 1])

    plot_spacetime(data, ax=ax_st)
    plot_density(data, ax=ax_rho)
    plot_current(data, ax=ax_j)

    fig.suptitle(
        f'1D TASEP  L={data["L"]}  n_steps={data["n_steps"]}  '
        f'α={data["alpha"]}  β={data["beta"]}',
        fontsize=12,
    )

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f'Saved to {save_path}')

    return fig


# ---------------------------------------------------------------------------
# Network-aware plots (operate on either a NetworkResult or a NetworkSpec)
# ---------------------------------------------------------------------------

def _config_layout(result):
    """Pull (junctions, road_endpoints) dicts out of result.config['layout']."""
    layout = (result.config or {}).get('layout', {})
    j = {int(d['id']): (float(d['x']), float(d['y']))
         for d in layout.get('junctions', [])}
    r = {int(d['id']): (tuple(d['end_1']), tuple(d['end_2']))
         for d in layout.get('roads', [])}
    return j, r


def _draw_network(ax, road_endpoints, junction_xy,
                  road_lane_info=None,
                  road_colours=None, road_label=None,
                  open_endpoints=(), junction_annotations=None,
                  title=None):
    """Geometry helper used by both spec-preview and result-heatmap callers.

    Args:
        road_endpoints: dict[rid -> ((x1, y1), (x2, y2))]
        junction_xy:    dict[jid -> (x, y)]
        road_lane_info: optional dict[rid -> list[{"flow_direction": int, "colour": float|None}]].
                        If given, each road is drawn as N parallel offset lines
                        (one per lane) with an arrow indicating each lane's
                        flow direction.  Per-lane colour (0–1) is used for the
                        heatmap overlay.
        road_colours:   legacy per-road colour dict; ignored if road_lane_info supplies colours.
        road_label:     callable rid -> str, OR None to default to 'R{rid}'
        open_endpoints: iterable of (rid, end_index) (end_index 1 or 2) to mark with red squares
        junction_annotations: optional dict[jid -> str], placed slightly above each junction
        title:          axes title

    Zoom-aware labels: road labels, α/β annotations and junction-routing
    annotations all start hidden and become visible once the axes are zoomed
    in to roughly 70% of the initial extent (matplotlib navigation toolbar
    triggers the xlim_changed callback).

    Returns: (fig, ax)
    """
    fig = ax.figure
    rid_list = sorted(road_endpoints.keys())

    # Auto spacing between parallel lanes is a fraction of the diagram extent.
    all_pts = [pt for ends in road_endpoints.values() for pt in ends] + list(junction_xy.values())
    if all_pts:
        xs = [p[0] for p in all_pts]; ys = [p[1] for p in all_pts]
        span = max(max(xs) - min(xs), max(ys) - min(ys), 1.0)
    else:
        span = 1.0
    lane_spacing = 0.022 * span

    # ----- Road geometry: either multi-lane (offset parallels + arrows)
    #       or legacy single-line per road. -----
    if road_lane_info is not None:
        _draw_roads_multilane(ax, road_endpoints, road_lane_info, lane_spacing)
    else:
        segs = [list(road_endpoints[rid]) for rid in rid_list]
        if road_colours is not None:
            colours = np.array([road_colours.get(rid, 0.0) for rid in rid_list])
            lc = LineCollection(segs, cmap='viridis', linewidths=4)
            lc.set_array(colours)
            lc.set_clim(0, 1)
            ax.add_collection(lc)
            fig.colorbar(lc, ax=ax, label='occupancy', shrink=0.7)
        else:
            ax.add_collection(LineCollection(segs, colors='steelblue', linewidths=3))

    # ----- Road midpoint labels (zoom-aware). -----
    label_fn = road_label or (lambda rid: f'R{rid}')
    road_label_artists = []
    for rid in rid_list:
        e1, e2 = road_endpoints[rid]
        mx, my = 0.5 * (e1[0] + e2[0]), 0.5 * (e1[1] + e2[1])
        # Offset labels slightly perpendicular so they don't sit on the road.
        dx, dy = e2[0] - e1[0], e2[1] - e1[1]
        seg_len = max(1e-9, (dx * dx + dy * dy) ** 0.5)
        nx, ny = -dy / seg_len, dx / seg_len
        off = 0.03 * span
        t = ax.text(mx + off * nx, my + off * ny, label_fn(rid),
                    fontsize=7, color='dimgray',
                    ha='center', va='center', backgroundcolor='white')
        t.set_visible(False)
        road_label_artists.append(t)

    # ----- Junction nodes + zoom-aware routing annotations. -----
    junction_label_artists = []
    if junction_xy:
        ax.scatter([p[0] for p in junction_xy.values()],
                   [p[1] for p in junction_xy.values()],
                   s=160, c='black', zorder=3)
        for jid, (x, y) in junction_xy.items():
            ax.text(x, y, f'J{jid}', color='white', fontsize=8,
                    ha='center', va='center', zorder=4, fontweight='bold')
            if junction_annotations and jid in junction_annotations:
                ann = ax.text(x, y + 0.05 * span, junction_annotations[jid],
                              fontsize=6, color='darkblue', ha='center', va='bottom',
                              backgroundcolor='white', zorder=4)
                ann.set_visible(False)
                junction_label_artists.append(ann)

    # ----- Open-boundary markers. -----
    for rid, end in open_endpoints:
        if rid not in road_endpoints:
            continue
        e = road_endpoints[rid][end - 1]
        ax.plot(*e, marker='s', color='tomato', markersize=8, zorder=5)

    # ----- Extent. -----
    if all_pts:
        pad = 0.3 * span
        ax.set_xlim(min(xs) - pad, max(xs) + pad)
        ax.set_ylim(min(ys) - pad, max(ys) + pad)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    if title is not None:
        ax.set_title(title)

    # ----- Zoom-aware visibility. -----
    _install_zoom_label_visibility(
        ax,
        road_label_artists + junction_label_artists,
    )

    return fig, ax


def _draw_roads_multilane(ax, road_endpoints, road_lane_info, lane_spacing):
    """Render each road as N parallel offset lines with directional arrows."""
    fig = ax.figure

    segs_all = []
    colour_vals = []
    has_colours = False
    for rid, (p1, p2) in road_endpoints.items():
        lanes = road_lane_info.get(rid)
        if not lanes:
            lanes = [{"flow_direction": 1, "colour": None}]
        p1a = np.array(p1, dtype=float)
        p2a = np.array(p2, dtype=float)
        dvec = p2a - p1a
        length = float(np.hypot(*dvec))
        if length < 1e-9:
            continue
        perp = np.array([-dvec[1], dvec[0]]) / length
        n = len(lanes)
        offsets = (np.arange(n) - (n - 1) / 2.0) * lane_spacing
        for lane_info, off in zip(lanes, offsets):
            fd = int(lane_info.get("flow_direction", 1))
            col = lane_info.get("colour")
            lp1 = p1a + off * perp
            lp2 = p2a + off * perp
            # Make the segment direction match the lane's flow.
            if fd < 0:
                lp1, lp2 = lp2, lp1
            segs_all.append([lp1.tolist(), lp2.tolist()])
            if col is not None:
                colour_vals.append(float(col))
                has_colours = True
            else:
                colour_vals.append(0.0)

            # Direction arrow at the midpoint.
            mid = 0.5 * (lp1 + lp2)
            arrow_half = 0.04 * length
            tip  = mid + (arrow_half / length) * (lp2 - lp1)
            base = mid - (arrow_half / length) * (lp2 - lp1)
            ax.annotate(
                '', xy=tip, xytext=base,
                arrowprops=dict(arrowstyle='->', color='dimgray',
                                lw=1.4, shrinkA=0, shrinkB=0),
                zorder=2,
            )

    if has_colours:
        lc = LineCollection(segs_all, cmap='viridis', linewidths=3)
        lc.set_array(np.array(colour_vals))
        lc.set_clim(0, 1)
        ax.add_collection(lc)
        fig.colorbar(lc, ax=ax, label='occupancy', shrink=0.7)
    else:
        ax.add_collection(LineCollection(segs_all, colors='steelblue', linewidths=2.4))


def _install_zoom_label_visibility(ax, label_artists, ratio: float = 0.70):
    """Show ``label_artists`` only when the axes are zoomed in below ``ratio``
    times the initial span.  Wired via matplotlib's xlim/ylim callbacks.
    """
    if not label_artists:
        return
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    initial_span = max(abs(x1 - x0), abs(y1 - y0), 1e-9)

    def _refresh(_evt_ax=None):
        cx0, cx1 = ax.get_xlim()
        cy0, cy1 = ax.get_ylim()
        cur = max(abs(cx1 - cx0), abs(cy1 - cy0), 1e-9)
        show = (cur / initial_span) < ratio
        for art in label_artists:
            art.set_visible(show)
        ax.figure.canvas.draw_idle()

    ax.callbacks.connect('xlim_changed', _refresh)
    ax.callbacks.connect('ylim_changed', _refresh)


def _result_road_colours(result, t):
    """Per-road mean occupancy at timestep t, keyed by road_id."""
    out = {}
    lane_road = result.lane_road_id
    lane_len  = result.lane_length
    occ_t = result.occupancy[t]
    for rid in np.unique(lane_road):
        lanes_of_r = np.where(lane_road == rid)[0]
        total_cells = sum(int(lane_len[li]) for li in lanes_of_r)
        total_occ = sum(int(occ_t[li, :lane_len[li]].sum()) for li in lanes_of_r)
        out[int(rid)] = total_occ / max(1, total_cells)
    return out


def _config_road_flow_directions(result):
    """dict[rid -> list[int]] of flow_direction per lane from result.config['roads']."""
    out = {}
    cfg_roads = ((result.config or {}).get('roads', []) or [])
    for rd in cfg_roads:
        try:
            rid = int(rd['id'])
            out[rid] = [int(ln.get('flow_direction', 1)) for ln in rd.get('lanes', [])]
        except (KeyError, TypeError, ValueError):
            continue
    return out


def _result_lane_info(result, t):
    """Per-road list of lane info dicts: {flow_direction, colour}.

    flow_direction comes from ``result.config['roads'][*]['lanes']``.
    colour is per-lane occupancy at timestep ``t`` (0–1), or ``None`` if
    ``t is None`` (preview mode).
    """
    out = {}
    lane_road = np.asarray(result.lane_road_id)
    lane_len  = np.asarray(result.lane_length)
    fd_map = _config_road_flow_directions(result)
    occ_t = result.occupancy[t] if t is not None else None

    for rid in np.unique(lane_road):
        rid_int = int(rid)
        lanes_of_r = np.where(lane_road == rid)[0]
        fds = fd_map.get(rid_int, [1] * len(lanes_of_r))
        info_list = []
        for k, li in enumerate(lanes_of_r):
            fd = fds[k] if k < len(fds) else 1
            if occ_t is not None:
                cells = int(lane_len[li])
                occ = int(occ_t[li, :cells].sum())
                col = occ / max(1, cells)
            else:
                col = None
            info_list.append({"flow_direction": fd, "colour": col})
        out[rid_int] = info_list
    return out


def _spec_lane_info(spec):
    """Per-road lane info from a NetworkSpec (no occupancy data)."""
    return {
        r.id: [{"flow_direction": ln.flow_direction, "colour": None}
               for ln in r.lanes]
        for r in spec.roads
    }


def _result_open_endpoints(result):
    """List of (rid, end) tuples (end 1 or 2) where the road has an open boundary."""
    open_ids = set(int(r) for r in result.lane_road_id[result.lane_open_in.astype(bool)])
    open_ids |= set(int(r) for r in result.lane_road_id[result.lane_open_out.astype(bool)])
    out = []
    for rid in sorted(open_ids):
        if rid - 1 >= len(result.road_end_junction):
            continue
        ej = result.road_end_junction[rid - 1]
        if ej[0] == 0:
            out.append((rid, 1))
        if ej[1] == 0:
            out.append((rid, 2))
    return out


def plot_network_layout(result, ax=None, occupancy_t: Optional[int] = None,
                        title: Optional[str] = None):
    """Draw a NetworkResult's network: roads, junctions, open boundaries.

    Each road is rendered as N parallel lines (one per lane) with directional
    arrows.  If ``occupancy_t`` is given, each lane is shaded individually by
    its own occupancy at that timestep.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 6))
    j_xy, r_ends = _config_layout(result)
    return _draw_network(
        ax, road_endpoints=r_ends, junction_xy=j_xy,
        road_lane_info=_result_lane_info(result, occupancy_t),
        open_endpoints=_result_open_endpoints(result),
        title=title or f'Network layout  (t={occupancy_t if occupancy_t is not None else "—"})',
    )


def plot_network_spec(spec, layout, ax=None, alpha_beta_labels: bool = False,
                      title: Optional[str] = None):
    """Draw a NetworkSpec preview (no occupancy overlay).

    Used by the GUI before a run completes.  Optionally annotates each road
    with α/β values from its open lanes, and each 4-way junction with its
    L/S/R turn probabilities.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 5))

    # Open endpoints: roads with end_junction == 0 and at least one open lane on that end.
    open_eps = []
    for r in spec.roads:
        has_open = any(ln.open_in or ln.open_out for ln in r.lanes)
        if not has_open:
            continue
        if r.end_junction[0] == 0:
            open_eps.append((r.id, 1))
        if r.end_junction[1] == 0:
            open_eps.append((r.id, 2))

    # Road labels: append α/β if requested and the road has open boundaries.
    road_by_id = {r.id: r for r in spec.roads}
    def _label(rid):
        if not alpha_beta_labels:
            return f'R{rid}'
        r = road_by_id[rid]
        a_set = any(ln.open_in for ln in r.lanes)
        b_set = any(ln.open_out for ln in r.lanes)
        if not (a_set or b_set):
            return f'R{rid}'
        bits = []
        if a_set:
            a = next(ln.alpha for ln in r.lanes if ln.open_in)
            bits.append(f'α={a:.2f}')
        if b_set:
            b = next(ln.beta for ln in r.lanes if ln.open_out)
            bits.append(f'β={b:.2f}')
        return f'R{rid}\n' + ' '.join(bits)

    # Junction annotations: show routing summary for every junction with
    # 4-way routes (crossroads and town).  Visibility is zoom-controlled so
    # they don't pile up when zoomed out.
    annotations = {}
    if alpha_beta_labels:
        for j in spec.junctions:
            if j.routes and j.n_in == 4 and len(j.routes[0]) == 4:
                row = j.routes[0]
                annotations[j.id] = f"L={row[1]:.2f} S={row[2]:.2f} R={row[3]:.2f}"

    return _draw_network(
        ax,
        road_endpoints=layout.road_endpoints,
        junction_xy=layout.junctions,
        road_lane_info=_spec_lane_info(spec),
        road_label=_label,
        open_endpoints=open_eps,
        junction_annotations=annotations,
        title=title or 'Network preview',
    )


def plot_network_density(result, ax=None, title=None, road_ids=None):
    """Per-road density vs time.  One line per road.

    Parameters
    ----------
    road_ids : list[int] | None
        0-based road indices to plot.  ``None`` plots all roads.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 4))
    else:
        fig = ax.figure
    n_steps, n_roads = result.road_density.shape
    t = np.arange(1, n_steps + 1)
    indices = road_ids if road_ids is not None else list(range(n_roads))
    for r in indices:
        ax.plot(t, result.road_density[:, r], label=f'R{r + 1}', linewidth=1.0)
    ax.set_xlabel('Time step')
    ax.set_ylabel('Density ρ')
    ax.set_ylim(0, 1)
    ax.set_title(title or 'Per-road density vs time')
    n_shown = len(indices)
    if n_shown <= 12:
        ax.legend(fontsize=8, ncol=min(4, n_shown), loc='upper right')
    return fig, ax


def plot_network_currents(result, ax=None, title=None, road_ids=None):
    """Per-road flow J(t) as an MCMC-style trace — one line per road.

    Matches the visual style of ``plot_network_density``: a noisy line per road
    vs time-step, y-axis bounded between 0 and 1.  J(t) is the exits-per-step
    count divided by the number of lanes of that road, so it is comparable
    across roads and bounded above by 1 (one vehicle/lane/step at most).
    Legend labels show the time-averaged J̄.

    Parameters
    ----------
    road_ids : list[int] | None
        0-based road indices to plot.  ``None`` plots all roads.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 4))
    else:
        fig = ax.figure
    n_steps, n_roads = result.road_exits.shape
    t = np.arange(1, n_steps + 1)
    indices = road_ids if road_ids is not None else list(range(n_roads))

    # Lanes per road, for normalising exits → flow per lane.
    lane_road = np.asarray(result.lane_road_id)
    lanes_per_road = {int(rid): int(np.sum(lane_road == rid))
                      for rid in np.unique(lane_road)}

    for r in indices:
        exits = result.road_exits[:, r].astype(float)
        n_lanes = max(1, lanes_per_road.get(r + 1, 1))
        flow = exits / n_lanes
        mean_j = float(flow.mean())
        ax.plot(t, flow, linewidth=1.0, alpha=0.85,
                label=f'R{r + 1}  (J̄={mean_j:.3f})')
    ax.set_xlabel('Time step')
    ax.set_ylabel('Flow J (exits/lane/step)')
    ax.set_ylim(0, 1)
    ax.set_title(title or 'Per-road flow vs time')
    if len(indices) <= 12:
        ax.legend(fontsize=8, ncol=min(4, max(1, len(indices))), loc='upper right')
    return fig, ax


def plot_network_spacetime(result, road_id: int, ax=None, title=None,
                           lane_index: int = 0):
    """Space-time diagram for a single lane of a single road.

    Parameters
    ----------
    lane_index : int
        Which lane of the road to display (0-based index into the lanes
        belonging to this road).  Ignored when the road has only one lane.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4))
    else:
        fig = ax.figure

    lane_idxs = np.where(result.lane_road_id == road_id)[0]
    if len(lane_idxs) == 0:
        ax.set_title(f'road {road_id}: not found')
        return fig, ax

    lane_index = max(0, min(int(lane_index), len(lane_idxs) - 1))
    li  = lane_idxs[lane_index]
    Lk  = int(result.lane_length[li])
    n_steps = result.occupancy.shape[0]

    occ = result.occupancy[:, li, :Lk].T   # (cells, time)

    ax.imshow(occ, aspect='auto', origin='lower', cmap='binary',
              interpolation='nearest', vmin=0, vmax=1,
              extent=[1, n_steps, 0, Lk])
    ax.set_xlabel('Time step')
    ax.set_ylabel('Cell')
    lane_label = f' — Lane {lane_index + 1}' if len(lane_idxs) > 1 else ''
    ax.set_title(title or f'Road {road_id}{lane_label} space-time')
    return fig, ax


def plot_network_summary(result, save_path: Optional[str] = None):
    """Four-panel summary: layout, density, currents, space-time of road 1."""
    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)

    ax_layout = fig.add_subplot(gs[0, 0])
    ax_dens   = fig.add_subplot(gs[0, 1])
    ax_cur    = fig.add_subplot(gs[1, 0])
    ax_st     = fig.add_subplot(gs[1, 1])

    plot_network_layout(result, ax=ax_layout, occupancy_t=result.occupancy.shape[0] - 1,
                        title='Final occupancy')
    plot_network_density(result, ax=ax_dens)
    plot_network_currents(result, ax=ax_cur)
    plot_network_spacetime(result, road_id=int(result.lane_road_id[0]), ax=ax_st)

    fig.suptitle(
        f'Network simulation  n_steps={result.n_steps}  '
        f'roads={result.road_density.shape[1]}  '
        f'lanes={result.occupancy.shape[1]}',
        fontsize=12,
    )

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f'Saved to {save_path}')
    return fig
