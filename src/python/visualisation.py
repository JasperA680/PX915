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
                  road_colours=None, road_label=None,
                  open_endpoints=(), junction_annotations=None,
                  title=None):
    """Geometry helper used by both spec-preview and result-heatmap callers.

    Args:
        road_endpoints: dict[rid -> ((x1, y1), (x2, y2))]
        junction_xy:    dict[jid -> (x, y)]
        road_colours:   optional dict[rid -> float in [0,1]]; draws viridis overlay + colorbar
        road_label:     callable rid -> str, OR None to default to 'R{rid}'
        open_endpoints: iterable of (rid, end_index) (end_index 1 or 2) to mark with red squares
        junction_annotations: optional dict[jid -> str], placed slightly above each junction
        title:          axes title
    Returns: (fig, ax)
    """
    fig = ax.figure
    rid_list = sorted(road_endpoints.keys())
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

    # Road midpoint labels.
    label_fn = road_label or (lambda rid: f'R{rid}')
    for rid in rid_list:
        e1, e2 = road_endpoints[rid]
        mx, my = 0.5 * (e1[0] + e2[0]), 0.5 * (e1[1] + e2[1])
        ax.text(mx, my, label_fn(rid), fontsize=7, color='dimgray',
                ha='center', va='center', backgroundcolor='white')

    # Junction nodes.
    if junction_xy:
        ax.scatter([p[0] for p in junction_xy.values()],
                   [p[1] for p in junction_xy.values()],
                   s=160, c='black', zorder=3)
        for jid, (x, y) in junction_xy.items():
            ax.text(x, y, f'J{jid}', color='white', fontsize=8,
                    ha='center', va='center', zorder=4, fontweight='bold')
            if junction_annotations and jid in junction_annotations:
                ax.text(x, y + 0.1, junction_annotations[jid],
                        fontsize=6, color='darkblue', ha='center', va='bottom')

    # Open-boundary markers.
    for rid, end in open_endpoints:
        if rid not in road_endpoints:
            continue
        e = road_endpoints[rid][end - 1]
        ax.plot(*e, marker='s', color='tomato', markersize=8, zorder=5)

    # Extent.
    all_pts = [pt for ends in road_endpoints.values() for pt in ends] + list(junction_xy.values())
    if all_pts:
        xs = [p[0] for p in all_pts]; ys = [p[1] for p in all_pts]
        pad = 0.3 * max(max(xs) - min(xs), max(ys) - min(ys), 1.0)
        ax.set_xlim(min(xs) - pad, max(xs) + pad)
        ax.set_ylim(min(ys) - pad, max(ys) + pad)
    ax.set_aspect('equal')
    ax.set_xticks([])
    ax.set_yticks([])
    if title is not None:
        ax.set_title(title)
    return fig, ax


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

    If `occupancy_t` is given, colour each road by its mean lane occupancy
    at that timestep; otherwise the roads are drawn in a uniform colour.
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(7, 6))
    j_xy, r_ends = _config_layout(result)
    return _draw_network(
        ax, road_endpoints=r_ends, junction_xy=j_xy,
        road_colours=_result_road_colours(result, occupancy_t) if occupancy_t is not None else None,
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

    # Junction annotations: only when there's a single junction.  When there
    # are several with identical defaults the labels overlap each other.
    annotations = {}
    if alpha_beta_labels and len(spec.junctions) == 1:
        j = spec.junctions[0]
        if j.n_in == 4 and j.routes:
            row = j.routes[0]
            annotations[j.id] = f"L={row[1]:.2f} S={row[2]:.2f} R={row[3]:.2f}"

    return _draw_network(
        ax,
        road_endpoints=layout.road_endpoints,
        junction_xy=layout.junctions,
        road_label=_label,
        open_endpoints=open_eps,
        junction_annotations=annotations,
        title=title or 'Network preview',
    )


def plot_network_density(result, ax=None, title=None):
    """Per-road density vs time.  One line per road."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 4))
    else:
        fig = ax.figure
    n_steps, n_roads = result.road_density.shape
    t = np.arange(1, n_steps + 1)
    for r in range(n_roads):
        ax.plot(t, result.road_density[:, r], label=f'R{r + 1}', linewidth=1.0)
    ax.set_xlabel('Time step')
    ax.set_ylabel('Density ρ')
    ax.set_ylim(0, 1)
    ax.set_title(title or 'Per-road density vs time')
    if n_roads <= 12:
        ax.legend(fontsize=8, ncol=min(4, n_roads), loc='upper right')
    return fig, ax


def plot_network_currents(result, ax=None, title=None):
    """Per-road cumulative entries and exits."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 4))
    else:
        fig = ax.figure
    n_steps, n_roads = result.road_entries.shape
    cum_in  = np.cumsum(result.road_entries, axis=0)
    cum_out = np.cumsum(result.road_exits,   axis=0)
    t = np.arange(1, n_steps + 1)
    for r in range(n_roads):
        if cum_in[-1, r] == 0 and cum_out[-1, r] == 0:
            continue
        ax.plot(t, cum_in[:, r],  linestyle='-',  linewidth=1.0, label=f'R{r + 1} in')
        ax.plot(t, cum_out[:, r], linestyle='--', linewidth=1.0, label=f'R{r + 1} out')
    ax.set_xlabel('Time step')
    ax.set_ylabel('Cumulative count')
    ax.set_title(title or 'Per-road cumulative entries (—) and exits (--)')
    if n_roads <= 8:
        ax.legend(fontsize=7, ncol=2, loc='upper left')
    return fig, ax


def plot_network_spacetime(result, road_id: int, ax=None, title=None):
    """Concatenate all lanes of a single road as a space-time diagram."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 4))
    else:
        fig = ax.figure

    lane_idxs = np.where(result.lane_road_id == road_id)[0]
    if len(lane_idxs) == 0:
        ax.set_title(f'road {road_id}: not found')
        return fig, ax

    lane_lens = result.lane_length[lane_idxs]
    n_steps   = result.occupancy.shape[0]
    rows = []
    for li, Lk in zip(lane_idxs, lane_lens):
        rows.append(result.occupancy[:, li, :Lk])
        rows.append(np.full((n_steps, 1), 0.5))   # separator stripe
    if rows:
        rows = rows[:-1]
    stacked = np.concatenate(rows, axis=1).T   # (cell-across-lanes, time)

    ax.imshow(stacked, aspect='auto', origin='lower', cmap='binary',
              interpolation='nearest', vmin=0, vmax=1,
              extent=[1, n_steps, 0, stacked.shape[0]])
    ax.set_xlabel('Time step')
    ax.set_ylabel('Cell (lanes stacked)')
    ax.set_title(title or f'Road {road_id} space-time')
    return fig, ax


def plot_speed_limit_comparison(panel_datasets, save_path=None):
    """Four-panel figure comparing shock behaviour under different speed limits.

    Parameters
    ----------
    panel_datasets:
        List of (label, datasets) tuples, one per panel (e.g. 4 panels for
        Greenshields-LF, Greenshields-Godunov, Newell-LF, Newell-Godunov).
        Each `datasets` is a list of (v_limit, data_dict) tuples where
        data_dict is the return value of load_pde_netcdf.
    save_path:
        Optional path to save the figure.

    The figure shows, for each panel:
        • Left half: fundamental diagram curves colour-coded by v_limit
        • Right half (main area): final-time density profiles, one per v_limit
    """
    from pde_runner import (
        q_of_rho, q_newell, rho_critical, rho_critical_newell, NEWELL_W,
    )

    n_panels = len(panel_datasets)
    ncols = 2
    nrows = (n_panels + 1) // 2

    fig = plt.figure(figsize=(7 * ncols, 4.5 * nrows))
    gs_outer = fig.add_gridspec(nrows, ncols, hspace=0.52, wspace=0.38)

    for panel_idx, (label, datasets) in enumerate(panel_datasets):
        row = panel_idx // ncols
        col = panel_idx % ncols
        gs_inner = gs_outer[row, col].subgridspec(1, 2, wspace=0.45, width_ratios=[1, 1.6])
        ax_fd   = fig.add_subplot(gs_inner[0])
        ax_snap = fig.add_subplot(gs_inner[1])

        # attrs from the first dataset (shared params)
        attrs   = datasets[0][1]['attrs']
        v_max   = float(attrs.get('v_max',   1.0))
        rho_max = float(attrs.get('rho_max', 1.0))
        closure = str(attrs.get('flux_type', ''))
        is_newell = ('newell' in closure)

        n_vl   = len(datasets)
        colors = plt.cm.plasma(np.linspace(0.15, 0.85, n_vl))

        rho_plot = np.linspace(0, rho_max, 400)

        for (v_limit, data), color in zip(datasets, colors):
            vl_label = f'$v_{{lim}}={v_limit:.2f}$'

            # --- fundamental diagram ---
            if is_newell:
                q_plot = q_newell(rho_plot, rho_max, v_limit)
                rc = rho_critical_newell(rho_max, v_limit)
            else:
                q_plot = q_of_rho(rho_plot, v_max, rho_max, v_limit)
                rc = rho_critical(rho_max)
            ax_fd.plot(rho_plot, q_plot, color=color, linewidth=1.4, label=vl_label)
            ax_fd.axvline(rc, color=color, linewidth=0.7, linestyle=':')

            # --- final-time density snapshot ---
            density = data['density']
            x       = data['x']
            ax_snap.plot(x, density[-1], color=color, linewidth=1.4, label=vl_label)

        n_steps = attrs.get('n_steps', '?')
        ic      = attrs.get('ic_type', '?')

        ax_fd.set_xlabel('Density ρ', fontsize=9)
        ax_fd.set_ylabel('Flow q', fontsize=9)
        ax_fd.set_xlim(0, rho_max)
        ax_fd.set_ylim(0)
        ax_fd.set_title(f'{label}\nFundamental diagrams', fontsize=9, fontweight='bold')
        ax_fd.legend(fontsize=7, loc='upper right')

        ax_snap.set_xlabel('Position x', fontsize=9)
        ax_snap.set_ylabel('Density ρ', fontsize=9)
        ax_snap.set_xlim(float(datasets[0][1]['x'][0]), float(datasets[0][1]['x'][-1]))
        ax_snap.set_ylim(-0.02, rho_max * 1.08)
        ax_snap.set_title(f'Final density  (IC: {ic}, steps: {n_steps})', fontsize=9)
        ax_snap.legend(fontsize=7, loc='best')

    fig.suptitle('Speed limit comparison — shock behaviour', fontsize=13, y=1.01)

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f'Saved to {save_path}')

    return fig


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
