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
    ax.plot(rho_theory, rho_theory * (1 - rho_theory), 'k--',
            linewidth=1.5, label='J = ρ(1−ρ)  theory')
    ax.scatter(rho_vals, J_vals, s=18, color='steelblue',
               alpha=0.8, zorder=3, label='simulation')
    ax.set_xlabel('Density ρ')
    ax.set_ylabel('Current J')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 0.30)
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
# Network-aware plots (operate on io.NetworkResult)
# ---------------------------------------------------------------------------

def _config_layout(result):
    """Pull (junctions, road_endpoints) dicts out of result.config['layout']."""
    layout = (result.config or {}).get('layout', {})
    j = {int(d['id']): (float(d['x']), float(d['y']))
         for d in layout.get('junctions', [])}
    r = {int(d['id']): (tuple(d['end_1']), tuple(d['end_2']))
         for d in layout.get('roads', [])}
    return j, r


def plot_network_layout(result, ax=None, occupancy_t: Optional[int] = None,
                        title: Optional[str] = None):
    """Draw road segments and junction nodes.

    If `occupancy_t` is given, colour each road by its mean lane occupancy
    at that timestep.  Otherwise roads are drawn in a uniform colour.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 6))
    else:
        fig = ax.figure

    j_xy, r_ends = _config_layout(result)
    n_roads = len(r_ends)

    # Per-road colour value at time t (mean occupancy across all lanes of that road).
    colours = None
    if occupancy_t is not None:
        colours = np.zeros(n_roads)
        # Map road_id -> mean occupancy
        rid_list = sorted(r_ends.keys())
        lane_road = result.lane_road_id
        lane_len  = result.lane_length
        occ_t = result.occupancy[occupancy_t]   # (lane_index, cell)
        for i, rid in enumerate(rid_list):
            lanes_of_r = np.where(lane_road == rid)[0]
            if len(lanes_of_r) == 0:
                colours[i] = 0.0
                continue
            total_cells = 0
            total_occ = 0
            for li in lanes_of_r:
                Lk = lane_len[li]
                total_cells += Lk
                total_occ += int(occ_t[li, :Lk].sum())
            colours[i] = total_occ / max(1, total_cells)

    # Draw roads.
    segs = []
    rid_list = sorted(r_ends.keys())
    for rid in rid_list:
        e1, e2 = r_ends[rid]
        segs.append([e1, e2])
    if colours is not None:
        lc = LineCollection(segs, cmap='viridis', linewidths=4)
        lc.set_array(colours)
        lc.set_clim(0, 1)
        ax.add_collection(lc)
        fig.colorbar(lc, ax=ax, label='occupancy', shrink=0.7)
    else:
        lc = LineCollection(segs, colors='steelblue', linewidths=3)
        ax.add_collection(lc)

    # Road labels at midpoint.
    for rid in rid_list:
        e1, e2 = r_ends[rid]
        mx, my = 0.5 * (e1[0] + e2[0]), 0.5 * (e1[1] + e2[1])
        ax.text(mx, my, f'R{rid}', fontsize=8, color='gray',
                ha='center', va='center', backgroundcolor='white')

    # Draw junctions.
    if j_xy:
        jx = [p[0] for p in j_xy.values()]
        jy = [p[1] for p in j_xy.values()]
        ax.scatter(jx, jy, s=140, c='black', zorder=3)
        for jid, (x, y) in j_xy.items():
            ax.text(x, y, f'J{jid}', color='white', fontsize=8,
                    ha='center', va='center', zorder=4, fontweight='bold')

    # Open-boundary markers: square at the "open" endpoint of each road.
    open_in_roads  = result.lane_road_id[result.lane_open_in.astype(bool)]
    open_out_roads = result.lane_road_id[result.lane_open_out.astype(bool)]
    open_roads = set(int(r) for r in open_in_roads) | set(int(r) for r in open_out_roads)
    for rid in sorted(open_roads):
        if rid not in r_ends:
            continue
        e1, e2 = r_ends[rid]
        # The "open" end is whichever endpoint is at an open boundary
        # (end_junction == 0).  Look up from result.road_end_junction.
        if rid - 1 < len(result.road_end_junction):
            ej = result.road_end_junction[rid - 1]
            if ej[0] == 0:
                ax.plot(*e1, marker='s', color='tomato', markersize=8, zorder=5)
            if ej[1] == 0:
                ax.plot(*e2, marker='s', color='tomato', markersize=8, zorder=5)

    # Layout sizing.
    all_pts = [pt for ends in r_ends.values() for pt in ends] + list(j_xy.values())
    if all_pts:
        xs = [p[0] for p in all_pts]
        ys = [p[1] for p in all_pts]
        pad = 0.3 * max(max(xs) - min(xs), max(ys) - min(ys), 1.0)
        ax.set_xlim(min(xs) - pad, max(xs) + pad)
        ax.set_ylim(min(ys) - pad, max(ys) + pad)
    ax.set_aspect('equal')
    ax.set_title(title or f'Network layout  (t={occupancy_t if occupancy_t is not None else "—"})')
    ax.set_xticks([])
    ax.set_yticks([])
    return fig, ax


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
