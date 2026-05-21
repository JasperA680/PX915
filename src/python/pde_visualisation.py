'''Visualisation routines for the LWR PDE and multilane PDE models.

This module contains plotting functions for the continuum traffic flow models.
The single lane routines visualise density fields, density snapshots, boundary 
flow and fundamental diagrams. The multilane routines visualise per-lane 
density fields, total density, lane density time series and total mass
conservation.

The expected input to most routines is the dictionary returned by 
``load_pde_netcdf``. For single lane data, the density is expected to have shape 
``(time, x)``. For multilane data, ``density_per_lane`` is expected to have 
shape ``(time, lane, x)``.
'''

import numpy as np
import matplotlib.pyplot as plt


def plot_pde_spacetime(data, ax=None, title=None):
    """Plot a space-time heatmap of the single lane PDE density.

    This is the main visual diagnostic for the LWR PDE solver. It shows the 
    density field as a function of position and time.

    Parameters
    ----------
    data : dict
        Dictionary returned by ``load_pde_netCDF``. It must contain the keys 
        ``"density"``, ``"x"``, ``"time"``, and ``"attrs"``. The density array
        is expected to have shape ``(n_times, n_cells)``.
    ax : matplotlib.axes.Axes, optional
        Existing axes on which to draw the plot. If ``None``, a new figure and 
        axes are created.
    title : str, optional
        Custom plot title. If ``None``, a title is generated from the initial 
        condition and flux type stored in ``data["attrs"]``.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the plot.
    ax : matplotlib.axes.Axes
        The axes object on which the plot was drawn.
    """
    density = data['density']   
    x       = data['x']
    time    = data['time']
    attrs   = data['attrs']
    rho_max = float(attrs.get('rho_max', 1.0))

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 5))
    else:
        fig = ax.figure

    im = ax.imshow(
        density,
        aspect='auto',
        origin='lower',
        cmap='viridis',
        vmin=0,
        vmax=rho_max,
        extent=[float(x[0]), float(x[-1]), float(time[0]), float(time[-1])],
    )
    fig.colorbar(im, ax=ax, label='Density ρ')
    ax.set_xlabel('Position x')
    ax.set_ylabel('Time t')
    ic   = attrs.get('ic_type',   '')
    flux = attrs.get('flux_type', '')
    ax.set_title(title or f'Space-time density  (IC: {ic},  flux: {flux})')
    return fig, ax


def plot_pde_snapshots(data, n_snapshots=6, ax=None, title=None):
    """Plot density profiles at evenly spaced times.

    This function overlays several line plots of ``rho(x)`` so that the
    temporal evolution of the density profile can be inspected.

    Parameters
    ----------
    data : dict
        Dictionary returned by ``load_pde_netcdf``. It must contain ``"density"``,
        ``"x"``, ``"time"`` and ``"attrs"``.
    n_snapshots : int, default=6
        Number of time snapshots to plot.
    ax : matplotlib.axes.Axes, optional
        Existing axes on which to draw the plot. If ``None``, a new figure and
        axes are created.
    title : str, optional
        Custom plot title. If ``None``, a title is generated from the initial
        condition stored in ``data["attrs"]``.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the plot.
    ax : matplotlib.axes.Axes
        The axes object on which the plot was drawn.
    """
    density = data['density']  
    x       = data['x']
    time    = data['time']
    attrs   = data['attrs']
    rho_max = float(attrs.get('rho_max', 1.0))

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4))
    else:
        fig = ax.figure

    indices = np.linspace(0, len(time) - 1, n_snapshots, dtype=int)
    colors  = plt.cm.plasma(np.linspace(0.15, 0.85, n_snapshots))

    for idx, color in zip(indices, colors):
        ax.plot(x, density[idx], color=color, linewidth=1.5,
                label=f't = {float(time[idx]):.3f}')

    ax.set_xlabel('Position x')
    ax.set_ylabel('Density ρ')
    ax.set_ylim(-0.02, rho_max * 1.08)
    ax.legend(fontsize=8, loc='best')
    ic = attrs.get('ic_type', '')
    ax.set_title(title or f'Density snapshots  (IC: {ic})')
    return fig, ax


def plot_pde_flow(data, ax=None, title=None):
    """Plot right-boundary flow as a function of time.

    The flow history is compared with the theoretical maximum flow for the
    selected fundamental diagram. For the Greenshields closure,

    ``q_max = v_max * rho_max / 4``.

    For the Newell-Daganzo triangular closure,

    ``q_max = v_max * NEWELL_W * rho_max / (v_max + NEWELL_W)``.

    Parameters
    ----------

    data : dict
        Dictionary returned by ``load_pde_netcdf``. It must contain ``"flow"``,
        ``"time"`` and ``"attrs"``.
    ax : matplotlib.axes.Axes, optional
        Existing axes on which to draw the plot. If ``None``, a new figure and
        axes are created.
    title : str, optional
        Custom plot title. If ``None``, a default title is used.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the plot.
    ax : matplotlib.axes.Axes
        The axes object on which the plot was drawn.
    """
    flow  = data['flow']
    time  = data['time']
    attrs = data['attrs']
    v_max     = float(attrs.get('v_max',   1.0))
    rho_max   = float(attrs.get('rho_max', 1.0))
    flux_type = str(attrs.get('flux_type', 'lf'))
    if flux_type == 'newell':
        NEWELL_W = 0.5
        q_max = v_max * NEWELL_W * rho_max / (v_max + NEWELL_W)
    else:
        q_max = v_max * rho_max / 4.0

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 3))
    else:
        fig = ax.figure

    ax.plot(time, flow, color='steelblue', linewidth=1.2)
    ax.axhline(q_max, color='tomato', linestyle='--', linewidth=1,
               label=f'q_max = {q_max:.3f}')
    ax.set_xlabel('Time t')
    ax.set_ylabel('Flow q(ρ_M)')
    ax.set_ylim(0, q_max * 1.25)
    ax.set_title(title or 'Right-boundary flow vs time')
    ax.legend(fontsize=9)
    return fig, ax


def plot_pde_summary(data, save_path=None):
    """Create a three-panel summary figure for a single-lane PDE run.

    The summary contains a space-time density heatmap, several density
    snapshots, and the right-boundary flow history.

    Parameters
    ----------
    data : dict
        Dictionary returned by ``load_pde_netcdf``.
    save_path : str or pathlib.Path, optional
        If provided, the figure is saved to this location using
        ``bbox_inches="tight"``.

    Returns
    -------
    matplotlib.figure.Figure
        Figure containing the three-panel PDE summary.
    """
    attrs  = data['attrs']
    M      = attrs.get('M',        '?')
    n_steps = attrs.get('n_steps', '?')
    ic     = attrs.get('ic_type',   '?')
    flux   = attrs.get('flux_type', '?')
    bc     = attrs.get('bc_type',   '?')

    fig = plt.figure(figsize=(14, 9))
    gs  = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.38)

    ax_st   = fig.add_subplot(gs[0, :])
    ax_snap = fig.add_subplot(gs[1, 0])
    ax_flow = fig.add_subplot(gs[1, 1])

    plot_pde_spacetime(data, ax=ax_st)
    plot_pde_snapshots(data, ax=ax_snap)
    plot_pde_flow(data, ax=ax_flow)

    fig.suptitle(
        f'LWR PDE  M={M}  n_steps={n_steps}  IC={ic}  flux={flux}  BC={bc}',
        fontsize=12,
    )

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f'Saved to {save_path}')

    return fig


def plot_pde_fundamental_diagram(rho, q, v_max=1.0, rho_max=1.0, ax=None, title=None):
    """Plot a PDE fundamental diagram against the Greenshields curve.

    The supplied density and flow samples are plotted as scatter points and
    compared with the analytical Greenshields relation,

    ``q(rho) = v_max * rho * (1 - rho/rho_max)``.

    Parameters
    ----------

    rho : array-like
        Density values from a PDE fundamental diagram sweep.
    q : array-like
        Flow values corresponding to ``rho``.
    v_max : float, default=1.0
        Maximum/free-flow velocity used in the analytical Greenshields curve.
    rho_max : float, default=1.0
        Maximum/jam density used in the analytical Greenshields curve.
    ax : matplotlib.axes.Axes, optional
        Existing axes on which to draw the plot. If ``None``, a new figure and
        axes are created.
    title : str, optional
        Custom plot title. If ``None``, a default title is used.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the plot.
    ax : matplotlib.axes.Axes
        The axes object on which the plot was drawn.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))
    else:
        fig = ax.figure

    rho_theory = np.linspace(0, rho_max, 300)
    q_theory = v_max * rho_theory * (1.0 - rho_theory / rho_max)
    q_max = v_max * rho_max / 4.0
    rho_c = rho_max / 2.0

    ax.plot(rho_theory, q_theory, 'k--', linewidth=1.5,
            label='q(ρ) = v·ρ(1−ρ/ρ_max)  [Greenshields]')
    ax.scatter(rho, q, s=25, color='steelblue', alpha=0.9, zorder=3,
               label='PDE sweep')
    ax.axvline(rho_c, color='grey', linestyle=':', linewidth=1,
               label=f'ρ_c = {rho_c:.2f}')
    ax.axhline(q_max, color='grey', linestyle=':', linewidth=1,
               label=f'q_max = {q_max:.3f}')
    ax.set_xlabel('Density ρ')
    ax.set_ylabel('Flow q')
    ax.set_xlim(0, rho_max)
    ax.set_ylim(0, q_max * 1.3)
    ax.set_title(title or 'PDE fundamental diagram')
    ax.legend(fontsize=9)
    return fig, ax


def plot_space_time_per_lane(data, fig=None, title=None):
    """Plot one space-time density heatmap for each lane.

    This function creates stacked heatmaps of the lane-resolved density field.
    It is useful for checking how lane-changing redistributes density between
    lanes over time.

    Parameters
    ----------
    data : dict
        Dictionary returned by ``load_pde_netcdf`` for a multilane PDE run. It
        must contain ``"density_per_lane"``, ``"x"``, ``"time"``, ``"attrs"``
        and ``"n_lanes"``. The density array must have shape
        ``(n_times, n_lanes, n_cells)``.

    fig : matplotlib.figure.Figure, optional
        Existing figure to draw into. If ``None``, a new figure and axes are
        created.

    title : str, optional
        Custom figure title. If ``None``, a default title is used.

     Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the plot.
    ax : matplotlib.axes.Axes
        The axes object on which the plot was drawn.
    """
    density_pl = data["density_per_lane"]   
    n_lanes    = density_pl.shape[1]
    x          = data["x"]
    time       = data["time"]
    attrs      = data["attrs"]
    rho_max    = float(attrs.get("rho_max", 1.0))

    if fig is None:
        fig, axes = plt.subplots(n_lanes, 1, figsize=(10, 3 * n_lanes),
                                 sharex=True, sharey=True)
        if n_lanes == 1:
            axes = [axes]
    else:
        axes = fig.get_axes()

    for lane_idx in range(n_lanes):
        ax = axes[lane_idx]
        im = ax.imshow(
            density_pl[:, lane_idx, :],
            aspect='auto', origin='lower', cmap='viridis',
            vmin=0, vmax=rho_max,
            extent=[float(x[0]), float(x[-1]), float(time[0]), float(time[-1])],
        )
        fig.colorbar(im, ax=ax, label='ρ')
        ax.set_ylabel(f'Time t\n(lane {lane_idx + 1})')

    axes[-1].set_xlabel('Position x')
    fig.suptitle(title or 'Space-time density per lane', y=1.01)
    fig.tight_layout()
    return fig, axes


def plot_space_time_total(data, ax=None, title=None):
    """Plot a space-time heatmap of total density over all lanes.

    The total density is defined as

    ``rho_tot(x, t) = sum_lanes rho_lane(x, t)``.

    Parameters
    ----------
    data : dict
        Dictionary returned by ``load_pde_netcdf`` for a multilane PDE run. It
        must contain ``"density_per_lane"``, ``"x"``, ``"time"``, ``"attrs"``
        and ``"n_lanes"``.
    ax : matplotlib.axes.Axes, optional
        Existing axes on which to draw the plot. If ``None``, a new figure and
        axes are created.
    title : str, optional
        Custom plot title. If ``None``, a title is generated from the number of
        lanes.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the plot.
    ax : matplotlib.axes.Axes
        The axes object on which the plot was drawn.
    """
    from analysis import compute_total_density
    rho_tot = compute_total_density(data)   
    x       = data["x"]
    time    = data["time"]
    attrs   = data["attrs"]
    n_lanes = data["n_lanes"]
    rho_max = float(attrs.get("rho_max", 1.0)) * n_lanes

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 5))
    else:
        fig = ax.figure

    im = ax.imshow(
        rho_tot, aspect='auto', origin='lower', cmap='viridis',
        vmin=0, vmax=rho_max,
        extent=[float(x[0]), float(x[-1]), float(time[0]), float(time[-1])],
    )
    fig.colorbar(im, ax=ax, label='Total density ρ_tot')
    ax.set_xlabel('Position x')
    ax.set_ylabel('Time t')
    ax.set_title(title or f'Total space-time density  ({n_lanes} lanes)')
    return fig, ax


def plot_lane_densities(data, x_pos=0.5, ax=None, title=None):
    """Plot lane-resolved density against time at a fixed position.

    The spatial cell closest to ``x_pos`` is selected, then the density in each
    lane at that cell is plotted as a function of time.

    Parameters
    ----------
    data : dict
        Dictionary returned by ``load_pde_netcdf`` for a multilane PDE run. It
        must contain ``"density_per_lane"``, ``"x"``, ``"time"`` and
        ``"n_lanes"``.
    x_pos : float, default=0.5
        Physical position along the road at which to sample the density. The
        nearest grid cell is used.
    ax : matplotlib.axes.Axes, optional
        Existing axes on which to draw the plot. If ``None``, a new figure and
        axes are created.
    title : str, optional
        Custom plot title. If ``None``, a title is generated from the selected
        cell position.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the plot.
    ax : matplotlib.axes.Axes
        The axes object on which the plot was drawn.
    """
    density_pl = data["density_per_lane"]   
    x          = data["x"]
    time       = data["time"]
    n_lanes    = data["n_lanes"]

    i_cell = int(np.argmin(np.abs(x - x_pos)))

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 4))
    else:
        fig = ax.figure

    colors = plt.cm.tab10(np.linspace(0, 0.9, n_lanes))
    for lane_idx in range(n_lanes):
        ax.plot(time, density_pl[:, lane_idx, i_cell],
                color=colors[lane_idx], linewidth=1.4,
                label=f'Lane {lane_idx + 1}')

    ax.set_xlabel('Time t')
    ax.set_ylabel('Density ρ')
    ax.set_title(title or f'Lane densities at x ≈ {float(x[i_cell]):.3f}')
    ax.legend(fontsize=9)
    return fig, ax


def plot_total_mass(data, ax=None, title=None):
    """Plot the deviation of total mass from its initial value.

    This is a conservation diagnostic for the multilane PDE solver. The plotted
    quantity is

    ``Delta mass(t) = mass(t) - mass(0)``.

    Plotting the deviation rather than the absolute mass avoids Matplotlib's
    offset notation and makes conservation errors easier to see.

    Parameters
    ----------
    data : dict
        Dictionary returned by ``load_pde_netcdf`` for a PDE run. It must contain
        the data required by ``analysis.compute_total_mass`` and the key
        ``"time"``.
    ax : matplotlib.axes.Axes, optional
        Existing axes on which to draw the plot. If ``None``, a new figure and
        axes are created.
    title : str, optional
        Custom plot title. If ``None``, a default conservation-check title is
        used.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the plot.
    ax : matplotlib.axes.Axes
        The axes object on which the plot was drawn.
    """
    from analysis import compute_total_mass
    mass      = compute_total_mass(data)
    deviation = mass - mass[0]
    time      = data["time"]

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 3))
    else:
        fig = ax.figure

    variation = float(deviation.max() - deviation.min())
    ax.plot(time, deviation, color='steelblue', linewidth=1.2)
    ax.axhline(0.0, color='tomato', linestyle='--', linewidth=1,
               label=f'zero (initial mass = {float(mass[0]):.4f}),  range = {variation:.2e}')
    ax.set_xlabel('Time t')
    ax.set_ylabel('Δ Total mass')
    ax.set_title(title or 'Mass conservation check  (deviation from initial mass)')
    ax.legend(fontsize=9)
    return fig, ax
