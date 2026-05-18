"""Visualisation functions for the LWR PDE (continuum) model.

All PDE-specific plots live here; CA/network plots are in visualisation.py.
"""

import numpy as np
import matplotlib.pyplot as plt


def plot_pde_spacetime(data, ax=None, title=None):
    """Heatmap of ρ(x, t) — the primary visual diagnostic for the LWR PDE.

    data must be the dict returned by load_pde_netcdf:
        density (n_steps+1, M), x (M,), time (n_steps+1,), attrs dict.
    """
    density = data['density']   # (time, x)
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
    """Line plots of ρ(x) at n_snapshots evenly-spaced times."""
    density = data['density']   # (time, x)
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
    """Right-boundary flow q(ρ_M) vs time."""
    flow  = data['flow']
    time  = data['time']
    attrs = data['attrs']
    v_max     = float(attrs.get('v_max',   1.0))
    rho_max   = float(attrs.get('rho_max', 1.0))
    flux_type = str(attrs.get('flux_type', 'lf'))
    if flux_type == 'newell':
        # Newell: q_max = v_f * rho_c = v_f * w * rho_max / (v_f + w)
        # NEWELL_W must match the Fortran constant in pde_flux.f90
        NEWELL_W = 0.5
        q_max = v_max * NEWELL_W * rho_max / (v_max + NEWELL_W)
    else:
        # Greenshields: q_max = v_max * rho_max / 4
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
    """Three-panel PDE summary: space-time heatmap, density snapshots, boundary flow."""
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
    """Flow q vs density ρ with the Greenshields analytical parabola overlaid.

    Parameters
    ----------
    rho, q : array-like
        Sweep points from ``pde_fundamental_diagram()``.
    v_max, rho_max : float
        Greenshields parameters for the analytical curve.
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


# ---------------------------------------------------------------------------
# Multi-lane PDE visualisation
# ---------------------------------------------------------------------------

def plot_space_time_per_lane(data, fig=None, title=None):
    """N stacked heatmaps of ρ_lane(x, t), one subplot per lane.

    Parameters
    ----------
    data : dict returned by load_pde_netcdf
    """
    density_pl = data["density_per_lane"]   # (time, lane, x)
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
    """Heatmap of total density ρ_tot(x, t) = Σ_lane ρ_lane(x, t)."""
    from analysis import compute_total_density
    rho_tot = compute_total_density(data)   # (time, x)
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
    """Per-lane density vs time at a fixed spatial position x_pos.

    Parameters
    ----------
    x_pos : float
        Position along the road (in domain units). Nearest cell is used.
    """
    density_pl = data["density_per_lane"]   # (time, lane, x)
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
    """Mass deviation Δmass(t) = mass(t) − mass(0) vs time.

    Plotting the deviation from initial mass rather than the absolute value
    avoids matplotlib's offset notation and makes conservation quality
    immediately readable: a perfect simulation gives a flat line at zero.
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
