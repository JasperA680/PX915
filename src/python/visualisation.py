import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


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
