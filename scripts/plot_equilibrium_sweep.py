"""Validate the multi-lane PDE equilibrium against the closed-form analytical solution.

For two lanes with periodic boundary conditions, equal jam densities
rho_max = 1, and fixed total mass M0 = rho1(0) + rho2(0), the lane-change
source term reaches a steady state when velocities equalise:

    v_1(rho_1_eq) = v_2(rho_2_eq)
    => 1·(1 - rho_1_eq) = r·(1 - rho_2_eq),       r = v_max_2 / v_max_1.

Combined with mass conservation rho_1 + rho_2 = M0 this yields:

    rho_2_eq = (r - 1 + M0) / (r + 1)
    rho_1_eq = M0 - rho_2_eq.

This script sweeps r over a range, runs the multi-lane solver to steady state
for each r, and overlays the simulated points on the analytical curves.

Usage
-----
    python scripts/plot_equilibrium_sweep.py
    python scripts/plot_equilibrium_sweep.py --save
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / "src" / "python"))

from pde_runner import run_pde, load_pde_netcdf

PLOTS_DIR = ROOT / "plots"


def analytical_equilibrium(r, M0=0.9):
    """Closed-form equilibrium densities for two-lane Greenshields with rho_max=1.

    Parameters
    ----------
    r : float
        Velocity ratio v_max_2 / v_max_1 (>= 1).
    M0 : float
        Total mass per cell pair (rho_1 + rho_2 at every cell).

    Returns
    -------
    (rho_1_eq, rho_2_eq) : tuple of floats
    """
    rho2 = (r - 1.0 + M0) / (r + 1.0)
    rho1 = M0 - rho2
    return rho1, rho2


def simulate_equilibrium(r, M0=0.9, n_steps=2000, M_cells=100, burnin_frac=0.75):
    """Run a two-lane simulation with v_max_lanes=[1.0, r] and return steady-state densities.

    Spatial structure is irrelevant (constant IC, periodic BCs), so M_cells can
    be small for speed.  The last (1 - burnin_frac) of time steps is averaged
    to get the steady-state estimate.
    """
    rho_init = M0 / 2.0
    with tempfile.TemporaryDirectory() as tmpdir:
        out = Path(tmpdir) / f"r_{r:.3f}.nc"
        run_pde(
            dict(
                n_lanes=2,
                lane_change_rate=0.5,
                v_max=1.0, rho_max=1.0,
                v_max_lanes=[1.0, r],
                rho_left_bc=rho_init, rho_right_bc=rho_init,
                ic_type="constant",
                flux_type="godunov",
                bc_type="periodic",
                M=M_cells, n_steps=n_steps,
            ),
            output_path=out,
        )
        data = load_pde_netcdf(out)

    rho_pl = data["density_per_lane"]                # (time, lane, x)
    burnin = int(burnin_frac * n_steps)
    rho1_ss = float(rho_pl[burnin:, 0, :].mean())
    rho2_ss = float(rho_pl[burnin:, 1, :].mean())
    return rho1_ss, rho2_ss


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--save", action="store_true",
                        help="Save figure to plots/multilane_equilibrium_sweep.png")
    parser.add_argument("--no-show", action="store_true",
                        help="Don't display the figure interactively")
    parser.add_argument("--n-steps", type=int, default=3000,
                        help="Time steps per simulation (default: 3000)")
    args = parser.parse_args()

    print("Building solver...")
    subprocess.run(["make", "pde"], check=True, cwd=ROOT, capture_output=True)

    # Sweep over velocity ratios
    r_values = np.array([1.0, 1.1, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0])
    rho1_sim = np.zeros_like(r_values)
    rho2_sim = np.zeros_like(r_values)

    print(f"\nRunning sweep over {len(r_values)} velocity ratios "
          f"(n_steps={args.n_steps} per run)...")
    print(f"  {'r':>6}  {'rho1_sim':>10}  {'rho2_sim':>10}  "
          f"{'rho1_th':>10}  {'rho2_th':>10}  {'|err|':>10}")
    for i, r in enumerate(r_values):
        rho1_sim[i], rho2_sim[i] = simulate_equilibrium(r, n_steps=args.n_steps)
        rho1_th, rho2_th = analytical_equilibrium(r)
        err = max(abs(rho1_sim[i] - rho1_th), abs(rho2_sim[i] - rho2_th))
        print(f"  {r:6.2f}  {rho1_sim[i]:10.4f}  {rho2_sim[i]:10.4f}  "
              f"{rho1_th:10.4f}  {rho2_th:10.4f}  {err:10.2e}")

    # Analytical curves (fine grid)
    r_fine = np.linspace(1.0, 3.0, 300)
    rho1_curve = np.array([analytical_equilibrium(r)[0] for r in r_fine])
    rho2_curve = np.array([analytical_equilibrium(r)[1] for r in r_fine])

    # Plot
    fig, ax = plt.subplots(figsize=(9, 6))

    ax.plot(r_fine, rho1_curve, color="steelblue", linewidth=2, alpha=0.7,
            label=r"$\rho_1^{\mathrm{eq}}$ (slow lane) — theory")
    ax.plot(r_fine, rho2_curve, color="firebrick", linewidth=2, alpha=0.7,
            label=r"$\rho_2^{\mathrm{eq}}$ (fast lane) — theory")
    ax.scatter(r_values, rho1_sim, s=80, color="steelblue", marker="o",
               edgecolors="black", linewidths=1, zorder=5,
               label=r"$\rho_1^{\mathrm{ss}}$ (slow lane) — simulation")
    ax.scatter(r_values, rho2_sim, s=80, color="firebrick", marker="s",
               edgecolors="black", linewidths=1, zorder=5,
               label=r"$\rho_2^{\mathrm{ss}}$ (fast lane) — simulation")
    ax.axhline(0.45, color="gray", linestyle=":", linewidth=1, alpha=0.6,
               label="initial density (both lanes)")

    ax.set_xlabel(r"Velocity ratio $r = v_{\max,2}\,/\,v_{\max,1}$", fontsize=12)
    ax.set_ylabel("Steady-state density", fontsize=12)
    ax.set_title(
        "Multi-lane PDE lane changing equilibrium: simulation vs analytical theory\n"
        r"Initial mass $M_0 = 0.9$,  $\rho_{\max} = 1$,  $k = 0.5$,  periodic BCs",
        fontsize=11,
    )
    ax.legend(loc="center right", fontsize=10, framealpha=0.95)
    ax.set_xlim(0.95, 3.05)
    ax.set_ylim(0, 1)
    ax.grid(alpha=0.3)

    fig.tight_layout()

    if args.save:
        PLOTS_DIR.mkdir(exist_ok=True)
        path = PLOTS_DIR / "multilane_equilibrium_sweep.png"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        print(f"\nSaved {path}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
