"""Multi-lane fundamental diagram: compare 1-lane vs 2-lane capacity.

For a road with N identical lanes (each Greenshields, v_max=1, rho_max=1),
the total flow at total density rho_tot is theoretically

    q_tot(rho_tot) = rho_tot · (1 − rho_tot/N)

obtained by letting each lane carry rho_tot/N at flow (rho_tot/N)(1 − rho_tot/N)
and summing.  Peak capacity is q_max = N/4 at rho_tot = N/2 (= 2 × single-lane
capacity for N = 2).

This script runs the existing single-lane fundamental-diagram sweep and the
multi-lane sweep, then overlays both against their analytical curves.
"""

import argparse
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / "src" / "python"))

from analysis import pde_fundamental_diagram, pde_multilane_fundamental_diagram

PLOTS_DIR = ROOT / "plots"

V_MAX   = 1.0
RHO_MAX = 1.0


def theory_curve(rho_tot, n_lanes, v_max=V_MAX, rho_max=RHO_MAX):
    """Greenshields total flow for N identical lanes, each carrying rho_tot/N."""
    rho_per_lane = rho_tot / n_lanes
    valid = rho_per_lane <= rho_max
    q = np.where(valid,
                 v_max * rho_tot * (1.0 - rho_per_lane / rho_max),
                 0.0)
    return q


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--save", action="store_true",
                        help="Save figure to plots/multilane_fundamental_diagram.png")
    parser.add_argument("--no-show", action="store_true",
                        help="Don't display the figure interactively")
    parser.add_argument("--n-points", type=int, default=14,
                        help="Sweep points per branch (default 14)")
    parser.add_argument("--n-steps", type=int, default=2000,
                        help="Time steps per simulation (default 2000)")
    args = parser.parse_args()

    print("Building solver...")
    subprocess.run(["make", "pde"], check=True, cwd=ROOT, capture_output=True)

    print(f"\nSweeping single-lane (n_points={args.n_points} per branch)...")
    rho1, q1 = pde_fundamental_diagram(
        n_points=args.n_points, n_steps=args.n_steps,
        v_max=V_MAX, rho_max=RHO_MAX,
    )

    print(f"\nSweeping two-lane k=0 (independent lanes)...")
    rho2_k0, q2_k0 = pde_multilane_fundamental_diagram(
        n_points=args.n_points, n_steps=args.n_steps,
        v_max=V_MAX, rho_max=RHO_MAX,
        n_lanes=2, lane_change_rate=0.0,
    )

    print(f"\nSweeping two-lane k=0.5 (coupled lanes)...")
    rho2_k5, q2_k5 = pde_multilane_fundamental_diagram(
        n_points=args.n_points, n_steps=args.n_steps,
        v_max=V_MAX, rho_max=RHO_MAX,
        n_lanes=2, lane_change_rate=0.5,
    )

    # Theory curves
    rho_fine_1 = np.linspace(0, RHO_MAX, 300)
    rho_fine_2 = np.linspace(0, 2 * RHO_MAX, 300)
    q_theory_1 = theory_curve(rho_fine_1, n_lanes=1)
    q_theory_2 = theory_curve(rho_fine_2, n_lanes=2)

    # Capacity points
    rho_c_1 = RHO_MAX / 2.0
    q_max_1 = V_MAX * RHO_MAX / 4.0
    rho_c_2 = RHO_MAX
    q_max_2 = V_MAX * RHO_MAX / 2.0

    # ---------------- Plot ----------------
    fig, ax = plt.subplots(figsize=(10, 6.5))

    # Theory curves
    ax.plot(rho_fine_1, q_theory_1, color="steelblue", linewidth=2.5, alpha=0.85,
            label=r"1-lane theory:  $q = \rho(1-\rho)$,  $q_{\max} = 0.25$")
    ax.plot(rho_fine_2, q_theory_2, color="firebrick", linewidth=2.5, alpha=0.85,
            label=r"2-lane theory:  $q_{\rm tot} = \rho_{\rm tot}(1 - \rho_{\rm tot}/2)$,  $q_{\max} = 0.50$")

    # Simulation points
    ax.scatter(rho1, q1, s=70, color="steelblue", marker="o",
               edgecolors="black", linewidths=1, zorder=4,
               label="1-lane simulation")
    ax.scatter(rho2_k0, q2_k0, s=70, color="firebrick", marker="s",
               edgecolors="black", linewidths=1, zorder=4,
               label="2-lane simulation  ($k = 0$)")
    ax.scatter(rho2_k5, q2_k5, s=55, color="goldenrod", marker="D",
               edgecolors="black", linewidths=1, zorder=4, alpha=0.85,
               label="2-lane simulation  ($k = 0.5$)")

    # Capacity markers
    ax.axhline(q_max_1, color="steelblue", linestyle=":", linewidth=1.2, alpha=0.6)
    ax.axhline(q_max_2, color="firebrick", linestyle=":", linewidth=1.2, alpha=0.6)
    ax.axvline(rho_c_1, color="steelblue", linestyle=":", linewidth=1.2, alpha=0.5)
    ax.axvline(rho_c_2, color="firebrick", linestyle=":", linewidth=1.2, alpha=0.5)

    ax.set_xlabel(r"Total density  $\rho_{\rm tot}$", fontsize=12)
    ax.set_ylabel(r"Total flow  $q_{\rm tot}$", fontsize=12)
    ax.set_title(r"Multi-lane fundamental diagram"
                 "\n"
                 r"Greenshields, $v_{\max} = 1$, $\rho_{\max} = 1$",
                 fontsize=11)
    ax.legend(loc="upper right", fontsize=9.5, framealpha=0.95)
    ax.set_xlim(0, 2 * RHO_MAX * 1.02)
    ax.set_ylim(0, q_max_2 * 1.25)
    ax.grid(alpha=0.3)

    fig.tight_layout()

    if args.save:
        PLOTS_DIR.mkdir(exist_ok=True)
        path = PLOTS_DIR / "multilane_fundamental_diagram.png"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        print(f"\nSaved {path}")

    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
