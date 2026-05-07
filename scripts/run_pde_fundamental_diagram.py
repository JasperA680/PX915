"""Generate the LWR PDE fundamental diagram (q vs rho) by sweeping boundary densities.

Two branches are swept, mirroring the TASEP fundamental_diagram() approach:

  Free-flow  (rho < rho_c): rho_left_bc varies from 0 to rho_c; rho_right_bc = 0.
  Congested  (rho > rho_c): rho_right_bc varies from rho_c to rho_max; rho_left_bc = rho_max.

Each run uses a Riemann IC so a genuine transient forms before measurement.

Usage:
    python scripts/run_pde_fundamental_diagram.py
    python scripts/run_pde_fundamental_diagram.py --save
    python scripts/run_pde_fundamental_diagram.py --points 30 --M 400
"""

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src' / 'python'))

from analysis import pde_fundamental_diagram
from visualisation import plot_pde_fundamental_diagram
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(
        description='LWR PDE fundamental diagram via boundary-density sweep',
    )
    parser.add_argument('--save', action='store_true',
                        help='save figure to plots/pde_fundamental.png')
    parser.add_argument('--points', type=int, default=20,
                        help='sweep points per branch (default 20, total = 2x)')
    parser.add_argument('--M', type=int, default=200,
                        help='spatial cells per run (default 200)')
    parser.add_argument('--steps', type=int, default=2000,
                        help='time steps per run (default 2000)')
    parser.add_argument('--burnin', type=float, default=0.5,
                        help='burn-in fraction discarded before averaging (default 0.5)')
    args = parser.parse_args()

    total_runs = 2 * args.points
    print(f'PDE fundamental diagram sweep: M={args.M}, {args.points} pts/branch '
          f'({total_runs} runs, {args.steps} steps each)')

    rho, q = pde_fundamental_diagram(
        n_points=args.points,
        n_steps=args.steps,
        M=args.M,
        burnin_frac=args.burnin,
    )

    print(f'Done. Plotting {len(rho)} points.')

    fig, ax = plot_pde_fundamental_diagram(rho, q)
    ax.set_title(
        f'PDE fundamental diagram  (M={args.M}, {args.points} pts/branch)'
    )

    if args.save:
        plots_dir = ROOT / 'plots'
        plots_dir.mkdir(exist_ok=True)
        path = plots_dir / 'pde_fundamental.png'
        fig.savefig(path, dpi=150, bbox_inches='tight')
        print(f'Saved to {path}')

    plt.show()


if __name__ == '__main__':
    main()
