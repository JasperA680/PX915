"""
Generate the TASEP fundamental diagram (J vs rho) by sweeping alpha and beta.

Usage:
    python scripts/run_fundamental_diagram.py
    python scripts/run_fundamental_diagram.py --save        # write plots/fundamental.png
    python scripts/run_fundamental_diagram.py --L 100       # larger lattice
    python scripts/run_fundamental_diagram.py --points 50   # finer sweep
"""

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))

from python.analysis import fundamental_diagram
from python.CA_visualisation import plot_fundamental_diagram
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--save',   action='store_true', help='save figure to plots/fundamental.png')
    parser.add_argument('--L',      type=int, default=50,  help='lattice size (default 50)')
    parser.add_argument('--points', type=int, default=30,  help='points per sweep (default 30)')
    args = parser.parse_args()

    print(f'Sweeping L={args.L}, {args.points} points per branch...')
    rho, J = fundamental_diagram(L=args.L, n_points=args.points)
    print(f'Done. Plotting...')

    fig, ax = plot_fundamental_diagram(rho, J)
    ax.set_title(f'Fundamental diagram  (L={args.L})')

    if args.save:
        plots_dir = ROOT / 'plots'
        plots_dir.mkdir(exist_ok=True)
        path = str(plots_dir / 'fundamental.png')
        fig.savefig(path, dpi=150, bbox_inches='tight')
        print(f'Saved to {path}')

    plt.show()


if __name__ == '__main__':
    main()
