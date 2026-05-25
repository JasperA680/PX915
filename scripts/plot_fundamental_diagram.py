"""
Plot the fundamental diagram from Fortran sweep output.

Usage:
    python scripts/plot_fundamental_diagram.py
    python scripts/plot_fundamental_diagram.py --save
    python scripts/plot_fundamental_diagram.py --input path/to/other.nc
"""

import argparse
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))

from python.CA_visualisation import plot_fundamental_diagram
import matplotlib.pyplot as plt


def load_netcdf(path):
    with nc.Dataset(path) as ds:
        rho = np.concatenate([ds.variables['alpha_rho'][:], ds.variables['beta_rho'][:]])
        J   = np.concatenate([ds.variables['alpha_J'][:],   ds.variables['beta_J'][:]])
        L        = int(ds.getncattr('L'))
        n_points = int(ds.getncattr('n_points'))
    return rho, J, L, n_points


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', default=str(ROOT / 'data/output/fundamental_diagram.nc'),
                        help='NetCDF from make run-fd (default: data/output/fundamental_diagram.nc)')
    parser.add_argument('--save', action='store_true',
                        help='save figure to plots/fundamental_fortran.png')
    args = parser.parse_args()

    rho, J, L, n_points = load_netcdf(args.input)
    print(f'Loaded {len(rho)} points from {args.input}  (L={L}, {n_points} per branch)')

    fig, ax = plot_fundamental_diagram(rho, J)
    ax.set_title(f'Fundamental diagram  (Fortran, L={L})')

    if args.save:
        plots_dir = ROOT / 'plots'
        plots_dir.mkdir(exist_ok=True)
        path = plots_dir / 'fundamental_fortran.png'
        fig.savefig(path, dpi=150, bbox_inches='tight')
        print(f'Saved to {path}')

    plt.show()


if __name__ == '__main__':
    main()
