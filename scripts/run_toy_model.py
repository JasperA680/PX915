"""
Run the compiled Fortran simulation then visualise the NetCDF output.

Usage:
    python scripts/run_toy_model.py              # uses default data/output/simulation.nc
    python scripts/run_toy_model.py --no-run     # skip Fortran, just plot existing file
    python scripts/run_toy_model.py --save       # also write plots/summary.png
"""

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))

from python.io import load_netcdf
from python.visualisation import plot_summary
import matplotlib.pyplot as plt

NC_PATH = ROOT / 'data' / 'output' / 'simulation.nc'
BINARY  = ROOT / 'build' / 'test_simulation'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--no-run', action='store_true',
                        help='skip Fortran build/run, just plot existing NetCDF')
    parser.add_argument('--save', action='store_true',
                        help='save summary figure to plots/summary.png')
    args = parser.parse_args()

    if not args.no_run:
        print('Building...')
        subprocess.run(['make'], cwd=ROOT, check=True)
        print('Running simulation...')
        subprocess.run([str(BINARY)], cwd=ROOT, check=True)

    print(f'Loading {NC_PATH}')
    data = load_netcdf(str(NC_PATH))

    save_path = None
    if args.save:
        plots_dir = ROOT / 'plots'
        plots_dir.mkdir(exist_ok=True)
        save_path = str(plots_dir / 'summary.png')

    plot_summary(data, save_path=save_path)
    plt.show()


if __name__ == '__main__':
    main()
