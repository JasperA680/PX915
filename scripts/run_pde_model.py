"""Run the LWR PDE solver and display a summary figure.

Usage
-----
    python scripts/run_pde_model.py                        # build, run, plot
    python scripts/run_pde_model.py --no-run               # skip build/run, replot existing .nc
    python scripts/run_pde_model.py --save                 # also save plots/pde_summary.png
    python scripts/run_pde_model.py --ic gaussian --flux godunov
    python scripts/run_pde_model.py --ic riemann --rho-left 0.8 --rho-right 0.2 --M 400
"""

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / 'src' / 'python'))

import matplotlib.pyplot as plt
from pde_runner import run_pde, load_pde_netcdf
from visualisation import plot_pde_summary

DEFAULT_NC   = ROOT / 'data' / 'output' / 'pde_simulation.nc'
DEFAULT_SAVE = ROOT / 'plots' / 'pde_summary.png'


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--no-run', action='store_true',
                        help='Skip build and run; just replot an existing NetCDF file')
    parser.add_argument('--save', action='store_true',
                        help='Save figure to plots/pde_summary.png')
    parser.add_argument('--M', type=int, default=200,
                        help='Number of spatial cells (default: 200)')
    parser.add_argument('--steps', type=int, default=500,
                        help='Number of time steps (default: 500)')
    parser.add_argument('--ic', default='riemann',
                        choices=['constant', 'riemann', 'gaussian', 'sine'],
                        help='Initial condition (default: riemann)')
    parser.add_argument('--flux', default='lf',
                        choices=['lf', 'godunov'],
                        help='Numerical flux scheme (default: lf)')
    parser.add_argument('--bc', default='open',
                        choices=['open', 'periodic'],
                        help='Boundary condition type (default: open)')
    parser.add_argument('--rho-left', type=float, default=0.1,
                        help='Left BC / left Riemann state (default: 0.1)')
    parser.add_argument('--rho-right', type=float, default=0.9,
                        help='Right BC / right Riemann state (default: 0.9)')
    parser.add_argument('--output', type=Path, default=DEFAULT_NC,
                        help='NetCDF output path')
    args = parser.parse_args()

    if not args.no_run:
        print('Building PDE solver...')
        subprocess.run(['make', 'pde'], check=True, cwd=ROOT)
        print('Running solver...')
        run_pde(
            params=dict(
                M=args.M,
                n_steps=args.steps,
                ic_type=args.ic,
                flux_type=args.flux,
                bc_type=args.bc,
                rho_left_bc=args.rho_left,
                rho_right_bc=args.rho_right,
            ),
            output_path=args.output,
        )

    print(f'Loading {args.output}')
    data = load_pde_netcdf(args.output)
    d    = data['density']
    print(f'  shape : {d.shape}  (time × space)')
    print(f'  ρ range: [{d.min():.4f}, {d.max():.4f}]')
    print(f'  attrs : { {k: v for k, v in data["attrs"].items()} }')

    save_path = DEFAULT_SAVE if args.save else None
    plot_pde_summary(data, save_path=save_path)
    plt.show()


if __name__ == '__main__':
    main()
