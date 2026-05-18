"""Run the LWR PDE solver and display a summary figure.

Usage
-----
    python scripts/run_pde_model.py                              # build, run, plot
    python scripts/run_pde_model.py --no-run                     # replot existing .nc
    python scripts/run_pde_model.py --save                       # save to plots/pde_summary.png
    python scripts/run_pde_model.py --save shock.png             # save to plots/shock.png
    python scripts/run_pde_model.py --ic gaussian --flux godunov
    python scripts/run_pde_model.py --ic riemann --rho-left 0.8 --rho-right 0.2 --M 400

Presets (set rho-left and rho-right automatically)
---------------------------------------------------
    python scripts/run_pde_model.py --preset shock
    python scripts/run_pde_model.py --preset rarefaction
    python scripts/run_pde_model.py --preset shock-critical
    python scripts/run_pde_model.py --preset sonic-rarefaction
    python scripts/run_pde_model.py --preset sonic-rarefaction --flux godunov --save
"""

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / 'src' / 'python'))

import matplotlib.pyplot as plt
from pde_runner import run_pde, load_pde_netcdf
from pde_visualisation import plot_pde_summary

DEFAULT_NC = ROOT / 'data' / 'output' / 'pde_simulation.nc'

# Canonical Riemann states and default save names for each wave type
PRESETS = {
    'shock': dict(
        rho_left=0.1, rho_right=0.7,
        description='rightward-moving shock  s = +0.2  (both sides of ρ_c)',
        default_save='shock.png',
    ),
    'rarefaction': dict(
        rho_left=0.7, rho_right=0.3,
        description='rarefaction fan  (both supercritical)',
        default_save='rarefaction.png',
    ),
    'shock-critical': dict(
        rho_left=0.2, rho_right=0.8,
        description='shock through critical  (ρ_L < ρ_c < ρ_R)',
        default_save='shock_critical.png',
    ),
    'sonic-rarefaction': dict(
        rho_left=0.8, rho_right=0.2,
        description='sonic rarefaction  (ρ_R < ρ_c < ρ_L, fan passes through ρ_c)',
        default_save='sonic_rarefaction.png',
    ),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--preset', choices=PRESETS.keys(), default=None,
                        help='Canonical wave-type preset — sets rho-left and rho-right '
                             'automatically and picks a matching default save name')
    parser.add_argument('--no-run', action='store_true',
                        help='Skip build and run; just replot an existing NetCDF file')
    parser.add_argument('--save', metavar='FILENAME', nargs='?', const='',
                        help='Save figure to plots/FILENAME. '
                             'With --preset, omitting FILENAME uses the preset default name; '
                             'without --preset it falls back to pde_summary.png')
    parser.add_argument('--M', type=int, default=200,
                        help='Number of spatial cells (default: 200)')
    parser.add_argument('--steps', type=int, default=500,
                        help='Number of time steps (default: 500)')
    parser.add_argument('--ic', default='riemann',
                        choices=['constant', 'riemann', 'gaussian', 'sine'],
                        help='Initial condition (default: riemann)')
    parser.add_argument('--flux', default='lf',
                        choices=['lf', 'godunov', 'newell'],
                        help='Numerical flux scheme (default: lf)')
    parser.add_argument('--bc', default='open',
                        choices=['open', 'periodic', 'sponge'],
                        help='Boundary condition type (default: open)')
    parser.add_argument('--rho-left', type=float, default=None,
                        help='Left BC / left Riemann state (default: 0.1, or preset value)')
    parser.add_argument('--rho-right', type=float, default=None,
                        help='Right BC / right Riemann state (default: 0.9, or preset value)')
    parser.add_argument('--output', type=Path, default=DEFAULT_NC,
                        help='NetCDF output path')
    args = parser.parse_args()

    # Apply preset, then let explicit --rho-left/--rho-right override it
    if args.preset:
        p = PRESETS[args.preset]
        rho_left  = args.rho_left  if args.rho_left  is not None else p['rho_left']
        rho_right = args.rho_right if args.rho_right is not None else p['rho_right']
        print(f'Preset: {args.preset} — {p["description"]}')
    else:
        rho_left  = args.rho_left  if args.rho_left  is not None else 0.1
        rho_right = args.rho_right if args.rho_right is not None else 0.9

    # Resolve save path
    if args.save is None:
        save_path = None
    elif args.save == '':
        # --save with no filename: use preset default or fallback
        name = PRESETS[args.preset]['default_save'] if args.preset else 'pde_summary.png'
        save_path = ROOT / 'plots' / name
    else:
        save_path = ROOT / 'plots' / args.save

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
                rho_left_bc=rho_left,
                rho_right_bc=rho_right,
            ),
            output_path=args.output,
        )

    print(f'Loading {args.output}')
    data = load_pde_netcdf(args.output)
    d = data['density']
    print(f'  shape : {d.shape}  (time × space)')
    print(f'  ρ range: [{d.min():.4f}, {d.max():.4f}]')

    plot_pde_summary(data, save_path=save_path)
    plt.show()


if __name__ == '__main__':
    main()
