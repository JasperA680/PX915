"""Compare shock behaviour under different speed limits.

Runs the LWR PDE solver once per (closure, flux, v_limit) combination and
produces a four-panel figure — one panel per (closure × flux) pair — with
fundamental diagrams and final-time density profiles overlaid for each speed
limit.

Usage
-----
    python scripts/plot_speed_limit_comparison.py
    python scripts/plot_speed_limit_comparison.py --v-limits 0.3 0.5 0.7 1.0
    python scripts/plot_speed_limit_comparison.py --preset rarefaction --save
    python scripts/plot_speed_limit_comparison.py --no-run --save comparison.png

Preset Riemann states (same as run_pde_model.py)
-------------------------------------------------
    shock            rho_left=0.1  rho_right=0.7
    rarefaction      rho_left=0.7  rho_right=0.3
    shock-critical   rho_left=0.2  rho_right=0.8
    sonic-rarefaction rho_left=0.8  rho_right=0.2
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / 'src' / 'python'))

import matplotlib.pyplot as plt
from pde_runner import run_pde, load_pde_netcdf
from visualisation import plot_speed_limit_comparison

PRESETS = {
    'shock':             dict(rho_left=0.1, rho_right=0.7),
    'rarefaction':       dict(rho_left=0.7, rho_right=0.3),
    'shock-critical':    dict(rho_left=0.2, rho_right=0.8),
    'sonic-rarefaction': dict(rho_left=0.8, rho_right=0.2),
}

# Four panels: (closure, numerical flux) combinations
PANELS = [
    ('Greenshields + LF',      'lf',      'lf'),
    ('Greenshields + Godunov', 'godunov', 'godunov'),
    ('Newell + LF',            'newell',  'lf'),
    ('Newell + Godunov',       'newell',  'godunov'),
]


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--preset', choices=PRESETS.keys(), default='shock',
                        help='Riemann initial condition preset (default: shock)')
    parser.add_argument('--rho-left',  type=float, default=None,
                        help='Override left density / left BC')
    parser.add_argument('--rho-right', type=float, default=None,
                        help='Override right density / right BC')
    parser.add_argument('--v-limits', type=float, nargs='+',
                        default=[0.4, 0.6, 0.8, 1.0],
                        help='Speed limits to compare (default: 0.4 0.6 0.8 1.0)')
    parser.add_argument('--M',     type=int, default=300,
                        help='Number of spatial cells (default: 300)')
    parser.add_argument('--steps', type=int, default=600,
                        help='Number of time steps (default: 600)')
    parser.add_argument('--bc',    default='open',
                        choices=['open', 'periodic', 'sponge'],
                        help='Boundary condition (default: open)')
    parser.add_argument('--n-lanes', type=int, default=1,
                        help='Number of lanes (default: 1). With n_lanes>1 each panel '
                             'shows total density summed over lanes.')
    parser.add_argument('--k', type=float, default=0.0,
                        help='Lane-change rate (default: 0). Only used when n_lanes>1.')
    parser.add_argument('--save', metavar='FILENAME', nargs='?', const='',
                        help='Save figure. Omit FILENAME to use default name.')
    parser.add_argument('--no-run', action='store_true',
                        help='Skip build and solver runs; expects temp files already present')
    args = parser.parse_args()

    preset = PRESETS[args.preset]
    rho_left  = args.rho_left  if args.rho_left  is not None else preset['rho_left']
    rho_right = args.rho_right if args.rho_right is not None else preset['rho_right']

    if args.save is None:
        save_path = None
    elif args.save == '':
        lanes_tag = f'_lanes{args.n_lanes}' if args.n_lanes > 1 else ''
        save_path = ROOT / 'plots' / f'speed_limit_{args.preset}{lanes_tag}.png'
    else:
        save_path = ROOT / 'plots' / args.save

    if not args.no_run:
        print('Building PDE solver...')
        subprocess.run(['make', 'pde'], check=True, cwd=ROOT)

    # Run solver for every (flux_type, v_limit) combination and collect results.
    # flux_type controls both the closure and the numerical scheme:
    #   'lf'      → Greenshields + Lax-Friedrichs
    #   'godunov' → Greenshields + Godunov
    #   'newell'  → Newell + Godunov  (handled inside Fortran dispatcher)
    # For Newell + LF we pass flux_type='newell' and the Fortran driver falls
    # through to the LF branch with the Newell q function evaluated there.
    # We encode that distinction purely at the Python level via the label.

    # Map panel spec → unique solver flux_type token understood by Fortran:
    #   Greenshields+LF     → flux_type='lf'
    #   Greenshields+Godunov→ flux_type='godunov'
    #   Newell+LF           → flux_type='newell_lf'  (new token, see note below)
    #   Newell+Godunov      → flux_type='newell'
    #
    # The Fortran pde_step dispatches on:
    #   flux_type == 'godunov' or 'newell' → Godunov branch
    #   anything else                       → Lax-Friedrichs branch
    # So 'newell_lf' correctly selects LF while q_dispatch uses Newell closure.

    with tempfile.TemporaryDirectory(dir=ROOT / 'data' / 'output') as tmpdir:
        tmp = Path(tmpdir)
        panel_datasets = []

        for label, closure_key, flux_key in PANELS:
            # Determine the flux_type string the Fortran driver understands:
            #   Newell+LF  → 'newell_lf'  (LF path, Newell q via q_dispatch)
            #   Newell+Godunov → 'newell'
            #   Greenshields+LF → 'lf'
            #   Greenshields+Godunov → 'godunov'
            if closure_key == 'newell' and flux_key == 'lf':
                fortran_flux = 'newell_lf'
            elif closure_key == 'newell':
                fortran_flux = 'newell'
            else:
                fortran_flux = flux_key

            datasets = []
            for v_limit in args.v_limits:
                nc_path = tmp / f'{fortran_flux}_vl{v_limit:.3f}.nc'
                if not args.no_run:
                    print(f'  Running {label}  v_limit={v_limit:.2f} ...')
                    run_pde(
                        params=dict(
                            M=args.M,
                            n_steps=args.steps,
                            ic_type='riemann',
                            flux_type=fortran_flux,
                            bc_type=args.bc,
                            rho_left_bc=rho_left,
                            rho_right_bc=rho_right,
                            v_limit=v_limit,
                            n_lanes=args.n_lanes,
                            lane_change_rate=args.k,
                        ),
                        output_path=nc_path,
                    )
                data = load_pde_netcdf(nc_path)
                datasets.append((v_limit, data))

            panel_datasets.append((label, datasets))

    print('Plotting...')
    plot_speed_limit_comparison(panel_datasets, save_path=save_path)
    plt.show()


if __name__ == '__main__':
    main()
