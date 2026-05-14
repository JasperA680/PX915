"""
Plot wallclock scaling of the serial 1D TASEP from benchmark_tasep output.

Usage:
    python scripts/plot_benchmark.py
    python scripts/plot_benchmark.py --input path/to/benchmark.nc
"""

import argparse
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent


def load_netcdf(path):
    with nc.Dataset(path) as ds:
        L_values  = ds.variables['L_values'][:]
        N_values  = ds.variables['N_values'][:]
        wallclock = ds.variables['wallclock'][:]
        rho       = ds.variables['steady_density'][:]
        J         = ds.variables['steady_current'][:]
        alpha     = float(ds.getncattr('alpha'))
        beta      = float(ds.getncattr('beta'))
    # NetCDF (C-order) reverses Fortran's (n_L, n_N, n_rep) into (n_rep, n_N, n_L);
    # transpose back to (n_L, n_N, n_rep) for natural indexing.
    wc = np.asarray(wallclock, dtype=float).transpose(2, 1, 0)
    wc = np.where(wc < 0, np.nan, wc)
    return (np.asarray(L_values), np.asarray(N_values),
            wc, np.asarray(rho), np.asarray(J), alpha, beta)


def fit_log_slope(x, y):
    mask = np.isfinite(y) & (y > 0)
    if mask.sum() < 2:
        return np.nan
    slope, _ = np.polyfit(np.log10(x[mask]), np.log10(y[mask]), 1)
    return slope


def plot_vs_L(L_values, N_values, mean, std, out_path):
    fig, ax = plt.subplots(figsize=(7, 5))
    for j, N in enumerate(N_values):
        y = mean[:, j]
        yerr = std[:, j]
        mask = np.isfinite(y)
        if mask.sum() == 0:
            continue
        ax.errorbar(L_values[mask], y[mask], yerr=yerr[mask],
                    marker='o', label=f'N = {N:,}')
        slope = fit_log_slope(L_values, y)
        print(f'  vs L  (N={N:>8}): slope = {slope:.3f}')

    ref_L = np.array([L_values.min(), L_values.max()], dtype=float)
    ref = ref_L * (mean[mean > 0].min() / ref_L[0]) * 0.3
    ax.plot(ref_L, ref, 'k--', alpha=0.4, label='slope = 1 (reference)')

    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('Lattice length L')
    ax.set_ylabel('Wallclock per run (s)')
    ax.set_title('TASEP serial wallclock vs L')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved {out_path}')


def plot_vs_N(L_values, N_values, mean, std, out_path):
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, L in enumerate(L_values):
        y = mean[i, :]
        yerr = std[i, :]
        mask = np.isfinite(y)
        if mask.sum() == 0:
            continue
        ax.errorbar(N_values[mask], y[mask], yerr=yerr[mask],
                    marker='o', label=f'L = {L}')
        slope = fit_log_slope(N_values, y)
        print(f'  vs N  (L={L:>5}): slope = {slope:.3f}')

    ref_N = np.array([N_values.min(), N_values.max()], dtype=float)
    ref = ref_N * (mean[mean > 0].min() / ref_N[0]) * 0.3
    ax.plot(ref_N, ref, 'k--', alpha=0.4, label='slope = 1 (reference)')

    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('Number of time steps N')
    ax.set_ylabel('Wallclock per run (s)')
    ax.set_title('TASEP serial wallclock vs n_steps')
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved {out_path}')


def plot_heatmap(L_values, N_values, mean, out_path):
    fig, ax = plt.subplots(figsize=(7, 5))
    Z = np.log10(np.where(mean > 0, mean, np.nan))
    mesh = ax.pcolormesh(N_values, L_values, Z, shading='auto', cmap='viridis')
    cb = fig.colorbar(mesh, ax=ax)
    cb.set_label('log10(mean wallclock / s)')
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('n_steps'); ax.set_ylabel('L')
    ax.set_title('TASEP serial wallclock — log10 seconds')
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'  saved {out_path}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', default=str(ROOT / 'data/output/benchmark.nc'),
                        help='NetCDF from make benchmark (default: data/output/benchmark.nc)')
    args = parser.parse_args()

    L_values, N_values, wc, rho, J, alpha, beta = load_netcdf(args.input)
    print(f'Loaded {args.input}  (alpha={alpha}, beta={beta})')
    print(f'  L_values = {L_values.tolist()}')
    print(f'  N_values = {N_values.tolist()}')
    print(f'  wallclock shape = {wc.shape}  (n_L, n_N, n_rep)')

    mean = np.nanmean(wc, axis=2)
    std  = np.nanstd(wc, axis=2)

    print('\nSteady-state baseline per L:')
    for L_, rho_, J_ in zip(L_values, rho, J):
        print(f'  L={L_:>5}  rho={rho_:.4f}  J={J_:.4f}')

    plots_dir = ROOT / 'plots'
    plots_dir.mkdir(exist_ok=True)

    print('\nLog-log slope fits:')
    plot_vs_L(L_values, N_values, mean, std, plots_dir / 'benchmark_vs_L.png')
    plot_vs_N(L_values, N_values, mean, std, plots_dir / 'benchmark_vs_N.png')
    plot_heatmap(L_values, N_values, mean, plots_dir / 'benchmark_heatmap.png')


if __name__ == '__main__':
    main()
