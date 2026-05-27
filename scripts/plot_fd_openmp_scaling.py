"""
Measure and plot OpenMP scaling of the Fortran fundamental-diagram sweep.

Runs ``build/fd_sweep TASEP`` at a fixed problem size for a range of
``OMP_NUM_THREADS`` values, records wall time, and plots both runtime and
speedup vs thread count against the ideal linear curve. Also verifies that
the NetCDF output is byte-identical across thread counts (per-iteration
RNG reseeding inside ``fundamental_diagram_mod`` guarantees this) and
annotates the figure with the parallel efficiency at the largest thread
count.

Usage:
    python scripts/plot_fd_openmp_scaling.py
    python scripts/plot_fd_openmp_scaling.py --save
    python scripts/plot_fd_openmp_scaling.py --L 400 --points 30 --max-threads 8
    python scripts/plot_fd_openmp_scaling.py --model NS --L 200
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np


ROOT = Path(__file__).resolve().parent.parent
FD_SWEEP = ROOT / "build" / "fd_sweep"


def run_sweep(model: str, L: int, n_points: int, n_steps: int,
              v_max: int, p_slow: float, seed: int, n_threads: int,
              out_path: Path) -> float:
    """Run one sweep with the requested thread count; return wall time (s)."""
    env = {**os.environ, "OMP_NUM_THREADS": str(n_threads)}
    env.setdefault("SDKROOT", "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk")
    cmd = [
        str(FD_SWEEP), model,
        str(L), str(n_points), str(n_steps),
        str(v_max), str(p_slow),
        str(out_path), str(seed),
    ]
    t0 = time.perf_counter()
    res = subprocess.run(cmd, env=env, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    if res.returncode != 0:
        raise RuntimeError(
            f"fd_sweep failed (exit {res.returncode}):\n"
            f"stdout:\n{res.stdout}\nstderr:\n{res.stderr}"
        )
    return elapsed


def load(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with nc.Dataset(path, "r") as ds:
        return (np.array(ds.variables["rho"][:]),
                np.array(ds.variables["J"][:]))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", choices=["TASEP", "NS"], default="TASEP",
                        help="sweep model (default TASEP — has 2 branches so "
                             "more parallel work)")
    parser.add_argument("--L", type=int, default=300, help="lattice size (default 300)")
    parser.add_argument("--points", type=int, default=20,
                        help="points per branch (default 20)")
    parser.add_argument("--n-steps", type=int, default=3000,
                        help="measurement window (default 3000)")
    parser.add_argument("--v-max", type=int, default=5, help="NS v_max (default 5)")
    parser.add_argument("--p-slow", type=float, default=0.2,
                        help="NS p_slow (default 0.2)")
    parser.add_argument("--seed", type=int, default=42,
                        help="RNG base seed (default 42) — kept fixed across "
                             "thread counts so output is bit-identical")
    parser.add_argument("--max-threads", type=int, default=None,
                        help="largest OMP_NUM_THREADS to test (default: "
                             "os.cpu_count())")
    parser.add_argument("--reps", type=int, default=3,
                        help="repeats per thread count, best-of taken (default 3)")
    parser.add_argument("--save", action="store_true",
                        help="save figure to plots/fd_openmp_scaling.png")
    args = parser.parse_args()

    if not FD_SWEEP.exists():
        sys.exit(f"fd_sweep binary not found at {FD_SWEEP}; run `make fd_sweep` first")

    max_threads = args.max_threads or (os.cpu_count() or 1)
    # 1, 2, 4, 8, 16, … up to the cap; clip the final entry to the cap exactly.
    thread_counts = []
    n = 1
    while n < max_threads:
        thread_counts.append(n)
        n *= 2
    if not thread_counts or thread_counts[-1] != max_threads:
        thread_counts.append(max_threads)

    tmp_dir = Path("/tmp") / "fd_scaling"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    print(f"Scaling sweep: model={args.model} L={args.L} points={args.points} "
          f"n_steps={args.n_steps} seed={args.seed}")
    print(f"Thread counts: {thread_counts}  (reps={args.reps}, best-of taken)")
    print()
    print(f"{'threads':>8}  {'time (s)':>10}  {'speedup':>8}  {'ideal':>6}  "
          f"{'efficiency':>10}")

    # Warm-up: pay the binary-load cost once so the first timed run isn't biased.
    run_sweep(args.model, args.L, args.points, args.n_steps,
              args.v_max, args.p_slow, args.seed, 1,
              tmp_dir / f"warm.nc")

    times = []
    nc_paths = []
    for nthr in thread_counts:
        out_path = tmp_dir / f"fd_{args.model}_{nthr}.nc"
        runs = []
        for r in range(args.reps):
            t = run_sweep(args.model, args.L, args.points, args.n_steps,
                          args.v_max, args.p_slow, args.seed, nthr, out_path)
            runs.append(t)
        best = min(runs)
        times.append(best)
        nc_paths.append(out_path)
        speedup = times[0] / best
        efficiency = speedup / nthr
        print(f"{nthr:>8d}  {best:>10.3f}  {speedup:>7.2f}x  {nthr:>5d}x  "
              f"{efficiency * 100:>9.1f}%")

    # Verify bit-identical output: per-iteration RNG seeding means thread
    # scheduling cannot affect the result.  This is a correctness check, not
    # a performance one.
    print()
    rho_ref, J_ref = load(nc_paths[0])
    all_match = True
    for nthr, path in zip(thread_counts[1:], nc_paths[1:]):
        rho, J = load(path)
        if not (np.array_equal(rho_ref, rho) and np.array_equal(J_ref, J)):
            print(f"  WARNING: output differs from serial at OMP_NUM_THREADS={nthr}")
            all_match = False
    if all_match:
        print("  All thread counts produced byte-identical NetCDF output ✓")

    # ----- Plot -----
    fig, (ax_time, ax_speed) = plt.subplots(1, 2, figsize=(11, 4.5))

    threads_arr = np.array(thread_counts)
    times_arr   = np.array(times)
    speedup_arr = times_arr[0] / times_arr
    ideal       = threads_arr.astype(float)

    ax_time.plot(threads_arr, times_arr, "o-", color="steelblue", linewidth=2,
                 markersize=7, label="measured")
    ax_time.plot(threads_arr, times_arr[0] / threads_arr, "--",
                 color="grey", linewidth=1.2, label="ideal (T₁ / N)")
    ax_time.set_xscale("log", base=2)
    ax_time.set_yscale("log")
    ax_time.set_xlabel("OMP_NUM_THREADS")
    ax_time.set_ylabel("wall time (s)")
    ax_time.set_title("Runtime")
    ax_time.set_xticks(threads_arr)
    ax_time.set_xticklabels(threads_arr)
    ax_time.grid(True, which="both", alpha=0.3)
    ax_time.legend()

    ax_speed.plot(threads_arr, speedup_arr, "o-", color="firebrick",
                  linewidth=2, markersize=7, label="measured")
    ax_speed.plot(threads_arr, ideal, "--", color="grey", linewidth=1.2,
                  label="ideal (= N)")
    ax_speed.set_xlabel("OMP_NUM_THREADS")
    ax_speed.set_ylabel("speedup  S = T₁ / Tₙ")
    ax_speed.set_title("Speedup")
    ax_speed.set_xticks(threads_arr)
    ax_speed.set_xticklabels(threads_arr)
    ax_speed.grid(True, alpha=0.3)
    ax_speed.legend()

    # Annotate efficiency at the largest thread count.
    final_eff = speedup_arr[-1] / threads_arr[-1] * 100
    ax_speed.annotate(
        f"efficiency @ {threads_arr[-1]} threads: {final_eff:.0f}%",
        xy=(threads_arr[-1], speedup_arr[-1]),
        xytext=(threads_arr[-1] * 0.6, speedup_arr[-1] * 0.65),
        fontsize=9,
        arrowprops=dict(arrowstyle="->", color="firebrick", alpha=0.6),
    )

    fig.suptitle(
        f"fd_sweep OpenMP scaling — {args.model}, L={args.L}, "
        f"{args.points} points/branch, n_steps={args.n_steps}",
        fontsize=11,
    )
    fig.tight_layout()

    if args.save:
        plots_dir = ROOT / "plots"
        plots_dir.mkdir(exist_ok=True)
        path = plots_dir / "fd_openmp_scaling.png"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        print(f"\nSaved to {path}")

    plt.show()


if __name__ == "__main__":
    main()
