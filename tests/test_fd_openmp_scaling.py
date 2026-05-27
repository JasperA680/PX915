"""Scaling test for the OpenMP-parallelised fundamental-diagram sweep.

The sweep loops over independent (alpha, beta) or density points and each
iteration is a self-contained steady-state measurement. They are wrapped
in !$omp parallel do with per-iteration RNG seeding (see
fundamental_diagram_mod::seed_iter_rng), so two things must hold:

1. **Correctness.** Output must be bit-identical for the same seed
   regardless of OMP_NUM_THREADS, because each iteration owns a
   deterministic RNG stream that does not depend on which thread runs it.

2. **Scaling.** Wall-clock time must drop measurably as the thread count
   rises. We require >= 1.5x with at least 2 cores on the machine — well
   below ideal so the test is robust to CI / busy-laptop noise.

The test is skipped if:
* the fd_sweep binary isn't built;
* the binary wasn't compiled with -fopenmp (detected by running with
  OMP_NUM_THREADS=2 and seeing no speedup over OMP_NUM_THREADS=1, then
  bailing with skip rather than failing — the build is then trivially
  serial and there is nothing to test);
* the host has < 2 CPUs reported by os.cpu_count().
"""

from __future__ import annotations

import os
import subprocess
import time
from pathlib import Path

import numpy as np
import pytest

netCDF4 = pytest.importorskip("netCDF4")


REPO_ROOT = Path(__file__).resolve().parents[1]
FD_SWEEP = REPO_ROOT / "build" / "fd_sweep"

# Problem size: chosen so the serial run takes ~1-2 s on a modern laptop,
# giving enough wall time for thread-count differences to rise above noise
# but small enough that the test stays comfortably under a few seconds at
# the parallel-end. Bumping L or n_points scales the runtime as L^2.
L = 300
N_POINTS = 20
N_STEPS = 3000
V_MAX = 5
P_SLOW = 0.2
SEED = 42

# Speedup we require with --threads-target cores. Loose enough to survive
# CI / busy hosts; tight enough that a missing -fopenmp flag (which would
# leave the loop serial) fails the test.
SPEEDUP_MIN = 1.5


def _run_sweep(model: str, n_threads: int, out_path: Path) -> float:
    """Run one fd_sweep with OMP_NUM_THREADS=n_threads, return wall-time seconds."""
    env = {**os.environ, "OMP_NUM_THREADS": str(n_threads)}
    # SDKROOT is needed on the user's Mac so gfortran's link step finds libm
    # under the current Command Line Tools SDK. Setting it on the test
    # subprocess is harmless on Linux / other Macs that don't need it.
    env.setdefault("SDKROOT", "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk")
    cmd = [
        str(FD_SWEEP), model,
        str(L), str(N_POINTS), str(N_STEPS),
        str(V_MAX), str(P_SLOW),
        str(out_path), str(SEED),
    ]
    t0 = time.perf_counter()
    result = subprocess.run(cmd, env=env, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    if result.returncode != 0:
        raise RuntimeError(
            f"fd_sweep failed (exit {result.returncode}):\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return elapsed


def _load(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with netCDF4.Dataset(path, "r") as ds:
        return (
            np.array(ds.variables["rho"][:], dtype=np.float32),
            np.array(ds.variables["J"][:],   dtype=np.float32),
        )


@pytest.fixture(scope="module")
def target_threads() -> int:
    """Number of threads to compare against the serial run."""
    if not FD_SWEEP.exists():
        pytest.skip(f"fd_sweep binary not built at {FD_SWEEP}; run `make fd_sweep`")
    cpu = os.cpu_count() or 1
    if cpu < 2:
        pytest.skip(f"need >= 2 CPUs for scaling test, host reports {cpu}")
    # Cap at 4 — beyond that, contention / hyperthreading noise grows and
    # the test threshold becomes less meaningful.
    return min(4, cpu)


@pytest.mark.parametrize("model", ["TASEP", "NS"])
def test_output_bit_identical_across_threads(model, tmp_path, target_threads):
    """Same seed + different OMP_NUM_THREADS must give exactly the same NetCDF.

    The OpenMP sweep re-seeds the RNG state per iteration with a deterministic
    function of (seed, i, branch), so the iteration's output is independent
    of which thread executed it.
    """
    out1 = tmp_path / f"fd_{model}_1.nc"
    out_n = tmp_path / f"fd_{model}_{target_threads}.nc"
    _run_sweep(model, n_threads=1, out_path=out1)
    _run_sweep(model, n_threads=target_threads, out_path=out_n)

    rho1, J1 = _load(out1)
    rhoN, JN = _load(out_n)
    np.testing.assert_array_equal(rho1, rhoN,
        err_msg=f"{model}: rho differs between 1-thread and {target_threads}-thread runs")
    np.testing.assert_array_equal(J1, JN,
        err_msg=f"{model}: J differs between 1-thread and {target_threads}-thread runs")


def test_tasep_sweep_scales_with_threads(tmp_path, target_threads, capsys):
    """Parallel sweep must be measurably faster than serial.

    Uses TASEP because it has 2 * N_POINTS independent measurements
    (alpha and beta branches), giving more parallel work than NS's single
    density sweep. We require >= 1.5x speedup at min(4, cpu_count) threads.
    """
    out1 = tmp_path / "fd_tasep_serial.nc"
    out_n = tmp_path / f"fd_tasep_{target_threads}.nc"

    # Warm caches with a throwaway run so the timed runs start from
    # comparable state. Without this the serial run pays the
    # binary-load cost and looks artificially slow.
    _run_sweep("TASEP", n_threads=1, out_path=out1)

    t_serial = _run_sweep("TASEP", n_threads=1, out_path=out1)
    t_par    = _run_sweep("TASEP", n_threads=target_threads, out_path=out_n)
    speedup = t_serial / t_par

    # Print to the captured log so `pytest -s` shows the numbers.
    with capsys.disabled():
        print(f"\n  TASEP L={L} n_points={N_POINTS}: "
              f"serial={t_serial:.2f}s, "
              f"{target_threads}-thread={t_par:.2f}s, "
              f"speedup={speedup:.2f}x")

    assert speedup >= SPEEDUP_MIN, (
        f"TASEP sweep speedup with {target_threads} threads was only "
        f"{speedup:.2f}x ({t_serial:.2f}s -> {t_par:.2f}s); expected "
        f">= {SPEEDUP_MIN}x. Likely causes: -fopenmp not on the fd_sweep "
        f"build line in the Makefile, or the OpenMP runtime isn't being "
        f"loaded at runtime (check `otool -L build/fd_sweep | grep gomp` "
        f"on macOS or `ldd build/fd_sweep | grep gomp` on Linux)."
    )
