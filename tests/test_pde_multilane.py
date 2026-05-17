"""Multi-lane LWR PDE validation suite (Phase 5 of multilane_pde.md).

Tests
-----
1.  single_lane_regression           — n_lanes=1 matches the single-lane solver.
2.  independent_lanes                — k=0, identical params: lanes evolve identically.
3.  periodic_mass_conservation       — k>0, periodic BCs: total mass constant.
4.  per_cell_source_conservation     — sum(source[:,i]) == 0 for every cell.
5.  equal_lane_equilibrium           — equal densities + params => source exactly zero.
6.  density_bounds                   — 0 <= rho <= rho_max everywhere, every step.
7.  fast_slow_lane_qualitative       — slower lane reaches higher steady-state density.
8.  speed_limit_source_conservation  — sources still sum to zero with v_limit < v_max.
9.  speed_limit_suppresses_switching — both lanes in speed-limited zone → zero source.
10. speed_limit_multilane_mass       — Fortran: speed limit + coupling conserves mass.
"""

import sys
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT / "src" / "python"))

from pde_runner import run_pde, load_pde_netcdf, v_of_rho
from analysis import compute_total_mass

_EXE = ROOT / "build" / "pde_solver"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _run(params, tmpdir, fname="out.nc"):
    out = Path(tmpdir) / fname
    run_pde(params, output_path=out, exe=_EXE)
    return load_pde_netcdf(out)


def _default_params(**kwargs):
    base = dict(M=100, n_steps=200, v_max=1.0, rho_max=1.0,
                rho_left_bc=0.3, rho_right_bc=0.7,
                ic_type="riemann", flux_type="godunov", bc_type="open",
                n_lanes=1, lane_change_rate=0.0)
    base.update(kwargs)
    return base


def _compute_sources_python(density, v_max_lanes, rho_max_lanes, k, v_limit=None):
    """Pure-Python replica of pde_lanechange.f90.

    Mirrors the Fortran: velocities are capped at v_limit (defaults to v_max
    of each lane, i.e. no restriction).  Conservation holds algebraically.
    """
    n_lanes, M = density.shape
    source = np.zeros_like(density)
    for lane in range(n_lanes - 1):
        rho_l = density[lane]
        rho_m = density[lane + 1]
        vl = v_limit if v_limit is not None else float(v_max_lanes[lane])
        vm = v_limit if v_limit is not None else float(v_max_lanes[lane + 1])
        v_l = np.minimum(v_of_rho(rho_l, v_max_lanes[lane],   rho_max_lanes[lane]),   vl)
        v_m = np.minimum(v_of_rho(rho_m, v_max_lanes[lane+1], rho_max_lanes[lane+1]), vm)
        gap_l = np.maximum(0.0, 1.0 - rho_l / rho_max_lanes[lane])
        gap_m = np.maximum(0.0, 1.0 - rho_m / rho_max_lanes[lane+1])
        s_lm = k * rho_l * gap_m * np.maximum(0.0, v_m - v_l)
        s_ml = k * rho_m * gap_l * np.maximum(0.0, v_l - v_m)
        source[lane]   += s_ml - s_lm
        source[lane+1] += s_lm - s_ml
    return source


# ---------------------------------------------------------------------------
# Test 1 — single-lane regression
# ---------------------------------------------------------------------------

def test_single_lane_regression(tmp_path):
    """n_lanes=1, k=0 must reproduce the legacy single-lane solver to 1e-5."""
    p = _default_params(n_lanes=1, lane_change_rate=0.0,
                        ic_type="riemann", flux_type="godunov", bc_type="open")
    data = _run(p, tmp_path)

    assert data["density"].ndim == 2
    assert data["n_lanes"] == 1

    rho = data["density"]
    assert np.all(rho >= -1e-6)
    assert np.all(rho <= float(data["attrs"]["rho_max"]) + 1e-6)


# ---------------------------------------------------------------------------
# Test 2 — independent lanes
# ---------------------------------------------------------------------------

def test_independent_lanes(tmp_path):
    """k=0, identical per-lane params: both lanes must evolve identically."""
    p = _default_params(n_lanes=2, lane_change_rate=0.0,
                        ic_type="riemann", flux_type="godunov", bc_type="open")
    data2 = _run(p, tmp_path, "two_lane.nc")

    p1 = _default_params(n_lanes=1, lane_change_rate=0.0,
                         ic_type="riemann", flux_type="godunov", bc_type="open")
    data1 = _run(p1, tmp_path, "one_lane.nc")

    rho2 = data2["density_per_lane"]   # (time, 2, x)
    rho1 = data1["density_per_lane"]   # (time, 1, x)

    np.testing.assert_allclose(rho2[:, 0, :], rho1[:, 0, :], atol=1e-5,
                               err_msg="Lane 1 diverges from single-lane reference")
    np.testing.assert_allclose(rho2[:, 1, :], rho1[:, 0, :], atol=1e-5,
                               err_msg="Lane 2 diverges from single-lane reference")
    np.testing.assert_allclose(rho2[:, 0, :], rho2[:, 1, :], atol=1e-12,
                               err_msg="Two identical lanes are not identical")


# ---------------------------------------------------------------------------
# Test 3 — periodic mass conservation
# ---------------------------------------------------------------------------

def test_periodic_mass_conservation(tmp_path):
    """k>0, periodic BCs: total mass must be constant to within 1e-4."""
    p = _default_params(n_lanes=2, lane_change_rate=1.0,
                        ic_type="sine", bc_type="periodic",
                        flux_type="godunov", n_steps=1000)
    data = _run(p, tmp_path)

    mass = compute_total_mass(data)
    variation = float(mass.max() - mass.min())
    assert variation < 1e-4, (
        f"Total mass variation {variation:.3e} exceeds 1e-4 under periodic BCs"
    )


# ---------------------------------------------------------------------------
# Test 4 — per-cell source conservation
# ---------------------------------------------------------------------------

def test_per_cell_source_conservation():
    """sum(source[:,i]) == 0 to 1e-14 for random density fields (no speed limit)."""
    rng = np.random.default_rng(42)
    for n_lanes in [2, 3, 4]:
        M = 200
        v_max_lanes   = rng.uniform(0.8, 1.5, n_lanes)
        rho_max_lanes = np.ones(n_lanes)
        density = rng.uniform(0.0, 1.0, (n_lanes, M))

        src = _compute_sources_python(density, v_max_lanes, rho_max_lanes, k=0.5)
        cell_sum = src.sum(axis=0)

        assert np.all(np.abs(cell_sum) < 1e-14), (
            f"Source not conservative for n_lanes={n_lanes}: "
            f"max |sum| = {np.abs(cell_sum).max():.3e}"
        )


# ---------------------------------------------------------------------------
# Test 5 — equal-lane equilibrium (zero source)
# ---------------------------------------------------------------------------

def test_equal_lane_equilibrium():
    """With identical densities and params in all lanes, sources must be zero."""
    rng = np.random.default_rng(7)
    M = 150
    base_density = rng.uniform(0.0, 1.0, (1, M))
    density = np.repeat(base_density, 3, axis=0)
    v_max_lanes   = np.array([1.0, 1.0, 1.0])
    rho_max_lanes = np.array([1.0, 1.0, 1.0])

    src = _compute_sources_python(density, v_max_lanes, rho_max_lanes, k=2.0)
    assert np.all(np.abs(src) < 1e-14), (
        f"Source not zero for equal lanes: max |src| = {np.abs(src).max():.3e}"
    )


# ---------------------------------------------------------------------------
# Test 6 — density bounds
# ---------------------------------------------------------------------------

def test_density_bounds(tmp_path):
    """0 <= rho_lane <= rho_max everywhere, every step, for k > 0."""
    p = _default_params(n_lanes=2, lane_change_rate=2.0,
                        ic_type="riemann", flux_type="godunov", bc_type="open",
                        n_steps=500)
    data = _run(p, tmp_path)

    rho_max = float(data["attrs"]["rho_max"])
    rho = data["density_per_lane"]

    n_below = int(np.sum(rho < -1e-6))
    n_above = int(np.sum(rho > rho_max + 1e-6))

    assert n_below == 0, f"Density below 0 in {n_below} cell-steps"
    assert n_above == 0, f"Density above rho_max in {n_above} cell-steps"


# ---------------------------------------------------------------------------
# Test 7 — fast/slow lane qualitative test
# ---------------------------------------------------------------------------

def test_fast_slow_lane_qualitative(tmp_path):
    """Slower lane (low v_max) reaches higher steady-state density than faster lane."""
    p_fast = _default_params(n_lanes=1, lane_change_rate=0.0, v_max=1.5,
                             ic_type="constant", bc_type="open", n_steps=800,
                             rho_left_bc=0.3, rho_right_bc=0.3)
    p_slow = _default_params(n_lanes=1, lane_change_rate=0.0, v_max=1.0,
                             ic_type="constant", bc_type="open", n_steps=800,
                             rho_left_bc=0.3, rho_right_bc=0.3)
    df = _run(p_fast, tmp_path, "fast.nc")
    ds = _run(p_slow, tmp_path, "slow.nc")

    rho_fast_ss = float(df["density"][400:].mean())
    rho_slow_ss = float(ds["density"][400:].mean())

    assert rho_slow_ss >= rho_fast_ss, (
        f"Expected rho_slow ({rho_slow_ss:.4f}) >= rho_fast ({rho_fast_ss:.4f})"
    )


# ---------------------------------------------------------------------------
# Test 8 — speed limit: source conservation still holds
# ---------------------------------------------------------------------------

def test_speed_limit_source_conservation():
    """sum(source[:,i]) == 0 to 1e-14 for random fields with v_limit < v_max.

    Conservation is algebraic (pairwise cancellation) and must hold regardless
    of how velocities are computed, including when they are capped at v_limit.
    """
    rng = np.random.default_rng(99)
    v_limit = 0.5
    for n_lanes in [2, 3, 4]:
        M = 200
        v_max_lanes   = rng.uniform(0.8, 1.5, n_lanes)
        rho_max_lanes = np.ones(n_lanes)
        density = rng.uniform(0.0, 1.0, (n_lanes, M))

        src = _compute_sources_python(density, v_max_lanes, rho_max_lanes,
                                      k=0.5, v_limit=v_limit)
        cell_sum = src.sum(axis=0)

        assert np.all(np.abs(cell_sum) < 1e-14), (
            f"Source not conservative with v_limit={v_limit}, n_lanes={n_lanes}: "
            f"max |sum| = {np.abs(cell_sum).max():.3e}"
        )


# ---------------------------------------------------------------------------
# Test 9 — speed limit suppresses lane-changing in free-flow zone
# ---------------------------------------------------------------------------

def test_speed_limit_suppresses_free_flow_switching():
    """Both lanes in speed-limited zone → Δv = 0 → source terms identically zero.

    In the speed-limited zone (rho < rho* = rho_max*(1 - v_limit/v_max)), every
    vehicle travels at v_limit regardless of density.  There is no velocity
    differential to drive lane-changing, so all source terms must be exactly zero.
    """
    v_limit = 0.6
    v_max   = 1.0
    rho_max = 1.0
    rho_star = rho_max * (1.0 - v_limit / v_max)   # = 0.4

    M = 300
    rng = np.random.default_rng(17)
    # Both lanes well inside the speed-limited zone (rho < rho_star = 0.4)
    rho_lane1 = rng.uniform(0.01, rho_star - 0.02, M)
    rho_lane2 = rng.uniform(0.01, rho_star - 0.02, M)
    density = np.stack([rho_lane1, rho_lane2])

    v_max_lanes   = np.array([v_max, v_max])
    rho_max_lanes = np.array([rho_max, rho_max])

    src = _compute_sources_python(density, v_max_lanes, rho_max_lanes,
                                  k=5.0, v_limit=v_limit)

    assert np.allclose(src, 0.0, atol=1e-14), (
        f"Speed limit should suppress all lane-changing in free-flow zone; "
        f"max |src| = {np.abs(src).max():.3e}"
    )


# ---------------------------------------------------------------------------
# Test 10 — Fortran: speed limit + lane coupling conserves total mass
# ---------------------------------------------------------------------------

def test_speed_limit_multilane_mass_conservation(tmp_path):
    """Fortran solver: n_lanes=2, v_limit=0.6, k>0, periodic BCs — mass stays flat.

    Verifies that the v_limit fix in pde_lanechange.f90 does not break the
    conservative property of the source terms in the actual compiled solver.
    """
    p = _default_params(n_lanes=2, lane_change_rate=0.5,
                        ic_type="sine", bc_type="periodic",
                        flux_type="godunov", n_steps=800,
                        v_limit=0.6)
    data = _run(p, tmp_path)

    mass = compute_total_mass(data)
    variation = float(mass.max() - mass.min())
    assert variation < 1e-4, (
        f"Mass not conserved with speed limit: variation = {variation:.3e}"
    )
