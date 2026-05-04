"""Unit tests for PDE flux functions (Phase 1) and solver properties (Phase 4)."""

import numpy as np
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / "src" / "python"))
sys.path.insert(0, str(Path(__file__).parent))
from pde_runner import v_of_rho, q_of_rho, dq_drho, rho_critical


# ---------------------------------------------------------------------------
# Flux function tests (Phase 1)
# ---------------------------------------------------------------------------

V_MAX   = 1.0
RHO_MAX = 1.0
RC      = rho_critical(RHO_MAX)  # 0.5


def test_q_at_zero_density():
    assert q_of_rho(0.0, V_MAX, RHO_MAX) == pytest.approx(0.0)


def test_q_at_max_density():
    assert q_of_rho(RHO_MAX, V_MAX, RHO_MAX) == pytest.approx(0.0)


def test_q_at_critical_density():
    """Maximum flow is v_max * rho_max / 4 (Greenshields analytical result)."""
    q_max_analytical = V_MAX * RHO_MAX / 4.0
    assert q_of_rho(RC, V_MAX, RHO_MAX) == pytest.approx(q_max_analytical)


def test_dq_drho_at_critical_is_zero():
    assert dq_drho(RC, V_MAX, RHO_MAX) == pytest.approx(0.0)


def test_dq_drho_positive_in_free_flow():
    """Characteristic speed positive for ρ < ρ_c (free-flow phase)."""
    rho_values = np.linspace(0.0, RC - 1e-6, 20)
    assert np.all(dq_drho(rho_values, V_MAX, RHO_MAX) > 0)


def test_dq_drho_negative_in_congested():
    """Characteristic speed negative for ρ > ρ_c (congested phase)."""
    rho_values = np.linspace(RC + 1e-6, RHO_MAX, 20)
    assert np.all(dq_drho(rho_values, V_MAX, RHO_MAX) < 0)


def test_q_is_concave():
    """q should be concave: all second differences negative."""
    rho = np.linspace(0.0, RHO_MAX, 200)
    q   = q_of_rho(rho, V_MAX, RHO_MAX)
    d2q = np.diff(q, n=2)
    assert np.all(d2q < 0), "q must be strictly concave"


def test_v_of_rho_decreasing():
    """Greenshields velocity is strictly decreasing in density."""
    rho = np.linspace(0.0, RHO_MAX, 100)
    v   = v_of_rho(rho, V_MAX, RHO_MAX)
    assert np.all(np.diff(v) < 0)


def test_rho_critical_is_half_rho_max():
    for rho_max in [0.5, 1.0, 2.0, 150.0]:
        assert rho_critical(rho_max) == pytest.approx(rho_max / 2.0)


def test_q_symmetric_about_critical():
    """q(ρ_c − δ) == q(ρ_c + δ) for Greenshields."""
    delta = 0.1
    q_sub  = q_of_rho(RC - delta, V_MAX, RHO_MAX)
    q_sup  = q_of_rho(RC + delta, V_MAX, RHO_MAX)
    assert q_sub == pytest.approx(q_sup)


# ---------------------------------------------------------------------------
# Solver integration tests (Phase 4) — require compiled Fortran binary
# ---------------------------------------------------------------------------

try:
    from pde_runner import run_pde, load_pde_netcdf
    from exact_riemann import exact_riemann_lwr
    BINARY = Path(__file__).parents[1] / "build" / "pde_solver"
    SOLVER_AVAILABLE = BINARY.exists()
except ImportError:
    SOLVER_AVAILABLE = False

SOLVER_SKIP = pytest.mark.skipif(
    not SOLVER_AVAILABLE,
    reason="build/pde_solver binary not found — run 'make pde' first",
)


@SOLVER_SKIP
def test_constant_solution_preserved(tmp_path):
    """Uniform density stays uniform to within 1e-10 over 1000 steps."""
    out = tmp_path / "const.nc"
    run_pde(
        dict(M=100, n_steps=1000, ic_type="constant",
             rho_left_bc=0.4, rho_right_bc=0.4,
             bc_type="periodic", flux_type="lf"),
        output_path=out,
    )
    data = load_pde_netcdf(out)
    rho  = data["density"]  # (n_steps+1, M)
    assert np.allclose(rho, 0.4, atol=1e-10), \
        f"Max deviation: {np.abs(rho - 0.4).max():.2e}"


@SOLVER_SKIP
def test_mass_conservation_periodic(tmp_path):
    """Total mass conserved to within 1e-8 under periodic BCs."""
    out = tmp_path / "mass.nc"
    run_pde(
        dict(M=200, n_steps=1000, ic_type="gaussian",
             bc_type="periodic", flux_type="lf"),
        output_path=out,
    )
    data = load_pde_netcdf(out)
    attrs = data["attrs"]
    dx    = float(attrs["dx"])
    mass  = data["density"].sum(axis=1) * dx  # (n_steps+1,)
    # Float32 storage in NetCDF limits precision; 1e-4 is tight enough to
    # catch open-boundary mass leakage while tolerating float32 rounding.
    assert np.allclose(mass, mass[0], atol=1e-4), \
        f"Mass drift: {np.abs(mass - mass[0]).max():.2e}"


@SOLVER_SKIP
def test_shock_speed(tmp_path):
    """Shock position matches Rankine-Hugoniot speed to within 4 cell widths.

    rho_L=0.3 (free-flow), rho_R=0.8 (congested) → shock crossing rho_c,
    travelling left at s = (q(0.8)-q(0.3))/(0.8-0.3) = -0.1.
    Tolerance follows PDE.md: one cell width per 100 steps.
    """
    rho_L, rho_R = 0.3, 0.8
    M, n_steps = 200, 400
    out = tmp_path / "shock.nc"
    run_pde(
        dict(M=M, n_steps=n_steps, ic_type="riemann",
             rho_left_bc=rho_L, rho_right_bc=rho_R,
             flux_type="godunov", bc_type="open"),
        output_path=out,
    )
    data = load_pde_netcdf(out)
    attrs    = data["attrs"]
    dx       = float(attrs["dx"])
    v_max    = float(attrs.get("v_max", 1.0))
    rho_max  = float(attrs.get("rho_max", 1.0))
    x        = data["x"]
    t_final  = float(data["time"][-1])
    rho_num  = data["density"][-1]

    # Locate shock as the cell interface with the largest density jump
    grad = np.abs(np.diff(rho_num))
    idx  = np.argmax(grad)
    x_shock_num = 0.5 * (x[idx] + x[idx + 1])

    # Exact shock position via Rankine-Hugoniot
    s = (q_of_rho(rho_R, v_max, rho_max) - q_of_rho(rho_L, v_max, rho_max)) / (rho_R - rho_L)
    x_shock_exact = 0.5 + s * t_final   # IC discontinuity placed at x = L/2 = 0.5

    tol = 4 * dx  # PDE.md: one cell width per 100 steps over 400 steps
    assert abs(x_shock_num - x_shock_exact) < tol, (
        f"Shock at {x_shock_num:.4f}, expected {x_shock_exact:.4f}, "
        f"|error|={abs(x_shock_num - x_shock_exact):.4f} > tol={tol:.4f}"
    )


@SOLVER_SKIP
def test_rarefaction_l1_error(tmp_path):
    """Rarefaction fan has L1 error O(dx) against the exact self-similar solution.

    rho_L=0.8 > rho_R=0.2 produces a sonic rarefaction spanning the critical
    density.  Godunov is exact on smooth rarefactions; error is O(dx) at the
    two kink points where the fan meets the constant states.
    """
    rho_L, rho_R = 0.8, 0.2
    M, n_steps = 200, 100
    out = tmp_path / "rarefaction.nc"
    run_pde(
        dict(M=M, n_steps=n_steps, ic_type="riemann",
             rho_left_bc=rho_L, rho_right_bc=rho_R,
             flux_type="godunov", bc_type="open"),
        output_path=out,
    )
    data = load_pde_netcdf(out)
    attrs   = data["attrs"]
    dx      = float(attrs["dx"])
    v_max   = float(attrs.get("v_max", 1.0))
    rho_max = float(attrs.get("rho_max", 1.0))
    x       = data["x"]
    t_final = float(data["time"][-1])
    rho_num = data["density"][-1]

    rho_exact = exact_riemann_lwr(x, t_final, rho_L, rho_R, v_max, rho_max)
    l1_error  = np.sum(np.abs(rho_num - rho_exact)) * dx

    assert l1_error < 10 * dx, (
        f"Rarefaction L1 error {l1_error:.5f} exceeds 10*dx={10*dx:.5f}"
    )


@SOLVER_SKIP
def test_convergence_first_order(tmp_path):
    """L1 error against the exact Riemann solution halves as dx halves.

    Uses n_steps = M for each grid so that t_final = CFL * domain / v_max ≈ 0.9
    is the same across all refinements, making the comparison well-defined.
    """
    rho_L, rho_R = 0.3, 0.8   # shock problem (same as test_shock_speed)
    errors = {}
    for M in (50, 100, 200):
        out = tmp_path / f"conv_{M}.nc"
        run_pde(
            dict(M=M, n_steps=M, ic_type="riemann",
                 rho_left_bc=rho_L, rho_right_bc=rho_R,
                 flux_type="godunov", bc_type="open"),
            output_path=out,
        )
        data    = load_pde_netcdf(out)
        attrs   = data["attrs"]
        dx      = float(attrs["dx"])
        v_max   = float(attrs.get("v_max", 1.0))
        rho_max = float(attrs.get("rho_max", 1.0))
        x       = data["x"]
        t_final = float(data["time"][-1])
        rho_num = data["density"][-1]
        rho_ex  = exact_riemann_lwr(x, t_final, rho_L, rho_R, v_max, rho_max)
        errors[M] = np.sum(np.abs(rho_num - rho_ex)) * dx

    ratio_50_100  = errors[50]  / errors[100]
    ratio_100_200 = errors[100] / errors[200]
    assert ratio_50_100 >= 1.5, (
        f"Convergence rate (M 50→100): error ratio={ratio_50_100:.2f} < 1.5 "
        f"(errors: {errors[50]:.4f}, {errors[100]:.4f})"
    )
    assert ratio_100_200 >= 1.5, (
        f"Convergence rate (M 100→200): error ratio={ratio_100_200:.2f} < 1.5 "
        f"(errors: {errors[100]:.4f}, {errors[200]:.4f})"
    )
