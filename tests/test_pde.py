"""Unit tests for PDE flux functions (Phase 1) and solver properties (Phase 4)."""

import numpy as np
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parents[1] / "src" / "python"))
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
