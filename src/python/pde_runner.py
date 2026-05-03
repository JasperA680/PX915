"""Python interface and flux helpers for the LWR PDE solver."""

import subprocess
import numpy as np
from pathlib import Path
from typing import Union
import netCDF4 as nc


# ---------------------------------------------------------------------------
# Greenshields flux — Python mirrors of the Fortran elemental functions.
# These are used directly by tests and by the analysis layer.
# ---------------------------------------------------------------------------

def v_of_rho(rho, v_max, rho_max):
    """Greenshields velocity: v(ρ) = v_max · (1 − ρ/ρ_max)."""
    return v_max * (1.0 - rho / rho_max)


def q_of_rho(rho, v_max, rho_max):
    """Greenshields flow: q(ρ) = v_max · ρ · (1 − ρ/ρ_max)."""
    return v_max * rho * (1.0 - rho / rho_max)


def dq_drho(rho, v_max, rho_max):
    """Characteristic speed: dq/dρ = v_max · (1 − 2ρ/ρ_max)."""
    return v_max * (1.0 - 2.0 * rho / rho_max)


def rho_critical(rho_max):
    """Argmax of q for Greenshields: ρ_c = ρ_max / 2."""
    return rho_max / 2.0


# ---------------------------------------------------------------------------
# Fortran runner
# ---------------------------------------------------------------------------

_DEFAULT_EXE = Path(__file__).parents[2] / "build" / "pde_solver"
_DEFAULT_OUT = Path(__file__).parents[2] / "data" / "output" / "pde_simulation.nc"


def run_pde(params: dict, output_path: Union[str, Path] = _DEFAULT_OUT,
            exe: Union[str, Path] = _DEFAULT_EXE) -> None:
    """Run the Fortran LWR solver.

    Parameters
    ----------
    params:
        Dictionary with any subset of the solver keys:
        M, n_steps, v_max, rho_max, rho_left_bc, rho_right_bc,
        ic_type, flux_type.
    output_path:
        NetCDF file to write.
    exe:
        Path to compiled pde_solver binary.
    """
    defaults = dict(
        M=200, n_steps=500, v_max=1.0, rho_max=1.0,
        rho_left_bc=0.1, rho_right_bc=0.9,
        ic_type="riemann", flux_type="lf", bc_type="open",
    )
    p = {**defaults, **params}
    cmd = [
        str(exe),
        str(p["M"]),
        str(p["n_steps"]),
        str(p["v_max"]),
        str(p["rho_max"]),
        str(p["rho_left_bc"]),
        str(p["rho_right_bc"]),
        str(p["ic_type"]),
        str(p["flux_type"]),
        str(p["bc_type"]),
        str(output_path),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"pde_solver failed (exit {result.returncode}):\n{result.stderr}"
        )


def load_pde_netcdf(path: Union[str, Path]) -> dict:
    """Load a PDE simulation NetCDF file.

    Returns a dict with keys:
        density  — np.ndarray (n_steps+1, M)  [time × space]
        flow     — np.ndarray (n_steps+1,)
        x        — np.ndarray (M,)
        time     — np.ndarray (n_steps+1,)
        attrs    — dict of global attributes
    """
    with nc.Dataset(path, "r") as ds:
        # Fortran dims [x, time] → C/Python reads as (n_steps+1, M) already
        density = np.array(ds.variables["density"])
        flow    = np.array(ds.variables["flow"])
        x       = np.array(ds.variables["x"])
        time    = np.array(ds.variables["time"])
        attrs   = {k: ds.getncattr(k) for k in ds.ncattrs()}
    return dict(density=density, flow=flow, x=x, time=time, attrs=attrs)
