"""Python interface and flux helpers for the LWR PDE solver."""

import subprocess
import numpy as np
from pathlib import Path
from typing import Union
import netCDF4 as nc


# ---------------------------------------------------------------------------
# Greenshields flux — Python mirrors of the Fortran elemental functions.
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
        ic_type, flux_type, bc_type, n_lanes, lane_change_rate.
    output_path:
        NetCDF file to write.
    exe:
        Path to compiled pde_solver binary.
    """
    defaults = dict(
        M=200, n_steps=500, v_max=1.0, rho_max=1.0,
        rho_left_bc=0.1, rho_right_bc=0.9,
        ic_type="riemann", flux_type="lf", bc_type="open",
        n_lanes=1, lane_change_rate=0.0,
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
        str(p["n_lanes"]),
        str(p["lane_change_rate"]),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"pde_solver failed (exit {result.returncode}):\n{result.stderr}"
        )


def load_pde_netcdf(path: Union[str, Path]) -> dict:
    """Load a PDE simulation NetCDF file (single-lane or multi-lane).

    Returns a dict with keys:
        density      — (n_steps+1, M)          single-lane, or sum over lanes for multi
        density_per_lane — (n_steps+1, n_lanes, M)  always present
        flow         — (n_steps+1,)             total flow (sum over lanes)
        flow_per_lane — (n_steps+1, n_lanes)    per-lane right-boundary flow
        x            — (M,)
        time         — (n_steps+1,)
        n_lanes      — int
        attrs        — dict of global attributes

    NetCDF dimension layout written by the Fortran solver:
      density:  [lane, x, time]  (Fortran column-major)
      Python reads reversed:  (time, x, lane)
      This function transposes to (time, lane, x) for density_per_lane,
      then sums over lane axis to produce the backward-compat density field.

    Old single-lane files without a 'lane' dimension are also handled; they
    return n_lanes=1 and density_per_lane with shape (n_steps+1, 1, M).
    """
    with nc.Dataset(path, "r") as ds:
        x    = np.array(ds.variables["x"])
        time = np.array(ds.variables["time"])
        attrs = {k: ds.getncattr(k) for k in ds.ncattrs()}

        if "lane" in ds.dimensions:
            # New multi-lane file.
            # density written as Fortran [lane, x, time] -> Python reads (time, x, lane)
            raw = np.array(ds.variables["density"])          # (time, x, lane)
            density_per_lane = raw.transpose(0, 2, 1)        # (time, lane, x)

            # flow written as [lane, time] -> Python reads (time, lane)
            flow_per_lane = np.array(ds.variables["flow"])   # (time, lane)
            flow_total = np.array(ds.variables["flow_total"])  # (time,)

            n_lanes = density_per_lane.shape[1]
        else:
            # Legacy single-lane file (no 'lane' dim).
            raw = np.array(ds.variables["density"])          # (time, x)
            density_per_lane = raw[:, np.newaxis, :]         # (time, 1, x)
            flow_1d = np.array(ds.variables["flow"])         # (time,)
            flow_per_lane = flow_1d[:, np.newaxis]           # (time, 1)
            flow_total = flow_1d
            n_lanes = 1

    # Backward-compat density: squeeze for single lane, sum for multi-lane
    if n_lanes == 1:
        density = density_per_lane[:, 0, :]   # (time, x)
    else:
        density = density_per_lane.sum(axis=1)  # (time, x) total

    return dict(
        density=density,
        density_per_lane=density_per_lane,
        flow=flow_total,
        flow_per_lane=flow_per_lane,
        x=x,
        time=time,
        n_lanes=n_lanes,
        attrs=attrs,
    )
