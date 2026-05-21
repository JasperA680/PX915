"""Python interface and flux helpers for the LWR PDE solver."""

import subprocess
import numpy as np
from pathlib import Path
from typing import Union
import netCDF4 as nc


# ---------------------------------------------------------------------------
# Greenshields flux with speed limit — Python mirrors of the Fortran functions.
# v_limit=None means no restriction (equivalent to v_limit=v_max).
# ---------------------------------------------------------------------------

def v_of_rho(rho, v_max, rho_max, v_limit=None):
    """Greenshields velocity capped at v_limit: v(ρ) = min(v_max·(1−ρ/ρ_max), v_limit)."""
    vl = v_max if v_limit is None else v_limit
    return np.minimum(v_max * (1.0 - rho / rho_max), vl)


def q_of_rho(rho, v_max, rho_max, v_limit=None):
    """Greenshields flow with speed limit: q(ρ) = v(ρ)·ρ."""
    return v_of_rho(rho, v_max, rho_max, v_limit) * rho


def dq_drho(rho, v_max, rho_max, v_limit=None):
    """Characteristic speed: v_limit in speed-limited zone, v_max·(1−2ρ/ρ_max) elsewhere."""
    vl = v_max if v_limit is None else v_limit
    rho_star = rho_max * (1.0 - vl / v_max)
    return np.where(rho < rho_star, vl, v_max * (1.0 - 2.0 * rho / rho_max))


def rho_critical(rho_max):
    """Argmax of q for Greenshields (unchanged by speed limit): ρ_c = ρ_max / 2."""
    return rho_max / 2.0


# ---------------------------------------------------------------------------
# Newell-Daganzo (triangular) fundamental diagram with speed limit
# ---------------------------------------------------------------------------
NEWELL_W = 0.5  # must match Fortran NEWELL_W parameter


def q_newell(rho, rho_max, v_limit):
    """Newell triangular flow with speed limit: q(ρ) = min(v_limit·ρ, w·(ρ_max−ρ))."""
    return np.minimum(v_limit * rho, NEWELL_W * (rho_max - rho))


def dq_drho_newell(rho, rho_max, v_limit):
    """Piecewise constant: +v_limit free-flow, -w congested."""
    rc = rho_critical_newell(rho_max, v_limit)
    return np.where(rho < rc, v_limit, np.where(rho > rc, -NEWELL_W, 0.0))


def rho_critical_newell(rho_max, v_limit):
    """Newell critical density with speed limit: ρ_c = w·ρ_max / (v_limit + w)."""
    return NEWELL_W * rho_max / (v_limit + NEWELL_W)


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
        ic_type, flux_type, bc_type, v_limit,
        n_lanes, lane_change_rate, v_max_lanes, rho_max_lanes.

        v_limit defaults to v_max (no speed restriction).
        v_max_lanes / rho_max_lanes accept either a scalar (broadcast
        to all lanes) or a list of length n_lanes, e.g. [1.0, 1.5].
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
        v_max_lanes=None, rho_max_lanes=None,
    )
    p = {**defaults, **params}
    # v_limit defaults to v_max if not provided
    if "v_limit" not in p:
        p["v_limit"] = p["v_max"]

    def _to_csv(val, n):
        """Return comma-separated string for a scalar or list parameter."""
        if val is None:
            return None
        if hasattr(val, '__len__'):
            return ",".join(str(v) for v in val)
        return ",".join(str(val) for _ in range(n))

    # Positional args mirror the Fortran driver:
    # 1-9: M n_steps v_max rho_max rho_left rho_right ic flux bc
    # 10: v_limit   11: output   12: n_lanes   13: lane_change_rate
    # 14: v_max_lanes (optional)   15: rho_max_lanes (optional)
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
        str(p["v_limit"]),
        str(output_path),
        str(p["n_lanes"]),
        str(p["lane_change_rate"]),
    ]

    # Only append per-lane args if they differ from the broadcast scalar
    v_csv   = _to_csv(p["v_max_lanes"],   p["n_lanes"])
    rho_csv = _to_csv(p["rho_max_lanes"], p["n_lanes"])
    if v_csv is not None:
        cmd.append(v_csv)
        if rho_csv is not None:
            cmd.append(rho_csv)

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"pde_solver failed (exit {result.returncode}):\n{result.stderr}"
        )


def load_pde_netcdf(path: Union[str, Path]) -> dict:
    """Load a PDE simulation NetCDF file (single-lane or multi-lane).

    Returns a dict with keys:
        density          — (n_steps+1, M)          single-lane, or sum over lanes
        density_per_lane — (n_steps+1, n_lanes, M)  always present
        flow             — (n_steps+1,)             total flow (sum over lanes)
        flow_per_lane    — (n_steps+1, n_lanes)     per-lane right-boundary flow
        x                — (M,)
        time             — (n_steps+1,)
        n_lanes          — int
        attrs            — dict of global attributes

    NetCDF layout (Fortran driver):
      density:  [lane, x, time]  -> Python reads (time, x, lane)
      Transposed here to (time, lane, x) for density_per_lane.
      Old single-lane files without a 'lane' dimension are also handled.
    """
    with nc.Dataset(path, "r") as ds:
        x    = np.array(ds.variables["x"])
        time = np.array(ds.variables["time"])
        attrs = {k: ds.getncattr(k) for k in ds.ncattrs()}

        if "lane" in ds.dimensions:
            # Multi-lane file: density written as [lane, x, time] -> (time, x, lane)
            raw = np.array(ds.variables["density"])          # (time, x, lane)
            density_per_lane = raw.transpose(0, 2, 1)        # (time, lane, x)
            flow_per_lane = np.array(ds.variables["flow"])   # (time, lane)
            flow_total = np.array(ds.variables["flow_total"])  # (time,)
            n_lanes = density_per_lane.shape[1]
        else:
            # Legacy single-lane file (no 'lane' dim)
            raw = np.array(ds.variables["density"])          # (time, x)
            density_per_lane = raw[:, np.newaxis, :]         # (time, 1, x)
            flow_1d = np.array(ds.variables["flow"])         # (time,)
            flow_per_lane = flow_1d[:, np.newaxis]           # (time, 1)
            flow_total = flow_1d
            n_lanes = 1

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
