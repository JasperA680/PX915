import numpy as np


# The 1D TASEP and Nagel-Schreckenberg fundamental-diagram sweeps live in
# the Fortran driver ``build/fd_sweep`` and are invoked from Python via
# ``python.fd_runner``.  Earlier pure-Python implementations of
# ``tasep_step``, ``run_tasep``, ``fundamental_diagram``,
# ``ns_step_periodic``, ``run_ns_density``, and ``fundamental_diagram_ns``
# duplicated physics already in the Fortran codebase and have been removed.


def pde_fundamental_diagram(
    n_points=20, n_steps=2000, M=200,
    v_max=1.0, rho_max=1.0,
    burnin_frac=0.5,
):
    """Fundamental diagram for the LWR PDE by sweeping inflow/outflow density.

    Mirrors the TASEP ``fundamental_diagram()`` function.  Two branches are swept:

    Free-flow branch (rho < rho_c):
        rho_left_bc varies from 0 to rho_c; rho_right_bc = 0 (free outflow).
        Left-travelling characteristics dominate; interior settles to rho_left_bc.

    Congested branch (rho > rho_c):
        rho_right_bc varies from rho_c to rho_max; rho_left_bc = rho_max.
        Right-travelling characteristics dominate; interior settles to rho_right_bc.

    Each run uses IC="riemann" so a genuine transient forms before measurement.
    The post-burnin spatial-mean bulk density and right-boundary flow are returned.

    Parameters
    ----------
    n_points : int
        Sweep points per branch (total = 2 * n_points).
    n_steps : int
        Time steps per solver run.
    M : int
        Spatial cells.
    v_max, rho_max : float
        Greenshields parameters.
    burnin_frac : float
        Fraction of steps discarded as burn-in before averaging.

    Returns
    -------
    rho : np.ndarray, shape (2 * n_points,)
    q   : np.ndarray, shape (2 * n_points,)
    """
    import sys
    import tempfile
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent))
    from pde_runner import run_pde, load_pde_netcdf

    rho_c = rho_max / 2.0
    eps = rho_max * 0.02          # keep slightly away from 0, rho_c, rho_max
    burnin_step = int(n_steps * burnin_frac)
    bulk = slice(M // 4, 3 * M // 4)   # middle half, avoids boundary layers

    free_flow = np.linspace(eps, rho_c - eps, n_points)
    congested = np.linspace(rho_c + eps, rho_max - eps, n_points)

    # Build run list: (rho_left_bc, rho_right_bc, filename)
    runs = (
        [(rho, 0.0,     f"ff_{i}.nc") for i, rho in enumerate(free_flow)]
        + [(rho_max, rho, f"cg_{i}.nc") for i, rho in enumerate(congested)]
    )
    total = len(runs)
    rho_list, q_list = [], []

    with tempfile.TemporaryDirectory() as tmpdir:
        for k, (rho_l, rho_r, fname) in enumerate(runs, 1):
            print(f'\r  [{k}/{total}] rho_left={rho_l:.3f}  rho_right={rho_r:.3f}',
                  end='', flush=True)
            out = Path(tmpdir) / fname
            run_pde(
                dict(M=M, n_steps=n_steps, v_max=v_max, rho_max=rho_max,
                     rho_left_bc=rho_l, rho_right_bc=rho_r,
                     ic_type="riemann", flux_type="godunov", bc_type="open"),
                output_path=out,
            )
            data = load_pde_netcdf(out)
            rho_post = data["density"][burnin_step:, bulk]
            rho_list.append(float(rho_post.mean()))
            q_list.append(float(data["flow"][burnin_step:].mean()))
    print()

    return np.array(rho_list), np.array(q_list)


# ---------------------------------------------------------------------------
# Multi-lane PDE diagnostics
# ---------------------------------------------------------------------------

def compute_total_density(data):
    """Total density ρ_tot(x, t) = Σ_lane ρ_lane(x, t).

    Parameters
    ----------
    data : dict returned by load_pde_netcdf

    Returns
    -------
    np.ndarray, shape (n_steps+1, M)
    """
    return data["density_per_lane"].sum(axis=1)


def compute_total_flow(data):
    """Total right-boundary flow q_tot(t) = Σ_lane q_lane(t).

    Returns
    -------
    np.ndarray, shape (n_steps+1,)
    """
    return data["flow_per_lane"].sum(axis=1)


def compute_total_mass(data):
    """Total vehicle mass Σ_{lane,i} ρ_{lane,i} · Δx vs time.

    For periodic BCs this should be constant; use as a conservation diagnostic.

    Returns
    -------
    np.ndarray, shape (n_steps+1,)
    """
    dx = float(data["attrs"].get("dx", data["x"][1] - data["x"][0]))
    return data["density_per_lane"].sum(axis=(1, 2)) * dx


def compute_multilane_fundamental_diagram(data, burnin_frac=0.5):
    """Scatter of (ρ_tot_avg, q_tot_avg) over the post-burnin period.

    Suitable for overlay against single-lane and TASEP diagrams.

    Parameters
    ----------
    data : dict returned by load_pde_netcdf
    burnin_frac : float
        Fraction of time steps treated as burn-in.

    Returns
    -------
    rho_mean : float  — time-averaged total spatial-mean density
    q_mean   : float  — time-averaged total right-boundary flow
    """
    n_times = len(data["time"])
    burnin = int(n_times * burnin_frac)

    rho_tot = compute_total_density(data)     # (time, x)
    q_tot   = compute_total_flow(data)        # (time,)

    rho_mean = float(rho_tot[burnin:].mean())
    q_mean   = float(q_tot[burnin:].mean())
    return rho_mean, q_mean


def pde_multilane_fundamental_diagram(
    n_points=20, n_steps=2000, M=200,
    v_max=1.0, rho_max=1.0,
    n_lanes=2, lane_change_rate=0.0,
    burnin_frac=0.5,
):
    """Fundamental diagram for the multi-lane LWR PDE by sweeping inflow density.

    Mirrors ``pde_fundamental_diagram()`` but runs with n_lanes and returns
    total (summed across lanes) density and flow.

    Returns
    -------
    rho : np.ndarray, shape (2 * n_points,)
    q   : np.ndarray, shape (2 * n_points,)
    """
    import sys
    import tempfile
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).parent))
    from pde_runner import run_pde, load_pde_netcdf

    rho_c = rho_max / 2.0
    eps = rho_max * 0.02
    bulk = slice(M // 4, 3 * M // 4)

    free_flow = np.linspace(eps, rho_c - eps, n_points)
    congested = np.linspace(rho_c + eps, rho_max - eps, n_points)
    runs = (
        [(rho, 0.0,     f"ff_{i}.nc") for i, rho in enumerate(free_flow)]
        + [(rho_max, rho, f"cg_{i}.nc") for i, rho in enumerate(congested)]
    )
    total = len(runs)
    rho_list, q_list = [], []

    with tempfile.TemporaryDirectory() as tmpdir:
        for k, (rho_l, rho_r, fname) in enumerate(runs, 1):
            print(f'\r  [{k}/{total}] rho_left={rho_l:.3f}  rho_right={rho_r:.3f}',
                  end='', flush=True)
            out = Path(tmpdir) / fname
            run_pde(
                dict(M=M, n_steps=n_steps, v_max=v_max, rho_max=rho_max,
                     rho_left_bc=rho_l, rho_right_bc=rho_r,
                     ic_type="riemann", flux_type="godunov", bc_type="open",
                     n_lanes=n_lanes, lane_change_rate=lane_change_rate),
                output_path=out,
            )
            data = load_pde_netcdf(out)
            burnin_step = int(n_steps * burnin_frac)
            rho_tot = compute_total_density(data)
            rho_list.append(float(rho_tot[burnin_step:, bulk].mean()))
            q_list.append(float(compute_total_flow(data)[burnin_step:].mean()))
    print()

    return np.array(rho_list), np.array(q_list)
