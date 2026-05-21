import numpy as np


def tasep_step(state, alpha, beta, rng):
    """One parallel-update TASEP step, matching the Fortran implementation exactly.

    Order of operations (same as tasep.f90):
      1. Exit:  particle at site L exits with prob beta (based on old state).
      2. Bulk:  particles hop right if next site empty (based on old state).
      3. Entry: particle enters site 1 with prob alpha if site 1 is empty in new state.
    """
    old = state
    new = old.copy()
    exits = 0

    if old[-1] == 1 and rng.random() < beta:
        new[-1] = 0
        exits = 1

    can_hop = (old[:-1] == 1) & (old[1:] == 0)
    new[:-1] = np.where(can_hop, 0, new[:-1])
    new[1:]  = np.where(can_hop, 1, new[1:])

    if new[0] == 0 and rng.random() < alpha:
        new[0] = 1

    return new, exits


def run_tasep(L, n_steps, alpha, beta, burnin=0, seed=None, bulk_slice=None):
    """Run TASEP and return (density, current) arrays for the post-burnin period.

    bulk_slice: slice of sites used for density measurement.  Defaults to all
    sites.  Pass slice(L//4, 3*L//4) to measure only the bulk and avoid
    boundary-layer bias when comparing against J = rho(1-rho).
    """
    if bulk_slice is None:
        bulk_slice = slice(0, L)
    bulk_L = len(range(*bulk_slice.indices(L)))

    rng = np.random.default_rng(seed)
    state = np.zeros(L, dtype=np.int8)
    density = np.empty(n_steps)
    current = np.empty(n_steps, dtype=int)

    for t in range(burnin + n_steps):
        state, exits = tasep_step(state, alpha, beta, rng)
        if t >= burnin:
            i = t - burnin
            density[i] = state[bulk_slice].sum() / bulk_L
            current[i] = exits

    return density, current


def fundamental_diagram(L=50, n_steps=3000, burnin=None, n_points=30, seed=None):
    """Return (rho, J) arrays by sweeping alpha and beta across the phase diagram.

    Alpha sweep (beta=1): traces the low-density branch and maximum-current peak.
    Beta sweep  (alpha=1): traces the high-density branch and maximum-current peak.
    Together these sample the full J = rho(1-rho) parabola.

    burnin defaults to 2*L**2 so that domain-wall fluctuations near the
    max-current phase boundary have time to equilibrate before statistics are recorded.
    """
    if burnin is None:
        burnin = 2 * L * L

    param_vals = np.linspace(0.02, 0.98, n_points)
    bulk = slice(3 * L // 8, 5 * L // 8)
    rho_list, J_list = [], []

    for alpha in param_vals:
        density, current = run_tasep(L, n_steps, alpha, beta=0.5, burnin=burnin, seed=seed, bulk_slice=bulk)
        rho_list.append(density.mean())
        J_list.append(current.mean())

    for beta in param_vals:
        density, current = run_tasep(L, n_steps, alpha=0.5, beta=beta, burnin=burnin, seed=seed, bulk_slice=bulk)
        rho_list.append(density.mean())
        J_list.append(current.mean())

    return np.array(rho_list), np.array(J_list)


# ---------------------------------------------------------------------------
# Nagel–Schreckenberg fundamental diagram (pure-Python, periodic ring)
# ---------------------------------------------------------------------------

def ns_step_periodic(positions, velocities, L, v_max, p_slow, rng):
    """One parallel-update Nagel–Schreckenberg step on a periodic L-cell ring.

    Parameters
    ----------
    positions   : sorted 1D int array of occupied cell indices, 0..L-1
    velocities  : 1D int array of velocities, aligned with ``positions``
    L           : ring length (number of cells)
    v_max       : maximum vehicle speed (cells/step)
    p_slow      : random dawdle probability per occupied cell per step
    rng         : ``numpy.random.Generator``

    Returns
    -------
    (new_positions, new_velocities, moves_this_step)
        new_positions are sorted ascending; ``moves_this_step`` is the sum of
        new velocities (used directly as the flow contribution).
    """
    n = len(positions)
    if n == 0:
        return positions, velocities, 0

    # 1. Accelerate.
    new_v = np.minimum(velocities + 1, v_max)

    # 2. Brake to gap (distance to the next car, periodic).
    next_pos = np.roll(positions, -1)
    gaps = (next_pos - positions - 1) % L
    new_v = np.minimum(new_v, gaps)

    # 3. Random dawdling.
    randoms = rng.random(n)
    slow_mask = (new_v > 0) & (randoms < p_slow)
    new_v = np.where(slow_mask, new_v - 1, new_v)

    # 4. Move.
    new_pos = (positions + new_v) % L
    moves = int(new_v.sum())

    # Maintain ascending order even after wrap-around so the gap calculation
    # on the next step stays correct.
    order = np.argsort(new_pos, kind='stable')
    return new_pos[order], new_v[order], moves


def run_ns_density(L, n_steps, n_vehicles, v_max=5, p_slow=0.2,
                   burnin=None, seed=None):
    """Run NS on a periodic L-cell ring with ``n_vehicles`` cars.

    Density is constant (no inflow/outflow), so this returns the (ρ, J̄) point
    where J̄ is the time-averaged flow over the post-burn-in window.
    """
    if burnin is None:
        burnin = max(500, 5 * L)
    if n_vehicles <= 0:
        return 0.0, 0.0
    if n_vehicles >= L:
        return 1.0, 0.0

    rng = np.random.default_rng(seed)
    positions = np.sort(rng.choice(L, n_vehicles, replace=False))
    velocities = np.zeros(n_vehicles, dtype=int)

    for _ in range(burnin):
        positions, velocities, _ = ns_step_periodic(
            positions, velocities, L, v_max, p_slow, rng)

    total_moves = 0
    for _ in range(n_steps):
        positions, velocities, moves = ns_step_periodic(
            positions, velocities, L, v_max, p_slow, rng)
        total_moves += moves

    rho = n_vehicles / L
    J = total_moves / (L * n_steps)
    return rho, J


def fundamental_diagram_ns(L=100, n_steps=1500, n_points=30,
                           v_max=5, p_slow=0.2, seed=None):
    """Return (rho, J) sweep for NS on a periodic ring.

    Unlike the TASEP function, this sweeps the global density directly
    (periodic boundary, n_vehicles fixed per point) rather than (α, β).
    With p_slow > 0 the resulting curve is the Nagel–Schreckenberg
    inverse-λ shape: linear rise up to ρ_c ≈ 1/(v_max+1), then a near-linear
    fall toward zero at jam density.
    """
    rho_targets = np.linspace(0.02, 0.98, n_points)
    rho_list, J_list = [], []
    for rho in rho_targets:
        n = max(1, int(round(rho * L)))
        rho_meas, J = run_ns_density(
            L, n_steps, n, v_max=v_max, p_slow=p_slow, seed=seed)
        rho_list.append(rho_meas)
        J_list.append(J)
    return np.array(rho_list), np.array(J_list)


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
