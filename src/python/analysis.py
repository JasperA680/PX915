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
