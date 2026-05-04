"""Exact LWR Riemann solution for the Greenshields fundamental diagram.

Consumed by the Phase 4 solver validation tests in test_pde.py.
"""

import numpy as np


def exact_riemann_lwr(x, t, rho_L, rho_R, v_max=1.0, rho_max=1.0, x0=0.5):
    """Exact solution to the LWR Riemann problem with Greenshields flux.

    Initial condition: rho_L for x < x0, rho_R for x >= x0.

    Parameters
    ----------
    x : array-like
        Spatial positions to evaluate (e.g. cell centres).
    t : float
        Evaluation time. t <= 0 returns the initial condition.
    rho_L, rho_R : float
        Left and right densities.
    v_max, rho_max : float
        Greenshields model parameters.
    x0 : float
        Initial discontinuity location (default 0.5 = domain midpoint).

    Returns
    -------
    np.ndarray, shape matching x
        Exact density at each position.

    Notes
    -----
    For a *concave* flux (Greenshields):

    rho_L < rho_R — shock at Rankine-Hugoniot speed
        s = (q(rho_R) - q(rho_L)) / (rho_R - rho_L)

    rho_L > rho_R — centred rarefaction fan; in the fan the self-similar
        variable xi = (x - x0) / t satisfies xi = dq/drho(rho), giving
        rho = (rho_max / 2) * (1 - xi / v_max).
    """
    x = np.asarray(x, dtype=float)

    if t <= 0:
        return np.where(x < x0, rho_L, rho_R).astype(float)

    xi = (x - x0) / t

    def q(r):
        return v_max * r * (1.0 - r / rho_max)

    def char_speed(r):
        return v_max * (1.0 - 2.0 * r / rho_max)

    if rho_L < rho_R:
        # Shock: Rankine-Hugoniot speed
        s = (q(rho_R) - q(rho_L)) / (rho_R - rho_L)
        return np.where(xi < s, rho_L, rho_R)

    if rho_L > rho_R:
        # Rarefaction fan bounded by the two characteristic speeds
        xi_L = char_speed(rho_L)   # left (slower) edge of fan
        xi_R = char_speed(rho_R)   # right (faster) edge of fan
        rho_fan = (rho_max / 2.0) * (1.0 - xi / v_max)
        return np.where(xi <= xi_L, rho_L,
               np.where(xi >= xi_R, rho_R, rho_fan))

    # rho_L == rho_R: trivial constant solution
    return np.full_like(x, rho_L)
