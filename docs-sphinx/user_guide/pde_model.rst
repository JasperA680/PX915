PDE traffic-flow model
======================

The PDE solver implements the Lighthill–Whitham–Richards (LWR) model: a
macroscopic, continuum description of traffic flow on a single lane.  Vehicles
are not tracked individually; instead the solver evolves the local traffic
density ρ(x, t) forward in time using a conservation law

.. math::

   \frac{\partial \rho}{\partial t} + \frac{\partial q}{\partial x} = 0

where q = ρ v(ρ) is the traffic flux and v(ρ) is the Greenshields velocity
relation

.. math::

   v(\rho) = v_\text{max} \left(1 - \frac{\rho}{\rho_\text{max}}\right).

The solver uses a finite-volume scheme (Lax–Friedrichs or Godunov) on a
uniform spatial grid and advances time with an explicit Euler step subject to
a CFL stability condition.


Running the solver
------------------

All arguments are positional.  A minimal run with default settings:

.. code-block:: bash

   build/pde_solver

To customise the run, pass arguments in order:

.. code-block:: bash

   build/pde_solver  <M>  <N_STEPS>  <V_MAX>  <RHO_MAX>  <RHO_LEFT>  <RHO_RIGHT> \
                     <IC_TYPE>  <FLUX_TYPE>  <BC_TYPE>  <V_LIMIT>  <OUTPUT>

Example — 400-cell Riemann problem, Godunov flux, 1000 steps:

.. code-block:: bash

   build/pde_solver 400 1000 1.0 1.0 0.2 0.8 riemann godunov open 1.0 data/output/run.nc

Output is written to a NetCDF file (``data/output/pde_simulation.nc`` by
default).


Parameters
----------

.. list-table::
   :header-rows: 1
   :widths: 20 45 15 20

   * - Parameter
     - Description
     - Default
     - Typical range
   * - ``M``
     - Number of spatial cells
     - 200
     - 100 – 1000
   * - ``N_STEPS``
     - Number of time steps
     - 500
     - 200 – 5000
   * - ``V_MAX``
     - Free-flow speed (dimensionless)
     - 1.0
     - 0.5 – 2.0
   * - ``RHO_MAX``
     - Jam density (dimensionless)
     - 1.0
     - 1.0 (fixed by convention)
   * - ``RHO_LEFT``
     - Left boundary density for Riemann IC
     - 0.1
     - 0.0 – 1.0
   * - ``RHO_RIGHT``
     - Right boundary density for Riemann IC
     - 0.9
     - 0.0 – 1.0
   * - ``V_LIMIT``
     - Speed cap applied to all vehicles (≤ V_MAX)
     - V_MAX (no cap)
     - 0.0 – V_MAX


Initial conditions
------------------

Set with ``IC_TYPE``:

``riemann``
    Step discontinuity: left half at ``RHO_LEFT``, right half at ``RHO_RIGHT``.
    Useful for studying shock waves and rarefaction fans.

``constant``
    Uniform density across the domain.  ``RHO_LEFT`` sets the value.

``sine``
    Sinusoidal perturbation around a mean density.  Good for testing numerical
    diffusion.

``periodic``
    Smooth periodic profile.  Use with ``BC_TYPE = periodic``.


Boundary conditions
-------------------

Set with ``BC_TYPE``:

``open``
    Inflow at the left boundary at density ``RHO_LEFT``; outflow at the right.
    Mimics a road segment embedded in a larger network.

``periodic``
    The domain wraps around: vehicles leaving the right re-enter at the left.
    Useful for studying bulk behaviour without boundary effects.

``sponge``
    A damping layer at each end gradually relaxes the density to the boundary
    values, reducing numerical reflections from strong shocks.


Flux schemes
------------

Set with ``FLUX_TYPE``:

``lf``
    Lax–Friedrichs: adds artificial diffusion proportional to the wave speed.
    Robust but slightly smears sharp fronts.

``godunov``
    Godunov (exact Riemann solver for the LWR model): sharper shocks, no extra
    diffusion.  Preferred for production runs.


Speed limit
-----------

Setting ``V_LIMIT`` below ``V_MAX`` caps the effective vehicle speed:

.. math::

   v_\text{eff}(\rho) = \min\!\left(v_\text{max}\!\left(1 -
   \frac{\rho}{\rho_\text{max}}\right),\; v_\text{limit}\right).

This creates a critical density

.. math::

   \rho^* = \rho_\text{max}\!\left(1 - \frac{v_\text{limit}}{v_\text{max}}\right)

below which all vehicles travel at exactly ``V_LIMIT``, and the flux–density
curve becomes linear for ρ < ρ*.

Example — speed limit at 60 % of free-flow speed:

.. code-block:: bash

   build/pde_solver 200 500 1.0 1.0 0.2 0.8 riemann lf open 0.6 data/output/speed_limit.nc


Visualising the output
----------------------

Use the provided Python script to plot a space–time diagram and fundamental
diagram from the NetCDF output:

.. code-block:: bash

   python scripts/run_pde_model.py --scenario riemann_shock --save

Plots are saved to ``plots/``.


Performance tuning
------------------

**Choosing M and N_STEPS**

Higher ``M`` gives sharper shock fronts but increases runtime linearly.  For
quick exploratory runs ``M = 100``, ``N_STEPS = 300`` is sufficient.  For
production quality results — particularly fundamental diagram sweeps or
speed-limit comparisons — use ``M = 200`` and ``N_STEPS ≥ 500``.

The CFL condition links spatial and temporal resolution: doubling ``M`` roughly
doubles the number of time steps needed to reach the same physical end time,
so total cost scales as M².  Keep this in mind when increasing resolution.

**Choosing the flux scheme**

``godunov`` is the recommended default.  It resolves shocks without numerical
smearing and costs roughly 1.5× more per step than ``lf``.

``lf`` (Lax–Friedrichs) adds artificial diffusion proportional to the wave
speed.  It is useful when you want a fast qualitative sweep or when the initial
condition is very smooth, but it noticeably smears sharp shock fronts at typical
resolutions.

**Fundamental diagram sweeps**

The Python sweep in ``analysis.pde_fundamental_diagram`` runs one solver
invocation per density point.  For a sweep with ``n_points = 20``, expect
roughly 40 solver calls.  Use ``n_steps = 2000`` with ``burnin_frac = 0.5``
to ensure the solver reaches steady state before measurement.  Reducing
``n_points`` to 10 halves the sweep time with only a small loss of resolution
in the diagram.

**Multilane runs**

Each additional lane scales the Fortran solver cost approximately linearly.
For two-lane runs the overhead is modest; beyond four lanes, consider reducing
``M`` to keep interactive response times acceptable.


Limitations
-----------

**Dimensionless units**

All quantities (density, velocity, flux) are dimensionless and normalised so
that ``RHO_MAX = 1`` and ``V_MAX = 1`` by default.  Mapping to physical units
requires independent calibration of free-flow speed and jam density for the
road of interest.

**Homogeneous traffic only**

The LWR model tracks a single vehicle class.  Heterogeneous traffic (e.g. a
mix of cars and trucks with different ``v_max`` values) is not supported; all
vehicles share the same Greenshields or Newell closure.

**No on-ramps or off-ramps in the 1D solver**

The single-lane solver treats the road as a closed segment with open or
periodic boundaries.  Merging or diverging flows require the network model
(``build/run_network``), which handles junctions through probabilistic routing
matrices.

**Open boundary artefacts**

With ``BC_TYPE = open``, the prescribed boundary densities are held fixed for
all time.  If the initial interior state is far from the boundary values, a
transient boundary layer forms at startup.  Discard an appropriate burn-in
period before measuring steady-state quantities.

**Speed limit as a hard velocity cap**

The ``V_LIMIT`` parameter applies a hard pointwise cap on the Greenshields
velocity.  This is a macroscopic approximation; it does not model variable
speed limits, enforcement zones, or reaction-time effects.

**Numerical diffusion at coarse resolution**

Even with the Godunov scheme, sharp shocks are spread over 2–3 cells.
Increasing ``M`` reduces this smearing but does not eliminate it entirely — the
scheme is first-order accurate in space.
