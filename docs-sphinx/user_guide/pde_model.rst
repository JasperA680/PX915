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
