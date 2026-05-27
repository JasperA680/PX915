Multilane PDE model
===================

The multilane extension runs one LWR solver per lane and couples them through
conservative lane-change source terms.  Vehicles move between lanes in
proportion to the velocity difference between adjacent lanes, so faster lanes
attract density from slower ones until the velocities equalise.

The density in lane *l* evolves as

.. math::

   \frac{\partial \rho_l}{\partial t}
   + \frac{\partial q_l}{\partial x}
   = S_l(\rho_1, \dots, \rho_L)

where the source term :math:`S_l` transfers mass from lane *l* to a faster
neighbour (or receives mass from a slower one).  The transfer rate between two
adjacent lanes is

.. math::

   S = k \, \rho_\text{source}
       \left(1 - \frac{\rho_\text{target}}{\rho_{\max}}\right)
       \max(0,\, \Delta v)

where *k* is the lane-change rate, the gap factor
:math:`(1 - \rho_\text{target}/\rho_\text{max})` prevents vehicles moving into
a full lane, and :math:`\Delta v` is the velocity advantage of the target lane.
The source terms conserve total density across all lanes at every time step.


Running the multilane solver
-----------------------------

Pass ``N_LANES`` as argument 12 and ``LANE_CHANGE_RATE`` as argument 13.
Arguments 1–11 are the same as the single-lane solver (see :doc:`pde_model`).

.. code-block:: bash

   build/pde_solver <M> <N_STEPS> <V_MAX> <RHO_MAX> <RHO_LEFT> <RHO_RIGHT> \
                    <IC_TYPE> <FLUX_TYPE> <BC_TYPE> <V_LIMIT> <OUTPUT>      \
                    <N_LANES> <LANE_CHANGE_RATE>

Two-lane example — slow and fast lane, moderate lane-change rate:

.. code-block:: bash

   build/pde_solver 200 500 1.0 1.0 0.5 0.5 constant lf open 1.0        \
       data/output/multilane.nc 2 0.5 "1.0,1.5" "1.0,1.0"

Here both lanes start at density 0.5 but lane 2 has a higher free-flow speed
(1.5 vs 1.0), so vehicles gradually migrate to lane 2 until velocities match.


Per-lane parameters
--------------------

Arguments 14 and 15 override the scalar ``V_MAX`` and ``RHO_MAX`` with
lane-specific values supplied as a comma-separated string.

.. code-block:: bash

   # Arg 14: per-lane v_max
   "1.0,1.5"          # lane 1 speed limit 1.0, lane 2 speed limit 1.5

   # Arg 15: per-lane rho_max
   "1.0,1.0"          # same jam density in both lanes

The number of values in each string must equal ``N_LANES`` exactly, or the
solver will exit with an error.


Parameters
----------

.. list-table::
   :header-rows: 1
   :widths: 25 45 15 15

   * - Parameter
     - Description
     - Default
     - Notes
   * - ``N_LANES``
     - Number of lanes
     - 1
     - Arg 12
   * - ``LANE_CHANGE_RATE`` (*k*)
     - Controls how quickly vehicles switch lanes
     - 0.0
     - Arg 13; k = 0 disables lane changing
   * - ``V_MAX_LANES``
     - Free-flow speed for each lane (comma-separated)
     - All set to ``V_MAX``
     - Arg 14; optional
   * - ``RHO_MAX_LANES``
     - Jam density for each lane (comma-separated)
     - All set to ``RHO_MAX``
     - Arg 15; optional


Effect of the lane-change rate *k*
------------------------------------

*k* = 0
    Lanes are completely independent.  An initial density difference between
    lanes persists for the entire simulation.

*k* small (≈ 0.1)
    Slow equilibration.  Useful for studying transient behaviour when lanes
    start far from equilibrium.

*k* moderate (≈ 0.5)
    Lanes reach a common velocity within a few hundred time steps.  This is the
    physically realistic regime for motorway traffic.

*k* large (> 2)
    Very rapid equilibration.  Can cause density overshoot and oscillation;
    reduce the time step or check the clip-count warning in the solver output.


Speed limit with multiple lanes
---------------------------------

``V_LIMIT`` applies a global speed cap across all lanes.  This is distinct from
setting different ``V_MAX_LANES`` values:

- **Different** ``V_MAX_LANES``: each lane has a different fundamental diagram.
  Lane changing still occurs because the velocity functions differ.
- **Global** ``V_LIMIT``: all lanes share the same cap.  At low densities
  (ρ < ρ*) every lane travels at ``V_LIMIT`` regardless of its ``V_MAX``, which
  suppresses the velocity difference driving lane changes.

Example with speed limit applied:

.. code-block:: bash

   build/pde_solver 200 800 1.0 1.0 0.2 0.2 constant lf open 0.6        \
       data/output/multilane_limit.nc 2 0.5 "1.0,1.5" "1.0,1.0"


Visualising multilane output
-----------------------------

.. code-block:: bash

   python scripts/run_multilane_pde.py --save

Omit ``--scenario`` to run all scenarios, or pass a letter prefix to run one.
The space–time panel shows density in each lane as a function of time at the
domain mid-point, making it straightforward to observe equilibration.
