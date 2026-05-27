Cellular automaton traffic models
==================================

Cellular automaton (CA) models take a different approach to traffic flow from
the LWR PDE.  Instead of treating traffic as a continuous fluid, they represent
individual vehicles as discrete particles on a lattice, updated in parallel at
each integer time step.

Two CA models are implemented: the **TASEP** (Totally Asymmetric Simple
Exclusion Process) and the **Nagel–Schreckenberg (NS)** model.  Both are driven
by the same Fortran network simulator (``build/run_network``) and share the same
road-network infrastructure: multi-lane roads, junction routing, and the optional
lane-change sub-step.


The TASEP model
---------------

The TASEP is the minimal open-boundary lattice gas.  The road is divided into
:math:`L` sites, each either occupied (:math:`1`) or empty (:math:`0`).  At
every time step all sites are updated simultaneously in three stages:

1. **Exit:** if site :math:`L` is occupied, the vehicle leaves with probability
   :math:`\beta`.
2. **Hop:** each vehicle at sites :math:`1, \ldots, L-1` moves one step forward
   if the next site is empty.
3. **Entry:** if site :math:`1` is empty, a new vehicle enters with probability
   :math:`\alpha`.

These three rules produce three macroscopically distinct steady states depending
on the inflow rate :math:`\alpha` and outflow rate :math:`\beta`:

.. list-table::
   :header-rows: 1
   :widths: 20 20 20 40

   * - Phase
     - Condition
     - Bulk density
     - Character
   * - Low density (LD)
     - :math:`\alpha < \beta`,\ :math:`\alpha < \tfrac{1}{2}`
     - :math:`\rho \approx \alpha`
     - Sparse diagonal streaks; inflow is the bottleneck
   * - Maximum current (MC)
     - :math:`\alpha > \tfrac{1}{2}`,\ :math:`\beta > \tfrac{1}{2}`
     - :math:`\rho \approx 0.5`
     - Dense, fluctuating; domain walls diffuse across the lattice
   * - High density (HD)
     - :math:`\beta < \alpha`,\ :math:`\beta < \tfrac{1}{2}`
     - :math:`\rho \approx 1 - \beta`
     - Dense jams; outflow is the bottleneck

All three phases share the same fundamental diagram

.. math::

   J(\rho) = \min(\rho,\; 1 - \rho),

a triangle whose peak :math:`J_\text{max} = 0.5` is reached at
:math:`\rho = 0.5`.  The two straight-line branches arise from the
deterministic bulk parallel update — at half-filling the lattice settles into
an alternating "X O X O" pattern in which every car hops every step.


The Nagel–Schreckenberg model
------------------------------

The NS model extends TASEP by giving each vehicle an integer velocity
:math:`v \in \{0, 1, \ldots, v_\text{max}\}`.  Its four parallel update rules
model human driving behaviour more realistically:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Step
     - Rule
   * - **Accelerate**
     - :math:`v \leftarrow \min(v + 1,\ v_\text{max})`
   * - **Brake**
     - :math:`v \leftarrow \min(v,\ d - 1)`, where :math:`d` is the gap to the
       next vehicle
   * - **Randomise**
     - with probability :math:`p_\text{slow}`:
       :math:`v \leftarrow \max(v - 1,\ 0)`
   * - **Move**
     - vehicle advances :math:`v` sites

The randomisation step is the key difference from TASEP.  Even at low densities
(:math:`\rho \ll 1`), if :math:`p_\text{slow} > 0` a single random deceleration
propagates backward through following vehicles and amplifies into a self-sustaining
**phantom jam** — a density wave travelling upstream at roughly one site per step
with no external cause.

Phantom jams appear reliably above a critical density

.. math::

   \rho_c = \frac{1}{v_\text{max} + 1}.

For the defaults :math:`v_\text{max} = 5`,\ :math:`p_\text{slow} = 0.2`,
this gives :math:`\rho_c \approx 0.17`.

The NS fundamental diagram differs markedly from the TASEP triangle:

- **Free-flow branch** (:math:`\rho < \rho_c`): roughly linear rise — each
  additional vehicle contributes :math:`\approx v_\text{max}` moves per step.
- **Congested branch** (:math:`\rho > \rho_c`): steep decline with the
  characteristic backward-wave signature of phantom jams.
- **Peak flux** is lower than TASEP's :math:`J = 0.5` and the maximum shifts
  to :math:`\rho_c < 0.5`.

The NS model is run on a **periodic ring** (``is_periodic = True`` on the lane)
for the fundamental diagram sweep, so there is no boundary bias.  Use the
``single_lane_periodic`` preset to reproduce this geometry.


Lane changing
-------------

The lane-change module (``lane_change.f90``) applies a sub-step before each NS
update.  A vehicle considers switching to an adjacent lane only when all four
conditions hold:

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Condition
     - Description
   * - **T1** — blocked
     - Gap ahead in the current lane :math:`< v + 1` (cannot reach full speed)
   * - **T2** — target ahead clear
     - Gap ahead in the target lane :math:`> v + 1`
   * - **T3** — target behind clear
     - Gap behind in the target lane :math:`> v_\text{max}`
   * - **T4** — random acceptance
     - Uniform draw :math:`< p_\text{change}`

Two lane-change modes are available, selected via ``lc_model`` in
:class:`~python.road_network.SimParams`:

``lc_model = 0`` (symmetric)
    Both adjacent lanes are considered equally.  A vehicle changes to whichever
    valid target satisfies T1–T4.

``lc_model = 1`` (asymmetric)
    Vehicles prefer the rightmost lane and only move left if it satisfies T1–T4,
    returning to the right as soon as the right lane is clear.

``lc_model = -1`` (disabled)
    No lane changes; each lane evolves independently (the default).

With symmetric lane changing enabled, blocked vehicles escape jammed lanes by
switching, which reduces peak jam density in either lane and slightly increases
mean throughput.  The two lanes become correlated: a dense patch in one lane
produces a transient shadow in the other as vehicles simultaneously leave the
jammed lane and occupy the clear one.


Road-network simulation
-----------------------

The Fortran binary supports arbitrary road networks: multiple roads connected at
junctions, each with a routing probability table specifying where vehicles go on
exit, and right-of-way rules to resolve simultaneous crossing conflicts.

**Junction right-of-way** is evaluated geometrically: each vehicle crossing a
junction is assigned a perimeter position on the junction boundary.  Two movements
conflict if their chord segments on the perimeter intersect (a left turn across an
oncoming straight-through movement, for example).  Conflicting movements are
resolved stochastically — one is allowed to proceed and the other waits.

Seven network presets are provided in :mod:`python.road_network`:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Preset
     - Description
   * - ``single_lane``
     - One road, one open-boundary lane.  Baseline for TASEP / NS runs.
   * - ``single_lane_periodic``
     - One road, one closed-ring NS lane, seeded with a fixed vehicle count.
       Used for the fundamental-diagram NS sweep.
   * - ``two_lane``
     - One road, two same-direction lanes.  Lane-change baseline.
   * - ``t_junction``
     - Three roads meeting at a T: West arm, East arm, and South stem.
   * - ``crossroads``
     - Four-arm crossroads; richer conflict patterns than the T-junction.
   * - ``roundabout``
     - Four T-junctions arranged in a clockwise ring.
   * - ``town``
     - 2 × 2 grid of crossroads with eight external entry arms.

A **T-junction** example: with ``p_through = 0.6``, 60 % of westbound traffic
continues straight to the East arm and 40 % diverts down the South stem.
Southbound traffic splits evenly (50 % West, 50 % East).  Road 3 (stem) therefore
receives roughly 1.8× the load of a single arm and will show the highest mean
density in the density time-series panel.


Running the simulator
---------------------

The simulator is driven through the Python API.  Import the preset builders and
the :func:`~python.run_simulation.run_simulation` function:

.. code-block:: python

   from python.road_network import single_lane, SimParams
   from python.run_simulation import run_simulation
   from python.io import load_network_netcdf

   spec, layout = single_lane(L=100, alpha=0.3, beta=0.8)
   params = SimParams(n_steps=400, model='TASEP', v_max=1, p_slow=0.0, rng_seed=42)

   nc_path = run_simulation(spec, params, layout, output_dir='data/output/run1',
                            binary='build/run_network')
   result  = load_network_netcdf(nc_path)

For the NS phantom-jam demo on a periodic ring:

.. code-block:: python

   from python.road_network import single_lane_periodic, SimParams
   from python.run_simulation import run_simulation

   spec, layout = single_lane_periodic(L=200, n_vehicles=40)
   params = SimParams(n_steps=600, model='NS', v_max=5, p_slow=0.3, rng_seed=42)

   nc_path = run_simulation(spec, params, layout, output_dir='data/output/ring',
                            binary='build/run_network')

To run a fundamental-diagram sweep across the full density range:

.. code-block:: python

   from python.fd_runner import run_fd_sweep, load_fd_netcdf

   nc = run_fd_sweep(model='NS', L=100, n_points=25, n_steps=1500,
                     v_max=5, p_slow=0.2, seed=42,
                     output_path='data/output/fd_ns.nc')
   fd = load_fd_netcdf(nc)
   rho, J = fd['rho'], fd['J']

The TASEP sweep traces two deterministic-boundary cuts: :math:`\alpha` swept
with :math:`\beta = 1` (LD branch) then :math:`\beta` swept with :math:`\alpha
= 1` (HD branch).  The NS sweep drives a periodic ring at each target density
and measures the mean flow after a burn-in period.


Parameters
----------

:class:`~python.road_network.SimParams` controls the simulation:

.. list-table::
   :header-rows: 1
   :widths: 20 45 15 20

   * - Parameter
     - Description
     - Default
     - Typical range
   * - ``n_steps``
     - Number of time steps
     - 2000
     - 300 – 5000
   * - ``dt``
     - Nominal time-step size (informational only; CA time is discrete)
     - 1.0
     - —
   * - ``v_max``
     - Maximum integer velocity (NS only; set to 1 for TASEP-like behaviour)
     - 5
     - 1 – 5
   * - ``p_slow``
     - Randomisation probability (NS only; 0 disables phantom jams)
     - 0.2
     - 0.0 – 0.5
   * - ``rng_seed``
     - Random seed for reproducibility
     - 42
     - any integer
   * - ``model``
     - Update rule: ``'NS'`` or ``'TASEP'``
     - ``'NS'``
     - —
   * - ``lc_model``
     - Lane-change mode: ``-1`` disabled, ``0`` symmetric, ``1`` asymmetric
     - ``-1``
     - —
   * - ``lc_p_change``
     - Probability of accepting a valid lane-change opportunity
     - 1.0
     - 0.0 – 1.0


Visualising the output
----------------------

Load the NetCDF result with :func:`~python.io.load_network_netcdf`, then use the
visualisation helpers in :mod:`python.CA_visualisation`:

.. code-block:: python

   from python.io import load_network_netcdf
   from python.CA_visualisation import (
       plot_network_spacetime,
       plot_network_layout,
       plot_network_density,
       plot_fundamental_diagram,
   )

   result = load_network_netcdf('data/output/run1/result.nc')

   plot_network_spacetime(result, road_id=1)   # space-time diagram
   plot_network_layout(result)                 # road-network topology
   plot_network_density(result)                # mean density per road vs time

   # Fundamental diagram (from fd_runner output)
   plot_fundamental_diagram(rho, J, model='NS')

Space-time diagrams show site occupancy (black = occupied, white = empty) with
site on the :math:`y`-axis and time on the :math:`x`-axis.  Diagonal bands
slanting to the right are freely-moving vehicles; diagonal bands slanting to the
left are backward-propagating jam waves.


Limitations
-----------

**Discrete representation**

CA models resolve individual vehicles and produce stochastic, sample-path
outputs.  Ensemble averages over many independent seeds are needed to obtain
smooth density profiles or reliable fundamental-diagram estimates.  Single runs
at small :math:`L` exhibit strong finite-size fluctuations.

**TASEP does not support variable speeds**

The TASEP update rule is binary (move / do not move); ``v_max`` and ``p_slow``
are ignored when ``model = 'TASEP'``.  Lane-changing is also not available for
TASEP runs: the model requires ``model = 'NS'`` with at least two same-direction
lanes.

**Finite-size effects at the phase boundary**

Near the LD/MC and HD/MC phase boundaries the domain wall performs a slow random
walk and steady-state convergence is slow.  Runs at :math:`L = 100`,
:math:`N_\text{steps} = 400` may not have fully converged; increase both by a
factor of 5 near :math:`\alpha \approx 0.5` or :math:`\beta \approx 0.5`.

**Open boundary artefacts**

With open boundaries (``alpha``, ``beta`` < 1) the density profile is not
spatially uniform: boundary layers of width :math:`\sim\!\sqrt{L}` form near
each end.  The mean density reported by the simulator includes these layers;
for bulk properties use a periodic ring (``single_lane_periodic``) or discard
sites within :math:`\sim 2\sqrt{L}` of each boundary.

**Periodic ring requires fixed vehicle count**

The ``single_lane_periodic`` preset seeds vehicles evenly and runs a closed ring
at fixed :math:`N`; density is conserved exactly.  This is appropriate for the
fundamental-diagram sweep (which needs a fixed density) but does not model
inflow/outflow dynamics.

**No on-ramps in the 1D lane model**

The single-lane presets treat the road as one segment.  Merging or diverging
flows require a junction preset (``t_junction``, ``crossroads``, or similar)
with a routing probability table.
