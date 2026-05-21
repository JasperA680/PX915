Numerical methods
=================

This page describes the numerical methods used in the traffic-flow models. The
codebase contains two main modelling approaches:

* a cellular automaton model, based on the one-dimensional TASEP;
* a continuum PDE model, based on finite-volume methods for the LWR traffic-flow
  equation.

Cellular automaton update
-------------------------

The one-dimensional TASEP model represents the road as a lattice of ``L`` cells.
Each cell is either empty or occupied by a single vehicle. This exclusion rule
means that two vehicles cannot occupy the same site.

At each time step, the open-boundary TASEP update consists of three stages:

#. A vehicle at the right boundary exits with probability ``beta``.
#. Vehicles in the bulk move one site to the right if the next site is empty.
#. A vehicle enters at the left boundary with probability ``alpha`` if the first
   site is empty.

The density is computed as

.. math::

   \rho = \frac{N}{L},

where ``N`` is the number of occupied cells and ``L`` is the road length.

The current is measured from the number of vehicles leaving the right boundary.
Over a measurement window of ``T`` time steps,

.. math::

   J = \frac{1}{T}\sum_{t=1}^{T} n_{\mathrm{exit}}(t).

The steady-state measurement routine first discards a burn-in period and then
averages the density and current over a measurement window. The density is
measured in the central part of the lattice to reduce boundary-layer effects.

Nagel-Schreckenberg update
--------------------------

The code also contains a simple Nagel-Schreckenberg-style update. Unlike the
basic TASEP model, each vehicle has an integer velocity. A single update uses:

#. acceleration;
#. deceleration to avoid collisions;
#. random slowing;
#. vehicle movement.

The acceleration step is

.. math::

   v \leftarrow \min(v + 1, v_{\max}).

The deceleration step prevents vehicles from moving through the car in front:

.. math::

   v \leftarrow \min(v, g),

where ``g`` is the gap to the next vehicle. Random slowing then decreases the
velocity by one with probability ``P_SLOW`` if ``v > 0``.

Finite-volume PDE update
------------------------

The continuum PDE model is based on the scalar conservation law

.. math::

   \frac{\partial \rho}{\partial t}
   +
   \frac{\partial q(\rho)}{\partial x}
   =
   0,

where ``rho`` is the traffic density and ``q(rho)`` is the traffic flux.

The solver uses a finite-volume update. If ``rho_i^n`` is the cell-average
density in cell ``i`` at time step ``n``, then

.. math::

   \rho_i^{n+1}
   =
   \rho_i^n
   -
   \frac{\Delta t}{\Delta x}
   \left(
   F_{i+1/2} - F_{i-1/2}
   \right),

where ``F`` is a numerical flux at the cell interface.

Greenshields flux
-----------------

The default PDE closure is the Greenshields fundamental diagram. The velocity is

.. math::

   v(\rho)
   =
   v_{\max}
   \left(
   1 - \frac{\rho}{\rho_{\max}}
   \right),

and the corresponding flux is

.. math::

   q(\rho)
   =
   \rho v(\rho)
   =
   v_{\max}\rho
   \left(
   1 - \frac{\rho}{\rho_{\max}}
   \right).

This flux is zero at ``rho = 0`` and ``rho = rho_max``. It reaches its maximum at

.. math::

   \rho_c = \frac{\rho_{\max}}{2}.

The characteristic speed is

.. math::

   \frac{dq}{d\rho}
   =
   v_{\max}
   \left(
   1 - \frac{2\rho}{\rho_{\max}}
   \right).

This is positive in free-flow traffic, zero at the critical density, and negative
in congested traffic.

Newell-Daganzo triangular flux
------------------------------

The code also supports a Newell-Daganzo triangular fundamental diagram:

.. math::

   q(\rho)
   =
   \min
   \left(
   v_f \rho,
   w(\rho_{\max} - \rho)
   \right),

where ``v_f`` is the free-flow speed and ``w`` is the backward congestion wave
speed.

The critical density is

.. math::

   \rho_c
   =
   \frac{w\rho_{\max}}{v_f + w}.

This flux is useful because it gives a piecewise-linear approximation to the
fundamental diagram, with a free-flow branch and a congested branch.

Lax-Friedrichs flux
-------------------

The Lax-Friedrichs flux is a simple numerical flux used by the PDE solver:

.. math::

   F_{i+1/2}
   =
   \frac{1}{2}
   \left[
   q(\rho_L) + q(\rho_R)
   \right]
   -
   \frac{\Delta x}{2\Delta t}
   \left(
   \rho_R - \rho_L
   \right).

Here ``rho_L`` and ``rho_R`` are the left and right states at a cell interface.
This flux is robust and easy to implement, but it is relatively diffusive.

Godunov flux
------------

The Godunov flux solves the local Riemann problem at each cell interface. For a
concave flux such as the Greenshields flux, the implementation uses a closed-form
case distinction based on the critical density.

For a shock case, ``rho_L <= rho_R``, the flux is selected from the minimum flux
over the interval between the two states. For a rarefaction case,
``rho_L > rho_R``, the flux is selected from the maximum flux over the interval.
If the rarefaction crosses the critical density, the interface flux is the
maximum physical flux:

.. math::

   F_{i+1/2} = q(\rho_c).

The Godunov method is less diffusive than Lax-Friedrichs and is better suited for
capturing shocks and rarefactions in the traffic-density field.

CFL condition
-------------

The PDE solver uses a CFL condition to choose a stable time step. A typical
condition is

.. math::

   \Delta t
   \leq
   C_{\mathrm{CFL}}
   \frac{\Delta x}
   {\max |dq/d\rho|},

where ``C_CFL`` is the chosen CFL number. The adaptive time-step routine computes
the maximum characteristic speed over the current density field and updates the
time step accordingly.

Multilane source update
-----------------------

The multilane PDE model adds conservative source terms that transfer density
between neighbouring lanes. For lane ``l``, the update can be written as

.. math::

   \frac{\partial \rho_l}{\partial t}
   +
   \frac{\partial q_l(\rho_l)}{\partial x}
   =
   S_l(\rho_1, \rho_2, \ldots).

The source terms are constructed pairwise between adjacent lanes. For a pair of
lanes ``l`` and ``l+1``, traffic moves preferentially towards the faster lane,
provided that the receiving lane has available capacity.

The transfer rate from lane ``l`` to lane ``l+1`` is

.. math::

   S_{l \rightarrow l+1, i}
   =
   k \rho_{l,i}
   \max\left(
   0,
   1 - \frac{\rho_{l+1,i}}{\rho_{\max,l+1}}
   \right)
   \max\left(
   0,
   v_{l+1}(\rho_{l+1,i}) - v_l(\rho_{l,i})
   \right).

The reverse transfer is defined similarly. The update is conservative because
what one lane loses, the neighbouring lane gains. Therefore, for every spatial
cell,

.. math::

   \sum_l S_{l,i} = 0.

This means lane changing redistributes density between lanes but does not create
or destroy vehicles.

Multilane equilibrium validation
--------------------------------

For two lanes with equal jam density and periodic boundary conditions, the
lane-change source term should drive the system towards a state where the lane
velocities are equal:

.. math::

   v_1(\rho_{1,\mathrm{eq}})
   =
   v_2(\rho_{2,\mathrm{eq}}).

Using conservation of total density,

.. math::

   \rho_{1,\mathrm{eq}} + \rho_{2,\mathrm{eq}} = M_0,

and defining

.. math::

   r = \frac{v_{\max,2}}{v_{\max,1}},

the analytical equilibrium is

.. math::

   \rho_{2,\mathrm{eq}}
   =
   \frac{r - 1 + M_0}{r + 1},

and

.. math::

   \rho_{1,\mathrm{eq}}
   =
   M_0 - \rho_{2,\mathrm{eq}}.

This provides a useful validation test for the multilane PDE implementation:
the numerical steady state should approach the analytical prediction after a
sufficiently long simulation.