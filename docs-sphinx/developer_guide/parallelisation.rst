Parallelisation
===============

This page outlines how parallelisation could be added to the traffic-flow
simulation code as a future development.

At present, the code is primarily written as a serial implementation. This makes
the update rules easier to validate and keeps the model behaviour transparent.
Parallelisation would be an additional performance feature rather than a core
part of the current implementation.

Motivation
----------

The simulation cost increases with the number of cells, lanes, roads, time steps
and parameter sweeps. Parallelisation could reduce runtime for larger studies,
especially when running many independent simulations for fundamental diagrams,
validation tests or uncertainty studies.

The main aim of any future parallel implementation should be to improve runtime
without changing the physical model or the numerical results.

Cellular automaton models
-------------------------

The cellular automaton models update vehicles on a discrete lattice. The main
challenge is that neighbouring cells are coupled: whether a vehicle can move
depends on whether the next cell is empty.

A future parallel implementation would need to preserve the existing
parallel-update semantics. This means that all movement decisions should be made
from an old-state snapshot, while accepted moves are written to a new state.

For the one-dimensional TASEP model, possible parallelisation strategies include:

* splitting the lattice into chunks and handling chunk boundaries carefully;
* parallelising over independent simulation runs in a parameter sweep;
* parallelising over repeated stochastic realisations used for averaging.

Parallelising over independent runs is likely to be the safest first step,
because each simulation can be run without communication with the others.

PDE models
----------

The PDE solver uses a finite-volume update. Each cell update depends mainly on
local neighbouring fluxes, giving a stencil-like structure:

.. math::

   \rho_i^{n+1}
   =
   \rho_i^n
   -
   \frac{\Delta t}{\Delta x}
   \left(
   F_{i+1/2} - F_{i-1/2}
   \right).

This locality makes the PDE solver a natural candidate for shared-memory
parallelisation, for example using OpenMP. A future implementation could compute
fluxes or cell updates in parallel, provided that all reads are taken from a
consistent old state.

Multilane PDE models
--------------------

The multilane PDE model adds conservative lane-changing source terms. These
source terms couple neighbouring lanes as well as neighbouring spatial cells.

A future parallel implementation should probably split the update into separate
stages:

#. compute flux updates from the old density field;
#. compute lane-change source terms from the old density field;
#. combine the flux and source updates into the new density field.

This structure would help avoid race conditions and make mass-conservation checks
easier.

Performance measurements
------------------------

If parallelisation is added, performance should be measured using speedup:

.. math::

   S_p = \frac{T_1}{T_p},

where ``T_1`` is the serial runtime and ``T_p`` is the runtime using ``p``
parallel workers.

An ideal implementation would give:

.. math::

   S_p \approx p,

but real performance will be limited by serial sections, memory bandwidth,
communication overhead and load imbalance.

Validation requirements
-----------------------

Any future parallel version should be validated against the current serial code.
Useful checks include:

* matching density and current statistics for the 1D TASEP model;
* reproducing the same fundamental diagrams;
* checking mass conservation in the PDE and multilane PDE solvers;
* verifying that multilane equilibrium tests still agree with the analytical
  prediction;
* confirming that differences in stochastic simulations are consistent with
  random variation rather than implementation error.

Summary
-------

Parallelisation is a possible future extension of the codebase. The current
serial implementation provides a correctness reference, and any future
parallel implementation should be judged by both runtime improvement and
agreement with the validated serial results.