Parallelisation
===============

This page outlines how parallelisation could be added to the traffic-flow
simulation code as a future development. One subsystem — the
fundamental-diagram sweep — has already been parallelised with OpenMP
because it sits at the easy end of the spectrum (independent measurements
with no shared state); the rest of the code is still serial.

At present, the network simulator (``run_network``), the PDE solver, and
the cellular automaton update rules themselves are written as serial
implementations. This makes the update rules easier to validate and keeps
the model behaviour transparent. Parallelisation of those components would
be an additional performance feature rather than a core part of the
current implementation.

Motivation
----------

The simulation cost increases with the number of cells, lanes, roads, time steps
and parameter sweeps. Parallelisation could reduce runtime for larger studies,
especially when running many independent simulations for fundamental diagrams,
validation tests or uncertainty studies.

The main aim of any future parallel implementation should be to improve runtime
without changing the physical model or the numerical results.

Fundamental-diagram sweep (implemented)
---------------------------------------

The fundamental-diagram sweep driver (``build/fd_sweep``, source
``src/fortran/fundamental_diagram.f90``) is parallelised with OpenMP. The
sweep runs ``n_points`` (TASEP: ``2 * n_points``, alpha and beta branches)
independent steady-state measurements, each of which builds its own
single-lane ``road_network_t`` and calls ``measure_steady_state_tasep`` or
``measure_steady_state_ns``. The measurement subroutines touch no global
state — every allocation lives on the calling thread's local
``road_network_t`` — so the outer ``do i = 1, n_points`` loops parallelise
directly with ``!$omp parallel do``.

The only subtlety is the random number generator. The Fortran intrinsic
``random_number`` has per-thread state under gfortran's OpenMP runtime,
so the calls themselves are thread-safe, but the *output* of a sweep
depends on which RNG state each iteration sees. To make the result
bit-identical across thread counts, the driver re-seeds the calling
thread's RNG state at the top of every iteration with a deterministic
function of ``(seed, i, branch)`` (see
``fundamental_diagram_mod::seed_iter_rng``). Each iteration then runs
with its own private, reproducible random stream.

The build is one flag: ``-fopenmp`` on the ``fd_sweep`` recipe in the
top-level ``Makefile``. On macOS the system clang doesn't ship libomp,
but this project's gfortran (``hyperspy-bundle``, GCC 13.4) ships its own
OpenMP runtime, so no extra setup is required.

At the time of writing, on an 8-core M-series Mac, a single TASEP sweep
at ``L = 300, n_points = 20, n_steps = 3000`` runs in 1.79 s on one
thread and 0.31 s on eight (5.8x speedup; the high-end superlinearity
comes from better L1 / L2 fit per worker on small networks).

The scaling and bit-identity properties are checked by
``tests/test_fd_openmp_scaling.py``, which runs the sweep at
``OMP_NUM_THREADS = 1`` and ``OMP_NUM_THREADS = min(4, cpu_count)``,
asserts byte-equal NetCDF outputs for both TASEP and NS, and requires a
>= 1.5x wall-clock speedup on the parallel run.

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