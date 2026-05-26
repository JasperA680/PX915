Fortran API
===========

This page contains automatically generated API documentation for the Fortran
modules in the traffic-flow codebase.

The Fortran code is split into two main modelling approaches:

* cellular automata models, including the one-dimensional TASEP and
  Nagel-Schreckenberg update routines;
* continuum PDE models, including the single-lane and multilane LWR solvers.


Cellular automata models
------------------------

1D TASEP model
~~~~~~~~~~~~~~

The ``tasep_model`` module contains the one-dimensional open-boundary TASEP
cellular automaton model. It defines the update rules for one step of the 
model, and the boundary entrance/exit rules.

.. f:autosrcfile:: tasep.f90


Nagel-Schreckenberg model
~~~~~~~~~~~~~~~~~~~~~~~~~

The ``ns_model`` module contains the update logic for the implementation of a
Nagel-Schreckenberg type cellular automaton model, introducing variable speeds
and more realistic driver behaviour. It defines the update logic and a universal
speed limit for all cars. Each lane is updated under either open or periodic
boundary conditions, selected per-lane via the ``is_periodic`` flag on
``lane_t``: open lanes are driven by the per-lane ``alpha`` / ``beta``
inflow and outflow rates (or by a junction at the holding cell), while
periodic lanes wrap site 1 onto site ``L`` and conserve the vehicle count.

.. f:autosrcfile:: ns_model.f90

Cellular automata lane-change module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``lane_change_mod`` module introduces the ability to include lane-changing 
behaviour to any of the CA models, allowing for simulation of more realistic driver 
behaviour.

.. f:autosrcfile:: lane_change.f90



Road-network models
-------------------

Junction module
~~~~~~~~~~~~~~~

The ``junction_mod`` module evaluates movements through road-network junctions,
including destination sampling, right-of-way rules, physical blocking and
stochastic deadlock resolution.

.. f:autosrcfile:: junction.f90

Road network module
~~~~~~~~~~~~~~~~~~~

The ``road_network_mod`` module defines the core data types for the
road-network model, including the cell, lane, road, junction and network
container types. It also provides utility routines for snapshotting lane state
and counting occupied cells across the network. A lane may be marked as
periodic via the ``is_periodic`` flag on ``lane_t``; this is used by the
fundamental-diagram sweep to drive an NS chain as a closed ring at fixed
vehicle count.

.. f:autosrcfile:: road_network.f90


Vehicle module
~~~~~~~~~~~~~~

The ``vehicle_mod`` module defines the occupancy constants ``V_EMPTY`` and
``V_OCCUPIED`` used by the junction holding-cell detection logic, and provides
a simple occupancy predicate.

.. f:autosrcfile:: vehicle.f90


Network initialisation module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``network_init_mod`` module provides hard-coded initialisers for canonical
road-network topologies, including a four-arm crossroad and a T-junction. It
also supplies a general lane allocator and a network deallocator.

.. f:autosrcfile:: network_init.f90


Network builder module
~~~~~~~~~~~~~~~~~~~~~~

The ``network_builder_mod`` module converts a flat network specification into a
fully allocated road-network object ready for simulation.

.. f:autosrcfile:: network_builder.f90


Network simulation module
~~~~~~~~~~~~~~~~~~~~~~~~~

The ``network_simulation_mod`` module provides the top-level driver for a
single network timestep. It orchestrates the full parallel update sequence:
network snapshot, optional lane-changing sub-step, junction evaluation, and
per-lane cellular automaton update using either the TASEP or
Nagel-Schreckenberg rule. The longitudinal updates themselves live in
``tasep_model`` and ``NS_model``; this module only sequences them with the
junction and lane-change sub-steps.

.. f:autosrcfile:: network_simulation.f90


Network I/O module
~~~~~~~~~~~~~~~~~~

The ``network_io_mod`` module writes the full network simulation output to a
NetCDF file. The output includes per-step occupancy and velocity histories,
per-lane and per-road metadata, flattened junction tables, and global
simulation attributes.

.. f:autosrcfile:: network_io.f90


Network driver program
~~~~~~~~~~~~~~~~~~~~~~

The ``run_network`` program is the executable entry point for the road-network
simulator. It reads a configuration file, builds the network, seeds the random
number generator deterministically, runs the simulation loop, and writes the
result to a NetCDF file.

.. f:autosrcfile:: run_network.f90


Fundamental-diagram sweep
-------------------------

Fundamental-diagram module
~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``fundamental_diagram_mod`` module measures the steady-state density and
current of a CA chain across a parameter sweep, and writes the resulting
``(rho, J)`` arrays to a NetCDF file. There is no separate per-chain update
kernel: each measurement run builds a minimal one-lane ``road_network_t`` and
drives it through the same ``tasep_lane_step`` (open boundary, alpha / beta
sweep) and ``NS_model_step`` (periodic ring, fixed vehicle count) routines
used by the full network simulator. The companion ``fd_sweep`` program in the
same source file is the executable entry point used by the Python launcher
(``python.fd_runner.run_fd_sweep``); it parses the CLI, runs the chosen
alpha / beta or density branch, and calls ``write_fd_netcdf`` to persist the
result.

For TASEP the sweep traverses the open-boundary phase diagram along two
deterministic-boundary cuts: alpha is swept with ``beta = 1`` (LD branch
closing off at the ``(0.5, 0.5)`` max-current point) and beta is swept with
``alpha = 1`` (HD branch back to the same point). Holding the fixed boundary
deterministic removes its stochastic contribution to the bulk and recovers
the textbook ``J(rho) = min(rho, 1 - rho)`` curve; choosing it on the
LD/HD-MC phase-boundary line at ``0.5`` instead would put the swept-parameter
> 0.5 half of each branch onto that boundary, and at large L the diagram
would develop a visible gap around ``rho = 0.5``.

.. f:autosrcfile:: fundamental_diagram.f90


Continuum PDE models
--------------------

PDE flux module
~~~~~~~~~~~~~~~

The ``pde_flux`` module contains velocity, flux, characteristic-speed and
numerical-flux functions for the PDE traffic-flow solver.

.. f:autosrcfile:: pde_flux.f90


PDE solver module
~~~~~~~~~~~~~~~~~

The ``pde_solver`` module contains the main finite-volume PDE solver, including
parameter setup, state initialisation, time stepping, adaptive time-step
selection, finalisation and NetCDF output.

.. f:autosrcfile:: pde_module.f90


PDE lane-change module
~~~~~~~~~~~~~~~~~~~~~~

The ``pde_lanechange`` module contains conservative source terms for density
exchange between adjacent lanes in the multilane PDE model.

.. f:autosrcfile:: pde_lanechange.f90


PDE driver program
~~~~~~~~~~~~~~~~~~

The ``pde_driver`` program is the executable entry point for the PDE solver.
It parses positional command-line arguments, constructs the parameter set,
runs the time loop, and writes the result to a NetCDF file.

.. f:autosrcfile:: pde_driver.f90