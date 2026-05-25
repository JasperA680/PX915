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
speed limit for all cars.

.. f:autosrcfile:: ns_model.f90

Cellular automata lane-change module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``lane_change_mod`` module introduces the ability to include lane-changing 
behaviour to any of the CA models, allowing for simulation of more realistic driver 
behaviour.

.. f:autosrcfile:: lane_change.f90



Cellular automata simulation module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``simulation`` module contains higher-level routines for running
one-dimensional cellular automaton simulations. It provides simulation drivers,
time-history output, steady-state measurement routines, and an optional
Nagel-Schreckenberg update rule.

.. f:autosrcfile:: simulation.f90

Road-network models
-------------------

Junction module
~~~~~~~~~~~~~~~

The ``junction_mod`` module evaluates movements through road-network junctions,
including destination sampling, right-of-way rules, physical blocking and
stochastic deadlock resolution.

.. f:autosrcfile:: junction.f90

Network builder module
~~~~~~~~~~~~~~~~~~~~~~

The ``network_builder_mod`` module converts a flat network specification into a
fully allocated road-network object ready for simulation.

.. f:autosrcfile:: network_builder.f90

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