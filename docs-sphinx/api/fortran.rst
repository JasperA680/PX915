Fortran API
===========

PDE flux module
---------------

The ``pde_flux`` module contains velocity, flux, characteristic-speed and
numerical-flux functions for the PDE traffic-flow solver.

.. f:autosrcfile:: pde_flux.f90

PDE solver module
-----------------

The ``pde_solver`` module contains the main finite-volume PDE solver, including
parameter setup, state initialisation, time stepping, adaptive time-step
selection, finalisation and NetCDF output.

.. f:autosrcfile:: pde_module.f90

PDE lane-change module
----------------------

The ``pde_lanechange`` module contains conservative source terms for density
exchange between adjacent lanes in the multilane PDE model.

.. f:autosrcfile:: pde_lanechange.f90
