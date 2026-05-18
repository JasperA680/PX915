.. PX915 Traffic Flow Modelling documentation master file, created by
   sphinx-quickstart on Sun May 17 10:43:35 2026.

PX915 Traffic Flow Modelling documentation
==========================================

This documentation describes the traffic-flow simulation code developed for
the PX915 group project.

Fortran PDE flux module
-----------------------

The ``pde_flux`` module contains the velocity, flux, characteristic speed and
numerical flux functions used by the PDE traffic-flow solver.

.. f:autosrcfile:: pde_flux.f90

Fortran PDE solver module
-------------------------

The ``pde_solver`` module contains the main finite-volume PDE solver, including

parameter setup, state initialisation, time stepping, adaptive time-step

selection, finalisation and NetCDF output.

.. f:autosrcfile:: pde_module.f90


.. toctree::
   :maxdepth: 2
   :caption: Contents: