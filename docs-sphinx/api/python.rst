Python API
==========

Road network module
-------------------

The ``road_network`` module defines the data structures and preset builders used
to describe a road network before a simulation run.  It provides dataclasses for
lanes, roads, junctions, and layouts (``LaneSpec``, ``RoadSpec``,
``JunctionSpec``, ``LayoutSpec``, ``NetworkSpec``, ``SimParams``), a
``validate`` function that checks the spec for consistency, and seven preset
factory functions — ``single_lane``, ``single_lane_periodic`` (closed NS
ring), ``two_lane``, ``t_junction``, ``crossroads``, ``roundabout``, and
``town`` — that return ready-to-use ``(NetworkSpec, LayoutSpec)`` pairs.
The ``LaneSpec`` dataclass also carries the optional ``is_periodic`` flag
and ``n_vehicles`` initial-count field consumed by the Fortran builder.

.. automodule:: road_network
   :members:
   :undoc-members:
   :show-inheritance:

I/O module
----------

The ``io`` module handles reading and writing simulation data.  It can load
single-lane CA results from NetCDF (``load_netcdf``), write a full network
configuration to NetCDF (``write_config_netcdf``), and reload a completed
network run — including occupancy history and junction routing — as a
``NetworkResult`` dataclass (``load_network_netcdf``, ``load_network_for_restart``).

.. automodule:: python.io
   :members:
   :undoc-members:
   :show-inheritance:

Simulation runner module
------------------------

The ``run_simulation`` module is the Python-side launcher for the Fortran network
driver.  ``run_simulation`` serialises a ``NetworkSpec`` + ``SimParams`` pair to
a NetCDF config file, invokes the compiled ``run_network`` binary as a
subprocess, and returns a ``NetworkResult``.  ``write_config`` exposes the
config-writing step independently for inspection or pre-staging.

.. automodule:: python.run_simulation
   :members:
   :undoc-members:
   :show-inheritance:

Analysis module
---------------

The ``analysis`` module provides post-processing for PDE outputs:
``compute_total_density``, ``compute_total_flow``, ``compute_total_mass``,
and ``compute_multilane_fundamental_diagram`` /
``pde_multilane_fundamental_diagram``.  The CA / NS fundamental-diagram
sweeps that used to live here are now run by the Fortran driver
``build/fd_sweep`` and invoked from Python through ``fd_runner``.

.. automodule:: analysis
   :members:
   :undoc-members:
   :show-inheritance:

PDE runner module
-----------------

The ``pde_runner`` module is the Python interface to the LWR PDE solver.
``run_pde`` invokes the compiled PDE binary and writes output to a NetCDF file;
``load_pde_netcdf`` reads it back as a plain dictionary.  The module also
exposes the Greenshields and Newell–Daganzo flux helpers used for fundamental
diagram plots: ``q_of_rho``, ``dq_drho``, ``rho_critical``, ``q_newell``,
``dq_drho_newell``, and ``rho_critical_newell``.

.. automodule:: pde_runner
   :members:
   :undoc-members:
   :show-inheritance:

Fundamental-diagram runner module
---------------------------------

The ``fd_runner`` module is the Python interface to the Fortran
fundamental-diagram sweep driver (``build/fd_sweep``).  ``run_fd_sweep``
invokes the binary for the chosen model (open-boundary 1D TASEP with an
alpha / beta sweep, or periodic-ring Nagel–Schreckenberg with a density
sweep); ``load_fd_netcdf`` reads the resulting ``rho`` and ``J`` arrays
back as a plain dictionary alongside the sweep metadata.

.. automodule:: fd_runner
   :members:
   :undoc-members:
   :show-inheritance:

CA visualisation module
-----------------------

The ``CA_visualisation`` module contains plotting functions for the cellular
automaton (TASEP / Nagel–Schreckenberg) simulations.  It covers single-lane
diagnostics (space-time diagrams, density and current traces, fundamental
diagrams) as well as road-network visualisations (layout previews, per-road
density and flow traces, space-time diagrams per lane, and multi-panel summary
figures).

.. automodule:: CA_visualisation
   :members:
   :undoc-members:
   :show-inheritance:

PDE visualisation module
------------------------

The ``pde_visualisation`` module contains plotting functions for the single-lane
and multilane PDE solvers, including space-time density plots, density
snapshots, flow plots, fundamental diagrams, mass-conservation diagnostics, and
speed-limit comparison figures.

.. automodule:: pde_visualisation
   :members:
   :undoc-members:
   :show-inheritance:
