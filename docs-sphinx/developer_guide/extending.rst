Extending the code
==================

This page describes how to add new functionality to the three main extension
points in the codebase: PDE flux schemes / closures, PDE initial conditions,
and network presets.  Each section lists the files to edit, the convention to
follow, and a minimal worked example.


Adding a new PDE flux scheme or closure
---------------------------------------

The PDE solver dispatches on the string ``flux_type`` at two levels:

1. **Fortran** — the closure functions and Godunov solver in
   ``src/fortran/pde_flux.f90``.
2. **Python runner** — the argument is passed through unchanged by
   ``src/python/pde_runner.py``; no Python changes are needed for the solver
   itself.

**Step 1 — implement the flux function in** ``src/fortran/pde_flux.f90``

Add a pure function that computes ``q(rho)`` for the new closure, following
the pattern of the existing ``q_greenshields`` and ``q_newell`` functions:

.. code-block:: fortran

   pure function q_myflux(rho, v_max, rho_max) result(q)
     real(dp), intent(in) :: rho, v_max, rho_max
     real(dp) :: q
     ! ... your closure here ...
   end function q_myflux

Add an equivalent ``dq_drho_myflux`` function for the characteristic speed
(needed by the CFL condition).

**Step 2 — register the new name in the dispatch functions**

In ``pde_flux.f90``, the functions ``q_dispatch``, ``dq_drho_dispatch``, and
``godunov_dispatch`` each contain a string-comparison block (lines ~343–400).
Add a new branch:

.. code-block:: fortran

   ! inside q_dispatch
   else if (index(trim(flux_type), 'myflux') > 0) then
     q = q_myflux(rho, v_max, rho_max)

Repeat the pattern for ``dq_drho_dispatch`` and, if you want Godunov support,
for ``godunov_dispatch``.

**Step 3 — expose it in the Python runner**

``src/python/pde_runner.py`` contains a helper ``q_of_rho(rho, params)`` used
by the visualisation and analysis layers.  Add a branch to its dispatch block
so that ``plot_pde_fundamental_diagram`` and the speed-limit comparison can
overlay the analytical curve:

.. code-block:: python

   # inside q_of_rho in pde_runner.py
   elif params.get('flux_type', '') == 'myflux':
       return q_myflux_python(rho, v_max, rho_max)

**Step 4 — rebuild**

.. code-block:: bash

   make


Adding a new initial condition
-------------------------------

Initial conditions are applied in ``src/fortran/pde_module.f90`` inside the
``initialise_pde`` subroutine (around line 154).  The ``select case`` block
dispatches on ``params%ic_type``.

**Step 1 — add a new** ``case`` **branch**:

.. code-block:: fortran

   case ('myic')
     ! fill rho(1:M) with the desired profile
     do i = 1, params%M
       rho(i) = ...
     end do

The array ``rho`` has already been allocated to size ``params%M`` when this
block is reached.  Use ``params%rho_left_bc``, ``params%rho_right_bc``,
``params%v_max``, and ``params%rho_max`` as needed.

**Step 2 — document the new string** in ``docs-sphinx/user_guide/pde_model.rst``
under the *Initial conditions* section, following the style of the existing
entries.

**Step 3 — rebuild** with ``make``.


Adding a new network preset
----------------------------

Network presets are Python functions in ``src/python/road_network.py``.  Each
preset returns a ``(NetworkSpec, LayoutSpec)`` tuple and is registered in the
``PRESETS`` dict at the bottom of the file so the GUI and runner can find it by
name.

For the data model and wiring rules a preset builder needs to satisfy
(lane flow direction, junction legs, routing matrices, perimeter ports),
see :doc:`custom_networks`. Promoting a one-off custom topology to a
preset is mostly a matter of wrapping it in a function and registering
it in ``PRESETS``.

**Step 1 — write the preset function**

Follow the signature and return types of an existing preset such as
``t_junction`` (line ~203):

.. code-block:: python

   def my_network(
       L: int = 20, alpha: float = 0.4, beta: float = 0.5
   ) -> Tuple[NetworkSpec, LayoutSpec]:
       """One-line description.

       Parameters
       ----------
       L : int
           Road segment length (cells).
       alpha, beta : float
           Default entry/exit probabilities for open ends.

       Returns
       -------
       NetworkSpec, LayoutSpec
       """
       # Build roads, junctions, and routing matrices.
       # Call validate(spec) before returning.
       ...
       validate(spec)
       return spec, layout

Key rules:

* Every road end that connects to a junction must appear in that junction's
  routing matrix.
* Every open road end (no junction) must have ``alpha`` / ``beta`` set on the
  relevant lane.
* Call ``validate(spec)`` — it checks routing matrix rows sum to 1 and that
  lane indices are consistent.

**Step 2 — register in** ``PRESETS``:

.. code-block:: python

   PRESETS = {
       ...
       "my_network": my_network,
   }

The key is the string the GUI and ``run_simulation`` use to look up the preset.

**Step 3 — add a layout entry** to the GUI's network widget if you want it to
appear in the preset dropdown; the widget reads from ``PRESETS`` automatically,
so no extra GUI code is needed.

**Step 4 — verify** by running the GUI or calling the preset directly:

.. code-block:: python

   from python.road_network import my_network
   spec, layout = my_network()
   print(f"{len(spec['roads'])} roads, {len(spec['junctions'])} junctions")


Code style conventions
-----------------------

**Fortran**

* Use ``implicit none`` in every module and subroutine.
* Declare all variables with explicit kinds; use the ``dp`` parameter
  (``real64``) defined in ``pde_module.f90`` for floating-point quantities.
* Pure functions should be marked ``pure``; elemental functions ``elemental``.
* Write docstring-style comments above each subroutine/function using the
  ``!>`` prefix so that ``sphinx-fortran`` picks them up automatically.

**Python**

* Follow NumPy docstring style (Parameters / Returns sections) — the Sphinx
  build uses the Napoleon extension to render these.
* Keep plotting functions in ``src/python/pde_visualisation.py`` or
  ``src/python/CA_visualisation.py`` as appropriate; keep data-generation
  and analysis logic in ``src/python/analysis.py`` or ``src/python/pde_runner.py``.
* New scripts that wrap the Fortran binaries go in ``scripts/``; new Python
  modules that provide reusable functionality go in ``src/python/``.
