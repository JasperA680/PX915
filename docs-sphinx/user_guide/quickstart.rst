Quickstart
==========

This page gives a brief introduction to running the traffic-flow simulation
code for the first time.

The project contains two main types of model:

* a one-dimensional cellular automaton model, including TASEP and
  Nagel-Schreckenberg style updates;
* a continuum PDE model, including single-lane and multilane LWR traffic-flow
  solvers.

From the project root directory, activate the Python virtual environment:

.. code-block:: bash

   source .venv/bin/activate

If the dependencies have not already been installed, install them with:

.. code-block:: bash

   pip install -r requirements.txt

Building the Fortran code
-------------------------

The Fortran simulation code can be built from the project root using:

.. code-block:: bash

   make

This should compile the Fortran modules and produce the relevant executable
programs. If the build fails, check that a Fortran compiler such as ``gfortran``
is available.

Running a 1D cellular automaton simulation
------------------------------------------

The one-dimensional cellular automaton model represents the road as a lattice of
empty or occupied cells. The TASEP model uses open boundaries with an entry
probability ``alpha`` and an exit probability ``beta``.

A typical workflow is:

#. choose the road length;
#. choose the number of time steps;
#. choose the boundary rates ``alpha`` and ``beta``;
#. run the simulation;
#. analyse the density and current histories.

The main Fortran routines for this model are documented in the
:doc:`../api/fortran` page under the cellular automata section.

Running a PDE simulation
------------------------

The PDE model evolves the traffic density as a continuous field using a
finite-volume method. The default closure is the Greenshields fundamental
diagram.

Typical quantities to set are:

* number of spatial cells;
* number of time steps;
* maximum velocity;
* maximum density;
* initial condition;
* boundary condition;
* numerical flux type.

The PDE solver writes output that can be visualised using the Python plotting
routines documented in the :doc:`../api/python` page.

Running the Python GUI
----------------------

A PyQt5 GUI bundles both models behind a single window so simulations can be
configured, run, and inspected interactively. The Fortran binaries (``make``)
must already be built — the GUI launches them as subprocesses.

From the project root, with the virtual environment active, run:

.. code-block:: bash

   python scripts/run_gui.py

The window has two top-level tabs:

* **Cellular Automaton (CA)** — choose a network preset
  (``single_lane``, ``two_lane``, ``t_junction``, ``crossroads``, ``roundabout``,
  ``town``) and a model (``NS`` or ``TASEP``); the network preview shows the
  layout, and clicking a road opens an editor for its α, β and length, while
  clicking a junction opens its routing matrix. Use the matplotlib toolbar
  above the preview to pan and zoom, or pick a road from the **Focus on**
  dropdown to jump-zoom to it — α/β/length labels only appear once you zoom
  in close enough to read them. A separate **Run FD sweep** button populates
  the fundamental-diagram tab from a parameter sweep; this is only meaningful
  for the ``single_lane`` preset, and the button is greyed out elsewhere.
* **PDE Continuum Model** — pick a single-lane Riemann preset
  (shock, rarefaction, …) or one of the multilane scenarios, or edit any
  parameter by hand; the form syncs ``n_steps`` and total time ``T`` via the
  CFL estimate.

Each tab has its own plot panel that refreshes after every run, with
several complementary views of the result; tabs that don't apply to the
current run (per-lane plots, mass-conservation diagnostics, the
fundamental-diagram view) stay disabled until the relevant data is
available.  A corner indicator marks results as *up to date*, *stale*,
*running…*, or empty; editing any parameter after a run flags the
displayed plots as stale until the next **Run**.

Every parameter input — spin boxes, sliders, and dropdowns alike — carries
an explanatory tooltip. Hover the cursor over any control to see what it
does and any caveats (for example, that ``v_limit`` should be set equal to
``v_max`` for the unrestricted Greenshields flux).

Output NetCDF files are written under ``data/output/gui/`` (CA) and
``data/output/gui_pde/`` (PDE); the bottom-of-window log dock streams the
solver's progress and prints per-road steady-state diagnostics on completion.

Building the documentation
--------------------------

The Sphinx documentation can be built from the ``docs-sphinx`` directory:

.. code-block:: bash

   cd docs-sphinx
   make clean
   make html

The rendered documentation can then be opened with:

.. code-block:: bash

   open _build/html/index.html

On Linux, use:

.. code-block:: bash

   xdg-open _build/html/index.html

or start a local server:

.. code-block:: bash

   cd _build/html
   python3 -m http.server 8000

and open ``http://localhost:8000`` in a web browser.

Where to go next
----------------

For a user-facing explanation of the continuum traffic model, see
:doc:`pde_model`.

For the multilane extension, see :doc:`multilane_pde`.

For implementation details, see :doc:`../developer_guide/numerical_methods`.

For automatically generated API documentation, see :doc:`../api/index`.