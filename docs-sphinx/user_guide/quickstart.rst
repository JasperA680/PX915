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

The project includes a Python GUI for running and visualising simulations.

From the project root, run:

.. code-block:: bash

   python src/python/gui.py

If the GUI file has a different name, replace ``gui.py`` with the correct script
name in ``src/python/``.

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