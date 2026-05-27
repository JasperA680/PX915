Installation
============

Prerequisites
-------------

The following must be installed before building the project.

**Fortran compiler**

The code requires ``gfortran`` (GCC Fortran compiler):

.. code-block:: bash

   # macOS
   brew install gcc

   # Ubuntu / Debian
   sudo apt install gfortran

**NetCDF Fortran library**

The solver reads and writes NetCDF files.  Both the C library and the Fortran
wrapper are required.  The Makefile detects them automatically via ``nf-config``
and ``nc-config``:

.. code-block:: bash

   # macOS
   brew install netcdf

   # Ubuntu / Debian
   sudo apt install libnetcdf-dev libnetcdff-dev

.. note::

   ``libnetcdff-dev`` (double-f) is the Fortran wrapper — it is a separate
   package from the C library and is easy to miss.

**Python 3**

Python 3.9 or later is required.  Check your version with:

.. code-block:: bash

   python3 --version


Getting the code
----------------

Clone the repository:

.. code-block:: bash

   git clone https://github.com/JasperA680/PX915.git
   cd PX915


Installing Python dependencies
-------------------------------

Create and activate a virtual environment, then install the package and its
dependencies in one step:

.. code-block:: bash

   python3 -m venv .venv
   source .venv/bin/activate   # macOS / Linux
   # .venv\Scripts\activate    # Windows

   pip install -e .

To also install the optional documentation-building dependencies:

.. code-block:: bash

   pip install -e ".[docs]"

This installs:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Package
     - Purpose
   * - ``numpy``
     - Array operations and data processing
   * - ``matplotlib``
     - Plotting space–time diagrams and fundamental diagrams
   * - ``netCDF4``
     - Reading simulation output files written by the Fortran solver
   * - ``PyQt5``
     - Graphical user interface for running simulations interactively
   * - ``sphinx``, ``sphinx-fortran``, ``sphinx-rtd-theme``
     - Building this documentation (optional, installed with ``.[docs]``)


Building the Fortran solver
----------------------------

From the repository root:

.. code-block:: bash

   make

This compiles all executables and places them in ``build/``:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Binary
     - Description
   * - ``build/pde_solver``
     - LWR PDE traffic-flow solver (single-lane and multilane)
   * - ``build/test_simulation``
     - TASEP cellular automaton simulation driver
   * - ``build/fundamental_diagram``
     - Fundamental diagram sweep utility
   * - ``build/run_network``
     - Road network simulation driver

To build only the PDE solver:

.. code-block:: bash

   make pde


Verifying the installation
---------------------------

Run the PDE solver with default parameters:

.. code-block:: bash

   make run-pde

A successful run prints a parameter summary followed by per-step progress and
writes output to ``data/output/pde_simulation.nc``:

.. code-block:: text

   === LWR PDE Solver (multi-lane) ======================
   M                = 200
   n_steps          = 500
   ...
   Done.

If the build fails with an error about ``nf-config`` not found, the NetCDF
Fortran library is not installed or not on your ``PATH`` — revisit the
NetCDF step above.
