# PX915 — Traffic Flow Simulator

A dual-model traffic flow simulator implementing both a microscopic cellular automaton model (Nagel–Schreckenberg / TASEP) and a macroscopic continuum PDE model (Lighthill–Whitham–Richards). The project is aimed at research and coursework in computational physics.

The computational core is written in Fortran for performance, with a Python frontend providing visualisation, a PyQt5 GUI, and a Jupyter tutorial.

---
## Tutorial & Documentation

After successful installation find the `tutorial.ipynb` at `PX915/tutorial.ipynb`.

Docs can be opened with:
```bash
open docs-sphinx/_build/html/index.html
```

---

## Features

- **Cellular automaton models** — TASEP and Nagel–Schreckenberg with configurable vehicle density, braking probability, and speed limit
- **LWR PDE solver** — finite-volume scheme with Greenshields and Newell–Daganzo flux closures; Godunov and Lax–Friedrichs schemes
- **Multi-lane roads** — per-lane parameters and conservative lane-change source terms for both CA and PDE
- **Network topologies** — single lane, two-lane, T-junction, crossroads, roundabout, and 2×2 town grid
- **Fundamental diagram sweeps** — OpenMP-parallelised steady-state α/β sweeps for TASEP and NS
- **Interactive PyQt5 GUI** — real-time parameter editing, live simulation progress, and tabbed result visualisation
- **NetCDF I/O** — configs and results stored as self-describing `.nc` files
- **API documentation** - Detailed information for both users and developers

---

## Installation

It is highly reccomended to read through the README and to refer to the Documentation for a clean installation.

However, for the most impatient of users (and those with `gfortran`, `netcdf`, and `Python >= 3.9`), the easiest way to jump straight into the software package is as follows:

```bash
git clone https://github.com/JasperA680/PX915.git   # Copy the repo locally
cd PX915    
python3 -m venv .venv                               # Create a virtual environment
source .venv/bin/activate           # macOS / Linux
# .venv\Scripts\activate            # Windows
make                                                # Install Python deps, register Jupyter kernel, compile Fortran binaries
python3 ./scripts/run_gui.py                        # Open the GUI
```
---

## Repository Structure

```
PX915/
├── Makefile                  # Builds Fortran binaries
├── requirements.txt          # Python dependencies
├── src/
│   ├── fortran/              # 21 Fortran source files (CA, PDE, FD sweep)
│   └── python/               # Python modules + PyQt5 GUI
│       └── gui/              # GUI application
├── tutorial.ipynb            # Interactive tutorial
├── scripts/                  # Standalone runners and plotters
├── docs-sphinx/              # Sphinx documentation source
├── data/
│   ├── input/                # Example input files
│   └── output/               # Simulation results (NetCDF)
└── tests/                    # UNSTABLE - Unit and integration tests

```

---

## Prerequisites

**Fortran**
- `gfortran` with OpenMP support
- NetCDF-Fortran and NetCDF-C libraries (located via `nf-config` and `nc-config`)

**Python ≥ 3.9**

| Package | Minimum version |
|---------|----------------|
| numpy | 2.0 |
| matplotlib | 3.9 |
| netCDF4 | 1.7 |
| PyQt5 | 5.15 |
| ipykernel | 6.31.0 |

Optional (documentation only): `sphinx>=7.0`, `sphinx-fortran>=1.1`, `sphinx-rtd-theme>=2.0`

Further information about prerequisites and dealing with issues can be found in the docs, under `User Guide/Installation/Prerequisites`.

---

## Quick Start

There are 3 main ways of interacting with the software package:

### Interactive GUI

```bash
python scripts/run_gui.py
```

The window has two tabs:

- **CA tab** — choose a network preset (single lane through town grid), select the NS or TASEP model, adjust parameters, and run the simulation. Results appear as space-time diagrams, density traces, and a network layout heatmap.
- **PDE tab** — configure the LWR solver (grid size, steps, flux type, boundary/initial conditions, lane count), run, and inspect density fields and fundamental diagrams.

### Jupyter Tutorial

```bash
jupyter notebook tutorial.ipynb
```

Covers LWR theory, Riemann problems, fundamental diagrams, multi-lane extensions, and cellular automaton models with live plots.

### Python API

The Jupyter tutorial will give you an introduction to how the Python API is used. Further information for this can be found in the documentation, under `User Guide/Quickstart`.

Even more detailed information about the Python API can also be found in the documentation, under `API reference/Python API`.


---


## Output Files

All results are written to `data/output/` as NetCDF (`.nc`) files. Use the Python loading utilities to read them back:

```python
from python.io import load_network_netcdf      # CA network results
from python.pde_runner import load_pde_netcdf  # PDE results
from python.fd_runner import load_fd_netcdf    # Fundamental diagram sweeps
```

---

## Tests

The folder `PX915/tests/` contains a number of test files, in both Fortran and Python, and a Makefile to help build them. **These files should be considered unstable and not be used by regular users.** They are included as they may be helpful for developers with a thorough understanding of the code.
