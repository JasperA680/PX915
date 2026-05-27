# PX915 — Traffic Flow Simulator

A dual-model traffic flow simulator combining a **microscopic cellular automaton** (Nagel–Schreckenberg / TASEP) with a **macroscopic continuum PDE** (Lighthill–Whitham–Richards). The project is aimed at research and coursework in computational physics.

The computational core is written in Fortran for performance, with a Python frontend providing visualisation, a PyQt5 GUI, and a Jupyter tutorial.

---

## Features

- **Cellular automaton models** — TASEP and Nagel–Schreckenberg with configurable vehicle density, braking probability, and speed limit
- **LWR PDE solver** — finite-volume scheme with Greenshields and Newell–Daganzo flux closures; Godunov and Lax–Friedrichs schemes
- **Multi-lane roads** — per-lane parameters and conservative lane-change source terms for both CA and PDE
- **Network topologies** — single lane, two-lane, T-junction, crossroads, roundabout, and 2×2 town grid
- **Fundamental diagram sweeps** — OpenMP-parallelised steady-state α/β sweeps for TASEP and NS
- **Interactive PyQt5 GUI** — real-time parameter editing, live simulation progress, and tabbed result visualisation
- **NetCDF I/O** — configs and results stored as self-describing `.nc` files
- **Sphinx API documentation** and physics notes

---

## Repository Structure - MIGHT NEED CHANGING

```
PX915/
├── Makefile                  # Builds Fortran binaries
├── requirements.txt          # Python dependencies
├── src/
│   ├── fortran/              # 21 Fortran source files (CA, PDE, FD sweep)
│   └── python/               # Python modules + PyQt5 GUI
│       └── gui/              # GUI application
├── scripts/                  # Standalone runners and plotters
├── notebooks/
│   └── tutorial.ipynb        # Interactive tutorial
├── tests/                    # Unit and integration tests
├── docs-sphinx/              # Sphinx documentation source
├── docs/                     # Markdown physics notes
├── data/
│   ├── input/                # Example input files
│   └── output/               # Simulation results (NetCDF)
└── plots/                    # Reference figures
```

---

## Prerequisites

**Fortran**
- `gfortran` with OpenMP support
- NetCDF-Fortran and NetCDF-C libraries (located via `nf-config` and `nc-config`)

**Python ≥ 3.10**

| Package | Minimum version |
|---------|----------------|
| numpy | 2.0 |
| matplotlib | 3.9 |
| netCDF4 | 1.7 |
| PyQt5 | 5.15 |

Optional (documentation only): `sphinx>=7.0`, `sphinx-fortran>=1.1`, `sphinx-rtd-theme>=2.0`

---

## Documentation

Build the Sphinx HTML docs:

```bash
cd docs-sphinx
make clean && make html         # Removes any existing files and makes the docs
Open _build/html/index.html     # Opens the docs in a browser
```

Physics background (Markdown): !! IS THIS STILL NEEDED

| File | Content |
|------|---------|
| `docs/pde_model.md` | LWR continuum model and Greenshields flux |
| `docs/multilane_pde_model.md` | Multi-lane PDE with lane-change coupling |
| `docs/1d_toy_model.md` | TASEP and Nagel–Schreckenberg introduction |

---

## Installation

Please find a thorough installation guide in the documentation, under `User guide/Installation`.

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
jupyter notebook notebooks/tutorial.ipynb
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
