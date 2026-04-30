# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

Traffic flow simulation software package for modelling vehicle dynamics on road networks. PX915 group project at the University of Warwick (Tristan McCarthy, Jasper Allen, Stephan Gambert, Lucas Belz-Koeling). **Deadline: 28 May 2026.**

The core is a 1D TASEP (Totally Asymmetric Simple Exclusion Process) cellular automaton model, implemented in Fortran. Python scripts and notebooks are currently empty stubs.

### Feature roadmap (from SDP)

| Priority | Feature |
|----------|---------|
| Basic | 1D TASEP CA model (done), configurable α/β/L, output density/flow/travel time/fundamental diagrams |
| Specific | Road network bifurcations & merges, 2D multi-lane model, lane changing & overtaking |
| Extension | PDE continuum model, comparison of PDE vs CA outputs |

### 5-week work plan

| Week | Dates | Tasks |
|------|-------|-------|
| 1 | 28 Apr | 1D toy model, road network |
| 2 | 05 May | Bifurcations, NetCDF I/O, Python front end |
| 3 | 12 May | Road merges, sensitivity analysis, PDE reading |
| 4 | 19 May | PDE implementation, 2D multi-lane |
| 5 | 26 May | Lane changing, overtaking logic, presentation prep |

## Build & run (Fortran)

```bash
# Build and run the simulation test
make run

# Build only
make

# Clean build artifacts
make clean
```

The Makefile compiles with `gfortran -Wall -O2` and places the binary in `build/test_simulation`. There is no separate lint step.

## Fortran architecture

The Fortran code follows a layered module structure:

| File | Module | Role |
|------|--------|------|
| `src/fortran/tasep.f90` | `tasep_model` | Core TASEP physics: lattice init, one-step update, density/occupancy calculations |
| `src/fortran/simulation.f90` | `simulation` | Driver that loops `tasep_step` for `n_steps` and records full history arrays |
| `src/fortran/test_simulation.f90` | program | End-to-end test: runs `run_simulation` and prints per-step state, density, and current |
| `src/fortran/test_tasep.f90` | program | Unit-level test: exercises `tasep_model` directly |
| `src/fortran/utils.f90` | — | Empty stub |

**TASEP update rules** (parallel update, open boundaries):
1. Particle at site `L` exits with probability `beta`.
2. Bulk particles hop right if the next site is empty.
3. A new particle enters at site `1` with probability `alpha` if site `1` is empty.

**Key array convention**: `history(i, t)` — site index first, time second.

**Output quantities**: `density_history(t)` = ρ = N/L; `current_history(t)` = number of exits at step t; `total_exits` = cumulative exit count.

## Python layer (planned)

The Python layer provides the high-level interface for running simulations, analysis, and visualisation. Uses NumPy and Matplotlib; simulation data stored in **NetCDF** format.

Stub files exist under `src/python/` and `scripts/`. They are intended to:
- `run_simulation.py` — call the Fortran binary or wrap the model in Python
- `analysis.py` — density, flow, fundamental diagrams, travel time statistics
- `visualisation.py` — space-time diagrams, flow-density plots, density vs time
- `io.py` — read config files, write/read NetCDF output
- `road_network.py` — road geometry, connectivity, bifurcations and merges
- `scripts/run_toy_model.py`, `run_pde_model.py`, `run_network_sim.py` — experiment entry points

Python tests live in `tests/` (also empty stubs). No Python package config or `requirements.txt` exists yet.

## Architecture overview

Five logical modules (per SDP design):

| Module | Responsibility |
|--------|---------------|
| Simulation | TASEP/CA update rules, time evolution, boundary conditions |
| Road Network | Road geometry, cell arrays, bifurcation/merge connectivity |
| Analysis | Density, flow, fundamental diagrams, travel time |
| Input/Output | Config files, NetCDF storage, result reuse |
| Visualisation | Space-time diagrams, flow-density plots |

## Data

`data/input/` and `data/output/` are tracked as empty placeholders (`.gitkeep`). Simulation output should be written here in NetCDF format.
