# 1D Toy Model — Developer Documentation

This document covers everything currently implemented for the 1D toy model: the
TASEP physics, the Fortran simulation core, the NetCDF I/O layer, the Python
analysis and visualisation, and the entry-point scripts.

---

## 1. Physics: open-boundary 1D TASEP

The model is a **Totally Asymmetric Simple Exclusion Process** (TASEP) on a 1D
lattice of length `L` with open boundaries.

- Each site is either empty (`0`) or occupied (`1`).
- Time advances in discrete steps with **parallel update** (all sites updated
  simultaneously based on the state at the start of the step).

### Update rules (per step)

1. **Exit:** if site `L` is occupied, the particle leaves with probability β.
2. **Bulk:** every occupied site `i` (for `1 ≤ i < L`) hops to `i+1` if
   `i+1` is empty.
3. **Entry:** if site `1` is empty (after the bulk step), a new particle
   enters with probability α.

### Theoretical phases

| Phase | Condition | Bulk density ρ | Current J |
|-------|-----------|----------------|-----------|
| Low density (LD) | α < β, α < 0.5 | α | α |
| High density (HD) | β < α, β < 0.5 | 1−β | β |
| Maximum current (MC) | α ≥ 0.5, β ≥ 0.5 | 0.5 | 0.25 |

### Fundamental diagram

For the **parallel-update** TASEP (equivalently, deterministic cellular
automaton Rule 184), the current–density relation is the **triangular** curve:

> **J = min(ρ, 1 − ρ)**   (peaks at J = 0.5 when ρ = 0.5)

This is **not** the parabola J = ρ(1−ρ), which applies to the continuous-time
sequential TASEP. The distinction matters for comparing simulation output
against theory.

Note: in this implementation J is measured as *exits per step* (a raw count of
0 or 1 per step), not as flux per lattice site. In the LD phase the exit rate
equals the entry rate α, so J ≈ α ≈ ρ for small ρ, consistent with the
triangular branch.

---

## 2. Repository layout (1D toy model only)

```
src/fortran/
├── tasep.f90                   # core physics module
├── simulation.f90              # time-evolution driver + steady-state measurer
├── io.f90                      # NetCDF writers (simulation + fundamental diagram)
├── test_simulation.f90         # end-to-end simulation program
├── fundamental_diagram.f90     # α/β sweep program → fundamental_diagram.nc
└── test_tasep.f90              # unit-level test of tasep_model

src/python/
├── io.py                       # NetCDF reader
├── visualisation.py            # space-time, density, current, fundamental diagram
└── analysis.py                 # Python TASEP + parameter sweep

scripts/
├── run_toy_model.py            # build, run, plot summary figure
├── run_fundamental_diagram.py  # Python-based α/β sweep, plot J vs ρ
└── plot_fundamental_diagram.py # plot fundamental_diagram.nc from Fortran sweep

data/output/                    # NetCDF output written here
plots/                          # generated figures (gitignored)
Makefile                        # builds Fortran binaries
requirements.txt                # numpy, matplotlib, netCDF4
```

---

## 3. Fortran core

### `tasep.f90` — module `tasep_model`

Public routines:

| Routine | Purpose |
|---------|---------|
| `initialise_lattice(state, L)` | Sets all `L` sites to empty (`state = 0`). |
| `tasep_step(state, L, alpha, beta, exit_count)` | One parallel-update step. Returns updated `state` and `exit_count` ∈ {0, 1}. |
| `count_occupied(state, L)` | Returns total number of occupied sites. |
| `compute_density(state, L)` | Returns ρ = N/L as a real. |

The step routine builds an `old_state` snapshot first, then applies the three
rules in order using `old_state` for all decisions, ensuring a true parallel
update.

### `simulation.f90` — module `simulation`

Public routines:

```fortran
subroutine run_simulation(L, n_steps, alpha, beta, &
                          history, density_history, current_history, total_exits)
```

- Loops `tasep_step` for `n_steps` iterations.
- `history(i, t)` — site `i` at step `t` (Fortran column-major layout).
- `density_history(t)` — ρ at step `t`.
- `current_history(t)` — exits at step `t` (0 or 1).
- `total_exits` — cumulative count.
- All output arrays must be `allocatable` and allocated by the caller before the call.

```fortran
subroutine measure_steady_state(L, n_burnin, n_measure, alpha, beta, &
                                 mean_density, mean_current)
```

Designed for parameter sweeps where storing the full history is unnecessary:

- Runs `n_burnin` steps to reach steady state (discards all output).
- Runs `n_measure` steps accumulating mean density and mean current.
- Density is measured on the **central L/4 slice**
  (`sites 3L/8+1 … 5L/8`, matching Python `slice(3*L//8, 5*L//8)`) to avoid
  boundary-layer bias.
- `mean_current` is mean exits per step.

### `io.f90` — module `tasep_io`

#### `write_netcdf` — simulation output

```fortran
subroutine write_netcdf(filename, L, n_steps, alpha, beta, &
                        history, density_history, current_history)
```

NetCDF schema:

| Item | Kind | Notes |
|------|------|-------|
| `site` | dimension | length `L` |
| `time` | dimension | length `n_steps` |
| `history` | `int(site, time)` | lattice occupancy |
| `density` | `float(time)` | mean density per step |
| `current` | `int(time)` | exits per step |
| `L`, `n_steps`, `alpha`, `beta`, `model` | global attributes | |

#### `write_fundamental_diagram_netcdf` — sweep output

```fortran
subroutine write_fundamental_diagram_netcdf(filename, L, n_points, n_burnin, n_measure, &
                                             fixed_beta, fixed_alpha, &
                                             alpha_param, alpha_rho, alpha_J, &
                                             beta_param,  beta_rho,  beta_J)
```

NetCDF schema:

| Item | Kind | Notes |
|------|------|-------|
| `n_points` | dimension | points per sweep branch |
| `alpha_param` | `float(n_points)` | α values swept (β held at `fixed_beta`) |
| `alpha_rho` | `float(n_points)` | mean bulk density (α sweep) |
| `alpha_J` | `float(n_points)` | mean current (α sweep) |
| `beta_param` | `float(n_points)` | β values swept (α held at `fixed_alpha`) |
| `beta_rho` | `float(n_points)` | mean bulk density (β sweep) |
| `beta_J` | `float(n_points)` | mean current (β sweep) |
| `L`, `n_points`, `n_burnin`, `n_measure`, `fixed_beta`, `fixed_alpha` | global attributes | |

### `test_simulation.f90` — simulation driver program

Accepts command-line arguments (all optional, with defaults):

```bash
./build/test_simulation [L] [n_steps] [alpha] [beta]
# defaults: L=10  n_steps=100  alpha=0.5  beta=0.5
```

Allocates history arrays at runtime, calls `run_simulation`, prints per-step
state to stdout, and writes `data/output/simulation.nc`.

### `fundamental_diagram.f90` — sweep program

Accepts command-line arguments (all optional):

```bash
./build/fundamental_diagram [L] [n_points] [n_measure]
# defaults: L=100  n_points=30  n_measure=3000
# n_burnin is set automatically to 2*L^2
```

- Alpha sweep: varies α ∈ [0.02, 0.98] with β = 0.5 fixed → LD branch and MC peak.
- Beta sweep: varies β ∈ [0.02, 0.98] with α = 0.5 fixed → HD branch and MC peak.
- Writes `data/output/fundamental_diagram.nc`.

### Build

```bash
make              # build both binaries (test_simulation + fundamental_diagram)
make run          # build + run test_simulation with defaults
make run-fd       # build + run fundamental_diagram sweep with defaults

# Override simulation parameters (any subset):
make run L=50 N_STEPS=500 ALPHA=0.3 BETA=0.7
```

`gfortran -Wall -O2` with NetCDF Fortran bindings via `nf-config`.

---

## 4. Python layer

### `src/python/io.py`

```python
load_netcdf(filename) -> dict
```

Returns `{L, n_steps, alpha, beta, history, density, current}`.

**Important:** `history` is transposed on load. Fortran writes the array as
`(site, time)` in column-major order, but `netCDF4` returns it in C order, so
the raw read is shape `(time, site)`. The loader transposes back so downstream
code sees `(L, n_steps)` as expected.

### `src/python/visualisation.py`

| Function | Output |
|----------|--------|
| `plot_spacetime(data, ax, title)` | 2D occupancy grid: site (y) vs time (x), black = occupied. |
| `plot_density(data, ax, title)` | ρ vs time line plot with time-mean overlay. |
| `plot_current(data, ax, title)` | Exits-per-step bar chart with mean overlay. |
| `plot_fundamental_diagram(rho, J, ax, title)` | Scatter of (ρ, J) against the triangular theory curve J = min(ρ, 1−ρ). |
| `plot_summary(data, save_path)` | 3-panel composite: space-time on top, density and current below. |

### `src/python/analysis.py`

A standalone NumPy implementation of TASEP, separate from the Fortran binary,
used for parameter sweeps from Python. The update rules match the Fortran
exactly.

| Function | Purpose |
|----------|---------|
| `tasep_step(state, alpha, beta, rng)` | One parallel-update step (NumPy-vectorised bulk hops). |
| `run_tasep(L, n_steps, alpha, beta, burnin, seed, bulk_slice)` | Run for `burnin + n_steps` steps, return density and current arrays for the post-burnin period. `bulk_slice` lets you measure density only in a chosen region. |
| `fundamental_diagram(L, n_steps, burnin, n_points, seed)` | Sweeps α with β=0.5, then β with α=0.5; returns `(rho_vals, J_vals)`. Defaults: `burnin = 2*L²`, `bulk_slice = slice(3*L//8, 5*L//8)`. |

---

## 5. Entry-point scripts

### `scripts/run_toy_model.py`

Wraps the full pipeline: builds the Fortran binary → runs simulation →
loads NetCDF → plots summary figure.

```bash
python scripts/run_toy_model.py              # full run
python scripts/run_toy_model.py --no-run     # skip Fortran, replot existing .nc
python scripts/run_toy_model.py --save       # also write plots/summary.png
```

### `scripts/run_fundamental_diagram.py`

Runs the **Python** TASEP sweep (via `analysis.fundamental_diagram`) and
produces J vs ρ figure.

```bash
python scripts/run_fundamental_diagram.py
python scripts/run_fundamental_diagram.py --L 100 --points 40 --save
```

Flags: `--L` (lattice size, default 50), `--points` (per-branch resolution,
default 30), `--save` (write `plots/fundamental.png`).

### `scripts/plot_fundamental_diagram.py`

Reads `data/output/fundamental_diagram.nc` produced by the **Fortran** sweep
(`make run-fd`) and plots J vs ρ using `visualisation.plot_fundamental_diagram`.

```bash
python scripts/plot_fundamental_diagram.py
python scripts/plot_fundamental_diagram.py --save   # write plots/fundamental_fortran.png
python scripts/plot_fundamental_diagram.py --input path/to/other.nc
```

---

## 6. Pipeline overview

```
# Single simulation run
Fortran binary  ──►  data/output/simulation.nc  ──►  Python loader  ──►  plots/summary.png
   (tasep rules)       (history, ρ, J + attrs)        (io.py)            (run_toy_model.py)

# Fundamental diagram — Fortran sweep (faster)
make run-fd  ──►  data/output/fundamental_diagram.nc  ──►  plot_fundamental_diagram.py

# Fundamental diagram — Python sweep
run_fundamental_diagram.py  ──►  plots/fundamental.png
```

---

## 7. Conventions and gotchas

- **Array order:** Fortran is column-major, NumPy is row-major. The transpose
  in `io.py` is what makes `history` line up — don't remove it.
- **`history(i, t)`** in Fortran corresponds to `history[i, t]` in Python
  (shape `(L, n_steps)`).
- **Theory curve:** the parallel-update TASEP follows J = min(ρ, 1−ρ)
  (triangular), not the parabola J = ρ(1−ρ). The two agree only in the
  LD and HD phases for small α/β.
- **Burn-in:** equilibration time scales as ~L² near the max-current phase
  boundary. Discard ≥ L² steps before recording statistics.
  `measure_steady_state` and `fundamental_diagram` both default to
  `n_burnin = 2*L²`.
- **Boundary layers:** the spatial average of density differs from the bulk
  by an O(1/L) correction. Both `measure_steady_state` and the Python
  `fundamental_diagram` measure density over the central L/4 slice
  (`sites 3L/8+1 … 5L/8`) to avoid this bias.

---

## 8. Setup

```bash
# Fortran
make                                    # needs gfortran + netCDF Fortran (nf-config)

# Python
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# Single simulation run
make run
python scripts/run_toy_model.py --save

# Fundamental diagram (Fortran sweep, then plot)
make run-fd
python scripts/plot_fundamental_diagram.py --save

# Fundamental diagram (Python sweep)
python scripts/run_fundamental_diagram.py --save
```

Output figures land in `plots/`, NetCDF data in `data/output/`.
