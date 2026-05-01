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
| Low density (LD) | α < β, α < 0.5 | α | α(1−α) |
| High density (HD) | β < α, β < 0.5 | 1−β | β(1−β) |
| Maximum current (MC) | α ≥ 0.5, β ≥ 0.5 | 0.5 | 0.25 |

In all three phases the current lies on the parabola **J = ρ(1−ρ)**, which
is the basis of the fundamental diagram.

---

## 2. Repository layout (1D toy model only)

```
src/fortran/
├── tasep.f90              # core physics module
├── simulation.f90         # time-evolution driver
├── io.f90                 # NetCDF writer
├── test_simulation.f90    # end-to-end test program
└── test_tasep.f90         # unit-level test of tasep_model

src/python/
├── io.py                  # NetCDF reader
├── visualisation.py       # space-time, density, current, fundamental diagram
└── analysis.py            # Python TASEP + parameter sweep

scripts/
├── run_toy_model.py       # build, run, plot summary figure
└── run_fundamental_diagram.py  # sweep α and β, plot J vs ρ

data/output/               # NetCDF output written here
plots/                     # generated figures (gitignored)
Makefile                   # builds Fortran binary
requirements.txt           # numpy, matplotlib, netCDF4
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

Public routine:

```fortran
subroutine run_simulation(L, n_steps, alpha, beta, &
                         history, density_history, current_history, total_exits)
```

- Loops `tasep_step` for `n_steps` iterations.
- `history(i, t)` — site `i` at step `t` (Fortran column-major layout).
- `density_history(t)` — ρ at step `t`.
- `current_history(t)` — exits at step `t` (0 or 1).
- `total_exits` — cumulative count.

### `io.f90` — module `tasep_io`

Public routine:

```fortran
subroutine write_netcdf(filename, L, n_steps, alpha, beta, &
                        history, density_history, current_history)
```

NetCDF schema:

| Item | Kind | Notes |
|------|------|-------|
| `site` | dimension | length `L` |
| `time` | dimension | length `n_steps` |
| `history` | variable, `int(site, time)` | lattice occupancy |
| `density` | variable, `float(time)` | mean density per step |
| `current` | variable, `int(time)` | exits per step |
| `L`, `n_steps`, `alpha`, `beta`, `model` | global attributes | for round-trip |

### `test_simulation.f90` — entry program

Compile-time parameters (currently `L=10`, `n_steps=100`, `α=β=0.5`). Calls
`run_simulation`, prints per-step state to stdout, and writes
`data/output/simulation.nc`. Edit the `parameter` lines to change the run.

### Build

```bash
make run        # build + run
make            # build only
make clean      # remove build/
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
| `plot_fundamental_diagram(rho, J, ax, title)` | Scatter of (ρ, J) against the theoretical J = ρ(1−ρ) parabola. |
| `plot_summary(data, save_path)` | 3-panel composite: space-time on top, density and current below. |

### `src/python/analysis.py`

A standalone NumPy implementation of TASEP, separate from the Fortran binary,
used for parameter sweeps where recompiling each time would be impractical.
The update rules match the Fortran exactly.

| Function | Purpose |
|----------|---------|
| `tasep_step(state, alpha, beta, rng)` | One parallel-update step (NumPy-vectorised bulk hops). |
| `run_tasep(L, n_steps, alpha, beta, burnin, seed, bulk_slice)` | Run for `burnin + n_steps` steps, return density and current arrays for the post-burnin period. `bulk_slice` lets you measure density only in a chosen region. |
| `fundamental_diagram(L, n_steps, burnin, n_points, seed)` | Sweeps α with β=1, then β with α=1, returns `(rho_vals, J_vals)`. Defaults: `burnin = 2*L²` (auto), `bulk_slice = slice(L//4, 3*L//4)` to avoid boundary-layer bias. |

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

Runs the Python TASEP sweep and produces J vs ρ figure.

```bash
python scripts/run_fundamental_diagram.py
python scripts/run_fundamental_diagram.py --L 100 --points 40 --save
```

Flags: `--L` (lattice size, default 50), `--points` (per-branch sweep
resolution, default 30), `--save` (write `plots/fundamental.png`).

---

## 6. Pipeline overview

```
Fortran binary  ──►  data/output/simulation.nc  ──►  Python loader  ──►  plots
   (rules)            (history, ρ, J + attrs)         (transpose)         (figures)
```

For the fundamental diagram the Fortran path is bypassed — the sweep runs
entirely in Python so α and β can be varied without recompiling.

---

## 7. Conventions and gotchas

- **Array order:** Fortran is column-major, NumPy is row-major. The transpose
  in `io.py` is what makes everything line up — don't remove it.
- **`history(i, t)`** in Fortran corresponds to `history[i, t]` in Python
  (i.e. shape `(L, n_steps)`).
- **Compile-time parameters:** `L`, `n_steps`, `α`, `β` are baked into the
  Fortran binary. Changing them requires editing `test_simulation.f90` and
  rebuilding. This will be replaced with runtime input when the road network
  lands.
- **Burn-in:** equilibration time scales as ~L² near the max-current phase
  boundary. For statistics in the steady state, discard the first ~L² steps.
  `fundamental_diagram` does this automatically via `burnin = 2*L²`.
- **Boundary layers:** the spatial average of density differs from the bulk
  density by an O(1/L) correction. For comparing against J = ρ(1−ρ), measure
  ρ from the middle of the lattice (the default `bulk_slice` in
  `fundamental_diagram`).

---

## 8. Setup

```bash
# Fortran
make                                    # needs gfortran + netCDF Fortran (nf-config)

# Python
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# Run
make run
python scripts/run_toy_model.py --save
python scripts/run_fundamental_diagram.py --save
```

Output figures land in `plots/`, NetCDF data in `data/output/`.
