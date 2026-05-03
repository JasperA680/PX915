# PDE Model — Developer Documentation

This document covers the LWR (Lighthill–Whitham–Richards) PDE implementation:
the physics, the Fortran solver core, the NetCDF I/O layer, the Python
interface, and the test suite. It tracks **Phases 1 and 2** of `PDE.md`.

---

## 1. Physics: LWR scalar conservation law

The model is the **Lighthill–Whitham–Richards** traffic flow equation — a scalar
hyperbolic conservation law:

```
∂ρ/∂t + ∂q(ρ)/∂x = 0
```

where ρ(x, t) is vehicle density and q(ρ) is the traffic flow (vehicles per
unit time). The **Greenshields fundamental diagram** is used as the closure:

```
v(ρ)  = v_max · (1 − ρ/ρ_max)          (speed decreases linearly with density)
q(ρ)  = ρ · v(ρ) = v_max · ρ · (1 − ρ/ρ_max)    (concave parabola)
dq/dρ = v_max · (1 − 2ρ/ρ_max)         (characteristic speed; zero at ρ_c)
```

### Critical density

The flow q is concave and achieves its maximum at:

> **ρ_c = ρ_max / 2**   (critical density, argmax of q)

This divides the domain into two phases:

| Phase | Condition | dq/dρ | Behaviour |
|-------|-----------|-------|-----------|
| Free flow (subcritical) | ρ < ρ_c | > 0 | density waves travel forward |
| Congested (supercritical) | ρ > ρ_c | < 0 | density waves travel backward |

### Wave structures

For the Riemann problem (piecewise-constant initial data):

| ρ_L vs ρ_R | Wave type | Propagation |
|------------|-----------|-------------|
| ρ_L < ρ_R | **Shock** | speed s = (q(ρ_R) − q(ρ_L)) / (ρ_R − ρ_L) (Rankine–Hugoniot) |
| ρ_L > ρ_R | **Rarefaction** | fan of characteristics connecting ρ_L to ρ_R |
| ρ_L < ρ_c < ρ_R | Shock through critical | Godunov flux = min(q(ρ_L), q(ρ_R)) |
| ρ_R < ρ_c < ρ_L | Sonic rarefaction | Godunov flux = q(ρ_c) = v_max·ρ_max/4 |

---

## 2. Repository layout (PDE files)

```
src/fortran/
├── pde_flux.f90          # Greenshields flux functions + numerical flux schemes
└── pde_module.f90        # solver types, initialise/step/finalise, NetCDF output

src/fortran/              # (driver program)
└── pde_driver.f90        # command-line entry point; runs simulation + writes .nc

src/python/
└── pde_runner.py         # Python flux mirrors, run_pde(), load_pde_netcdf()

tests/
└── test_pde.py           # Phase 1 flux unit tests + Phase 4 solver tests (partial)

data/output/
└── pde_simulation.nc     # default output path
```

---

## 3. Fortran core

### `pde_flux.f90` — module `pde_flux`

All flux routines are **elemental** (work on scalars and arrays transparently).

| Routine | Signature | Purpose |
|---------|-----------|---------|
| `v_of_rho` | `(rho, v_max, rho_max)` | Greenshields speed v(ρ) |
| `q_of_rho` | `(rho, v_max, rho_max)` | Flow q(ρ) = ρ·v(ρ) |
| `dq_drho`  | `(rho, v_max, rho_max)` | Characteristic speed dq/dρ |
| `rho_critical` | `(rho_max)` | ρ_c = ρ_max/2 |

Two non-elemental **numerical flux** functions:

#### `lax_friedrichs_flux(rho_L, rho_R, v_max, rho_max, dx, dt)`

```
F_LF = [q(ρ_L) + q(ρ_R)] / 2  −  (Δx / 2Δt) · (ρ_R − ρ_L)
```

Stable for CFL ≤ 1 but adds numerical diffusion proportional to Δx/(2Δt).
Used as the default scheme in Phases 1–2 because it is simple and robust.

#### `godunov_flux(rho_L, rho_R, v_max, rho_max)`

Exact Godunov flux for a concave flux. Closed-form case table
(LeVeque, *Finite Volume Methods for Hyperbolic Problems*, §12.1):

```
If ρ_L ≤ ρ_R  (shock):
  both subcritical  (ρ_R ≤ ρ_c): F = q(ρ_L)
  both supercritical (ρ_L ≥ ρ_c): F = q(ρ_R)
  mixed (ρ_L < ρ_c < ρ_R):       F = min(q(ρ_L), q(ρ_R))

If ρ_L > ρ_R  (rarefaction):
  both subcritical  (ρ_L ≤ ρ_c): F = q(ρ_L)
  both supercritical (ρ_R ≥ ρ_c): F = q(ρ_R)
  sonic (ρ_R < ρ_c < ρ_L):       F = q(ρ_c)
```

Does not take Δt as an argument — this is a deliberate design choice: the
Godunov flux is independent of the time step, which makes it easier to use
with adaptive CFL.

---

### `pde_module.f90` — module `pde_solver`

#### Derived types

```fortran
type :: pde_params_t
  real    :: dx, dt, domain_length, v_max, rho_max
  real    :: rho_left_bc, rho_right_bc, cfl_number
  integer :: M, n_steps, C_checkpoint
  character(len=16) :: ic_type    ! "constant" | "riemann" | "gaussian" | "sine"
  character(len=16) :: bc_type    ! "open" | "periodic"
  character(len=16) :: flux_type  ! "lf" | "godunov"
  logical :: use_adaptive_dt
end type

type :: pde_state_t
  real, allocatable :: density(:)   ! physical cells 1:M
  real, allocatable :: rho_ext(:)   ! ghost cells 0:M+1  (scratch, pre-allocated)
  real, allocatable :: flux(:)      ! interface fluxes 0:M (scratch, pre-allocated)
  real    :: t_current
  integer :: step
end type
```

The scratch arrays `rho_ext` and `flux` are allocated once in `pde_initialise`
and reused at every step to avoid per-step allocation overhead.

#### `pde_initialise(state, params)`

Allocates all arrays and fills `state%density` from the chosen IC:

| `ic_type` | Initial condition |
|-----------|------------------|
| `"constant"` | Uniform density = `rho_left_bc` |
| `"riemann"` | `rho_left_bc` for x < L/2, `rho_right_bc` for x ≥ L/2 |
| `"gaussian"` | Base 0.2·ρ_max + Gaussian bump of amplitude 0.6·ρ_max, σ = 5% of domain |
| `"sine"` | Mean 0.5·ρ_max + sinusoidal perturbation of amplitude 0.15·ρ_max |

For the `"riemann"` IC, `rho_left_bc` and `rho_right_bc` double as the
left/right initial densities. This intentionally links the IC to the open
boundary values so that the inflow/outflow match the initial state.

#### `pde_step(state, params)`

One conservative finite-volume update:

```
ρ_i^{n+1} = ρ_i^n − (Δt/Δx) · [F_{i+1/2} − F_{i-1/2}]
```

where `F_{i+1/2}` is either `lax_friedrichs_flux` or `godunov_flux` depending
on `params%flux_type`.

**Boundary conditions** are applied via ghost cells before the flux loop:

| `bc_type` | Left ghost (`rho_ext(0)`) | Right ghost (`rho_ext(M+1)`) |
|-----------|--------------------------|------------------------------|
| `"open"` | `rho_left_bc` | `rho_right_bc` |
| `"periodic"` | `density(M)` | `density(1)` |

#### `compute_dt(state, params, dt_out)`

Adaptive CFL time step, **capped** at the fully conservative bound:

```
dt_adaptive = cfl · Δx / max_i |dq/dρ(ρ_i)|
dt_out      = min(dt_adaptive, cfl · Δx / v_max)
```

The cap is essential. Near the sonic point (ρ → ρ_c), dq/dρ → 0, so
`dt_adaptive` → ∞. The Lax–Friedrichs scheme's numerical diffusion term
is Δx/(2Δt); if Δt → ∞ this vanishes and LF degenerates to an unstable
centred scheme. The cap keeps Δt ≤ cfl·Δx/v_max, which is always
sufficient for stability because |dq/dρ| ≤ v_max for all ρ ∈ [0, ρ_max]
under the Greenshields model.

In practice, because |dq/dρ|_max = v_max (achieved at ρ = 0 or ρ = ρ_max),
the cap usually applies and the actual dt used is `cfl·Δx/v_max`.

#### `write_pde_netcdf(filename, params, density_history, flow_history)`

Writes all simulation output to a single NetCDF file.

**Schema:**

| Item | Kind | Notes |
|------|------|-------|
| `time` | unlimited dimension | length n_steps+1 (includes t=0) |
| `x` | dimension | length M |
| `time` | `float(time)` | physical time at each step |
| `x` | `float(x)` | cell-centre x-coordinates |
| `density` | `float(x, time)` | full density history |
| `flow` | `float(time)` | q(ρ_M) — flux at the right boundary |

Global attributes stored for exact restarts: `model`, `M`, `domain_length`,
`dx`, `dt`, `v_max`, `rho_max`, `rho_left_bc`, `rho_right_bc`, `ic_type`,
`bc_type`, `flux_type`, `n_steps`.

---

### `pde_driver.f90` — program `pde_driver`

Command-line entry point analogous to `test_simulation.f90`.

```bash
./build/pde_solver [M] [n_steps] [v_max] [rho_max] [rho_left] [rho_right] \
                   [ic_type] [flux_type] [bc_type] [output_file]
```

All arguments are positional and optional; defaults:

| Arg | Default | Notes |
|-----|---------|-------|
| M | 200 | number of cells |
| n_steps | 500 | time steps to run |
| v_max | 1.0 | free-flow speed |
| rho_max | 1.0 | jam density |
| rho_left | 0.1 | left BC / left Riemann state |
| rho_right | 0.9 | right BC / right Riemann state |
| ic_type | `riemann` | initial condition |
| flux_type | `lf` | numerical flux |
| bc_type | `open` | boundary condition type |
| output_file | `data/output/pde_simulation.nc` | NetCDF output path |

The driver:
1. Parses arguments and computes `dx = domain_length / M`.
2. Sets an initial `dt` from the CFL condition (`cfl · dx / v_max`).
3. Allocates `density_history(M, n_steps+1)` and `flow_history(n_steps+1)`.
4. Calls `pde_initialise`, then loops `pde_step`, updating `dt` adaptively
   each step via `compute_dt`.
5. Records the density field and right-boundary flow at every step.
6. Calls `write_pde_netcdf` and `pde_finalise`.

Progress is printed every 100 steps: step number, physical time, and mean density.

### Build

```bash
make pde              # build only build/pde_solver
make run-pde          # build + run with defaults

# Override any parameter:
make run-pde PDE_M=400 PDE_IC=riemann PDE_FLUX=godunov PDE_BC=open
make run-pde PDE_IC=gaussian PDE_BC=periodic PDE_STEPS=1000
```

Compiler: `gfortran -Wall -O2` plus NetCDF Fortran bindings via `nf-config`.
The three PDE source files must be compiled in order (pde_flux → pde_module →
pde_driver) because each `use`s the previous module.

---

## 4. Python layer

### `src/python/pde_runner.py`

#### Flux functions (Python mirrors)

```python
v_of_rho(rho, v_max, rho_max)   # Greenshields speed; accepts NumPy arrays
q_of_rho(rho, v_max, rho_max)   # flow q(ρ)
dq_drho(rho, v_max, rho_max)    # characteristic speed dq/dρ
rho_critical(rho_max)           # ρ_c = ρ_max / 2
```

These are pure NumPy and mirror the Fortran elemental functions exactly.
They are used directly by the test suite and by any Python-side analysis
that needs the flux without running the solver.

#### `run_pde(params, output_path, exe)`

Invokes the Fortran `pde_solver` binary as a subprocess.

```python
run_pde(
    params=dict(
        M=200, n_steps=500, v_max=1.0, rho_max=1.0,
        rho_left_bc=0.1, rho_right_bc=0.9,
        ic_type="riemann", flux_type="lf", bc_type="open",
    ),
    output_path="data/output/pde_simulation.nc",
)
```

All keys are optional; unspecified keys fall back to the defaults above.
Raises `RuntimeError` if the binary exits non-zero.

#### `load_pde_netcdf(path)`

Reads a PDE output file and returns:

```python
{
    "density": np.ndarray,  # shape (n_steps+1, M)  — time × space
    "flow":    np.ndarray,  # shape (n_steps+1,)
    "x":       np.ndarray,  # shape (M,)
    "time":    np.ndarray,  # shape (n_steps+1,)
    "attrs":   dict,        # all global attributes from the NetCDF header
}
```

**No transpose is applied to `density`.** The Fortran variable is defined
with dimensions `[x_dimid, time_dimid]` (Fortran order: x fastest), which
the C-convention NetCDF library stores as `(time, x)` in the file. Python's
`netCDF4` reads it directly as shape `(n_steps+1, M)`. Applying `.T` here
would give shape `(M, n_steps+1)` — a cell-indexed array — which is wrong
for time-series analysis. See §5 Gotchas for the contrast with the TASEP
`io.py` behaviour.

---

## 5. Tests — `tests/test_pde.py`

### Phase 1 — flux unit tests (always run, no binary required)

| Test | What it checks |
|------|----------------|
| `test_q_at_zero_density` | q(0) = 0 |
| `test_q_at_max_density` | q(ρ_max) = 0 |
| `test_q_at_critical_density` | q(ρ_c) = v_max · ρ_max / 4 (analytical maximum) |
| `test_dq_drho_at_critical_is_zero` | dq/dρ(ρ_c) = 0 |
| `test_dq_drho_positive_in_free_flow` | dq/dρ > 0 for ρ < ρ_c |
| `test_dq_drho_negative_in_congested` | dq/dρ < 0 for ρ > ρ_c |
| `test_q_is_concave` | all second differences of q are negative |
| `test_v_of_rho_decreasing` | v is strictly decreasing |
| `test_rho_critical_is_half_rho_max` | ρ_c = ρ_max/2 for several ρ_max values |
| `test_q_symmetric_about_critical` | q(ρ_c − δ) = q(ρ_c + δ) (parabola symmetry) |

### Phase 4 — solver integration tests (require `build/pde_solver`)

Tests are automatically skipped if the binary is absent. Run `make pde` first.

| Test | What it checks | Tolerance |
|------|----------------|-----------|
| `test_constant_solution_preserved` | Uniform IC stays uniform for 1000 steps with periodic BCs | 1 × 10⁻¹⁰ |
| `test_mass_conservation_periodic` | ∑ ρ_i · Δx constant for Gaussian IC with periodic BCs over 1000 steps | 1 × 10⁻⁴ (float32 limit) |

The mass conservation tolerance is 1×10⁻⁴ rather than the 1×10⁻¹⁰ quoted in
`PDE.md`. This is because the Fortran solver and NetCDF output both use
single precision (float32, ~7 significant digits). LF with periodic BCs is
exactly mass-conservative in exact arithmetic (the boundary flux telescopes),
so any deviation comes from float32 rounding in the stored NetCDF values.
If double precision is needed, change `NF90_FLOAT` to `NF90_DOUBLE` in
`write_pde_netcdf` and use `real(8)` throughout the Fortran.

```bash
.venv/bin/python -m pytest tests/test_pde.py -v
```

---

## 6. Pipeline overview

```
# Typical PDE run
make run-pde                         # build + run Fortran solver → data/output/pde_simulation.nc
python -c "
import sys; sys.path.insert(0,'src/python')
from pde_runner import load_pde_netcdf
d = load_pde_netcdf('data/output/pde_simulation.nc')
print(d['density'].shape, d['attrs']['ic_type'])
"

# Or drive from Python
from src.python.pde_runner import run_pde, load_pde_netcdf
run_pde(dict(M=400, n_steps=1000, ic_type='riemann', flux_type='godunov'))
data = load_pde_netcdf('data/output/pde_simulation.nc')
# data['density'] is (1001, 400) — (time, x)
```

```
pde_flux.f90       elemental flux functions (no I/O)
      │
pde_module.f90     solver types, step logic, NetCDF writer
      │
pde_driver.f90     CLI → runs loop → density_history → write_pde_netcdf
      │
pde_simulation.nc  density(time, x), flow(time), global attrs
      │
pde_runner.py      load_pde_netcdf → dict with (n_steps+1, M) density array
```

---

## 7. Conventions and gotchas

- **NetCDF dimension ordering.** The Fortran variable `density_history(M, n_steps+1)`
  is defined with `[x_dimid, time_dimid]` (Fortran: x fastest). The C-convention
  NetCDF library reverses this so the file stores `(time, x)`. Python's `netCDF4`
  reads it directly as shape `(n_steps+1, M)`. **Do not transpose** — unlike the
  TASEP `io.py`, `load_pde_netcdf` returns density already in `(time, x)` order.

- **TASEP vs PDE array shape.** The TASEP loader (`io.py`) does apply `.T` to
  `history`, yielding shape `(L, n_steps)` = `(site, time)`. The PDE loader does
  not. This is intentional: TASEP history is used for space-time plots where
  site-indexed layout is convenient; PDE density is used for time-series analysis
  where time-indexed layout is natural.

- **CFL and the sonic point.** The Lax–Friedrichs scheme becomes an unstable
  centred scheme when Δt → ∞. With Greenshields, dq/dρ → 0 as ρ → ρ_c, so
  a naïve adaptive CFL gives Δt → ∞ for near-uniform fields. `compute_dt`
  caps Δt at `cfl · Δx / v_max` to prevent this.

- **Godunov vs Lax–Friedrichs.** Godunov is the correct scheme for the
  validation suite; LF is the default for early development because it is
  simpler and never requires case analysis. Switch with `--flux godunov` on
  the command line or `flux_type="godunov"` in `run_pde`. Expect sharper
  shocks from Godunov and more smeared ones from LF.

- **Open BC and the Riemann IC.** For `ic_type="riemann"`, the left and right
  initial densities reuse `rho_left_bc` and `rho_right_bc`. This means the
  open boundary ghost cells are already consistent with the IC, so there is no
  spurious transient at t=0 from a mismatch between the IC and the BC.

- **Flow variable.** `flow(t)` in the NetCDF is q(ρ_M) — the flux computed at
  the rightmost physical cell. Under open BCs this approximates the outflow
  rate. It is analogous to the TASEP `current` variable (exits per step) but
  is a continuous quantity in [0, v_max·ρ_max/4].

---

## 8. Setup

```bash
# Build the PDE solver
make pde                            # produces build/pde_solver

# Run with defaults (Riemann IC, Lax-Friedrichs, open BCs)
make run-pde

# Run with Godunov flux and a Gaussian IC
make run-pde PDE_IC=gaussian PDE_FLUX=godunov PDE_BC=periodic PDE_STEPS=1000

# Python tests
source .venv/bin/activate
python -m pytest tests/test_pde.py -v
```

Fortran prerequisites: `gfortran`, NetCDF Fortran library (`nf-config`,
`nc-config`). Python prerequisites: `numpy`, `netCDF4`, `pytest` (see
`requirements.txt`).
