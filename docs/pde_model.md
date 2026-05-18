# PDE Model — Developer Documentation

This document covers the LWR (Lighthill–Whitham–Richards) PDE implementation:
the physics, the Fortran solver core, the NetCDF I/O layer, the Python
interface, and the test suite. It tracks **Phases 1 and 2** of `PDE.md`.

---

## 1. Physics: LWR scalar conservation law

The model is the **Lighthill–Whitham–Richards** traffic flow equation — a scalar
hyperbolic conservation law:

$$
\frac{\partial \rho}{\partial t} + \frac{\partial q(\rho)}{\partial x} = 0
$$

where $\rho(x, t)$ is vehicle density and $q(\rho)$ is the traffic flow (vehicles per
unit time). The **Greenshields fundamental diagram** is used as the closure:

$$
v(\rho) = v_{\max}\left(1 - \frac{\rho}{\rho_{\max}}\right)
$$

Speed decreases linearly with density.

$$
q(\rho) = \rho v(\rho)
= v_{\max}\rho\left(1 - \frac{\rho}{\rho_{\max}}\right)
$$

This gives a concave parabolic flow density relation.

$$
\frac{dq}{d\rho}= v_{\max}\left(1 - \frac{2\rho}{\rho_{\max}}\right)
$$

This is the characteristic speed, which becomes zero at 

$$
\rho_c = \frac{\rho_{\max}}{2}
$$

### Critical density

The flow $q$ is concave and achieves its maximum at:

> $ρ_c = ρ_{\max} / 2$ 

This divides the domain into two phases:

| Phase | Condition | dq/dρ | Behaviour |
|-------|-----------|-------|-----------|
| Free flow (subcritical) | $ρ < ρ_c$ | > 0 | density waves travel forward |
| Congested (supercritical) | $ρ > ρ_c$ | < 0 | density waves travel backward |

### Wave structures

The finite-volume solver evolves the density by computing fluxes across cell interfaces. At each interface, the numerical method locally sees a jump between a left density and a right density. This local jump problem is called a **Riemann problem**.

#### The Riemann problem

A Riemann problem is an initial-value problem where the density is piecewise constant:

$$
\rho(x,0) =
\begin{cases}
\rho_L, & x < 0, \\
\rho_R, & x > 0.
\end{cases}
$$

Here, $\rho_L$ is the density immediately to the left of an interface and \(\rho_R\) is the density immediately to the right. Even if the full simulation contains a smooth or complicated density profile, the finite-volume method approximates the solution as piecewise constant inside each cell. Therefore, every cell boundary can be treated as a small Riemann problem.

The solution of the Riemann problem determines how information travels away from the interface. For the LWR traffic equation,

$$
\frac{\partial \rho}{\partial t} + \frac{\partial q(\rho)}{\partial x} = 0,
$$

the relevant wave speed is the characteristic speed

$$
\frac{dq}{d\rho}.
$$

For the Greenshields flux,

$$
q(\rho) = v_{\max}\rho\left(1-\frac{\rho}{\rho_{\max}}\right),
$$

the characteristic speed is

$$
\frac{dq}{d\rho}=v_{\max}\left(1-\frac{2\rho}{\rho_{\max}}\right).
$$

This speed is positive in the free-flow regime, zero at the critical density, and negative in the congested regime. As a result, density disturbances can move either downstream or upstream depending on the local density.

#### Shock waves

A **shock** forms when characteristics move into each other and the density profile steepens into a discontinuity. For the Greenshields flux, which is concave, this occurs when

$$
\rho_L < \rho_R.
$$

This corresponds to a lower-density region feeding into a higher-density region. In traffic terms, vehicles are moving from a freer region into a more congested region, so the transition compresses into a sharp front.

The shock travels at the Rankine--Hugoniot speed

$$
s =\frac{q(\rho_R)-q(\rho_L)}{\rho_R-\rho_L}.
$$

This expression comes from conservation: the speed of the discontinuity is determined by the difference in flow across the jump divided by the difference in density. Depending on the values of $\rho_L$ and $\rho_R$, the shock may move downstream or upstream.

A particularly important case is when the shock connects a free-flow state to a congested state:

$$
\rho_L < \rho_c < \rho_R.
$$

Here the interface crosses the critical density. In the Godunov method, the flux is limited by the smaller of the two possible flows,

$$
F_G = \min(q(\rho_L), q(\rho_R)).
$$

This reflects the idea that the amount of traffic passing through the interface is constrained by the bottleneck between upstream demand and downstream supply.

#### Rarefaction waves

A **rarefaction** forms when characteristics spread apart. For the Greenshields flux, this occurs when

$$
\rho_L > \rho_R.
$$

This corresponds to a higher-density region opening into a lower-density region. In traffic terms, congestion is releasing into a freer part of the road, so the transition spreads out rather than forming a sharp jump.

Instead of a single discontinuity, the solution becomes a fan of characteristics connecting $\rho_L$ to $\rho_R$. Each density within the fan travels at its own characteristic speed

$$
\frac{dq}{d\rho}.
$$


The most important rarefaction case is when the left state is congested and the right state is free-flowing:

$$
\rho_R < \rho_c < \rho_L.
$$

This rarefaction crosses the critical density, where

$$
\frac{dq}{d\rho}=0.
$$

Because one characteristic has zero speed, this is sometimes called a **sonic rarefaction**. In this case, the Godunov flux takes the maximum possible value of the fundamental diagram:

$$
F_G = q(\rho_c)
= \frac{v_{\max}\rho_{\max}}{4}.
$$

Physically, this means the interface is operating at road capacity.

#### Summary of wave types

| $\rho_L$ vs $\rho_R$ | Wave type | Propagation |
|------------|-----------|-------------|
| $\rho_L < \rho_R$ | **Shock** | speed $s = \dfrac{q(\rho_R)-q(\rho_L)}{\rho_R-\rho_L}$ using the Rankine--Hugoniot condition |
| $\rho_L > \rho_R$ | **Rarefaction** | fan of characteristics connecting $\rho_L$ to $\rho_R$ |
| $\rho_L < \rho_c < \rho_R$ | Shock through critical | Godunov flux $= \min(q(\rho_L), q(\rho_R))$ |
| $\rho_R < \rho_c < \rho_L$ | Sonic rarefaction | Godunov flux $= q(\rho_c) = \dfrac{v_{\max}\rho_{\max}}{4}$ |

## 2. Repository layout (PDE files)

```
src/fortran/
├── pde_flux.f90          # Greenshields flux functions + numerical flux schemes
├── pde_module.f90        # solver types, initialise/step/finalise, NetCDF output
└── pde_driver.f90        # command-line entry point; runs simulation + writes .nc

src/python/
├── pde_runner.py         # Python flux mirrors, run_pde(), load_pde_netcdf()
└── visualisation.py      # (shared) — PDE plots added alongside existing TASEP plots

scripts/
└── run_pde_model.py      # build → run → load → plot summary figure

tests/
├── test_pde.py           # Phase 1 flux unit tests + Phase 4 solver tests (complete)
└── exact_riemann.py      # exact LWR Riemann solution helper (consumed by Phase 4 tests)

data/output/
└── pde_simulation.nc     # default output path

plots/
└── pde_summary.png       # written by run_pde_model.py --save
```

---

## 3. Fortran core

### `pde_flux.f90` — module `pde_flux`

All flux routines are **elemental** (work on scalars and arrays transparently).

| Routine | Signature | Purpose |
|---------|-----------|---------|
| `v_of_rho` | `(rho, v_max, rho_max)` | Greenshields speed $v(\rho)$ |
| `q_of_rho` | `(rho, v_max, rho_max)` | Flow $q(ρ) = \rho v(ρ)$ |
| `dq_drho`  | `(rho, v_max, rho_max)` | Characteristic speed $dq/d\rho$|
| `rho_critical` | `(rho_max)` | $\rho_c = \rho_\max/2$ |

Two non-elemental **numerical flux** functions:

#### `lax_friedrichs_flux(rho_L, rho_R, v_max, rho_max, dx, dt)`

$$
F_{LF} = \frac{q(ρ_L) + q(ρ_R)}{2}  −  \frac{\Delta x}{2 \Delta t} (\rho_R − \rho_L)
$$

Stable for CFL ≤ 1 but adds numerical diffusion proportional to $\Delta x/(2\Delta t)$.
Used as the default scheme in Phases 1–2 because it is simple and robust.

#### `godunov_flux(rho_L, rho_R, v_max, rho_max)`

Exact Godunov flux for a concave flux. Closed-form case table
(LeVeque, *Finite Volume Methods for Hyperbolic Problems*, §12.1):

$$
F_G(\rho_L,\rho_R)=\begin{cases}
q(\rho_L),
& \rho_L \leq \rho_R \leq \rho_c
\quad \text{shock, both subcritical},
\\[4pt]
q(\rho_R),
& \rho_c \leq \rho_L \leq \rho_R
\quad \text{shock, both supercritical},
\\[4pt]
\min\!\left(q(\rho_L),q(\rho_R)\right),
& \rho_L < \rho_c < \rho_R
\quad \text{shock through critical},
\\[8pt]
q(\rho_L),
& \rho_R < \rho_L \leq \rho_c
\quad \text{rarefaction, both subcritical},
\\[4pt]
q(\rho_R),
& \rho_c \leq \rho_R < \rho_L
\quad \text{rarefaction, both supercritical},
\\[4pt]
q(\rho_c),
& \rho_R < \rho_c < \rho_L
\quad \text{sonic rarefaction}.
\end{cases}
$$

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
| `"riemann"` | `rho_left_bc` for $x < L/2$, `rho_right_bc` for $x ≥ L/2$ |
| `"gaussian"` | Base $0.2\rho_\max$ + Gaussian bump of amplitude $0.6\rho_\max$, $\sigma = 0.05L$ |
| `"sine"` | Mean $0.5\rho_\max$ + sinusoidal perturbation of amplitude $0.15\rho_\max$ |

For the `"riemann"` IC, `rho_left_bc` and `rho_right_bc` double as the
left/right initial densities. This intentionally links the IC to the open
boundary values so that the inflow/outflow match the initial state.

#### `pde_step(state, params)`

One conservative finite-volume update:

$$
\rho_i^{n+1} = \rho_i^n − \frac{\Delta t}{\Delta x} [F_{i+1/2} − F_{i-1/2}]
$$

where $F_{i+1/2}$ is either `lax_friedrichs_flux` or `godunov_flux` depending
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

### `src/python/visualisation.py` — PDE plotting functions

Four functions added to the shared visualisation module (alongside the
existing TASEP functions). All accept the dict returned by `load_pde_netcdf`.

#### `plot_pde_spacetime(data, ax=None, title=None)`

`viridis` heatmap of ρ(x, t) with x on the horizontal axis and time
increasing upward. This is the primary visual diagnostic — shocks appear
as straight diagonal lines whose slope gives the wave speed, and rarefaction
fans appear as triangular colour gradients.

#### `plot_pde_snapshots(data, n_snapshots=6, ax=None, title=None)`

Line plots of ρ(x) at `n_snapshots` evenly-spaced times, coloured from dark
to light using the `plasma` colourmap. Useful for reading off shock position
at each time and checking the profile shape against the analytical solution.

#### `plot_pde_flow(data, ax=None, title=None)`

Time series of q(ρ_M) — the flow at the rightmost physical cell — with a
dashed reference line at q_max = v_max·ρ_max/4. Under open BCs this tracks
the outflow rate; it should remain constant while the shock is away from the
boundary and then change when the shock exits.

#### `plot_pde_summary(data, save_path=None)`

Three-panel composite figure: space-time heatmap (top, full width), density
snapshots (bottom-left), and boundary flow (bottom-right). The suptitle
includes M, n_steps, IC, flux, and BC type from the NetCDF attributes.
Pass a path to `save_path` to write a PNG at 150 dpi.

---

## 5. Entry-point script — `scripts/run_pde_model.py`

Wraps the full pipeline: builds the Fortran binary → runs the solver →
loads the NetCDF → displays the summary figure. Mirrors `run_toy_model.py`
in style.

```bash
python scripts/run_pde_model.py                        # build, run, plot
python scripts/run_pde_model.py --no-run               # replot existing .nc
python scripts/run_pde_model.py --save                 # write plots/pde_summary.png
```

**Flags:**

| Flag | Default | Notes |
|------|---------|-------|
| `--M` | 200 | spatial cells |
| `--steps` | 500 | time steps |
| `--ic` | `riemann` | `constant` / `riemann` / `gaussian` / `sine` |
| `--flux` | `lf` | `lf` (Lax–Friedrichs) / `godunov` |
| `--bc` | `open` | `open` / `periodic` |
| `--rho-left` | 0.1 | left BC and left Riemann state |
| `--rho-right` | 0.9 | right BC and right Riemann state |
| `--output` | `data/output/pde_simulation.nc` | NetCDF path |
| `--no-run` | — | skip build and solver; load existing file |
| `--save` | — | save figure to `plots/pde_summary.png` |

**Example commands for different wave structures:**

```bash
# Stationary shock (default — symmetric states, s=0)
python scripts/run_pde_model.py

# Rightward-moving shock  s = (q(0.7)−q(0.1))/0.6 = +0.2
python scripts/run_pde_model.py --rho-left 0.1 --rho-right 0.7

# Leftward-moving shock   s = (q(0.8)−q(0.3))/0.5 = −0.1
python scripts/run_pde_model.py --rho-left 0.3 --rho-right 0.8

# Sonic rarefaction (dense left, sparse right, passes through ρ_c)
python scripts/run_pde_model.py --rho-left 0.8 --rho-right 0.1

# Sharper shocks with Godunov flux
python scripts/run_pde_model.py --rho-left 0.1 --rho-right 0.7 --flux godunov --save

# Gaussian bump on periodic domain
python scripts/run_pde_model.py --ic gaussian --bc periodic --M 400 --steps 1000
```

---

## 6. Tests — `tests/test_pde.py`

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
The exact Riemann solution used by tests (c)–(e) lives in `tests/exact_riemann.py`.

| Test | What it checks | Tolerance |
|------|----------------|-----------|
| `test_constant_solution_preserved` | Uniform IC stays uniform for 1000 steps with periodic BCs | 1 × 10⁻¹⁰ |
| `test_mass_conservation_periodic` | ∑ ρ_i · Δx constant for Gaussian IC with periodic BCs over 1000 steps | 1 × 10⁻⁴ (float32 limit) |
| `test_shock_speed` | Shock position (ρ_L=0.3, ρ_R=0.8) matches Rankine–Hugoniot speed s=−0.1 | 4 cell widths (1 per 100 steps) |
| `test_rarefaction_l1_error` | Sonic rarefaction (ρ_L=0.8, ρ_R=0.2) L1 error vs exact fan | < 10 · Δx |
| `test_convergence_first_order` | L1 error ratio ≥ 1.5 as M doubles (50→100→200); n_steps=M keeps t_final≈0.9 | — |

The mass conservation tolerance is 1×10⁻⁴ rather than the 1×10⁻¹⁰ quoted in
`PDE.md`. This is because the Fortran solver and NetCDF output both use
single precision (float32, ~7 significant digits). LF with periodic BCs is
exactly mass-conservative in exact arithmetic (the boundary flux telescopes),
so any deviation comes from float32 rounding in the stored NetCDF values.
If double precision is needed, change `NF90_FLOAT` to `NF90_DOUBLE` in
`write_pde_netcdf` and use `real(8)` throughout the Fortran.

#### `tests/exact_riemann.py` — `exact_riemann_lwr(x, t, rho_L, rho_R, ...)`

Returns the exact density field for an LWR Riemann problem with Greenshields
flux, evaluated at positions `x` and time `t`:

- **ρ_L < ρ_R** (shock): discontinuity at `x0 + s·t` where `s = (q(ρ_R)−q(ρ_L))/(ρ_R−ρ_L)`.
- **ρ_L > ρ_R** (rarefaction): centred fan; in the fan `ρ = (ρ_max/2)·(1 − (x−x0)/(v_max·t))`.

The function uses the self-similar variable ξ = (x−x0)/t, so it correctly
handles sonic rarefactions (fan crossing ρ_c) without any special-casing.

```bash
.venv/bin/python -m pytest tests/test_pde.py -v
```

---

## 7. Pipeline overview

```
# Full pipeline via script
python scripts/run_pde_model.py --ic riemann --flux godunov --save

# Or step by step
make run-pde                  # Fortran solver → data/output/pde_simulation.nc
python -c "
import sys; sys.path.insert(0,'src/python')
from pde_runner import load_pde_netcdf
from visualisation import plot_pde_summary
import matplotlib.pyplot as plt
data = load_pde_netcdf('data/output/pde_simulation.nc')
plot_pde_summary(data, save_path='plots/pde_summary.png')
plt.show()
"
```

```
pde_flux.f90         elemental flux functions (no I/O)
       │
pde_module.f90       solver types, step logic, NetCDF writer
       │
pde_driver.f90       CLI → runs loop → density_history → write_pde_netcdf
       │
pde_simulation.nc    density(time, x), flow(time), global attrs
       │
pde_runner.py        load_pde_netcdf → dict with (n_steps+1, M) density array
       │
visualisation.py     plot_pde_summary → 3-panel figure
       │
run_pde_model.py     end-to-end CLI: build → run → load → show/save
```

---

## 8. Conventions and gotchas

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

## 9. Setup

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
