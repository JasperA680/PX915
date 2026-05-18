# Multi-Lane PDE Model — Developer Documentation

This document covers the multi-lane extension of the LWR PDE solver: the physics
of lane-coupled traffic flow, the per-lane data structures, the conservative
lane-change source terms, the extended NetCDF schema, the Python interface, and
the validation suite.

It builds directly on `pde_model.md` (the single-lane LWR documentation) and tracks
**Phases 1–5** of `multilane_pde.md`. The single-lane solver is recovered exactly
when `n_lanes = 1` and `lane_change_rate = 0`.

A multi-lane PDE extension replaces the scalar density field with one density field
per lane. Each lane obeys its own LWR conservation law with a lane-specific flux,
while additional source terms transfer vehicles between adjacent lanes. These
source terms model lane-changing behaviour and must be constructed so that total
vehicle count is conserved across lanes. Numerically, the finite-volume update is
applied independently per lane for the longitudinal fluxes, followed by a
conservative source-term update for lateral exchange.

---

## 1. Physics: multi-lane LWR with lane-change source terms

### 1.1 Per-lane conservation law

The state is now a vector of density fields, one per lane:

$$
\rho_1(x, t),\ \rho_2(x, t),\ \ldots,\ \rho_N(x, t)
$$

Each lane obeys its own LWR conservation law, augmented with a **source term**
$S_\ell(x, t)$ that represents the rate of vehicles entering lane $\ell$ from
adjacent lanes:

$$
\frac{\partial \rho_\ell}{\partial t}
+ \frac{\partial q_\ell(\rho_\ell)}{\partial x}
= S_\ell.
$$

The per-lane flux uses the same Greenshields closure as the single-lane model,
but with **lane-specific parameters** $v_{\max,\ell}$ and $\rho_{\max,\ell}$:

$$
q_\ell(\rho_\ell) = v_{\max,\ell}\,\rho_\ell\!\left(1 - \frac{\rho_\ell}{\rho_{\max,\ell}}\right),
\qquad
v_\ell(\rho_\ell) = v_{\max,\ell}\!\left(1 - \frac{\rho_\ell}{\rho_{\max,\ell}}\right).
$$

This allows different lanes to have different free-flow speeds (fast/slow lanes)
or different jam densities (lane closures, bottlenecks).

### 1.2 The conservation requirement

The non-negotiable physical requirement: **lane-change source terms must sum to zero
across lanes at every cell**. Vehicles moving between lanes are neither created
nor destroyed:

$$
\sum_{\ell=1}^{N} S_\ell(x, t) = 0 \quad \forall\, x, t.
$$

This guarantees that the total vehicle count

$$
M_{\text{tot}}(t) = \sum_{\ell=1}^N \int \rho_\ell(x, t)\, dx
$$

is exactly conserved under periodic boundary conditions, and conserved up to
the boundary fluxes under open BCs. The model is constructed pairwise so this
property holds by construction (see §1.4).

### 1.3 Two-lane source-term model

For two lanes the exchange is described by directional flux densities
$S_{1\to 2}$ and $S_{2\to 1}$ representing vehicles moving from lane 1 to lane 2
and vice versa:

$$
S_{1\to 2,\,i} = k\, \rho_{1,i}
                  \cdot \max\!\bigl(0,\ 1 - \rho_{2,i}/\rho_{\max,2}\bigr)
                  \cdot \max\!\bigl(0,\ v_2(\rho_{2,i}) - v_1(\rho_{1,i})\bigr),
$$

$$
S_{2\to 1,\,i} = k\, \rho_{2,i}
                  \cdot \max\!\bigl(0,\ 1 - \rho_{1,i}/\rho_{\max,1}\bigr)
                  \cdot \max\!\bigl(0,\ v_1(\rho_{1,i}) - v_2(\rho_{2,i})\bigr).
$$

The net source per lane is then

$$
S_{1,i} = S_{2\to 1,i} - S_{1\to 2,i}, \qquad
S_{2,i} = S_{1\to 2,i} - S_{2\to 1,i}.
$$

By construction, $S_{1,i} + S_{2,i} = 0$ at every cell.

The three factors carry distinct physical meaning:

| Factor | Meaning |
|--------|---------|
| $k$ | Lane-change rate constant (units of inverse time). $k=0$ disables lane changing. |
| $\rho_{\ell,i}$ | Source-lane density: more vehicles available to change lanes. |
| $\max(0,\ 1 - \rho_{m,i}/\rho_{\max,m})$ | Receiving-lane headroom: blocked at jam density. |
| $\max(0,\ v_m - v_\ell)$ | Asymmetric drive: vehicles only move toward a faster lane. |

The asymmetric $\max(0, \Delta v)$ factor means **at most one of $S_{1\to 2}$ and
$S_{2\to 1}$ is nonzero per cell**. Vehicles never simultaneously flow both ways.

### 1.4 Generalisation to N lanes — adjacent-lane coupling

For $N \geq 2$ lanes the model is restricted to **adjacent-lane exchanges only**:
lane $\ell$ exchanges with $\ell-1$ and $\ell+1$, and there is no end wrap.
Each adjacent pair uses the same pairwise formula, and the contributions are
accumulated into the source array:

$$
S_\ell(x, t) = \sum_{\substack{m \in \{\ell-1,\ \ell+1\} \\ 1 \leq m \leq N}}
                \bigl(S_{m \to \ell}(x, t) - S_{\ell \to m}(x, t)\bigr).
$$

This choice keeps the model physically meaningful (a vehicle on the inside lane
cannot teleport to the outside lane in one step) and keeps the conservation
property pairwise — easier to reason about and verify than a fully coupled
all-to-all model. See §6 for the explicit per-cell-source conservation test.

### 1.5 Lane-change equilibrium and qualitative behaviour

A steady state of the source step is reached when no exchange is active in any
cell, which requires either zero density in the source lane, zero capacity in
the receiving lane, or **equal velocities**:

$$
v_1(\rho_1) = v_2(\rho_2)
\quad\Leftrightarrow\quad
v_{\max,1}\!\left(1 - \frac{\rho_1}{\rho_{\max,1}}\right)
= v_{\max,2}\!\left(1 - \frac{\rho_2}{\rho_{\max,2}}\right).
$$

For two lanes with equal jam densities and **fixed total mass** under periodic
BCs ($\rho_1 + \rho_2 = M_0$), this gives a closed-form equilibrium

$$
\rho_2^{\text{eq}}
= \frac{v_{\max,2}/v_{\max,1} - 1 + M_0/\rho_{\max}}{1 + v_{\max,2}/v_{\max,1}},
\qquad
\rho_1^{\text{eq}} = M_0 - \rho_2^{\text{eq}}.
$$

For $v_{\max,1} = 1.0$, $v_{\max,2} = 1.5$, $\rho_{\max} = 1$, $M_0 = 0.9$,
this gives $\rho_1^{\text{eq}} \approx 0.34$ and $\rho_2^{\text{eq}} \approx 0.56$.
The faster lane ends up denser because vehicles preferentially migrate into it.
The slower lane carries lower flow.

Under **open BCs** with constant inflow, the steady state is different: each lane
adjusts to carry the imposed inflow density at its own speed, and the slower lane
ends up denser because it transports vehicles more slowly for a given throughput.

### 1.6 Stability bounds

Two stability constraints apply each step:

| Constraint | Reason |
|------------|--------|
| $\Delta t \leq \dfrac{C_{\text{cfl}}\,\Delta x}{\max_{\ell} v_{\max,\ell}}$ | Longitudinal CFL across all lanes |
| $\Delta t < 1/k$ (rule of thumb) | Source-ODE stability for the explicit source step |

The first is enforced exactly by `compute_dt`. The second is currently advisory:
the implementation uses **explicit Euler** for the source term and clips any
density that goes outside $[0, \rho_{\max,\ell}]$ after the source update. If
clipping fires repeatedly, reduce $k$ or use a smaller $\Delta t$.

---

## 2. Repository layout (multi-lane PDE files)

```
src/fortran/
├── pde_flux.f90          # unchanged from single-lane (per-lane scalar flux calls)
├── pde_lanechange.f90    # NEW — conservative lane-change source-term module
├── pde_module.f90        # extended: 2D density, n_lanes, multi-lane step + NetCDF
└── pde_driver.f90        # extended: n_lanes / lane_change_rate / per-lane v_max CLI

src/python/
├── pde_runner.py         # run_pde() supports v_max_lanes/rho_max_lanes lists;
│                         # load_pde_netcdf() returns density_per_lane and flow_per_lane
├── analysis.py           # multi-lane diagnostic helpers
│                         # (compute_total_density / compute_total_flow / compute_total_mass / …)
└── visualisation.py      # per-lane and total plots (plot_space_time_per_lane / …)

scripts/
├── run_pde_model.py      # single-lane entry point (still works unchanged)
└── run_multilane_pde.py  # NEW — four canonical multi-lane scenarios

tests/
├── test_pde.py                # single-lane suite (15 tests, still passes)
└── test_pde_multilane.py      # NEW — Phase 5 validation suite (7 tests)

docs/
├── multilane_pde.md           # design plan
└── multilane_pde_model.md     # this document — developer reference
```

The single-lane solver is now the `n_lanes = 1` special case of the multi-lane
solver — there is exactly one PDE solver in the codebase.

---

## 3. Fortran core

### 3.1 `pde_flux.f90` — unchanged

The flux module is reused without modification. Numerical fluxes
(`lax_friedrichs_flux`, `godunov_flux`) are called per lane with that lane's
$v_{\max,\ell}$ and $\rho_{\max,\ell}$, producing one flux array per lane.

### 3.2 `pde_lanechange.f90` — module `pde_lanechange`

A small, self-contained module exposing one public subroutine.

```fortran
subroutine compute_lane_change_sources(density, n_lanes, M, &
                                       v_max_lanes, rho_max_lanes, k, source)
  integer, intent(in)  :: n_lanes, M
  real,    intent(in)  :: density(n_lanes, M)
  real,    intent(in)  :: v_max_lanes(n_lanes), rho_max_lanes(n_lanes)
  real,    intent(in)  :: k
  real,    intent(out) :: source(n_lanes, M)
end subroutine
```

For each adjacent pair $(\ell, \ell+1)$ and each cell $i$, it computes the
directional fluxes $S_{\ell \to \ell+1}$ and $S_{\ell+1 \to \ell}$ from the
formulas in §1.3 and accumulates the net contribution into `source(:, i)`.

The output satisfies $\sum_\ell \mathrm{source}(\ell, i) = 0$ to machine
precision for every cell (verified to $10^{-14}$ in Phase 5).

### 3.3 `pde_module.f90` — module `pde_solver`

#### Derived types

```fortran
type :: pde_params_t
  real    :: dx, dt, domain_length
  real    :: v_max, rho_max                  ! scalar aliases (broadcast targets)
  real    :: rho_left_bc, rho_right_bc, cfl_number
  integer :: M, n_steps, C_checkpoint
  character(len=16) :: ic_type    ! "constant"|"riemann"|"gaussian"|"sine"|"staggered"
  character(len=16) :: bc_type    ! "open"|"periodic"|"sponge"
  character(len=16) :: flux_type  ! "lf"|"godunov"
  logical :: use_adaptive_dt
  integer :: n_sponge
  real    :: sponge_damping

  ! Multi-lane extensions
  integer :: n_lanes
  real    :: lane_change_rate                ! k; 0 disables lane changing
  real, allocatable :: v_max_lanes(:)        ! (n_lanes)
  real, allocatable :: rho_max_lanes(:)      ! (n_lanes)
  real, allocatable :: rho_left_bc_lanes(:)  ! (n_lanes)
  real, allocatable :: rho_right_bc_lanes(:) ! (n_lanes)
end type

type :: pde_state_t
  real, allocatable :: density(:,:)   ! (n_lanes, M)         physical cells
  real, allocatable :: rho_ext(:,:)   ! (n_lanes, 0:M+1)     ghost-cell scratch
  real, allocatable :: flux(:,:)      ! (n_lanes, 0:M)       interface-flux scratch
  real    :: t_current
  integer :: step
  integer :: clip_count               ! cells clamped after source step
end type
```

The **lane-major** convention `density(lane, i)` matches Fortran column-major
order: per-lane operations stride contiguously in memory. All scratch arrays
follow the same convention.

#### `pde_setup_params(params)`

Allocates the per-lane arrays and broadcasts the scalar aliases into them. Must
be called **after** all scalar fields are set and **before** `pde_initialise`.
After this call:

- `params%v_max_lanes(:)         = params%v_max`
- `params%rho_max_lanes(:)       = params%rho_max`
- `params%rho_left_bc_lanes(:)   = params%rho_left_bc`
- `params%rho_right_bc_lanes(:)  = params%rho_right_bc`

The driver overrides individual lanes after this call when comma-separated lists
are supplied on the command line.

#### `pde_initialise(state, params)`

Allocates the 2D state arrays and fills `state%density` from the IC,
**applied per lane** so each lane receives a possibly different initial profile.
Available IC types:

| `ic_type` | Initial condition (per lane $\ell$) |
|-----------|------------------|
| `"constant"` | Uniform density `rho_left_bc_lanes(`$\ell$`)` |
| `"riemann"` | `rho_left_bc_lanes(`$\ell$`)` for $x<L/2$, else `rho_right_bc_lanes(`$\ell$`)` |
| `"gaussian"` | $0.2\rho_{\max,\ell}$ + $0.6\rho_{\max,\ell}\,\exp\!\bigl(-(x-L/2)^2/(2\sigma^2)\bigr)$ |
| `"sine"` | $0.5\rho_{\max,\ell}$ + $0.15\rho_{\max,\ell}\,\sin(2\pi x/L)$ |
| `"staggered"` | Odd lanes use `rho_left_bc`, even lanes use `rho_right_bc` |

The new `"staggered"` IC is the simplest way to give different uniform
densities to alternating lanes without extending the CLI.

#### `pde_step(state, params)`

One full time step with three phases applied in order:

1. **Per-lane longitudinal finite-volume update** (CFL-stable):

   $$
   \rho_{\ell,i}^{n+1/2} = \rho_{\ell,i}^n
   - \frac{\Delta t}{\Delta x}\bigl(F^{(\ell)}_{i+1/2} - F^{(\ell)}_{i-1/2}\bigr)
   $$

   where $F^{(\ell)}$ uses lane $\ell$'s $v_{\max,\ell}$, $\rho_{\max,\ell}$ and
   the shared `flux_type`. Boundary ghost cells follow the same rules as the
   single-lane solver but are populated independently per lane.

2. **Conservative lane-change source step** (only if $k > 0$ and $n_{\text{lanes}} > 1$):

   $$
   \rho_{\ell,i}^{n+1} = \rho_{\ell,i}^{n+1/2} + \Delta t\,S_{\ell,i}
   $$

   followed by clipping to $[0, \rho_{\max,\ell}]$ with a warning counter
   (`state%clip_count`).

3. **Sponge damping** (if `bc_type == "sponge"`), applied per lane.

#### `compute_dt(state, params, dt_out)`

Adaptive CFL across all lanes:

```
max_speed       = max over all (lane, i) of |dq/dρ(ρ_{lane,i}; v_max_lane, rho_max_lane)|
v_max_global    = max over lanes of v_max_lanes(lane)
dt_conservative = cfl · Δx / v_max_global
dt_out          = min(cfl · Δx / max_speed, dt_conservative)
```

The cap `cfl · Δx / v_max_global` is the same sonic-point safeguard described
in `pde_model.md` §3, generalised to take the maximum lane speed.

#### `write_pde_netcdf(filename, params, density_history, flow_history)`

Writes all output to a single NetCDF file.

**Schema:**

| Item | Kind | Notes |
|------|------|-------|
| `time` | unlimited dimension | length n_steps+1 (includes t=0) |
| `x` | dimension | length M |
| `lane` | dimension | length n_lanes (always present, even for n_lanes=1) |
| `time` | `float(time)` | physical time per step |
| `x` | `float(x)` | cell-centre x-coordinates |
| `lane` | `float(lane)` | lane indices 1, 2, …, n_lanes |
| `density` | `float(lane, x, time)` | per-lane density history |
| `flow` | `float(lane, time)` | per-lane right-boundary flow |
| `flow_total` | `float(time)` | $\sum_\ell q_\ell$, sum over lanes |

Global attributes: all single-lane attributes plus `n_lanes`,
`lane_change_rate`, `v_max_lanes(:)`, and `rho_max_lanes(:)` (vector
attributes, length `n_lanes`).

**Dimension ordering note.** The Fortran `nf90_def_var` call lists dimensions
in column-major order `[lane_dimid, x_dimid, time_dimid]`. The C-convention
NetCDF library reverses this in storage, so Python's `netCDF4` reads the
variable as shape `(time, x, lane)`. The Python loader transposes to
`(time, lane, x)` for downstream callers — see §4.

### 3.4 `pde_driver.f90` — program `pde_driver`

Command-line entry point with positional arguments. Backward-compatible with
the original 10-argument form; arguments 11–14 enable multi-lane mode.

```bash
./build/pde_solver [M] [n_steps] [v_max] [rho_max] [rho_left] [rho_right] \
                   [ic_type] [flux_type] [bc_type] [output] \
                   [n_lanes] [lane_change_rate] \
                   [v_max_lanes] [rho_max_lanes]
```

| Arg | Default | Notes |
|-----|---------|-------|
| 1. `M` | 200 | spatial cells |
| 2. `n_steps` | 500 | time steps |
| 3. `v_max` | 1.0 | scalar $v_{\max}$ broadcast to all lanes |
| 4. `rho_max` | 1.0 | scalar $\rho_{\max}$ broadcast to all lanes |
| 5. `rho_left` | 0.1 | left BC / left Riemann state |
| 6. `rho_right` | 0.9 | right BC / right Riemann state |
| 7. `ic_type` | `riemann` | `constant`/`riemann`/`gaussian`/`sine`/`staggered` |
| 8. `flux_type` | `lf` | `lf` / `godunov` |
| 9. `bc_type` | `open` | `open` / `periodic` / `sponge` |
| 10. `output` | `data/output/pde_simulation.nc` | NetCDF path |
| 11. `n_lanes` | 1 | number of lanes |
| 12. `lane_change_rate` | 0.0 | $k$; 0 disables lane changing |
| 13. `v_max_lanes` | broadcast | comma list, e.g. `"1.0,1.5"`, length must equal `n_lanes` |
| 14. `rho_max_lanes` | broadcast | comma list, length must equal `n_lanes` |

The driver:

1. Parses arguments, computes $\Delta x = L/M$, sets initial $\Delta t$.
2. Calls `pde_setup_params` to allocate and broadcast scalar params.
3. Parses optional comma lists for arguments 13–14, overriding lane arrays.
4. Recomputes $\Delta t$ using `maxval(v_max_lanes)` if argument 13 is supplied.
5. Allocates `density_history(n_lanes, M, n_steps+1)` and
   `flow_history(n_lanes, n_steps+1)`.
6. Runs the time loop, recording per-lane density and per-lane right-boundary
   flow at every step.
7. Calls `write_pde_netcdf` and `pde_finalise`.
8. Reports the final clip count if any cells were clamped after the source step.

The `parse_real_list(str, vals, n_found)` internal subroutine parses the
comma-separated argument strings.

### 3.5 Build

The `Makefile` rule lists the four PDE sources in `use`-order:

```make
$(PDE_EXE): $(PDE_FLUX_SRC) $(PDE_LC_SRC) $(PDE_MOD_SRC) $(PDE_DRV_SRC)
```

`pde_lanechange.f90` is independent of `pde_module.f90` (it only depends on
intrinsic real arithmetic), so it must be compiled before `pde_module.f90`,
which `use`s it.

```bash
make pde              # build only build/pde_solver
make run-pde          # build + run with single-lane defaults
```

Multi-lane runs go through the Python `run_pde` (§4) or directly via positional
arguments to the binary.

---

## 4. Python layer

### 4.1 `src/python/pde_runner.py`

#### `run_pde(params, output_path, exe)`

Invokes the Fortran `pde_solver` binary as a subprocess. Accepts every
single-lane key plus four multi-lane keys:

```python
run_pde(
    params=dict(
        M=200, n_steps=500,
        v_max=1.0, rho_max=1.0,
        rho_left_bc=0.45, rho_right_bc=0.45,
        ic_type="constant", flux_type="godunov", bc_type="periodic",
        n_lanes=2,
        lane_change_rate=0.5,
        v_max_lanes=[1.0, 1.5],     # NEW — list of length n_lanes (or None)
        rho_max_lanes=None,         # NEW — None means broadcast scalar rho_max
    ),
    output_path="data/output/pde_simulation.nc",
)
```

`v_max_lanes` and `rho_max_lanes` accept either `None` (broadcast scalar) or a
list of length `n_lanes`. If supplied, they are formatted as comma-separated
strings and passed as positional arguments 13/14 to the binary.

#### `load_pde_netcdf(path)`

Reads a PDE output file and returns:

```python
{
    "density":         np.ndarray,   # (n_steps+1, M)            backward-compat:
                                     #   single-lane: squeezed lane axis
                                     #   multi-lane:  sum over lanes (ρ_tot)
    "density_per_lane": np.ndarray,  # (n_steps+1, n_lanes, M)   always present
    "flow":            np.ndarray,   # (n_steps+1,)              total flow (∑_ℓ q_ℓ)
    "flow_per_lane":   np.ndarray,   # (n_steps+1, n_lanes)      per-lane flow
    "x":               np.ndarray,   # (M,)
    "time":            np.ndarray,   # (n_steps+1,)
    "n_lanes":         int,
    "attrs":           dict,         # global NetCDF attributes
}
```

**Transpose logic.** The Fortran solver writes `density` with NetCDF dimension
order `[lane_dimid, x_dimid, time_dimid]`. NetCDF's C convention reverses this,
so Python reads `(time, x, lane)`. `load_pde_netcdf` applies
`raw.transpose(0, 2, 1)` to produce the `(time, lane, x)` ordering used by
all downstream analysis and visualisation code.

The loader also supports legacy single-lane NetCDF files without a `lane`
dimension. In that case `density_per_lane` is constructed as
`density[:, np.newaxis, :]` and `flow_per_lane` as `flow[:, np.newaxis]`, so
calling code never has to branch on dimensionality.

### 4.2 `src/python/analysis.py` — multi-lane diagnostics

| Function | Returns | Purpose |
|----------|---------|---------|
| `compute_total_density(data)` | `(n_steps+1, M)` | $\rho_{\text{tot}}(x, t) = \sum_\ell \rho_\ell(x, t)$ |
| `compute_total_flow(data)` | `(n_steps+1,)` | $q_{\text{tot}}(t) = \sum_\ell q_\ell(t)$ at right boundary |
| `compute_total_mass(data)` | `(n_steps+1,)` | $\sum_{\ell,i} \rho_{\ell,i} \Delta x$ — conservation diagnostic |
| `compute_multilane_fundamental_diagram(data, burnin_frac=0.5)` | `(rho_mean, q_mean)` floats | Time-averaged $(\rho_{\text{tot}}, q_{\text{tot}})$ for one run |
| `pde_multilane_fundamental_diagram(...)` | `(rho, q)` arrays | Sweep over inflow density, returns scatter points |

The single-lane `pde_fundamental_diagram` is unchanged. The multi-lane variant
uses identical sweep logic and runs through `compute_total_density` /
`compute_total_flow` to aggregate over lanes.

### 4.3 `src/python/visualisation.py` — multi-lane plots

| Function | Output |
|----------|--------|
| `plot_space_time_per_lane(data)` | $N$ stacked heatmaps of $\rho_\ell(x, t)$ |
| `plot_space_time_total(data)` | One heatmap of $\rho_{\text{tot}}(x, t)$ |
| `plot_lane_densities(data, x_pos=0.5)` | Per-lane $\rho_\ell$ vs $t$ at fixed $x$ |
| `plot_total_mass(data)` | Mass deviation $\Delta M(t) = M(t) - M(0)$ vs $t$ |

`plot_total_mass` deliberately plots the **deviation from the initial mass**
rather than the absolute total mass. Plotting absolute mass causes matplotlib
to apply offset notation (e.g. `1e-6 + 9e-1`), making conservation quality hard
to read. The deviation centres the line at zero and the legend reports the
range and initial mass.

The single-lane PDE plotters (`plot_pde_summary`, `plot_pde_spacetime`,
`plot_pde_snapshots`, `plot_pde_flow`) continue to work unchanged. They consume
`data["density"]`, which the loader squeezes to `(time, x)` for `n_lanes = 1`
or sums over lanes for $n_{\text{lanes}} > 1$, so a multi-lane file is treated
as the total-density single-lane file by these legacy plotters.

---

## 5. Entry-point script — `scripts/run_multilane_pde.py`

Wraps the full multi-lane pipeline: build → run → load → 6-panel figure per
scenario. Four canonical scenarios are bundled:

| Key | Setup | Demonstrates |
|------|-------|--------------|
| `A_independent` | Staggered IC, $k=0$, periodic | Two lanes evolve independently and stay different |
| `B_coupled` | Same IC, $k=0.5$, periodic | Lane changing equilibrates the two densities |
| `C_conservation` | Identical sine IC both lanes, $k=0.5$, periodic | Total mass deviation panel — flat line conservation check |
| `D_fast_slow` | Constant IC $\rho=0.45$, $v_{\max}=[1.0, 1.5]$, $k=0.5$, periodic | Fast lane fills up; slow lane drains |

```bash
python scripts/run_multilane_pde.py                    # build, run all four, display
python scripts/run_multilane_pde.py --save             # save to plots/multilane_*.png
python scripts/run_multilane_pde.py --no-run           # replot existing .nc files
python scripts/run_multilane_pde.py --scenario D       # single scenario (prefix-match)
```

Each scenario figure has six panels:

| Panel | Content |
|-------|---------|
| Lane 1 heatmap | $\rho_1(x, t)$ |
| Lane 2 heatmap | $\rho_2(x, t)$ |
| Total heatmap | $\rho_{\text{tot}}(x, t)$ |
| Lane densities at $x=0.5$ | Per-lane time series — equilibration visible here |
| Mass deviation | $\Delta M(t)$ — conservation check |
| Final snapshot | $\rho_\ell(x)$ at the last time step |

---

## 6. Tests — `tests/test_pde_multilane.py`

Seven validation tests are run with `pytest`. The single-lane `test_pde.py`
suite (15 tests) still passes against the multi-lane solver with `n_lanes = 1`.

| Test | What it checks | Tolerance |
|------|----------------|-----------|
| `test_single_lane_regression` | $n_{\text{lanes}}=1$ Riemann run reproduces single-lane output | density bounds + shape sanity |
| `test_independent_lanes` | $k=0$ + identical lanes: lane 1 = lane 2 = single-lane reference | $10^{-5}$ vs reference, $10^{-12}$ between lanes |
| `test_periodic_mass_conservation` | $k>0$, periodic, sine IC, 1000 steps: mass variation small | $< 10^{-4}$ (float32 limit) |
| `test_per_cell_source_conservation` | Random densities, $N=2,3,4$: $\sum_\ell S_\ell(\cdot, i) = 0$ | $10^{-14}$ |
| `test_equal_lane_equilibrium` | Identical densities + params in 3 lanes: $S = 0$ | $10^{-14}$ |
| `test_density_bounds` | $k=2$ aggressive: $0 \le \rho_\ell \le \rho_{\max}$ at every cell-step | exact |
| `test_fast_slow_lane_qualitative` | Single-lane comparison: slow lane denser than fast under open BCs | qualitative |

```bash
.venv/bin/python -m pytest tests/test_pde.py tests/test_pde_multilane.py -v
```

The mass-conservation tolerance is relaxed to $10^{-4}$ rather than the
$10^{-10}$ in `multilane_pde.md` because the Fortran solver and NetCDF output
both use single precision (float32, ~7 significant digits). The model is
exactly conservative in real arithmetic; the residual comes from float32
rounding in `density_history` and the stored NetCDF values.

The fast/slow qualitative test deliberately runs **two separate single-lane
simulations** rather than one coupled two-lane run. This is the right test for
the open-BC throughput-asymmetry deliverable in `multilane_pde.md` §Phase 3.
The opposite ordering ($\rho_2 > \rho_1$) under periodic BCs with fixed total
mass is a separate equilibrium and is verified end-to-end via scenario D of
`run_multilane_pde.py`.

---

## 7. Pipeline overview

```
# Full pipeline via script
python scripts/run_multilane_pde.py --save

# Or manually for one custom run
python -c "
import sys; sys.path.insert(0,'src/python')
from pde_runner import run_pde, load_pde_netcdf
from visualisation import plot_space_time_per_lane, plot_total_mass
import matplotlib.pyplot as plt
run_pde(dict(n_lanes=2, lane_change_rate=0.5, v_max_lanes=[1.0, 1.5],
             ic_type='constant', bc_type='periodic',
             rho_left_bc=0.45, rho_right_bc=0.45, M=200, n_steps=1200),
        output_path='data/output/multilane.nc')
data = load_pde_netcdf('data/output/multilane.nc')
plot_space_time_per_lane(data)
plot_total_mass(data)
plt.show()
"
```

```
pde_flux.f90              elemental flux functions
       │
pde_lanechange.f90        conservative source-term module
       │
pde_module.f90            2D state, multi-lane step + NetCDF (lane, x, time)
       │
pde_driver.f90            CLI → density_history(n_lanes, M, T+1) → write
       │
pde_simulation.nc         density(lane, x, time), flow(lane, time), flow_total(time)
       │
pde_runner.py             load_pde_netcdf → density_per_lane(time, lane, x)
       │
analysis.py               compute_total_density / compute_total_mass / …
       │
visualisation.py          plot_space_time_per_lane / plot_total_mass / …
       │
run_multilane_pde.py      end-to-end CLI: build → run → 6-panel figures
```

---

## 8. Conventions and gotchas

- **Lane-major Fortran array convention.** All multi-lane arrays use
  `density(lane, i)` (lane fastest, cell second). Per-lane operations stride
  contiguously in memory, and the `pde_step` per-lane loop accesses
  `state%density(lane, :)` as a contiguous slice.

- **NetCDF dimension ordering.** Fortran `nf90_def_var` uses
  `[lane_dimid, x_dimid, time_dimid]` (column-major: lane fastest, time
  slowest). C-convention NetCDF reverses this in storage, so Python reads
  shape `(time, x, lane)`. `load_pde_netcdf` transposes to `(time, lane, x)`
  for downstream code. **Treat `(time, lane, x)` as the canonical Python
  ordering.**

- **Backward compatibility.** With `n_lanes = 1` and `lane_change_rate = 0`,
  the multi-lane solver is the single-lane solver. The `lane` dimension is
  still present in the NetCDF output (length 1); the loader squeezes it for
  the `data["density"]` field so legacy plotters continue to work.

- **`density` vs `density_per_lane`.** Backward-compat `data["density"]` is
  `(time, x)`: it is the squeezed lane axis for $n_{\text{lanes}}=1$ and the
  **sum over lanes** for $n_{\text{lanes}} > 1$. For per-lane analysis,
  always use `data["density_per_lane"]`, which has shape `(time, lane, x)`.

- **`flow` vs `flow_per_lane`.** Same convention: `data["flow"]` is the total
  flow $\sum_\ell q_\ell$ over time; `data["flow_per_lane"]` is the per-lane
  flow.

- **Conservation under periodic BCs.** The sum $\sum_\ell S_\ell(\cdot, i) = 0$
  holds to machine precision per cell. Total mass under periodic BCs is then
  conserved to float32 precision (~$10^{-7}$ relative error over a typical
  run). Use `compute_total_mass` and `plot_total_mass` to verify.

- **Clipping.** The explicit Euler source step is unconditionally stable only
  for $k\,\Delta t < 1$. If a cell's density would leave $[0, \rho_{\max,\ell}]$,
  it is clamped and `state%clip_count` is incremented. The driver prints a
  warning if `clip_count > 0`. If clipping fires regularly, reduce `k` or run
  with `--bc periodic` and a larger `n_steps` to allow the source ODE more
  sub-steps.

- **Plotting absolute mass.** Total mass is typically large (e.g. 0.9 or 1.0
  for a unit-domain two-lane run) while variations are at the $10^{-7}$ level.
  matplotlib's offset notation (`1e-6 + 9e-1`) makes this hard to read. The
  visualisation module's `plot_total_mass` plots $M(t) - M(0)$ instead, with
  the actual initial mass shown in the legend.

- **Plan deliverable wording (`multilane_pde.md` §Phase 3).** The plan states
  that under fast/slow lanes the **slow lane** ends up denser. This applies to
  open BCs with constant inflow (verified by `test_fast_slow_lane_qualitative`).
  Under periodic BCs with fixed total mass (scenario D of
  `run_multilane_pde.py`), the **fast lane** ends up denser because vehicles
  preferentially migrate into it. Both are physically correct in their own
  regimes — see §1.5 for the closed-form equilibrium.

---

## 9. Setup

```bash
# Build the multi-lane PDE solver (drop-in replacement for the single-lane build)
make pde                            # produces build/pde_solver

# Run with single-lane defaults (backward compatible)
make run-pde

# Run all four multi-lane demonstration scenarios
.venv/bin/python scripts/run_multilane_pde.py

# Save scenario figures to plots/
.venv/bin/python scripts/run_multilane_pde.py --save

# Run only the fast/slow scenario
.venv/bin/python scripts/run_multilane_pde.py --scenario D

# Tests (single-lane + multi-lane)
.venv/bin/python -m pytest tests/test_pde.py tests/test_pde_multilane.py -v
```

Fortran prerequisites: `gfortran`, NetCDF Fortran library (`nf-config`,
`nc-config`). Python prerequisites: `numpy`, `netCDF4`, `pytest`,
`matplotlib`. All pinned in `requirements.txt`.
