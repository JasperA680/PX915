# PDE Traffic Flow Model — Implementation Plan

This document specifies the implementation of the continuum PDE-based traffic flow model for the PX915 Traffic Flow Simulation project. It is intended as a working brief for Claude Code.

## Context

The project already contains a working 1D TASEP cellular automaton (Toy Model) implemented in Fortran with a Python interface, NetCDF-based I/O, and a modular architecture (see `src/fortran/`, `src/python/`). This plan adds the **PDE solver** as a parallel modelling approach. The PDE module must integrate with the existing I/O, analysis, and visualisation modules so that PDE outputs are interchangeable with TASEP outputs for downstream comparison.

Target model: the **Lighthill–Whitham–Richards (LWR)** scalar conservation law:

```
∂ρ/∂t + ∂q(ρ)/∂x = 0
```

with a **Greenshields fundamental diagram** as the default closure:

```
v(ρ) = v_max · (1 − ρ/ρ_max)
q(ρ) = ρ · v(ρ) = v_max · ρ · (1 − ρ/ρ_max)
```

## Architectural Constraints (must follow)

- **Language split**: Core solver in **Fortran 2008**; orchestration, analysis, and plotting in **Python 3**.
- **I/O**: All simulation output goes through **NetCDF**, using the same schema conventions already established by the TASEP code. Parameters live in the NetCDF header and must be sufficient for exact restarts.
- **Module layout**: Mirror the existing structure. Add files only where indicated below; do not restructure existing TASEP code.
- **Checkpoint/restart**: Every C time steps (user-configurable), write the full state to a `/restart` group in the NetCDF file. PDE restart fields per SDP §2.7: `density(1:M)`, `t_current`, and parameters `Δx, Δt, domain_length, boundary_values, C`.
- **No new external dependencies** beyond NumPy, Matplotlib, and the existing NetCDF library already used by the TASEP code.

## Files to Create

```
src/fortran/pde_module.f90        # core solver
src/fortran/pde_flux.f90          # fundamental diagram + numerical flux
src/python/pde_runner.py          # Python entry point analogous to the TASEP runner
scripts/run_pde.py                # command-line entry point
tests/test_pde.py                 # validation suite
notebooks/pde_demo.ipynb          # worked example for documentation
```

Update (do not rewrite) `src/python/io.py`, `src/python/analysis.py`, and `src/python/visualisation.py` only as needed to consume PDE output. The output schema must allow these modules to operate without branching on model type.

---

## Phase 1 — Mathematical Setup and Module Skeleton

**Goal**: lay the foundations so subsequent phases plug in cleanly.

1. Create `src/fortran/pde_flux.f90` exposing:
   - `function v_of_rho(rho, v_max, rho_max) result(v)` — Greenshields velocity
   - `function q_of_rho(rho, v_max, rho_max) result(q)` — flux
   - `function dq_drho(rho, v_max, rho_max) result(dq)` — characteristic speed (used for CFL)
   - `function rho_critical(rho_max) result(rc)` — argmax of q (= rho_max/2 for Greenshields)
   - All routines must accept scalar or array input via `elemental`.

2. Create `src/fortran/pde_module.f90` with the following public interface (skeletons only at this phase):
   - `subroutine pde_initialise(state, params)` — allocate arrays, set initial condition
   - `subroutine pde_step(state, params)` — advance one time step
   - `subroutine pde_finalise(state)` — deallocate
   - Define a `pde_state_t` derived type holding `density(:)`, `t_current`, `step`.
   - Define a `pde_params_t` derived type holding `dx, dt, domain_length, M, v_max, rho_max, rho_left_bc, rho_right_bc, C_checkpoint`.

3. Add unit tests in `tests/test_pde.py` for the flux functions:
   - `q(0) == 0`, `q(rho_max) == 0`
   - `q(rho_critical) == v_max * rho_max / 4` (analytical maximum)
   - `dq_drho(rho_critical) == 0`
   - Signs of `dq_drho` correct on either side of `rho_critical`

**Deliverable**: code compiles, flux tests pass, no solver logic yet.

---

## Phase 2 — Numerical Scheme

**Goal**: working finite volume solver, validated on simple cases.

4. Implement the **Godunov numerical flux** for a concave flux function in `pde_flux.f90`:

   ```
   F(ρ_L, ρ_R) =
     if ρ_L ≤ ρ_R:
       min over ρ ∈ [ρ_L, ρ_R] of q(ρ)
     else:
       max over ρ ∈ [ρ_R, ρ_L] of q(ρ)
   ```

   For concave q with maximum at `rho_critical`, this has a closed form:
   - If `ρ_L ≤ ρ_R`: `F = q(ρ_L) if ρ_L ≥ rho_critical else (q(ρ_R) if ρ_R ≤ rho_critical else q(ρ_L))` — implement carefully, the standard cases are: both subcritical → q(ρ_L); both supercritical → q(ρ_R); ρ_L sub, ρ_R super → min(q(ρ_L), q(ρ_R)); ρ_L super, ρ_R sub → q(rho_critical).
   - If `ρ_L > ρ_R`: same case structure but using max.

   Reference any standard text (LeVeque, *Finite Volume Methods for Hyperbolic Problems*) and add a docstring citing the case table.

5. Implement `pde_step` in `pde_module.f90` using the conservative update:

   ```
   ρ_i^(n+1) = ρ_i^n − (Δt/Δx) · [F(ρ_i, ρ_{i+1}) − F(ρ_{i−1}, ρ_i)]
   ```

   Use ghost cells for boundary handling (one ghost on each side is sufficient for first-order Godunov).

6. Implement an **adaptive CFL-controlled Δt** as a helper subroutine `compute_dt(state, params, cfl_number)`:

   ```
   Δt = cfl_number · Δx / max_i |dq_drho(ρ_i)|
   ```

   Default `cfl_number = 0.9`. The user may supply a fixed Δt instead, in which case validate it against the CFL bound and abort with a clear error if violated.

7. Optional but recommended: implement **Lax–Friedrichs** as a second flux choice in `pde_flux.f90` selectable via a parameter. Useful as a baseline because it is simpler and more diffusive — handy for debugging Godunov behaviour.

**Deliverable**: solver advances a constant initial condition without drift (machine-precision conservation under periodic BCs).

---

## Phase 3 — Initial and Boundary Conditions

**Goal**: enough flexibility to set up validation problems and realistic traffic scenarios.

8. Implement initial condition presets in `pde_module.f90`, selectable by a string parameter:
   - `"constant"` — uniform density `ρ_0`
   - `"riemann"` — step function: `ρ_L` for `x < x_split`, `ρ_R` otherwise
   - `"gaussian"` — smooth bump for visualising shock formation
   - `"sine"` — sinusoidal perturbation around a mean density

9. Implement boundary conditions:
   - **Open boundaries** (default, matching TASEP setup): fix `ρ_left_bc` at the inflow ghost cell, `ρ_right_bc` at the outflow ghost cell. This is consistent with the α/β interpretation in the TASEP model and is essential for the eventual PDE-vs-CA comparison.
   - **Periodic boundaries**: copy first physical cell into right ghost, last physical cell into left ghost. Used for conservation tests.
   - Boundary type selectable via parameter `bc_type ∈ {"open", "periodic"}`.

**Deliverable**: all four ICs and both BCs work; selecting them does not require recompilation.

---

## Phase 4 — Validation

**Goal**: prove the solver is correct before integrating downstream.

10. Add the following tests to `tests/test_pde.py`. Each test runs the Fortran solver via the Python runner and checks the output:

    a. **Constant solution preservation**: uniform ρ stays uniform to within 1e-12 over 1000 steps.

    b. **Mass conservation (periodic BCs)**: total mass `Σ ρ_i · Δx` is constant to within 1e-10 over 1000 steps for a Gaussian IC.

    c. **Shock speed (Riemann problem)**: with `ρ_L = 0.3·ρ_max`, `ρ_R = 0.8·ρ_max` (both on the congested side relative to the mixing — pick values that produce a shock under Greenshields), the discontinuity propagates at the Rankine–Hugoniot speed `s = [q]/[ρ] = (q(ρ_R) − q(ρ_L)) / (ρ_R − ρ_L)`. Allow tolerance of one cell width per 100 steps to account for numerical smearing.

    d. **Rarefaction shape**: with `ρ_L > ρ_R` chosen to produce a rarefaction, compare the numerical solution to the exact self-similar fan at a fixed time. L1 error should be O(Δx).

    e. **Convergence**: refine Δx by factors of 2 and verify L1 error against an exact Riemann solution decreases at the expected first-order rate.

11. Provide a Python helper `tests/exact_riemann.py` that returns the exact LWR Riemann solution for a Greenshields flux at given (x, t). Tests (c)–(e) consume this.

**Deliverable**: full test suite passes; CI (if configured) green.

---

## Phase 5 — I/O Integration

**Goal**: PDE outputs are first-class citizens of the existing analysis pipeline.

12. Extend the NetCDF schema to include PDE output. The schema should:
    - Store `density(time, x)` as the primary state variable (analogous to the TASEP `lattice(time, x)`).
    - Store `flow(time)` computed as `q(ρ)` evaluated at the right boundary, for direct comparison with TASEP exit flow.
    - Store all parameters in the file header (group attributes): `dx, dt, domain_length, M, v_max, rho_max, ic_type, bc_type, rho_left_bc, rho_right_bc, C_checkpoint, model="LWR-Greenshields"`.
    - Use the same dimension and variable naming conventions as the TASEP output where semantically equivalent (e.g. `time` dimension shared).

13. Implement checkpoint writing in `pde_module.f90`:
    - Every `C_checkpoint` steps, write `density(1:M)`, `t_current`, `step` to a `/restart` NetCDF group.
    - Overwrite the previous checkpoint (single rolling checkpoint to bound disk usage); make this configurable later if needed.

14. Implement restart in `pde_runner.py`:
    - Add a `--restart <path>` CLI flag.
    - On restart, read parameters from the NetCDF header, state from `/restart`, and resume from the saved step.
    - Validate that supplied CLI parameters do not conflict with the checkpointed parameters; abort with a clear error if they do.

**Deliverable**: a PDE run can be killed mid-execution and resumed exactly.

---

## Phase 6 — Python Interface, Analysis, and Visualisation

**Goal**: a user can run, analyse, and plot a PDE simulation with a single command, and compare it against a TASEP run with minimal code.

15. Implement `src/python/pde_runner.py`:
    - Function `run_pde(params: dict, output_path: str) -> None` that writes the parameter NetCDF header, invokes the Fortran executable, and returns when complete.
    - Mirror the function signature and conventions of the existing TASEP runner.

16. Implement `scripts/run_pde.py` as a thin CLI wrapper around `run_pde`, supporting:
    - `--ic`, `--bc`, `--rho-max`, `--v-max`, `--domain-length`, `--M`, `--cfl`, `--n-steps`, `--checkpoint-interval`, `--output`, `--restart`.

17. Extend `src/python/analysis.py` with PDE-specific helpers (only where the TASEP equivalents do not naturally apply):
    - `compute_density_profile(nc_file, t)` — return ρ(x) at time t.
    - `compute_flow_timeseries(nc_file)` — return q(ρ) at right boundary vs t.
    - `compute_fundamental_diagram(nc_file)` — scatter of (ρ_avg, q_avg) over time, for direct overlay against TASEP fundamental diagrams.

18. Extend `src/python/visualisation.py` with:
    - `plot_space_time(nc_file)` — heatmap of ρ(x, t). This is the primary visual diagnostic.
    - `plot_density_snapshots(nc_file, times)` — line plots at selected times.
    - The fundamental diagram plotter should already work via the analysis output if the schema is consistent.

19. Create `notebooks/pde_demo.ipynb` demonstrating:
    - A Riemann problem producing a shock (with the analytical shock speed overlaid).
    - A Riemann problem producing a rarefaction.
    - A Greenshields fundamental diagram derived from a parameter sweep over inflow density.

**Deliverable**: a new user can run `python scripts/run_pde.py --ic riemann --output run.nc` and then open the demo notebook to visualise the result.

---

## Out of Scope for This Plan

The following are explicitly **not** part of the initial PDE implementation. They are noted here so they are not accidentally included:

- **Bifurcations and merges** for the PDE — flagged in SDP §2.1.3 as future work; will be added in a separate phase once the single-road solver is validated.
- **PDE-vs-CA quantitative comparison** — Stephan owns this per SDP §3.1; this plan only ensures the I/O schema makes it possible.
- **UQ and parameter sensitivity analysis** — Jasper's contingent task per SDP §3.1.
- **Higher-order schemes** (MUSCL, WENO) — Godunov is sufficient for the project scope; revisit only if shock smearing turns out to be a problem for the model comparison.
- **Multi-class / heterogeneous traffic in the PDE** — explicitly listed as a microscopic-only extension in the SDP.

---

## Definition of Done

The PDE implementation is complete when:

1. All five validation tests in Phase 4 pass.
2. A PDE simulation can be run, checkpointed, killed, and restarted exactly.
3. A space-time plot of a Riemann problem visibly shows the correct shock or rarefaction structure.
4. The output NetCDF file can be loaded by the same `analysis.py` functions that handle TASEP output, with the fundamental diagram plotter producing a meaningful curve.
5. The demo notebook runs end-to-end without manual intervention.
6. Code is documented with Fortran docstrings and the Python public API has NumPy-style docstrings suitable for ReadTheDocs.

---

## Suggested First Working Session

If picking a single tractable starting point, do Phases 1 and 2 together: get `pde_flux.f90`, the `pde_module.f90` skeleton, and a working **Lax–Friedrichs** step (not yet Godunov) producing visible shock motion on a Riemann problem, dumped to NetCDF, plotted in Python. Once that end-to-end pipeline works, replacing Lax–Friedrichs with Godunov and adding the validation suite is incremental rather than structural.
