# Multi-Lane PDE Traffic Flow Model — Implementation Plan

This document specifies the multi-lane extension of the existing single-lane LWR PDE solver for the PX915 Traffic Flow Simulation project. It is intended as a working brief for Claude Code.

## Context

The single-lane PDE solver is already implemented per `PDE.md`. Before making any changes, **inspect the actual repository state** (`src/fortran/pde_*.f90`, `src/python/pde_runner.py`, `tests/test_pde.py`, the NetCDF schema in `src/python/io.py`) and confirm the current module layout, derived type names, and parameter conventions. This plan assumes those names but defers to the repository where they differ.

The single-lane model solves:

```
∂ρ/∂t + ∂q(ρ)/∂x = 0
```

with Greenshields closure, Lax–Friedrichs or Godunov fluxes, and open or periodic BCs. This plan extends it to one density field per lane:

```
ρ_1(x,t), ρ_2(x,t), …, ρ_N(x,t)
```

with a longitudinal flux per lane and conservative lateral source terms representing lane changes. For a two-lane case:

```
∂ρ_1/∂t + ∂q_1(ρ_1)/∂x = S_{2→1} − S_{1→2}
∂ρ_2/∂t + ∂q_2(ρ_2)/∂x = S_{1→2} − S_{2→1}
```

The non-negotiable conservation requirement: lane-change source terms must sum to zero across lanes at every cell, so total mass over all lanes is exactly conserved (to machine precision under periodic BCs).

## Architectural Constraints (must follow)

- **Backward compatibility**: with `n_lanes = 1`, the new solver must reproduce the existing single-lane solver bit-for-bit (or to within 1e-12 floating-point tolerance). This is enforced as a regression test in Phase 5.
- **No restructuring of unrelated code**: only modify TASEP, Nagel–Schreckenberg, or non-PDE Python modules where strictly necessary to accommodate the new NetCDF dimension.
- **Schema evolution**: the multi-lane NetCDF schema must remain readable by the existing `analysis.py` and `visualisation.py` functions in the single-lane case. For the multi-lane case, add a `lane` dimension; existing functions should treat absence of this dimension as `n_lanes = 1`.
- **Fortran array convention**: use `density(1:n_lanes, 1:M)` (lane-major). Document this in a header comment in `pde_module.f90` and stick to it everywhere — no transposition mid-pipeline.
- **No new external dependencies**.

## Files to Modify or Create

Inspect the repo first to confirm exact filenames. Expected changes:

```
src/fortran/pde_module.f90        # extend state type, finite-volume step
src/fortran/pde_flux.f90          # no API changes; called per-lane
src/fortran/pde_lanechange.f90    # NEW — lane-change source terms
src/python/pde_runner.py          # add multi-lane CLI parameters
src/python/io.py                  # add `lane` dimension to NetCDF schema
src/python/analysis.py            # add per-lane and total diagnostics
src/python/visualisation.py       # add per-lane and total plots
scripts/run_pde.py                # expose new CLI flags
tests/test_pde_multilane.py       # NEW — multi-lane validation suite
notebooks/pde_multilane_demo.ipynb # NEW — worked two-lane example
```

Do **not** create a parallel `pde_multilane_module.f90`. The single-lane solver becomes the `n_lanes = 1` special case of the multi-lane solver — there should be exactly one PDE solver in the codebase at the end.

---

## Phase 1 — Independent Lanes (no lane changing)

**Goal**: extend the data structures and finite-volume step to support N lanes with `S_{ℓ→m} = 0`. This is a pure refactor — physics is unchanged per lane.

1. Extend `pde_state_t` in `pde_module.f90`:
   - `density(:,:)` → shape `(n_lanes, M)`.
   - Add `n_lanes` as a state field.
   - Update allocation, deallocation, and any direct indexing.

2. Extend `pde_params_t`:
   - Add `n_lanes` (integer, default 1).
   - Add `v_max(:)` and `rho_max(:)` as allocatable arrays of length `n_lanes`. For Phase 1, populate all entries identically from the existing scalar parameters.
   - Keep the existing scalar parameter names as deprecated aliases that fill the arrays — this preserves CLI backward compatibility.

3. Update `pde_step` to loop over lanes for the longitudinal update:
   ```fortran
   do lane = 1, params%n_lanes
     call apply_boundary_conditions(state%density(lane,:), params, lane)
     do i = 0, M
       flux(lane, i) = numerical_flux( &
         rho_ext(lane, i), rho_ext(lane, i+1), &
         params%v_max(lane), params%rho_max(lane), params%scheme)
     end do
     do i = 1, M
       state%density(lane, i) = state%density(lane, i) &
         - params%dt/params%dx * (flux(lane, i) - flux(lane, i-1))
     end do
   end do
   ```

4. Update the CFL helper to take the maximum wave speed over all lanes:
   ```
   Δt ≤ CFL · Δx / max_{ℓ,i} |q'_ℓ(ρ_{ℓ,i})|
   ```
   For Greenshields the conservative bound is `Δt ≤ CFL · Δx / max_ℓ v_max(ℓ)`.

5. Update boundary condition handling:
   - **Periodic**: `ρ_{ℓ,0} = ρ_{ℓ,M}`, `ρ_{ℓ,M+1} = ρ_{ℓ,1}` for every lane.
   - **Open**: `ρ_{ℓ,0} = rho_left_bc(ℓ)`, `ρ_{ℓ,M+1} = rho_right_bc(ℓ)`. For Phase 1, broadcast the existing scalar boundary values to all lanes; lane-specific values come in Phase 3.

**Deliverable**: with `n_lanes = 2` and identical lane parameters, the two lanes evolve independently and identically. The `n_lanes = 1` case reproduces the existing solver bit-for-bit.

---

## Phase 2 — Conservative Lane-Change Source Terms

**Goal**: add lane-changing physics while guaranteeing total-mass conservation across lanes.

6. Create `src/fortran/pde_lanechange.f90` exposing:
   - `subroutine compute_lane_change_sources(density, params, source)`
     - Inputs: `density(1:n_lanes, 1:M)`, `params`.
     - Output: `source(1:n_lanes, 1:M)` such that `sum(source(:, i)) == 0` for every i.

7. For two lanes, implement the asymmetric exchange model:
   ```
   S_{1→2,i} = k · ρ_{1,i} · max(0, 1 − ρ_{2,i}/ρ_max(2)) · max(0, v_2(ρ_{2,i}) − v_1(ρ_{1,i}))
   S_{2→1,i} = k · ρ_{2,i} · max(0, 1 − ρ_{1,i}/ρ_max(1)) · max(0, v_1(ρ_{1,i}) − v_2(ρ_{2,i}))
   source(1,i) = S_{2→1,i} − S_{1→2,i}
   source(2,i) = S_{1→2,i} − S_{2→1,i}
   ```
   Note: at most one of `S_{1→2}` and `S_{2→1}` is nonzero per cell because of the asymmetric `max(0, Δv)` factors.

8. **Generalise to N ≥ 2 lanes** by restricting exchanges to **adjacent lanes only** (lane ℓ exchanges with ℓ−1 and ℓ+1, no end wrap). Apply the same pairwise exchange formula to each adjacent pair, and accumulate into `source`. This keeps the model physically meaningful (a car on the inside lane cannot teleport to the outside lane in one step) and keeps the conservation property pairwise — easier to verify than a fully-coupled all-to-all model.

9. Extend `pde_step` to apply the source after the longitudinal update:
   ```fortran
   call compute_lane_change_sources(state%density, params, source)
   do lane = 1, params%n_lanes
     do i = 1, M
       state%density(lane, i) = state%density(lane, i) + params%dt * source(lane, i)
     end do
   end do
   ```

10. Add a positivity safeguard after the source step. Two acceptable approaches; pick one and document it:
    - **Clipping**: `density = max(0, min(density, rho_max(lane)))` — simple, breaks strict conservation by a tiny amount, log a warning if clipping fires.
    - **Adaptive sub-stepping for the source**: split the `dt` source update into smaller sub-steps if any cell would go out of bounds. More work but conservation-preserving.

    For the initial implementation use clipping with a warning counter; flag adaptive sub-stepping as a TODO if clipping fires more than rarely in test scenarios.

11. Add a stability comment in `pde_lanechange.f90` noting that large `k` can violate stability of the source ODE; recommend `k · dt < 1` as a rule of thumb and consider adding it to the CFL check as `Δt ≤ min(CFL · Δx / max v_max, 1/k)`.

**Deliverable**: with `k > 0` and unequal initial densities, vehicles flow between lanes; total mass across lanes is conserved to within 1e-10 over 1000 steps under periodic BCs.

---

## Phase 3 — Lane-Dependent Parameters

**Goal**: allow lanes to have different `v_max` and `rho_max` so the model can represent fast/slow lanes, lane closures, and bottlenecks.

12. Wire `v_max(:)` and `rho_max(:)` through to the flux calls (already done structurally in Phase 1; this phase ensures they're actually different per lane rather than broadcast from a scalar).

13. Add lane-specific boundary conditions: `rho_left_bc(:)` and `rho_right_bc(:)`. If the user supplies a scalar via the CLI, broadcast it to all lanes. If they supply a comma-separated list of length `n_lanes`, use it directly.

14. Test that a fast lane (high `v_max`) draws vehicles from a slow lane (low `v_max`) under the lane-change model — this is the qualitative correctness check that the asymmetric `max(0, Δv)` factor does what it should.

**Deliverable**: a two-lane simulation with `v_max = [1.0, 1.5]` and an initially uniform total density develops a higher steady-state density in the slower lane and a higher flux in the faster lane.

---

## Phase 4 — I/O and Diagnostics

**Goal**: multi-lane outputs are first-class in the NetCDF schema and analysis pipeline.

15. Extend the NetCDF schema:
    - Add a `lane` dimension of length `n_lanes`.
    - Store `density(time, lane, x)` (Python view; the underlying Fortran storage is `(lane, x, time)` — document the transpose explicitly in `io.py`).
    - Store `flow(time, lane)` per lane plus `flow_total(time)` as the sum.
    - Header attributes: `n_lanes`, `v_max(:)`, `rho_max(:)`, `lane_change_rate` (k), and the existing single-lane attributes.
    - For `n_lanes = 1`, the `lane` dimension should still be present (length 1) — this keeps the analysis code uniform rather than branching on dimensionality.

16. Update checkpoint/restart (per SDP §2.7) to store the full `density(1:n_lanes, 1:M)` array and the `n_lanes`, `v_max(:)`, `rho_max(:)` parameter arrays in the `/restart` group.

17. Extend `src/python/analysis.py` with multi-lane helpers:
    - `compute_total_density(nc_file)` — `ρ_tot(x, t) = Σ_ℓ ρ_ℓ(x, t)`.
    - `compute_total_flow(nc_file)` — `q_tot(t) = Σ_ℓ q_ℓ(t)` at the right boundary.
    - `compute_total_mass(nc_file)` — `Σ_{ℓ,i} ρ_{ℓ,i} · Δx` vs t (for conservation diagnostics).
    - `compute_multilane_fundamental_diagram(nc_file)` — scatter of `(ρ_tot_avg, q_tot_avg)` over time, suitable for overlay against single-lane and TASEP diagrams.
    - Optional: `compute_lane_change_rate(nc_file)` — net flux between lanes vs time, if source histories are stored.

18. Extend `src/python/visualisation.py`:
    - `plot_space_time_per_lane(nc_file)` — N stacked heatmaps of `ρ_ℓ(x, t)`.
    - `plot_space_time_total(nc_file)` — single heatmap of `ρ_tot(x, t)`.
    - `plot_lane_densities(nc_file, x)` — per-lane density vs t at a fixed location.
    - `plot_total_mass(nc_file)` — total mass vs t (should be flat for periodic BCs; the visual check on conservation).

**Deliverable**: a two-lane run produces a NetCDF file that the existing fundamental-diagram plotter handles correctly via the `_total` aggregations, and that the new per-lane plotters can decompose for inspection.

---

## Phase 5 — Validation

**Goal**: prove the multi-lane solver is correct and backward compatible.

Add the following to `tests/test_pde_multilane.py`:

19. **Single-lane regression**: with `n_lanes = 1`, a Riemann problem produces output identical to the existing single-lane solver to within 1e-12. Run the existing single-lane test suite against the new solver as part of CI.

20. **Independent lanes**: with `k = 0` and identical per-lane parameters, two lanes initialised with different densities evolve independently — each lane matches a single-lane reference run with the same IC.

21. **Periodic mass conservation**: with periodic BCs and `k > 0`, total mass `Σ_{ℓ,i} ρ_{ℓ,i} · Δx` is constant to within 1e-10 over 1000 steps.

22. **Per-cell source conservation**: for randomly generated `(ρ_1, ρ_2, …, ρ_N)` over many cells, `sum(source(:, i)) == 0` to within 1e-14 for every i.

23. **Equal-lane equilibrium**: with identical densities and parameters in all lanes, lane-change source terms are exactly zero.

24. **Density bounds**: across all tests, assert `0 ≤ ρ_ℓ ≤ rho_max(ℓ)` at every cell, every step. If clipping ever fires, the test logs the count.

25. **Fast/slow lane qualitative test** (from Phase 3): the slower lane reaches higher steady-state density than the faster lane.

**Deliverable**: full multi-lane test suite passes; existing single-lane tests still pass; CI green.

---

## Phase 6 — CLI, Notebook, and Documentation

26. Extend `scripts/run_pde.py` with new flags:
    - `--n-lanes <int>` (default 1)
    - `--lane-change-rate <float>` (default 0)
    - `--v-max <float | comma-list>` — scalar broadcasts; list must have length `n_lanes`.
    - `--rho-max <float | comma-list>` — same convention.
    - `--rho-left-bc <float | comma-list>`, `--rho-right-bc <float | comma-list>`.

    Existing single-lane invocations must continue to work without modification.

27. Create `notebooks/pde_multilane_demo.ipynb` demonstrating:
    - A two-lane simulation with equal parameters and `k = 0` showing independent evolution.
    - The same setup with `k > 0` showing equilibration.
    - A fast/slow lane comparison reproducing the Phase 3 qualitative result.
    - A multi-lane fundamental diagram from a parameter sweep over total inflow density.

28. Update `docs/` (or wherever the PDE module is documented) to cover:
    - How the multi-lane model differs from the single-lane LWR.
    - The longitudinal flux vs lateral source split.
    - Why source terms must sum to zero across lanes.
    - The `(lane, x, time)` Fortran storage vs `(time, lane, x)` Python view convention.
    - How to run a basic two-lane simulation.
    - How to generate a multi-lane fundamental diagram.

    Include this short framing paragraph at the top of the docs:

    > A multi-lane PDE extension replaces the scalar density field with one density field per lane. Each lane obeys its own LWR conservation law with a lane-specific flux, while additional source terms transfer vehicles between adjacent lanes. These source terms model lane-changing behaviour and must be constructed so that total vehicle count is conserved across lanes. Numerically, the finite-volume update is applied independently per lane for the longitudinal fluxes, followed by a conservative source-term update for lateral exchange.

---

## Out of Scope

Explicitly **not** included in this plan:

- **Lane-changing in the TASEP or Nagel–Schreckenberg models** — Lucas owns these per SDP §3.1.
- **Bifurcations and merges in the multi-lane PDE** — handled separately once the single-road multi-lane solver is validated; junction conditions for multi-lane PDEs require their own design.
- **Higher-order schemes** — Godunov is sufficient.
- **All-to-all lane coupling** — adjacent-lane-only is the model. If a future user needs cross-lane coupling, the source-term module is the right place to extend.
- **Driver heterogeneity within a lane** — the PDE is a continuum model; per-driver behaviour belongs in the Nagel–Schreckenberg extension.
- **Adaptive sub-stepping for the source term** — clipping with a warning counter is the initial approach; only escalate if clipping fires in normal use.

---

## Definition of Done

The multi-lane PDE implementation is complete when:

1. The existing single-lane test suite passes against the new solver with `n_lanes = 1`.
2. All seven multi-lane validation tests in Phase 5 pass.
3. A two-lane Riemann problem can be run, checkpointed, killed, and restarted exactly.
4. The fast/slow lane qualitative test produces the expected steady-state asymmetry.
5. Multi-lane NetCDF output is consumed correctly by the existing fundamental-diagram plotter via the `_total` aggregations.
6. The demo notebook runs end-to-end without manual intervention.
7. Code has Fortran docstrings; new Python public API has NumPy-style docstrings; the docs page is updated.

---

## Suggested First Working Session

Do Phase 1 only. Get `n_lanes = 2` running with `k = 0` and identical lane parameters, output to NetCDF with the new `lane` dimension, and a per-lane space-time plot showing two identical evolutions. Once that works — and the `n_lanes = 1` regression test still passes — Phase 2 is incremental. The riskiest part of the whole plan is the array reshape in Phase 1; everything after that is additive.
