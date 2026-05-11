# PX915 — Road Network Junction Model: Progress Summary

> **Group**: Tristan McCarthy, Jasper Allen, Stephan Gambert, Lucas Belz-Koeling
> **Deadline**: 28 May 2026 — University of Warwick

---

## Table of Contents

[TOC]

---

## Overview

The project implements a traffic flow simulation package built on a **1D TASEP** (Totally Asymmetric Simple Exclusion Process) cellular automaton. Vehicles hop right with probability 1, enter at site 1 with probability α, and exit at site L with probability β.

The network extension connects multiple TASEP lanes at junctions with realistic UK right-of-way rules. The original single-lane model (`tasep.f90`, `simulation.f90`, `io.f90`) is untouched; all network code lives in separate modules.

---

## Current Architecture

Five modules, strictly ordered (no circular dependencies):

```
vehicle_mod          → road_network_mod → network_init_mod
                                       → junction_mod
                                       → network_simulation_mod
```

| File | Module | Role |
|------|--------|------|
| `src/fortran/vehicle.f90` | `vehicle_mod` | Cell state codes: `V_EMPTY=0`, `V_OCCUPIED=1` |
| `src/fortran/road_network.f90` | `road_network_mod` | All data types + snapshot/count utilities |
| `src/fortran/network_init.f90` | `network_init_mod` | Topology builders: `init_crossroad`, `init_t_junction` |
| `src/fortran/junction.f90` | `junction_mod` | Right-of-way logic, yield rules, deadlock breaking |
| `src/fortran/network_simulation.f90` | `network_simulation_mod` | Timestep driver and multi-step runner |

### Data Types

```fortran
lane_t          ! unidirectional 1D TASEP lane (cells, old, alpha, beta, open_in/out)
road_t          ! bundle of lanes; end_junction(2) records which junction at each end
leg_route_t     ! real, allocatable :: prob(:)  — outbound probability distribution
junction_t      ! flat in_road/in_lane/out_road/out_lane lists + in_routes + perimeter ports
road_network_t  ! allocatable roads(:) and junctions(:)
```

### Lane Orientation

`flow_direction = +1` means site 1 is at end_1 (entry), site L is at end_2 (holding cell toward junction).
`flow_direction = -1` means site 1 is at end_2 (entry), site L is at end_1.

The holding cell (site L of an inbound lane) is always where a vehicle waits before the junction evaluates it. The destination cell (site 1 of an outbound lane) is where the junction places it.

### Routing — No Intent on the Cell

Vehicles carry **no turning intent**. Every cell is either `V_EMPTY` (0) or `V_OCCUPIED` (1). Routing decisions live entirely on the junction:

- `junction_t%in_routes(k)%prob(:)` — length `n_out`, sums to 1.
- When a vehicle arrives at holding cell of inbound leg `k`, the junction samples `dst_idx` from `prob` to decide which outbound lane it targets.

For the crossroad with `p_left` / `p_right` the distributions are built as:

| Turn | Outbound index | Probability |
|------|---------------|-------------|
| LEFT (UK easy swing) | `mod(k, 4)+1` | `p_left` |
| STRAIGHT | `mod(k+1, 4)+1` | `1 − p_left − p_right` |
| RIGHT (cuts across) | `mod(k+2, 4)+1` | `p_right` |
| U-TURN | `k` | `0` |

### Move Category

Inside the junction, categories are derived purely from the leg indices — no cell code involved:

```fortran
d = modulo(out_idx - in_idx, n)
! d=0 → UTURN, d=1 → LEFT, d=n/2 → STRAIGHT, d=n-1 → RIGHT
```

### Junction Flat Lists

The junction stores inbound and outbound legs as parallel arrays in **clockwise order**:

```
in_road(:), in_lane(:)    ! length n_in, clockwise-ordered inbound legs
out_road(:), out_lane(:)  ! length n_out, clockwise-ordered outbound legs
```

No `connected_road_ids`, `end_at`, or `inbound/outbound_lane_idx` anywhere — those have been removed.

---

## Right-of-Way Rules

### Symmetric junctions (n_in = n_out = 4, crossroad)

Three rules build the `yields_to(i,j)` matrix:

| Rule | Condition | Action |
|------|-----------|--------|
| **R1** yield-to-right | `i`'s path conflicts with `right_of_i`, OR same destination | `yields_to(i, right_of_i) = true` |
| **R2** right-turn yields oncoming | `cat_i == RIGHT` and opposite leg is LEFT or STRAIGHT | `yields_to(i, opp) = true` |
| **R3** opposite-leg scan | (STRAIGHT, LEFT) pair: LEFT yields; (RIGHT, RIGHT): mutual | `yields_to` entries set accordingly |

Path conflict (`paths_conflict_sym`): two LEFT turns from adjacent or opposite legs never conflict; all other combinations do.

### Asymmetric junctions (T-junction, n_in ≠ 4 or n_out ≠ 4)

Only **R1** applies. Conflict is determined by the **chord-crossing predicate**:

- Each leg is assigned a port position on a cyclic perimeter (0-indexed, clockwise).
- Two vehicle paths conflict if their chords intersect, i.e. exactly one of `{in_b, out_b}` lies in the open clockwise arc from `in_a` to `out_a`.
- Same destination also triggers R1.

### Approval and deadlock breaking

After building the yield matrix:

1. **Approval pass** — a vehicle is approved if it has a clear destination and no live yield target.
2. **Deadlock breaker** — a vehicle is *in deadlock* if it yields to at least one non-approved, non-deadlocked vehicle **and** yields to no already-approved vehicle. One member of each deadlock set is chosen stochastically; the rest wait this timestep.

:::warning
**Critical invariant**: the two-clause deadlock condition must be preserved verbatim. Removing the second clause (`.and. .not. any(... .and. approved)`) causes mass non-conservation — vehicles that lose an R2/R3 interaction get incorrectly added to the deadlock set and approved, writing to an already-occupied destination.
:::

### Parallel-update ordering

Every timestep runs in this exact order:
```
snapshot_network      ← freeze old state for ALL lanes
evaluate_junctions    ← read old, write cells
lane_internal_step    ← read old, write cells (bulk hops + open boundaries)
```

A vehicle placed by the junction at outbound site 1 this step is absent in `old`, so it will not move again until the next step.

---

## Supported Topologies

### 4-arm Crossroad

```
        Road 1 (N)
           ↕
Road 4 (W) ✛ Road 2 (E)
           ↕
        Road 3 (S)
```

`init_crossroad(net, lane_length, alpha, beta, p_left, p_right)`

- 4 roads, each bidirectional (2 lanes), connecting at end_1.
- `n_in = 4`, `n_out = 4`. Full R1/R2/R3 rules active.
- Open alpha/beta boundaries at end_2 of every road.

### T-Junction

```
  [open]  ←── road 1 out ─────────── [junction]
  [open]  ──→ road 1 in  ──────────↗
                              ↑
                          road 2 in (stem)
                              ↑
                           [open]
                                   ──→ road 3 out ──→ [open]
```

`init_t_junction(net, lane_length, alpha, beta, route_in_west, route_in_stem)`

- Road 1 (west arm): bidirectional, end_2 at junction.
- Road 2 (south stem): inbound only, end_2 at junction.
- Road 3 (east arm): outbound only, end_1 at junction.
- `n_in = 2`, `n_out = 2`. Chord-crossing conflict predicate.
- Perimeter ports (CW, 0-indexed): `0=east_out, 1=stem_in, 2=west_in, 3=west_out`.

---

## Test Suite

All tests build with zero warnings. Run with `make run_test_<name>`.

| Executable | Tests | What it validates |
|-----------|-------|-------------------|
| `test_types` | 4 | `V_EMPTY`/`V_OCCUPIED`, lane-at-end helpers, `count_occupied_network`, `snapshot_network` |
| `test_init_crossroad` | ~20 assertions | Topology, lane counts, flow directions, boundary flags, `in_routes` probability distributions |
| `test_junction` | 9 | All crossroad junction behaviours (see table below) |
| `test_network_run` | 2 | 5000-step mass conservation + steady-state symmetry (< 20% spread) |
| `test_t_junction` | 6 | T-junction init, merge conflict, diverge, chord non-conflict, chord conflict (cars yield), 2000-step mass conservation |

### Crossroad junction test cases

| # | Test | Setup | Key assertion |
|---|------|-------|--------------|
| 1 | `test_dest_math` | Exercise `move_category` directly | `move_category(1,2,4)==LEFT`, `(1,3,4)==STRAIGHT`, `(1,4,4)==RIGHT`, etc. |
| 2 | `test_single_move` | Leg 1 → STRAIGHT | Holding clears; road 3 site 1 = `V_OCCUPIED` |
| 3 | `test_physical_block` | Leg 1 → STRAIGHT; road 3 pre-filled | Neither cell changes |
| 4 | `test_yield_right` | Legs 1, 4 both STRAIGHT | Leg 4 advances; **leg 1 holding = `V_OCCUPIED`** (waits) |
| 5 | `test_right_turn_yields` | Leg 1 RIGHT, leg 3 STRAIGHT | Leg 3 advances; **leg 1 holding = `V_OCCUPIED`** (waits) |
| 6 | `test_parallel_opposite_straight` | Legs 1, 3 both STRAIGHT | Both advance simultaneously |
| 7 | `test_parallel_left_left` | Legs 1, 3 both LEFT | Both advance simultaneously |
| 8 | `test_mutual_right_deadlock` | Legs 1, 3 both RIGHT | Exactly one advances (stochastic tie-break) |
| 9 | `test_4way_straight_deadlock` | All 4 legs STRAIGHT | Exactly one advances (ring deadlock) |

### T-junction test cases

| # | Test | Setup | Key assertion |
|---|------|-------|--------------|
| 1 | `test_t_init` | `init_t_junction` | n_in=2, n_out=2, perimeter ports correct, route probs sum to 1 |
| 2 | `test_t_merge_yield` | Both → east_out (same destination) | Exactly one advances |
| 3 | `test_t_diverge` | Stem → west_out only | Stem clears; west_out site 1 = `V_OCCUPIED` |
| 4 | `test_t_chord_no_conflict` | west_in→west_out, stem_in→east_out (non-crossing) | Both advance simultaneously |
| 5 | `test_t_chord_conflict` | west_in→east_out, stem_in→west_out (crossing chords) | Exactly one advances; **blocked vehicle = `V_OCCUPIED`** (waits) |
| 6 | `test_t_mass_conservation` | α=0.4, β=0.5, 2000 steps | `count(t+1) == count(t) + entries − exits` at every step |

---

## Build Reference

```bash
make clean                   # remove build/ and *.mod files
make all                     # build all six executables

make run_test_types
make run_test_init_crossroad
make run_test_junction
make run_test_network_run
make run_test_t_junction
```

---

## What Remains

- [ ] **NetCDF output** — write per-step density/flow to `.nc` using the existing `tasep_io` pattern
- [ ] **Python front end** — analysis scripts, visualisation, NetCDF reader
- [ ] **Sensitivity analysis** — sweep α, β, p_left, p_right; fundamental diagrams, phase boundaries
- [ ] **2D multi-lane model** — Week 4 milestone (SDP)
- [ ] **PDE continuum model** — Week 4/5 milestone (SDP)
- [ ] **Multi-lane crossroad** — architecture supports it (set `n_in > 4`, populate `in_perim`/`out_perim`); needs an initialiser and tests

---

## Key Design Decisions

| Decision | Rationale |
|----------|-----------|
| Lane-as-primitive, road-as-grouping | Enables variable lane counts per arm without restructuring junction logic |
| Routing on the junction, not the cell | Cells are just `V_EMPTY`/`V_OCCUPIED`; turning intent doesn't travel down the lane — only matters at the junction |
| Flat `(road, lane)` leg lists | Any topology expressible without special-casing; T-junction and crossroad share identical evaluation code |
| Chord-crossing perimeter predicate | Generalises conflict detection to arbitrary `(n_in, n_out)` without enumerating turn combinations |
| `p_left`/`p_right` preserved in `init_crossroad` signature | Backward-compatible with all existing callers; converted to `in_routes` internally |
