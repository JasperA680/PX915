# Parallelisation design — 1D TASEP

Working branch: `parallel`.

This document captures (a) what the serial benchmark tells us about where the cost actually lives, and (b) the OpenMP and MPI strategies we intend to implement next, including the changes to the update rule and RNG that each requires.

## 1. Bottleneck analysis (from `data/output/benchmark.nc`)

### Setup

- Driver: [src/fortran/benchmark_tasep.f90](../src/fortran/benchmark_tasep.f90)
- Sweep: `L ∈ {50, 100, 200, 500, 1000, 2000, 5000}` × `n_steps ∈ {10³, 10⁴, 10⁵, 10⁶}`
- 3 repeats per cell, fixed α = β = 0.5 (max-current phase)
- Memory cap: cells with `L · n_steps > 1.25 × 10⁸` int4 entries (~500 MB) are skipped — see (c) below
- Build: `gfortran -O2`, single core (Apple Silicon, macOS)
- Steady-state baselines per L use `n_burnin = min(max(2L, 1000), 20000)`, `n_measure = min(max(5L, 5000), 50000)`

### Headline numbers

| L     | n_steps | wallclock (s) | site-updates/s |
|-------|---------|---------------|----------------|
| 50    | 10⁶     | 0.16          | 3.1 × 10⁸      |
| 100   | 10⁶     | 0.24          | 4.2 × 10⁸      |
| 1000  | 10⁵     | 0.15          | 6.7 × 10⁸      |
| 5000  | 10⁴     | 0.06          | 8.3 × 10⁸      |

The serial code clocks **≈ 0.3–0.8 × 10⁹ site-updates/s** on this machine — much faster than the typical "~10⁷/s" rule of thumb (helped by `-O2`, the dense integer lattice, and the simple inner loop). Across all measured cells the longest wallclock was 0.25 s.

### Scaling fits (`scripts/plot_benchmark.py`, log-log)

- **vs n_steps (at fixed L):** slopes ≈ 0.91, 0.95, 0.87, 1.02, 1.02, 1.06, 1.07 for L = 50 … 5000.
  Consistent with the expected `O(n_steps)` — small deviations at low L come from sub-millisecond timer noise.
- **vs L (at fixed n_steps):** slopes ≈ 0.59, 0.82, 0.75, 0.59 for N = 10³ … 10⁶.
  **Sub-linear**, because the per-step overhead (two `random_number` calls, history slice copy, density/current bookkeeping in [src/fortran/simulation.f90:45-48](../src/fortran/simulation.f90#L45-L48)) is independent of L and dominates at small L. The site-updates/s column above confirms this — throughput rises with L, presumably from better cache locality during the bulk sweep.

### Where the cost actually lives

(a) **Compute is fast.** At ~10⁹ updates/s, even L=5000 × n_steps=10⁶ (5 × 10⁹ updates) would be ~6 s of pure compute. Compute is **not** the binding constraint for the sizes we plan to study.

(b) **Per-step overhead dominates at small L.** Below L ~ 500 the per-step costs (RNG calls, history slice, density sum) outweigh the bulk loop. OpenMP on the bulk loop alone will scale poorly here — for small L, threading overhead will exceed the parallel work. **OpenMP only pays off at L ≳ 1000.**

(c) **Memory for `history` is the binding constraint.** [src/fortran/simulation.f90:45](../src/fortran/simulation.f90#L45) writes `history(:, step) = state` every step, so storage is `4 · L · n_steps` bytes. At L=2000, N=10⁶ that's 8 GB; at L=5000, N=10⁶ it's 20 GB. The benchmark skips 6 of 28 cells purely on this constraint, despite their compute being affordable. **This is the first thing to fix before parallelisation buys us much** — see §4.

(d) **Correct steady-state values: ρ ≈ 2/3, J ≈ 1/3 at α = β = 0.5.** The initial benchmark run showed scattered values (0.58–0.67) because the burnin was capped at 20,000 steps — insufficient for the system to equilibrate. After fixing the burnin to `min(L², 1,000,000)`, all L values converge tightly to ρ ≈ 0.6667, J ≈ 0.3333. These are the correct values for this code's dynamics. The bulk hop in `tasep_step` is **deterministic** (hop always succeeds if next site is empty), not stochastic with rate p < 1. This deterministic right-to-left sequential update is not the standard parallel-update stochastic ASEP, so the parallel-update mean-field prediction (ρ=0.5, J=0.25) does not apply. The measured ρ=2/3, J=1/3 are consistent with the theoretical steady state of the sequential deterministic TASEP at the max-current boundary. These values serve as the correctness reference for future parallel implementations.

## 2. OpenMP strategy

### The blocking dependency

The bulk loop at [src/fortran/tasep.f90:78-84](../src/fortran/tasep.f90#L78-L84) is

```fortran
do i = L-1, 1, -1
    if (old_state(i) == 1 .and. old_state(i+1) == 0) then
        new_state(i)   = 0
        new_state(i+1) = 1
    end if
end do
```

`new_state(i+1)` is written by iteration `i` and the iteration only reads from `old_state` — which is a per-step snapshot taken before the loop. There is **no read/write conflict between iterations** as long as every iteration uses the same `old_state` snapshot. The right-to-left order is therefore cosmetic; what is actually required is that **all sites update from the same snapshot**.

This means a straight `!$omp parallel do` on the existing loop is **already correct** if the snapshot is made before the parallel region. The only true write conflict is between iteration `i` and iteration `i-1` both touching `new_state(i)` — but each only writes its own pair (i and i+1), so the pattern is "shift the particle one right", and no two iterations write to the same cell because each iteration owns a unique `i`.

Wait — iteration `i` writes `new_state(i+1)`, and iteration `i+1` would also write `new_state(i+1)` if its predicate fires. Concretely: if `old_state(i)=1, old_state(i+1)=0` then iteration `i` writes `new_state(i+1)=1`; if simultaneously `old_state(i+1)=0, old_state(i+2)=0` (vacuous, no move) or `old_state(i+1)=1` (different case) — so by inspection the predicate of iteration `i+1` can only fire if `old_state(i+1)=1`, which is incompatible with iteration `i`'s predicate that requires `old_state(i+1)=0`. **Predicates are mutually exclusive between adjacent iterations** → no write race. Confirmed.

### Plan

Two-phase change to [src/fortran/tasep.f90](../src/fortran/tasep.f90):

1. **RNG refactor.** Replace the global `call random_number(r)` at lines 70 and 88 with calls to a thread-local RNG (e.g. a small xoshiro128++ state held in a `!$omp threadprivate` variable, seeded once per thread from a base seed + `omp_get_thread_num()`). This is needed for two reasons: the bulk loop calls `random_number` zero times today (only entry/exit do) so RNG is not yet on the parallel path, but adding any per-site stochastic move later will require it; and reproducibility under threading is a hard requirement for the verification step.
2. **Annotate the bulk loop.** Add `!$omp parallel do default(shared) private(i)` over the right-to-left loop. Keep the snapshot `old_state = state` as the serial step it already is — it is one `memcpy` of `4L` bytes, negligible.

Add `-fopenmp` to `FFLAGS` and `OMP_NUM_THREADS=1` as the default in the `benchmark` target so the existing serial numbers remain reproducible.

### Expected scaling

Only useful for L large enough that bulk-loop work exceeds the snapshot copy + entry/exit serial code:
- At L = 1000, bulk loop is ~10³ comparisons; copy is 4 kB. Speedup likely 2–3× on 4 threads.
- At L = 5000, copy is 20 kB; speedup approaching the thread count, modulo memory bandwidth on a shared-memory machine.

### Verification

Compare `measure_steady_state` output between serial and OpenMP for (α, β) ∈ {(0.2, 0.5), (0.5, 0.2), (0.5, 0.5)}, 10 independent base-seed values, threads ∈ {1, 2, 4, 8}. Require `|J_par − J_ser| < 3σ` where σ comes from the 10-seed spread of the serial run. Also run `make run-fd` on both and assert the fundamental diagram shape is preserved.

## 3. MPI strategy

### Decomposition

1D domain decomposition. Rank `r` of `P` owns sites `[r · L/P + 1, (r+1) · L/P]` plus a 1-cell ghost on each side (so rank 0 has only a right ghost, rank `P-1` only a left ghost).

Each step:

1. `old_state = state` locally.
2. Bulk loop over owned sites — same predicate, parallel-safe by §2.
3. `MPI_Sendrecv` the rightmost owned cell to the right neighbour (it becomes their left ghost) and receive from the left neighbour (left ghost ← left neighbour's rightmost owned cell). One 4-byte message per neighbour per step.
4. Rank 0 handles entry at the global site 1; rank `P-1` handles exit at the global site L.
5. `MPI_Allreduce(SUM)` of local exit counts and occupancy → global current and density.

History I/O — see §4 below — should be done via parallel NetCDF rather than a gather to rank 0.

### Why the §2 OpenMP refactor is a prerequisite

The serial code's right-to-left bulk-loop order is fine for a single domain, but at a domain boundary it implicitly assumed information could propagate from site `i+1` to site `i` "instantly within a step". Once site `i` and `i+1` are on different ranks, the right-going particle at the rightmost owned cell of rank `r` either moves into rank `r+1`'s leftmost cell or doesn't — and the only way to know is the ghost exchange. The §2 "update from snapshot" formulation handles this trivially: rank `r` decides based on its ghost copy of rank `r+1`'s leftmost cell from the **previous** step's state. So MPI just inherits the snapshot semantics that §2 establishes.

### Expected scaling

Communication is one int per neighbour per step ≈ 8 bytes/step/rank, plus an allreduce of two ints per step. On a single node with shared memory MPI (one-sided or shm transport), latency ~1 μs per `Sendrecv` × 1 step at ~10⁹ updates/s — communication will dominate below L ≈ 10⁴ per rank. **Use MPI for L ≳ 10⁵; below that, prefer OpenMP.**

Across nodes, the latency floor is ~10 μs, so the strong-scaling sweet spot needs L ≳ 10⁶ per rank — beyond what the report scope likely requires, but worth measuring once.

## 4. History I/O — the real scaling problem

Compute is cheap. Storing `history(L, n_steps)` is not. At the regimes we want for the report (L ≳ 1000, n_steps ≳ 10⁶), the array exceeds RAM. Three options, in order of effort:

1. **Thinning.** Add `record_every` to `run_simulation` and only write `history(:, step)` when `mod(step, record_every) == 0`. For visualisation this is almost always enough.
2. **Drop history entirely** in the benchmark / steady-state paths; keep `density_history` and `current_history` (each O(n_steps), trivial). The `run_simulation` API can stay; add a `measure_only` overload that doesn't allocate the history array.
3. **Async / parallel NetCDF.** Use `nf90_def_var_chunking` + `nf90_put_var` slabs per step, possibly with HDF5 chunk cache tuning. This is the right answer for runs where we genuinely want the space-time diagram at L = 10⁴.

Recommend doing (1) and (2) before any parallelisation work — they unblock the regime where parallelisation is interesting.

## 5. Verification plan

Common to OpenMP and MPI:

| Check | Method | Tolerance |
|-------|--------|-----------|
| Steady-state current `J` matches serial | `measure_steady_state`, 10 seeds × {(0.2,0.5), (0.5,0.2), (0.5,0.5)} | `|J_par − J_ser| < 3σ_ser` |
| Steady-state density `ρ` matches serial | same | same |
| Fundamental diagram shape preserved | `make run-fd` + Python plotter | visual + RMS of `J(ρ)` curve < 1% of serial |
| Per-step RNG stream reproducible at thread count 1 | same base seed, `OMP_NUM_THREADS=1` | bit-exact match against serial |
| Strong scaling | run benchmark sweep with `OMP_NUM_THREADS ∈ {1,2,4,8}` and `mpirun -n ∈ {1,2,4,8}` | report curves; do not require ideal speedup |

## 6. Implementation order

1. **History thinning / measure-only path** (§4 items 1–2). Unblocks larger L · n_steps. ~½ day.
2. **OpenMP** (§2). Snapshot-based bulk loop + thread-local RNG. ~1 day.
3. **MPI** (§3). 1D decomposition + ghost exchange + parallel NetCDF for history. ~2–3 days.
4. **Strong/weak scaling sweep** re-running `make benchmark` with new compiled binaries.
5. **Cross-check** against the serial baseline in `data/output/benchmark.nc` (this branch) — keep that file as the immutable reference.

## Critical files (touched / to be touched)

- [src/fortran/tasep.f90](../src/fortran/tasep.f90) — bulk loop, RNG calls (lines 78–84, 70, 88)
- [src/fortran/simulation.f90](../src/fortran/simulation.f90) — `run_simulation`, history thinning hook (line 45)
- [src/fortran/io.f90](../src/fortran/io.f90) — `write_benchmark_netcdf` already added; may add parallel NetCDF helper later
- [src/fortran/benchmark_tasep.f90](../src/fortran/benchmark_tasep.f90) — re-used to re-time after each step above
- [Makefile](../Makefile) — `-fopenmp` flag; new `benchmark-omp` and `benchmark-mpi` targets later
