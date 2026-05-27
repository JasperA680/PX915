module fundamental_diagram_mod
    ! Steady-state density / current measurement and NetCDF I/O for the
    ! fundamental-diagram sweep driver.
    !
    ! Both measurement routines drive a minimal one-lane ``road_network_t``
    ! through the shared network primitives (``snapshot_network`` +
    ! ``tasep_lane_step`` for TASEP, ``NS_model_step`` for NS): there is no
    ! separate per-chain update kernel. The module also bundles the NetCDF
    ! writer used by the sweep program in the second program unit of this
    ! file (``program fd_sweep``).
    !
    ! The TASEP routine matches the convention used by the Python
    ! ``analysis.fundamental_diagram`` it replaces: density is averaged only
    ! on the central quarter of the chain (sites ``3L/8 + 1`` through
    ! ``5L/8``) to suppress boundary-layer bias near the open boundaries,
    ! and the current is the time-averaged number of right-boundary exits
    ! detected as ``old(L)%has_car .and. .not. cells(L)%has_car`` (the same
    ! pattern ``run_network`` uses for ``road_exits``).
    !
    ! The NS routine drives a periodic ring at fixed vehicle count. The
    ! density is constant by construction (no inflow / outflow), so the
    ! measured value is just ``n_vehicles / L``; the current is the
    ! per-cell time-averaged sum of vehicle velocities, equivalent to the
    ! flow definition used by the Python NS fundamental-diagram code.
    !
    ! The NetCDF writer stores a flat list of ``(rho, J)`` points. For
    ! TASEP the list is the concatenation of the alpha-branch and
    ! beta-branch sweeps (length ``2 * n_points``). For NS it is a single
    ! density sweep (length ``n_points``). Sweep metadata (model, L,
    ! burnin, v_max, p_slow) is stored as global attributes so the
    ! plotter can label the figure without consulting the original CLI
    ! arguments.

    use netcdf
    use road_network_mod
    use network_init_mod, only: init_lane, free_network, place_evenly
    use tasep_model,      only: tasep_lane_step
    use NS_model,         only: NS_model_step
    use network_io_mod,   only: check

    implicit none
    private

    public :: measure_steady_state_tasep
    public :: measure_steady_state_ns
    public :: write_fd_netcdf
    public :: seed_iter_rng

contains

    subroutine seed_iter_rng(base_seed, i, branch)
        ! Deterministically seed the current thread's RNG state for one
        ! independent sweep iteration. Called at the top of every iteration
        ! of the OpenMP-parallelised sweep loops so that the output is
        ! identical regardless of how the loop iterations are scheduled
        ! across threads. Within an iteration, ``random_number`` runs on
        ! the calling thread's private RNG state (gfortran's RNG is
        ! per-thread under OpenMP).
        integer, intent(in) :: base_seed, i, branch

        integer :: n_seed, k
        integer, allocatable :: seed_array(:)

        call random_seed(size=n_seed)
        allocate(seed_array(n_seed))
        do k = 1, n_seed
            seed_array(k) = base_seed + 1009 * i + 7919 * branch + 37 * k
        end do
        call random_seed(put=seed_array)
        deallocate(seed_array)
    end subroutine seed_iter_rng



    subroutine measure_steady_state_tasep(L, n_burnin, n_measure, alpha, beta, &
                                          mean_density, mean_current)
        ! Run open-boundary TASEP for ``n_burnin`` steps to equilibrate, then
        ! average density and exit current over the next ``n_measure`` steps.
        integer, intent(in) :: L, n_burnin, n_measure
        real, intent(in) :: alpha, beta
        real, intent(out) :: mean_density, mean_current

        type(road_network_t) :: net
        integer :: step, current_acc
        integer :: bulk_lo, bulk_hi, bulk_L
        real :: density_acc

        bulk_lo = 3 * L / 8 + 1
        bulk_hi = 5 * L / 8
        bulk_L  = bulk_hi - bulk_lo + 1

        call build_chain_network(net, L, alpha, beta, &
                                 open_in=.true., open_out=.true., is_periodic=.false.)

        do step = 1, n_burnin
            call snapshot_network(net)
            call tasep_lane_step(net)
        end do

        density_acc = 0.0
        current_acc = 0
        do step = 1, n_measure
            call snapshot_network(net)
            call tasep_lane_step(net)
            if (net%roads(1)%lane(1)%old(L)%has_car .and. &
                .not. net%roads(1)%lane(1)%cells(L)%has_car) then
                current_acc = current_acc + 1
            end if
            density_acc = density_acc + &
                real(count(net%roads(1)%lane(1)%cells(bulk_lo:bulk_hi)%has_car)) / real(bulk_L)
        end do

        mean_density = density_acc / real(n_measure)
        mean_current = real(current_acc) / real(n_measure)

        call free_network(net)
    end subroutine measure_steady_state_tasep


    subroutine measure_steady_state_ns(L, n_burnin, n_measure, n_vehicles, &
                                       v_max, p_slow, mean_density, mean_current)
        ! Run NS on a length-``L`` periodic ring with ``n_vehicles`` cars
        ! for ``n_burnin`` steps, then average the per-cell flow over
        ! ``n_measure`` steps.  Density is constant; both empty
        ! (``n_vehicles == 0``) and fully-jammed (``n_vehicles == L``)
        ! cases short-circuit to a zero flow.
        integer, intent(in) :: L, n_burnin, n_measure, n_vehicles, v_max
        real, intent(in) :: p_slow
        real, intent(out) :: mean_density, mean_current

        type(road_network_t) :: net
        integer :: step, moves_acc

        mean_density = real(n_vehicles) / real(L)
        if (n_vehicles <= 0 .or. n_vehicles >= L) then
            mean_current = 0.0
            return
        end if

        call build_chain_network(net, L, alpha=0.0, beta=0.0, &
                                 open_in=.false., open_out=.false., is_periodic=.true.)
        call place_evenly(net%roads(1)%lane(1), n_vehicles)

        do step = 1, n_burnin
            call snapshot_network(net)
            call NS_model_step(net, v_max, p_slow)
        end do

        moves_acc = 0
        do step = 1, n_measure
            call snapshot_network(net)
            call NS_model_step(net, v_max, p_slow)
            moves_acc = moves_acc + sum(net%roads(1)%lane(1)%cells%velocity)
        end do

        mean_current = real(moves_acc) / (real(L) * real(n_measure))

        call free_network(net)
    end subroutine measure_steady_state_ns


    subroutine write_fd_netcdf(filename, model, L, n_burnin, n_measure, &
                               v_max, p_slow, rho, J)
        ! Write a fundamental-diagram sweep to NetCDF.
        !
        ! ``rho`` and ``J`` are the per-point density and current arrays
        ! (length ``2 * n_points`` for TASEP, ``n_points`` for NS). Sweep
        ! metadata is written as global attributes.
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: model
        integer, intent(in) :: L, n_burnin, n_measure, v_max
        real, intent(in) :: p_slow
        real, intent(in) :: rho(:), J(:)

        integer :: ncid, dim_pts, vid_rho, vid_J, n

        if (size(rho) /= size(J)) then
            print *, 'write_fd_netcdf: rho and J have different lengths'
            stop 1
        end if
        n = size(rho)

        call check( nf90_create(filename, NF90_CLOBBER, ncid) )

        call check( nf90_def_dim(ncid, 'n_points', n, dim_pts) )

        call check( nf90_put_att(ncid, NF90_GLOBAL, 'model',     trim(model)) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'L',         L) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'n_burnin',  n_burnin) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'n_measure', n_measure) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'v_max',     v_max) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'p_slow',    p_slow) )

        call check( nf90_def_var(ncid, 'rho', NF90_FLOAT, [dim_pts], vid_rho) )
        call check( nf90_put_att(ncid, vid_rho, 'long_name', 'mean steady-state density') )

        call check( nf90_def_var(ncid, 'J', NF90_FLOAT, [dim_pts], vid_J) )
        call check( nf90_put_att(ncid, vid_J, 'long_name', 'mean steady-state current') )

        call check( nf90_enddef(ncid) )

        call check( nf90_put_var(ncid, vid_rho, rho) )
        call check( nf90_put_var(ncid, vid_J,   J) )

        call check( nf90_close(ncid) )
    end subroutine write_fd_netcdf


    subroutine build_chain_network(net, L, alpha, beta, open_in, open_out, is_periodic)
        ! Allocate a minimal one-road / one-lane network suitable for FD sweeps.
        ! Junctions are intentionally left unallocated; the TASEP and NS step
        ! routines only iterate over ``net%roads`` so this is safe, and
        ! ``free_network`` checks ``allocated(net%junctions)`` before touching it.
        type(road_network_t), intent(out) :: net
        integer, intent(in) :: L
        real, intent(in) :: alpha, beta
        logical, intent(in) :: open_in, open_out, is_periodic

        allocate(net%roads(1))
        net%roads(1)%id = 1
        net%roads(1)%end_junction = [0, 0]

        allocate(net%roads(1)%lane(1))
        call init_lane(net%roads(1)%lane(1), L, +1)

        net%roads(1)%lane(1)%alpha       = alpha
        net%roads(1)%lane(1)%beta        = beta
        net%roads(1)%lane(1)%open_in     = open_in
        net%roads(1)%lane(1)%open_out    = open_out
        net%roads(1)%lane(1)%is_periodic = is_periodic
    end subroutine build_chain_network


end module fundamental_diagram_mod


program fd_sweep
    ! Fundamental-diagram sweep driver.
    !
    ! Replaces the pure-Python ``analysis.fundamental_diagram`` and
    ! ``analysis.fundamental_diagram_ns`` routines.  Runs entirely in
    ! Fortran and writes a NetCDF file the Python plotter consumes.
    !
    ! Usage:
    !
    !     fd_sweep MODEL L N_POINTS N_STEPS V_MAX P_SLOW OUTPUT_PATH [SEED]
    !
    ! where MODEL is ``TASEP`` or ``NS``.  V_MAX and P_SLOW are ignored
    ! for TASEP but must be supplied for argument-position consistency
    ! with the NS branch.  SEED is optional; if omitted the Fortran RNG
    ! default seeding is used.
    !
    ! TASEP sweep:
    !   Alpha branch (beta = 1): alpha from 0.02 to 0.98 in N_POINTS steps.
    !   Beta  branch (alpha = 1): beta  from 0.02 to 0.98 in N_POINTS steps.
    !   Output array has 2 * N_POINTS rows.
    !
    !   The fixed branch parameter is held at 1 (deterministic boundary) so
    !   the only stochasticity comes from the swept parameter and the bulk
    !   sees an effectively noise-free reservoir on one side. With p = 1
    !   parallel-update bulk this recovers the textbook
    !   ``J(rho) = min(rho, 1 - rho)`` fundamental diagram: the alpha branch
    !   traces the LD line ``J = rho`` from (0, 0) to (0.5, 0.5), the beta
    !   branch traces the HD line ``J = 1 - rho`` from (1, 0) to (0.5, 0.5),
    !   and the two meet at the max-current point. Holding the boundary at
    !   0.5 puts the swept-parameter > 0.5 half of each branch on the
    !   LD/HD-MC phase-boundary line where the steady-state density is not
    !   ``rho = 0.5`` and the diagram develops a visible gap at large L.
    !
    ! NS sweep:
    !   Density branch: target density from 0.02 to 0.98 in N_POINTS steps,
    !   round to integer vehicle counts.  Output array has N_POINTS rows.

    use fundamental_diagram_mod

    implicit none

    character(len=32)  :: arg, model
    character(len=256) :: outfile
    integer :: L, n_points, n_steps, n_burnin, v_max
    integer :: i, n_total, n_vehicles, seed_val, nargs
    real    :: p_slow, alpha, beta, target_rho
    real    :: mean_rho, mean_J
    real, allocatable :: rho(:), J(:)

    nargs = command_argument_count()
    if (nargs < 7) then
        print '(A)', 'Usage: fd_sweep MODEL L N_POINTS N_STEPS V_MAX P_SLOW OUTPUT_PATH [SEED]'
        stop 1
    end if

    call get_command_argument(1, model)
    call get_command_argument(2, arg);  read(arg, *) L
    call get_command_argument(3, arg);  read(arg, *) n_points
    call get_command_argument(4, arg);  read(arg, *) n_steps
    call get_command_argument(5, arg);  read(arg, *) v_max
    call get_command_argument(6, arg);  read(arg, *) p_slow
    call get_command_argument(7, outfile)

    if (nargs >= 8) then
        call get_command_argument(8, arg);  read(arg, *) seed_val
    else
        ! Wall-clock fallback so an un-seeded sweep still varies between runs;
        ! per-iteration seeding below makes the output bit-reproducible for any
        ! fixed seed_val regardless of OpenMP thread count.
        call system_clock(count=seed_val)
    end if

    select case (trim(model))
    case ('TASEP')
        ! Burn-in scales like the autocorrelation time near the max-current
        ! line, ~ L^2. L^2 / 2 is enough headroom well inside LD/HD (the
        ! Python version used 2 L^2 which is safe but ~4x slower than needed).
        n_burnin = L * L / 2
        n_total = 2 * n_points
        allocate(rho(n_total), J(n_total))

        print '(A)',      'TASEP fundamental-diagram sweep'
        print '(A,I6)',   '  L         =', L
        print '(A,I6)',   '  n_points  =', n_points
        print '(A,I9)',   '  n_burnin  =', n_burnin
        print '(A,I6)',   '  n_measure =', n_steps

        ! Alpha branch: deterministic exit (beta = 1) so the right boundary
        ! adds no noise. Bulk density tracks alpha/(1 + alpha) and the LD
        ! branch closes off at (0.5, 0.5) when alpha = 1, recovering the
        ! textbook J = min(rho, 1 - rho) shape.
        !
        ! Each iteration is an independent steady-state measurement, so the
        ! sweep parallelises trivially. Per-iteration RNG seeding makes the
        ! output bit-identical to the serial run regardless of how the loop
        ! is scheduled.
        beta = 1.0
        !$omp parallel do default(shared) private(i, alpha, mean_rho, mean_J) &
        !$omp& schedule(dynamic)
        do i = 1, n_points
            call seed_iter_rng(seed_val, i, 1)
            alpha = 0.02 + (0.96 / real(n_points - 1)) * real(i - 1)
            call measure_steady_state_tasep(L, n_burnin, n_steps, alpha, beta, &
                                            mean_rho, mean_J)
            rho(i) = mean_rho
            J(i)   = mean_J
        end do
        !$omp end parallel do

        ! Beta branch: symmetric — deterministic injection (alpha = 1) so
        ! the left boundary adds no noise, tracing the HD line back from
        ! rho ~ 1 to the same (0.5, 0.5) max-current point.
        alpha = 1.0
        !$omp parallel do default(shared) private(i, beta, mean_rho, mean_J) &
        !$omp& schedule(dynamic)
        do i = 1, n_points
            call seed_iter_rng(seed_val, i, 2)
            beta = 0.02 + (0.96 / real(n_points - 1)) * real(i - 1)
            call measure_steady_state_tasep(L, n_burnin, n_steps, alpha, beta, &
                                            mean_rho, mean_J)
            rho(n_points + i) = mean_rho
            J(n_points + i)   = mean_J
        end do
        !$omp end parallel do

    case ('NS')
        ! NS relaxes much faster than TASEP near criticality; the Python
        ! version used max(500, 5 L) so we match.
        n_burnin = max(500, 5 * L)
        n_total = n_points
        allocate(rho(n_total), J(n_total))

        print '(A)',      'NS fundamental-diagram sweep (periodic ring)'
        print '(A,I6)',   '  L         =', L
        print '(A,I6)',   '  n_points  =', n_points
        print '(A,I6)',   '  v_max     =', v_max
        print '(A,F6.3)', '  p_slow    =', p_slow
        print '(A,I9)',   '  n_burnin  =', n_burnin
        print '(A,I6)',   '  n_measure =', n_steps

        !$omp parallel do default(shared) &
        !$omp& private(i, target_rho, n_vehicles, mean_rho, mean_J) &
        !$omp& schedule(dynamic)
        do i = 1, n_points
            call seed_iter_rng(seed_val, i, 3)
            target_rho = 0.02 + (0.96 / real(n_points - 1)) * real(i - 1)
            n_vehicles = max(1, nint(target_rho * real(L)))
            call measure_steady_state_ns(L, n_burnin, n_steps, n_vehicles, &
                                         v_max, p_slow, mean_rho, mean_J)
            rho(i) = mean_rho
            J(i)   = mean_J
        end do
        !$omp end parallel do

    case default
        print '(A,A)', 'fd_sweep: unknown model ', trim(model)
        stop 1
    end select

    call write_fd_netcdf(trim(outfile), trim(model), L, n_burnin, n_steps, &
                         v_max, p_slow, rho, J)
    print '(A,A)', '  done. wrote ', trim(outfile)
    deallocate(rho, J)

end program fd_sweep
