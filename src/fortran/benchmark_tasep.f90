program benchmark_tasep
    ! Wallclock sweep of the serial 1D TASEP across (L, n_steps).
    ! Uses run_simulation_thin so all cells run regardless of history size.
    ! Produces data/output/benchmark.nc and prints a per-run summary.
    use simulation
    use tasep_io
    implicit none

    integer, parameter :: n_L     = 7
    integer, parameter :: n_N     = 4
    integer, parameter :: n_rep   = 3
    integer, parameter :: L_values(n_L) = [50, 100, 200, 500, 1000, 2000, 5000]
    integer, parameter :: N_values(n_N) = [1000, 10000, 100000, 1000000]
    real,    parameter :: alpha   = 0.5
    real,    parameter :: beta    = 0.5

    ! Target history memory: ~50 MB = 12.5 M int4 entries.
    integer(kind=8), parameter :: max_sites = 12500000_8

    real    :: wallclock(n_L, n_N, n_rep)
    real    :: steady_density(n_L)
    real    :: steady_current(n_L)
    integer :: L_arr(n_L), N_arr(n_N)

    integer :: iL, iN, rep, L, n_steps, total_exits
    integer :: n_burnin, n_measure, record_every, n_recorded
    integer(kind=8) :: site_count
    integer, allocatable :: history(:,:)
    real,    allocatable :: density_history(:)
    integer, allocatable :: current_history(:)
    integer :: t0, t1, count_rate
    real    :: elapsed

    L_arr = L_values
    N_arr = N_values
    wallclock = -1.0

    call system_clock(count_rate = count_rate)

    ! --- Steady-state baselines ---
    ! n_burnin = min(L², 1 M) so the system is well-equilibrated before measuring.
    ! n_measure = max(10·L, 20000) for good statistics.
    print '(A)', '=== Steady-state baselines (alpha = beta = 0.5) ==='
    do iL = 1, n_L
        L         = L_arr(iL)
        n_burnin  = min(L*L, 1000000)
        n_measure = max(10*L, 20000)
        call measure_steady_state(L, n_burnin, n_measure, alpha, beta, &
                                  steady_density(iL), steady_current(iL))
        print '(A,I6,A,I8,A,I8,A,F8.5,A,F8.5)', &
            '  L=', L, '  burnin=', n_burnin, '  measure=', n_measure, &
            '  rho=', steady_density(iL), '  J=', steady_current(iL)
    end do

    ! --- Wallclock sweep ---
    ! record_every is chosen so history(L, n_recorded) stays within max_sites.
    ! record_every = 1 means full history (no thinning).
    print '(/,A)', '=== Wallclock sweep ==='
    do iL = 1, n_L
        L = L_arr(iL)
        do iN = 1, n_N
            n_steps    = N_arr(iN)
            site_count = int(L, kind=8) * int(n_steps, kind=8)
            record_every = max(1, int((site_count + max_sites - 1_8) / max_sites))
            n_recorded   = (n_steps + record_every - 1) / record_every

            if (record_every > 1) then
                print '(A,I6,A,I8,A,I4,A,I8,A)', &
                    '  L=', L, '  N=', n_steps, '  (thinning every ', record_every, &
                    ' steps -> history size ', n_recorded, ')'
            end if

            allocate(history(L, n_recorded))
            allocate(density_history(n_steps))
            allocate(current_history(n_steps))
            do rep = 1, n_rep
                call system_clock(t0)
                call run_simulation_thin(L, n_steps, record_every, n_recorded, &
                                         alpha, beta, history, density_history, &
                                         current_history, total_exits)
                call system_clock(t1)
                elapsed = real(t1 - t0) / real(count_rate)
                wallclock(iL, iN, rep) = elapsed
                print '(A,I6,A,I8,A,I2,A,F10.4,A)', '  L=', L, '  N=', n_steps, &
                    '  rep=', rep, '  t=', elapsed, ' s'
            end do
            deallocate(history, density_history, current_history)
        end do
    end do

    call write_benchmark_netcdf('data/output/benchmark.nc', n_L, n_N, n_rep, &
                                alpha, beta, L_arr, N_arr, &
                                wallclock, steady_density, steady_current)

    print '(/,A)', 'Wrote data/output/benchmark.nc'

end program benchmark_tasep
