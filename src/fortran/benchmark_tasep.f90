program benchmark_tasep
    ! Wallclock sweep of the serial 1D TASEP across (L, n_steps).
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

    ! Memory cap on history(L, n_steps) — skip cells exceeding ~500 MB of int4 storage.
    integer(kind=8), parameter :: max_sites = 125000000_8

    real    :: wallclock(n_L, n_N, n_rep)
    real    :: steady_density(n_L)
    real    :: steady_current(n_L)
    integer :: L_arr(n_L), N_arr(n_N)

    integer :: iL, iN, rep, L, n_steps, total_exits
    integer :: n_burnin, n_measure
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

    print '(A)', '=== Steady-state baselines (alpha = beta = 0.5) ==='
    do iL = 1, n_L
        L = L_arr(iL)
        n_burnin  = min(max(2*L, 1000),  20000)
        n_measure = min(max(5*L, 5000),  50000)
        call measure_steady_state(L, n_burnin, n_measure, alpha, beta, &
                                  steady_density(iL), steady_current(iL))
        print '(A,I6,A,F8.5,A,F8.5)', '  L=', L, &
            '  rho=', steady_density(iL), '  J=', steady_current(iL)
    end do

    print '(/,A)', '=== Wallclock sweep ==='
    do iL = 1, n_L
        L = L_arr(iL)
        do iN = 1, n_N
            n_steps    = N_arr(iN)
            site_count = int(L, kind=8) * int(n_steps, kind=8)
            if (site_count > max_sites) then
                print '(A,I6,A,I8,A,I0,A)', '  L=', L, '  N=', n_steps, &
                    '  SKIPPED (history would be ', site_count, ' ints)'
                cycle
            end if

            allocate(history(L, n_steps))
            allocate(density_history(n_steps))
            allocate(current_history(n_steps))
            do rep = 1, n_rep
                call system_clock(t0)
                call run_simulation(L, n_steps, alpha, beta, &
                                    history, density_history, current_history, total_exits)
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
