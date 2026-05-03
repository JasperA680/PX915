program test_network_run
    !--------------------------------------------------------------------
    ! Integration soak test: validates mass conservation at every step.
    !
    !   count(t+1) = count(t) + entries(t) - exits(t)
    !
    ! where entries = new vehicles at inbound site 1 (open alpha boundary)
    !       exits   = vehicles that leave at outbound site L (open beta boundary)
    !
    ! Any per-step violation stops the test with a diagnostic.
    !
    ! A steady-state symmetry check is also run: all four road legs should
    ! carry a similar outbound current after many steps.
    !--------------------------------------------------------------------
    use vehicle_mod
    use road_network_mod
    use network_init_mod
    use junction_mod
    use network_simulation_mod
    implicit none

    integer, parameter :: L       = 20
    integer, parameter :: N_STEPS = 5000
    integer, parameter :: BURNIN  = 500
    real,    parameter :: ALPHA   = 0.4
    real,    parameter :: BETA    = 0.5
    real,    parameter :: P_LEFT  = 0.25
    real,    parameter :: P_RIGHT = 0.25
    real,    parameter :: SYM_TOL = 0.20

    type(road_network_t) :: net
    integer :: cum_entries(4), cum_exits(4)
    integer :: n_before, n_after, step_entries, step_exits, r
    integer :: step, seed_size
    integer, allocatable :: seed(:)
    real    :: mean_exit(4)
    real    :: max_j, min_j

    call random_seed(size = seed_size)
    allocate(seed(seed_size))
    seed = 99999
    call random_seed(put = seed)
    deallocate(seed)

    call init_crossroad(net, L, ALPHA, BETA, P_LEFT, P_RIGHT)

    do step = 1, BURNIN
        call network_step(net)
    end do

    ! n_initial for production run is the network count AFTER burnin.
    n_before      = count_occupied_network(net)
    cum_entries   = 0
    cum_exits     = 0

    !---- Step-by-step mass conservation --------------------------------
    do step = 1, N_STEPS
        call network_step(net)

        step_entries = 0
        step_exits   = 0
        do r = 1, 4
            ! Entry: new vehicle appeared at inbound lane (lane 2) site 1.
            if (net%roads(r)%lane(2)%cells(1) /= V_EMPTY .and. &
                net%roads(r)%lane(2)%old(1)   == V_EMPTY) then
                step_entries        = step_entries + 1
                cum_entries(r) = cum_entries(r) + 1
            end if
            ! Exit: vehicle disappeared from outbound lane (lane 1) site L.
            if (net%roads(r)%lane(1)%cells(net%roads(r)%lane(1)%length) == V_EMPTY .and. &
                net%roads(r)%lane(1)%old(net%roads(r)%lane(1)%length)   /= V_EMPTY) then
                step_exits        = step_exits + 1
                cum_exits(r) = cum_exits(r) + 1
            end if
        end do

        n_after = count_occupied_network(net)

        if (n_after /= n_before + step_entries - step_exits) then
            print '(a,i6)', 'FAIL: mass conservation violated at step ', step
            print '(a,i6)', '  n_before     = ', n_before
            print '(a,i6)', '  step_entries = ', step_entries
            print '(a,i6)', '  step_exits   = ', step_exits
            print '(a,i6)', '  expected n   = ', n_before + step_entries - step_exits
            print '(a,i6)', '  actual n     = ', n_after
            stop 1
        end if

        n_before = n_after
    end do

    print *, 'PASS: mass conservation OK (all ', N_STEPS, ' steps)'

    !---- Symmetry check ------------------------------------------------
    print '(a)', ''
    print '(a)', '--- Per-road outbound currents ---'
    print '(a)', ' Road  Entries  Exits  J_in   J_out'
    do r = 1, 4
        mean_exit(r) = real(cum_exits(r)) / real(N_STEPS)
        print '(i5, 2i8, 2f8.4)', r, cum_entries(r), cum_exits(r), &
              real(cum_entries(r))/real(N_STEPS), mean_exit(r)
    end do

    max_j = maxval(mean_exit)
    min_j = minval(mean_exit)
    print '(a,f7.4,a,f7.4)', 'Exit current range: ', min_j, ' - ', max_j

    if (max_j < 1e-4) then
        print *, 'FAIL: no outflow detected'
        stop 1
    end if
    if ((max_j - min_j) / max_j > SYM_TOL) then
        print '(a,f6.3,a)', 'FAIL: asymmetry > ', SYM_TOL*100, '%'
        stop 1
    end if
    print *, 'PASS: per-road currents are non-zero and symmetric'

    call free_network(net)
    print *, ''
    print *, 'test_network_run: ALL OK'

end program test_network_run
