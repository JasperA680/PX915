program test_simulation
    use simulation
    use tasep_io
    use tasep_model, only: rng_init
    implicit none

    integer :: L, n_steps, nargs
    real    :: alpha, beta

    integer, allocatable :: history(:,:)
    real,    allocatable :: density_history(:)
    integer, allocatable :: current_history(:)
    integer :: total_exits, step
    character(len=32) :: arg

    ! Defaults
    L       = 10
    n_steps = 100
    alpha   = 0.5
    beta    = 0.5

    nargs = command_argument_count()
    if (nargs >= 1) then; call get_command_argument(1, arg); read(arg, *) L;       end if
    if (nargs >= 2) then; call get_command_argument(2, arg); read(arg, *) n_steps; end if
    if (nargs >= 3) then; call get_command_argument(3, arg); read(arg, *) alpha;   end if
    if (nargs >= 4) then; call get_command_argument(4, arg); read(arg, *) beta;    end if

    print '(A,I6)',   'L       =', L
    print '(A,I6)',   'n_steps =', n_steps
    print '(A,F6.3)', 'alpha   =', alpha
    print '(A,F6.3)', 'beta    =', beta

    allocate(history(L, n_steps))
    allocate(density_history(n_steps))
    allocate(current_history(n_steps))

    call rng_init(42)
    call run_simulation(L, n_steps, alpha, beta, history, density_history, current_history, total_exits)

    do step = 1, n_steps
        print *, "Step:", step
        print *, "State:", history(:, step)
        print *, "Density:", density_history(step)
        print *, "Current:", current_history(step)
        print *, ""
    end do

    print *, "Total exits:", total_exits
    print *, "Mean current:", real(total_exits)/real(n_steps)

    call write_netcdf('data/output/simulation.nc', L, n_steps, alpha, beta, &
                      history, density_history, current_history)
    print *, "Output written to data/output/simulation.nc"

    deallocate(history, density_history, current_history)

end program test_simulation
