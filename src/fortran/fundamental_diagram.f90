program fundamental_diagram_sweep
    use simulation
    use tasep_io
    implicit none

    integer :: L, n_points, n_measure, n_burnin
    integer :: i, nargs
    real    :: alpha, beta
    character(len=32)  :: arg
    character(len=256) :: outfile

    real, allocatable :: alpha_param(:), alpha_rho(:), alpha_J(:)
    real, allocatable :: beta_param(:),  beta_rho(:),  beta_J(:)

    ! Defaults
    L         = 100
    n_points  = 30
    n_measure = 3000
    outfile   = 'data/output/fundamental_diagram.nc'

    nargs = command_argument_count()
    if (nargs >= 1) then; call get_command_argument(1, arg); read(arg, *) L;         end if
    if (nargs >= 2) then; call get_command_argument(2, arg); read(arg, *) n_points;  end if
    if (nargs >= 3) then; call get_command_argument(3, arg); read(arg, *) n_measure; end if

    n_burnin = 2 * L * L

    print '(A)',      'Fundamental diagram sweep'
    print '(A,I6)',   '  L         =', L
    print '(A,I6)',   '  n_points  =', n_points
    print '(A,I8)',   '  n_burnin  =', n_burnin
    print '(A,I6)',   '  n_measure =', n_measure
    print '(A,A)',    '  output    = ', trim(outfile)

    allocate(alpha_param(n_points), alpha_rho(n_points), alpha_J(n_points))
    allocate(beta_param(n_points),  beta_rho(n_points),  beta_J(n_points))

    ! Alpha sweep: vary alpha in [0.02, 0.98] with beta=0.5 fixed.
    ! Traces the low-density branch and max-current peak.
    beta = 0.5
    do i = 1, n_points
        alpha = 0.02 + (0.96 / real(n_points - 1)) * real(i - 1)
        alpha_param(i) = alpha
        call measure_steady_state(L, n_burnin, n_measure, alpha, beta, alpha_rho(i), alpha_J(i))
        if (mod(i, max(1, n_points/5)) == 0) &
            print '(A,I3,A,I3)', '  alpha sweep: ', i, '/', n_points
    end do

    ! Beta sweep: vary beta in [0.02, 0.98] with alpha=0.5 fixed.
    ! Traces the high-density branch and max-current peak.
    alpha = 0.5
    do i = 1, n_points
        beta = 0.02 + (0.96 / real(n_points - 1)) * real(i - 1)
        beta_param(i) = beta
        call measure_steady_state(L, n_burnin, n_measure, alpha, beta, beta_rho(i), beta_J(i))
        if (mod(i, max(1, n_points/5)) == 0) &
            print '(A,I3,A,I3)', '  beta sweep:  ', i, '/', n_points
    end do

    call write_fundamental_diagram_netcdf(trim(outfile), L, n_points, n_burnin, n_measure, &
                                          0.5, 0.5, &
                                          alpha_param, alpha_rho, alpha_J, &
                                          beta_param,  beta_rho,  beta_J)

    print '(A,A)', '  Done. Results written to ', trim(outfile)

    deallocate(alpha_param, alpha_rho, alpha_J, beta_param, beta_rho, beta_J)

end program fundamental_diagram_sweep
