program test_simulation
    use simulation
    use tasep_io
    implicit none

    ! Select model: 'TASEP' or 'NS'
    character(len=10), parameter :: MODEL = 'TASEP'

    integer, parameter :: L       = 10
    integer, parameter :: n_steps = 100
    real,    parameter :: alpha   = 0.5, beta = 0.5

    integer :: history(L, n_steps)
    real    :: density_history(n_steps)
    integer :: current_history(n_steps)
    integer :: total_exits
    integer :: step

    call run_simulation(L, n_steps, alpha, beta, &
                        history, density_history, current_history, total_exits, &
                        MODEL)

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

end program test_simulation
