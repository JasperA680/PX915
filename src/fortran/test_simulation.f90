program test_simulation
    use simulation
    implicit none

    integer, parameter :: L = 10
    integer, parameter :: n_steps = 20

    integer :: history(L, n_steps)
    real    :: density_history(n_steps)
    integer :: current_history(n_steps)
    integer :: total_exits

    integer :: step

    call run_simulation(L, n_steps, 0.5, 0.5, history, density_history, current_history, total_exits)

    do step = 1, n_steps
        print *, "Step:", step
        print *, "State:", history(:, step)
        print *, "Density:", density_history(step)
        print *, "Current:", current_history(step)
        print *, ""
    end do

    print *, "Total exits:", total_exits
    print *, "Mean current:", real(total_exits)/real(n_steps)

end program test_simulation