module simulation 
    use tasep_model 
    implicit none 
    private 

    public :: run_simulation 

contains 

    subroutine run_simulation(L, n_steps, alpha, beta, history, density_history, current_history, total_exits)
        ! Run the 1D open boundary TASEP simulation for n_steps.
        ! 
        ! Inputs:
        !   L                - length of road
        !   n_steps          - number of time steps
        !   alpha            - entry probability at site 1
        !   beta             - exit probability at site L
        !
        ! Outputs:
        !   history          - lattice state at each timestep 
        !                      history(i, t) = state of site i at time t
        !   density_history  - density at each timestep
        !   current_history  - number of particles that exited at each timestep
        !   total_exits      - total number of particles that exited during the simulation

        integer, intent(in) :: L, n_steps
        real, intent(in) :: alpha, beta

        integer, intent(out) :: history(L, n_steps)
        real, intent(out) :: density_history(n_steps)
        integer, intent(out) :: current_history(n_steps)
        integer, intent(out) :: total_exits

        integer :: state(L)
        integer :: step, exit_count

        ! start with an empty lattice
        call initialise_lattice(state, L)
        total_exits = 0

        do step = 1, n_steps
            call tasep_step(state, L, alpha, beta, exit_count)

            history(:, step) = state
            density_history(step) = compute_density(state, L)
            current_history(step) = exit_count
            total_exits = total_exits + exit_count
        end do

    end subroutine run_simulation

end module simulation

