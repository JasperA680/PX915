module simulation
    use tasep_model
    implicit none
    private

    public :: run_simulation
    public :: measure_steady_state

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

    subroutine measure_steady_state(L, n_burnin, n_measure, alpha, beta, mean_density, mean_current)
        ! Run TASEP and return mean density and current over the measurement window.
        ! Density is measured on the central L/4 slice to avoid boundary-layer bias.
        integer, intent(in)  :: L, n_burnin, n_measure
        real,    intent(in)  :: alpha, beta
        real,    intent(out) :: mean_density, mean_current

        integer :: state(L)
        integer :: step, exit_count
        integer :: bulk_lo, bulk_hi, bulk_L
        real    :: density_acc
        integer :: current_acc

        ! Central L/4 slice (matches Python: slice(3*L//8, 5*L//8))
        bulk_lo = 3*L/8 + 1
        bulk_hi = 5*L/8
        bulk_L  = bulk_hi - bulk_lo + 1

        call initialise_lattice(state, L)

        do step = 1, n_burnin
            call tasep_step(state, L, alpha, beta, exit_count)
        end do

        density_acc = 0.0
        current_acc = 0
        do step = 1, n_measure
            call tasep_step(state, L, alpha, beta, exit_count)
            density_acc = density_acc + real(sum(state(bulk_lo:bulk_hi))) / real(bulk_L)
            current_acc = current_acc + exit_count
        end do

        mean_density = density_acc / real(n_measure)
        mean_current = real(current_acc) / real(n_measure)

    end subroutine measure_steady_state

end module simulation

