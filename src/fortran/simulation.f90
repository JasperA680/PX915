module simulation
    use road_network_mod
    use tasep_model
    implicit none
    private

    integer, parameter :: V_MAX_NS  = 5
    real,    parameter :: P_SLOW_NS = 0.2

    public :: run_simulation
    public :: measure_steady_state

contains

    subroutine run_simulation(L, n_steps, alpha, beta, &
                              history, density_history, current_history, total_exits, &
                              model)
        ! Run a 1D open-boundary simulation for n_steps.
        !
        ! Optional model argument selects the update rule:
        !   'TASEP' (default) -- single-cell parallel hop, open boundaries
        !   'NS'              -- Nagel-Schreckenberg with V_MAX=5, P_SLOW=0.2

        integer,          intent(in)  :: L, n_steps
        real,             intent(in)  :: alpha, beta
        integer,          intent(out) :: history(L, n_steps)
        real,             intent(out) :: density_history(n_steps)
        integer,          intent(out) :: current_history(n_steps)
        integer,          intent(out) :: total_exits
        character(len=*), intent(in), optional :: model

        character(len=10) :: use_model
        type(cell) :: state(L)
        integer :: step, exit_count

        use_model = 'TASEP'
        if (present(model)) use_model = trim(model)

        call initialise_lattice(state, L)
        total_exits = 0

        do step = 1, n_steps
            select case (trim(use_model))
            case ('NS')
                call ns_single_step(state, L, alpha, beta, exit_count)
            case default
                call tasep_step(state, L, alpha, beta, exit_count)
            end select

            history(:, step)    = merge(1, 0, state%has_car)
            density_history(step) = compute_density(state, L)
            current_history(step) = exit_count
            total_exits           = total_exits + exit_count
        end do

    end subroutine run_simulation

    subroutine measure_steady_state(L, n_burnin, n_measure, alpha, beta, mean_density, mean_current)
        ! Run TASEP and return mean density and current over the measurement window.
        ! Density is measured on the central L/4 slice to avoid boundary-layer bias.
        integer, intent(in)  :: L, n_burnin, n_measure
        real,    intent(in)  :: alpha, beta
        real,    intent(out) :: mean_density, mean_current

        type(cell) :: state(L)
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
            density_acc = density_acc + compute_density(state(bulk_lo:bulk_hi), bulk_L)
            current_acc = current_acc + exit_count
        end do

        mean_density = density_acc / real(n_measure)
        mean_current = real(current_acc) / real(n_measure)

    end subroutine measure_steady_state

    !-----------------------------------------------------------------
    ! One Nagel-Schreckenberg step on a standalone 1D open-boundary road.
    !
    ! Parallel update: all movements determined from old state.
    !   - Accelerate: v <- min(v+1, V_MAX)
    !   - Decelerate for gap: v <- min(v, gap)
    !   - Randomise: v <- v-1 with probability P_SLOW (if v > 0)
    !   - Move: new_pos = old_pos + v
    ! Open boundaries: exit at site L with beta, entry at site 1 with alpha
    ! (only if site 1 was empty at the start of the step).
    !-----------------------------------------------------------------
    subroutine ns_single_step(state, L, alpha, beta, exit_count)
        type(cell), intent(inout) :: state(L)
        integer,    intent(in)    :: L
        real,       intent(in)    :: alpha, beta
        integer,    intent(out)   :: exit_count

        type(cell) :: old_state(L), new_state(L)
        integer :: i, j, gap, v, new_pos
        real    :: rnd

        old_state  = state
        new_state  = old_state
        exit_count = 0

        ! Clear positions occupied at step start.
        do i = 1, L
            if (old_state(i)%has_car) then
                new_state(i)%has_car  = .false.
                new_state(i)%velocity = 0
            end if
        end do

        do i = 1, L
            if (.not. old_state(i)%has_car) cycle

            ! Exit at right boundary with probability beta.
            if (i == L) then
                call random_number(rnd)
                if (rnd < beta) then
                    exit_count = exit_count + 1
                else
                    new_state(L)%has_car  = .true.
                    new_state(L)%velocity = 0
                end if
                cycle
            end if

            ! Gap to next car (or to end of road).
            gap = L - i
            do j = i + 1, L
                if (old_state(j)%has_car) then
                    gap = j - i - 1
                    exit
                end if
            end do

            v = min(old_state(i)%velocity + 1, V_MAX_NS)
            v = min(v, gap)
            call random_number(rnd)
            if (rnd < P_SLOW_NS .and. v > 0) v = v - 1

            new_pos = i + v
            new_state(new_pos)%has_car  = .true.
            new_state(new_pos)%velocity = v
        end do

        ! Entry at site 1 only if it was empty at step start.
        if (.not. old_state(1)%has_car) then
            call random_number(rnd)
            if (rnd < alpha) then
                new_state(1)%has_car  = .true.
                new_state(1)%velocity = 0
            end if
        end if

        state = new_state
    end subroutine ns_single_step

end module simulation
