module simulation
    ! High-level simulation routines for one-dimensional cellular automaton models.
    !
    ! This module provides driver routines for running open-boundary one-dimensional
    ! traffic simulations. The default model is the TASEP cellular automaton, but
    ! the ``run_simulation`` routine can also use a simple Nagel-Schreckenberg
    ! update rule.
    !
    ! The open-boundary TASEP model uses two boundary probabilities:
    !
    ! .. math::
    !
    !    \alpha
    !
    ! for entry at the left boundary, and
    !
    ! .. math::
    !
    !    \beta
    !
    ! for exit at the right boundary.
    !
    ! The density history is computed as:
    !
    ! .. math::
    !
    !    \rho(t) = \frac{N(t)}{L},
    !
    ! where ``N(t)`` is the number of occupied sites at time step ``t`` and ``L``
    ! is the lattice length.
    !
    ! The instantaneous current is measured from the number of particles exiting
    ! the right boundary during each time step. The time-averaged current over a
    ! measurement window of length ``T`` is:
    !
    ! .. math::
    !
    !    J = \frac{1}{T}\sum_{t=1}^{T} n_{\mathrm{exit}}(t).

    use road_network_mod
    use tasep_model

    implicit none
    private

    integer, parameter :: V_MAX_NS  = 5
    ! Maximum vehicle speed used by the Nagel-Schreckenberg update.
    real, parameter :: P_SLOW_NS = 0.2
    ! Random slowing probability used by the Nagel-Schreckenberg update.

    public :: run_simulation
    public :: measure_steady_state

contains

    subroutine run_simulation(L, n_steps, alpha, beta, &
                              history, density_history, current_history, total_exits, &
                              model)
        ! Run a one-dimensional open-boundary cellular automaton simulation.
        !
        ! This routine initialises an empty one-dimensional lattice and advances it
        ! for ``n_steps`` time steps. At each step, it stores the lattice
        ! occupation, mean density, instantaneous exit current and cumulative
        ! number of exits.
        !
        ! The optional ``model`` argument selects the update rule:
        !
        ! * ``'TASEP'`` uses the single-site open-boundary TASEP update;
        ! * ``'NS'`` uses the Nagel-Schreckenberg update with fixed parameters
        !   ``V_MAX_NS`` and ``P_SLOW_NS``.
        !
        ! The occupation history is stored as integers:
        !
        ! .. math::
        !
        !    h_{i,t} =
        !    \begin{cases}
        !    1, & \text{if site } i \text{ is occupied at time } t, \\
        !    0, & \text{otherwise.}
        !    \end{cases}
        !
        ! The output ``current_history(step)`` is the number of vehicles exiting
        ! during that step. For the TASEP update this is either zero or one.
        integer, intent(in) :: L
        ! Number of sites in the one-dimensional lattice.
        integer, intent(in) :: n_steps
        ! Number of time steps to simulate.
        real, intent(in) :: alpha
        ! Entry probability at the left boundary.
        real, intent(in) :: beta
        ! Exit probability at the right boundary.
        integer, intent(out) :: history(L, n_steps)
        ! Occupation history. ``history(i, step)`` is 1 if site ``i`` is occupied
        ! at ``step`` and 0 otherwise.
        real, intent(out) :: density_history(n_steps)
        ! Mean density at each time step.
        integer, intent(out) :: current_history(n_steps)
        ! Number of vehicles exiting at each time step.
        integer, intent(out) :: total_exits
        ! Total number of vehicles exiting over the full simulation.
        character(len=*), intent(in), optional :: model
        ! Optional update rule selector: ``'TASEP'`` or ``'NS'``.

        character(len=10) :: use_model
        ! Update rule selected for this run.
        type(cell) :: state(L)
        ! Current road state.
        integer :: step
        ! Time-step index.
        integer :: exit_count
        ! Number of vehicles exiting during the current time step.

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

            history(:, step)      = merge(1, 0, state%has_car)
            density_history(step) = compute_density(state, L)
            current_history(step) = exit_count
            total_exits           = total_exits + exit_count
        end do

    end subroutine run_simulation


    subroutine measure_steady_state(L, n_burnin, n_measure, alpha, beta, mean_density, mean_current)
        ! Measure steady-state density and current for the TASEP model.
        !
        ! This routine first runs the TASEP model for ``n_burnin`` steps to reduce
        ! sensitivity to the initially empty lattice. It then measures the density
        ! and current over a further ``n_measure`` steps.
        !
        ! The density is measured only on the central quarter of the lattice to
        ! reduce boundary-layer bias from the open boundaries. The measurement
        ! region is:
        !
        ! .. math::
        !
        !    i \in \left[\frac{3L}{8}, \frac{5L}{8}\right].
        !
        ! The returned mean density is:
        !
        ! .. math::
        !
        !    \bar{\rho}
        !    =
        !    \frac{1}{T}
        !    \sum_{t=1}^{T}
        !    \rho_{\mathrm{bulk}}(t),
        !
        ! where ``T = n_measure``. The returned mean current is:
        !
        ! .. math::
        !
        !    \bar{J}
        !    =
        !    \frac{1}{T}
        !    \sum_{t=1}^{T}
        !    n_{\mathrm{exit}}(t).
        integer, intent(in) :: L
        ! Number of sites in the one-dimensional lattice.
        integer, intent(in) :: n_burnin
        ! Number of burn-in steps discarded before measurement.
        integer, intent(in) :: n_measure
        ! Number of time steps used for measurement.
        real, intent(in) :: alpha
        ! Entry probability at the left boundary.
        real, intent(in) :: beta
        ! Exit probability at the right boundary.
        real, intent(out) :: mean_density
        ! Time-averaged density measured on the central bulk region.
        real, intent(out) :: mean_current
        ! Time-averaged exit current over the measurement window.

        type(cell) :: state(L)
        ! Current road state.
        integer :: step
        ! Time-step index.
        integer :: exit_count
        ! Number of vehicles exiting during the current time step.
        integer :: bulk_lo
        ! First cell in the central bulk measurement region.
        integer :: bulk_hi
        ! Last cell in the central bulk measurement region.
        integer :: bulk_L
        ! Number of cells in the central bulk measurement region.
        real :: density_acc
        ! Accumulator for measured bulk density.
        integer :: current_acc
        ! Accumulator for measured exit current.

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


    subroutine ns_single_step(state, L, alpha, beta, exit_count)
        ! Advance the Nagel-Schreckenberg cellular automaton by one time step.
        !
        ! This is an internal helper routine used by ``run_simulation`` when the
        ! optional model selector is ``'NS'``. The Nagel-Schreckenberg model keeps
        ! an integer velocity for each vehicle. The update is parallel: all
        ! movements are determined from the old state before the new state is
        ! written.
        !
        ! The bulk update consists of four stages:
        !
        ! 1. Acceleration:
        !
        !    .. math::
        !
        !       v \leftarrow \min(v + 1, v_{\max}).
        !
        ! 2. Deceleration to avoid collisions:
        !
        !    .. math::
        !
        !       v \leftarrow \min(v, g),
        !
        !    where ``g`` is the gap to the next vehicle.
        !
        ! 3. Random slowing:
        !
        !    .. math::
        !
        !       v \leftarrow v - 1
        !
        !    with probability ``P_SLOW_NS`` if ``v > 0``.
        !
        ! 4. Movement:
        !
        !    .. math::
        !
        !       x \leftarrow x + v.
        !
        ! Open boundaries are treated using the same ``alpha`` and ``beta``
        ! probabilities as the TASEP driver.
        type(cell), intent(inout) :: state(L)
        ! Road state. On entry this is the current state; on return it is the
        ! updated state after one Nagel-Schreckenberg step.
        integer, intent(in) :: L
        ! Number of sites in the one-dimensional lattice.
        real, intent(in) :: alpha
        ! Entry probability at the left boundary.
        real, intent(in) :: beta
        ! Exit probability at the right boundary.
        integer, intent(out) :: exit_count
        ! Number of vehicles exiting during this time step.

        type(cell) :: old_state(L)
        ! Copy of the road state before the update.
        type(cell) :: new_state(L)
        ! Road state being constructed for the next time step.
        integer :: i
        ! Current lattice-site index.
        integer :: j
        ! Search index used to find the next vehicle ahead.
        integer :: gap
        ! Number of empty cells between the current vehicle and the next obstacle.
        integer :: v
        ! Updated vehicle velocity.
        integer :: new_pos
        ! New lattice position after movement.
        real :: rnd
        ! Uniform random number used for stochastic entry, exit and slowing.

        old_state  = state
        new_state  = old_state
        exit_count = 0

        ! Clear positions occupied at the start of the step.
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

            ! Gap to next car, or to the end of the road if no car is ahead.
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

        ! Entry at site 1 only if it was empty at the start of the step.
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