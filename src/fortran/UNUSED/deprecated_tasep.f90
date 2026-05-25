module tasep_model
    ! One-dimensional open-boundary TASEP cellular automaton model.
    !
    ! This module implements a simple totally asymmetric simple exclusion process
    ! (TASEP) on a one-dimensional lattice. Each lattice site is either empty or
    ! occupied by one vehicle. Vehicles move only to the right and obey exclusion:
    ! a vehicle can move into the next site only if that site is empty.
    !
    ! The model uses open boundaries. At each time step:
    !
    ! - a vehicle may leave the right boundary with probability ``beta``;
    ! - vehicles in the bulk move one site to the right if the next site is empty;
    ! - a vehicle may enter the left boundary with probability ``alpha`` if the
    !   first site is empty.
    !
    ! The density is defined as:
    !
    ! .. math::
    !
    !    \rho = \frac{N}{L},
    !
    ! where ``N`` is the number of occupied sites and ``L`` is the lattice length.
    !
    ! The measured exit current over a simulation can be estimated from:
    !
    ! .. math::
    !
    !    J = \frac{N_{\mathrm{exit}}}{N_{\mathrm{steps}}},
    !
    ! where ``N_exit`` is the total number of vehicles leaving the right boundary
    ! over ``N_steps`` update steps.

    use road_network_mod

    implicit none
    private

    public :: initialise_lattice
    public :: tasep_step
    public :: count_occupied
    public :: compute_density

contains

    subroutine initialise_lattice(state, L)
        ! Initialise an empty one-dimensional TASEP lattice.
        !
        ! This routine creates a road lattice of length ``L`` where all sites are
        ! initially empty. The ``velocity`` field is also reset to zero for every
        ! cell. In this TASEP implementation, the main dynamical variable is
        ! whether a cell contains a car; the velocity field is retained for
        ! compatibility with the wider road-network data structure.
        integer, intent(in) :: L
        ! Number of lattice sites.
        type(cell), intent(out) :: state(L)
        ! Output lattice state. On return, all sites are empty.

        state%has_car = .false.
        state%velocity = 0
    end subroutine initialise_lattice


    subroutine tasep_step(state, L, alpha, beta, exit_count)
        ! Advance the open-boundary TASEP model by one time step.
        !
        ! This routine performs one parallel-update TASEP step. The old lattice
        ! state is copied before any moves are applied, so all bulk movements are
        ! decided from the same previous-time configuration.
        !
        ! The update consists of three stages:
        !
        ! 1. If site ``L`` is occupied, the vehicle exits with probability
        !    ``beta``.
        ! 2. Each bulk vehicle attempts to hop one site to the right. The hop is
        !    accepted if the target site was empty in the old state.
        ! 3. If site ``1`` is empty after the bulk update, a new vehicle enters
        !    with probability ``alpha``.
        !
        ! The boundary rates satisfy:
        !
        ! .. math::
        !
        !    0 \leq \alpha \leq 1,
        !    \qquad
        !    0 \leq \beta \leq 1.
        !
        ! For a single step, ``exit_count`` is either zero or one. Over many steps,
        ! summing ``exit_count`` gives the total number of vehicles that have left
        ! the road, which can be used to estimate the current:
        !
        ! .. math::
        !
        !    J \approx \frac{1}{T}\sum_{t=1}^{T} n_{\mathrm{exit}}(t).
        integer, intent(in) :: L
        ! Number of lattice sites.
        type(cell), intent(inout) :: state(L)
        ! Lattice state. On entry, this is the current state. On return, it is the
        ! updated state after one TASEP step.
        real, intent(in) :: alpha
        ! Entry probability at the left boundary.
        real, intent(in) :: beta
        ! Exit probability at the right boundary.
        integer, intent(out) :: exit_count
        ! Number of vehicles that exited during this time step.

        type(cell) :: old_state(L)
        ! Copy of the lattice before the update.
        type(cell) :: new_state(L)
        ! Lattice state being constructed for the next time step.
        integer :: i
        ! Lattice-site index used for the bulk update.
        real :: r
        ! Uniform random number used for stochastic entry and exit events.

        old_state = state
        new_state = old_state
        exit_count = 0

        ! Exit at the right boundary.
        if (old_state(L)%has_car) then
            call random_number(r)
            if (r < beta) then
                new_state(L)%has_car = .false.
                new_state(L)%velocity = 0
                exit_count = 1
            end if
        end if

        ! Bulk motion.
        do i = L - 1, 1, -1
            if (old_state(i)%has_car .and. .not. old_state(i+1)%has_car) then
                new_state(i)%has_car = .false.
                new_state(i)%velocity = 0
                new_state(i+1)%has_car = .true.
                new_state(i+1)%velocity = 0
            end if
        end do

        ! Entry at the left boundary.
        if (.not. new_state(1)%has_car) then
            call random_number(r)
            if (r < alpha) then
                new_state(1)%has_car = .true.
                new_state(1)%velocity = 0
            end if
        end if

        state = new_state
    end subroutine tasep_step


    function count_occupied(state, L) result(nocc)
        ! Count the number of occupied sites in the lattice.
        !
        ! This routine returns:
        !
        ! .. math::
        !
        !    N = \sum_{i=1}^{L} \tau_i,
        !
        ! where ``tau_i = 1`` if site ``i`` contains a vehicle and ``tau_i = 0``
        ! otherwise.
        integer, intent(in) :: L
        ! Number of lattice sites.
        type(cell), intent(in) :: state(L)
        ! Current lattice state.
        integer :: nocc
        ! Number of occupied sites.

        nocc = count(state%has_car)
    end function count_occupied


    function compute_density(state, L) result(rho)
        ! Compute the mean lattice density.
        !
        ! The density is the fraction of lattice sites that are occupied:
        !
        ! .. math::
        !
        !    \rho = \frac{N}{L}
        !         = \frac{1}{L}\sum_{i=1}^{L} \tau_i.
        !
        ! Here ``tau_i`` is one for an occupied site and zero for an empty site.
        integer, intent(in) :: L
        ! Number of lattice sites.
        type(cell), intent(in) :: state(L)
        ! Current lattice state.
        real :: rho
        ! Mean density of the lattice.

        rho = real(count(state%has_car)) / real(L)
    end function compute_density

end module tasep_model