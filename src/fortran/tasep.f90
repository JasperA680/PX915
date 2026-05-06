module tasep_model
    use road_network_mod

    implicit none
    private

    !--------------------------------------------------------------------------
    ! Public routines exposed by this module
    !
    ! initialise_lattice  : Creates an empty 1D lattice of length L
    !                       representing the road (all cells initially empty).
    !
    ! tasep_step          : Performs one timestep of the open-boundary
    !                       TASEP update. Handles particle entry (alpha),
    !                       bulk movement, and particle exit (beta).
    !
    ! count_occupied      : Returns the number of occupied sites in the lattice.
    !                       Useful for computing traffic density.
    !
    ! compute_density     : Computes the traffic density rho = N/L,
    !                       where N is the number of vehicles and
    !                       L is the road length.
    !--------------------------------------------------------------------------
    public :: initialise_lattice
    public :: tasep_step
    public :: count_occupied
    public :: compute_density

contains 

   subroutine initialise_lattice(state, L)
        ! Initialise an empty lattice of length L
        integer, intent(in) :: L 
        type(cell), intent(out) :: state(L)

        ! state%has_car = .false.  ! Set all cells to empty, should be done automatically i think?
    end subroutine initialise_lattice

    subroutine tasep_step(state, L, alpha, beta, exit_count)
        ! Perform one parallel update tasep step with open boundaries.
        !
        ! Rules:
        ! 1. A particle at site L exits with probability beta.
        ! 2. Particles in the bulk hop one site to the right if the next site is empty.
        ! 3. A particle enters at site 1 with probability alpha if site 1 is empty.
        !
        ! Input:
        !    state       - current lattice state (0 empty, 1 occupied)
        !    L           - length of the lattice
        !    alpha       - entry probability at site 1
        !    beta        - exit probability at site L
        !
        ! Output:
        !    state       - updated lattice state after one TASEP step
        !    exit_count  - number of particles that exited during this step

        integer, intent(in) :: L 
        type(cell), intent(inout) :: state(L)
        real, intent(in) :: alpha, beta
        integer, intent(out) :: exit_count

        type(cell) :: old_state(L), new_state(L)
        integer :: i
        real    :: r

        old_state = state 
        new_state = old_state
        exit_count = 0

        ! Exit at the right boundary
        if (old_state(L)%has_car) then                  ! If there is a car at the end, let it leave with probablity beta
            call random_number(r) 
            if (r < beta) then 
                new_state(L)%has_car = .false.  
                exit_count = 1
            end if
        end if

        ! Bulk motion 
        do i = L-1, 1, -1
            if (old_state(i)%has_car .and. .not. old_state(i+1)%has_car) then 
                ! Move particle right by one site
                new_state(i)%has_car = .false.
                new_state(i+1)%has_car = .true.
            end if   
        end do

        ! Entry at the left boundary
        if (.not. new_state(1)%has_car) then 
            call random_number(r)
            if (r < alpha) then 
                new_state(1)%has_car = .true.
            end if
        end if

        state = new_state

    end subroutine tasep_step


    function count_occupied(state, L) result(nocc)
        ! Count the number of occupied states
        integer, intent(in) :: L 
        type(cell), intent(in) :: state(L)
        integer             :: nocc 

        nocc = count(state%has_car)
    end function count_occupied

    function compute_density(state, L) result(rho)
        ! compute the density = number of occupied states / L
        integer, intent(in) :: L 
        type(cell), intent(in) :: state(L)
        real                :: rho

        rho = real(count(state%has_car))/real(L)
    end function compute_density

end module tasep_model 

    