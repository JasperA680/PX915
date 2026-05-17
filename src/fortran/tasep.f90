module tasep_model
    use omp_lib
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
    !
    ! rng_init            : Seed the per-thread xoshiro128++ RNG. Call once
    !                       before any simulation. Each thread gets an
    !                       independent stream derived from base_seed.
    !--------------------------------------------------------------------------
    public :: initialise_lattice
    public :: tasep_step
    public :: count_occupied
    public :: compute_density
    public :: rng_init
    public :: tasep_step_omp

    ! Per-thread xoshiro128++ state. Default seed is non-zero so the RNG
    ! works even if rng_init is never called.
    integer(kind=4), save :: rng_xsh(4) = [123456789_4, -987654321_4, 362436069_4, 521288629_4]
    !$omp threadprivate(rng_xsh)

contains

    ! splitmix32 hash — used only by rng_init to derive independent per-thread seeds.
    ! Constants are the standard splitmix32 multipliers (Murmur / Stafford variant).
    function splitmix32(x) result(y)
        integer(kind=4), intent(in) :: x
        integer(kind=4) :: y
        y = x + int(z'9E3779B9', 4)
        y = ieor(y, ishft(y, -16))
        y = y * int(z'85EBCA6B', 4)
        y = ieor(y, ishft(y, -13))
        y = y * int(z'C2B2AE35', 4)
        y = ieor(y, ishft(y, -16))
    end function splitmix32

    subroutine rng_init(base_seed)
        ! Seed every OpenMP thread's xoshiro128++ state independently.
        ! Thread k gets a stream derived from base_seed + k * golden-ratio constant.
        ! Must be called from a serial region; spawns its own parallel region internally.
        integer, intent(in) :: base_seed
        integer :: tid
        integer(kind=4) :: s
        !$omp parallel private(tid, s)
        tid = omp_get_thread_num()
        s   = int(base_seed, 4) + int(tid, 4) * int(z'9E3779B9', 4)
        rng_xsh(1) = splitmix32(s)
        rng_xsh(2) = splitmix32(rng_xsh(1))
        rng_xsh(3) = splitmix32(rng_xsh(2))
        rng_xsh(4) = splitmix32(rng_xsh(3))
        !$omp end parallel
    end subroutine rng_init

    function rng_real() result(r)
        ! xoshiro128++ step. Returns a uniform real in [0, 1).
        ! Uses lower 24 bits of the result to stay within single-precision mantissa.
        real :: r
        integer(kind=4) :: res, t
        res        = ishftc(rng_xsh(1) + rng_xsh(4), 7) + rng_xsh(1)
        t          = ishft(rng_xsh(2), 9)
        rng_xsh(3) = ieor(rng_xsh(3), rng_xsh(1))
        rng_xsh(4) = ieor(rng_xsh(4), rng_xsh(2))
        rng_xsh(2) = ieor(rng_xsh(2), rng_xsh(3))
        rng_xsh(1) = ieor(rng_xsh(1), rng_xsh(4))
        rng_xsh(3) = ieor(rng_xsh(3), t)
        rng_xsh(4) = ishftc(rng_xsh(4), 11)
        r = real(iand(res, 16777215_4)) * (1.0 / 16777216.0)
    end function rng_real

   subroutine initialise_lattice(state, L)
        ! Initialise an empty lattice of length L
        integer, intent(in) :: L 
        integer, intent(out) :: state(L)

        state = 0  ! Set all sites to empty
    end subroutine initialise_lattice

    subroutine tasep_step(state, L, alpha, beta, exit_count)
        ! Perform one parallel update TASEP step with open boundaries.
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
        integer, intent(inout) :: state(L)
        real, intent(in) :: alpha, beta
        integer, intent(out) :: exit_count

        integer :: old_state(L), new_state(L)
        integer :: i
        real    :: r

        old_state = state 
        new_state = old_state
        exit_count = 0

        ! Exit at the right boundary
        if (old_state(L) == 1) then
            r = rng_real()
            if (r < beta) then
                new_state(L) = 0
                exit_count = 1
            end if
        end if

        ! Bulk motion — snapshot old_state is read-only; adjacent predicates are
        ! mutually exclusive, so no write race between iterations.
        !$omp parallel do default(shared) private(i)
        do i = L-1, 1, -1
            if (old_state(i) == 1 .and. old_state(i+1) == 0) then
                new_state(i)   = 0
                new_state(i+1) = 1
            end if
        end do
        !$omp end parallel do

        ! Entry at the left boundary
        if (new_state(1) == 0) then
            r = rng_real()
            if (r < alpha) then
                new_state(1) = 1
            end if
        end if

        state = new_state

    end subroutine tasep_step


    subroutine tasep_step_omp(state, old_state, new_state, L, alpha, beta, exit_count)
        ! Like tasep_step but designed to be called from within an !$omp parallel region.
        ! Uses !$omp single for serial boundary work and !$omp do for the bulk loop so
        ! the caller's persistent thread team is reused across steps — no per-step
        ! thread-team creation cost. old_state and new_state must be declared in the
        ! calling scope (shared) and passed as workspace.
        integer, intent(inout) :: state(L), old_state(L), new_state(L)
        integer, intent(in)    :: L
        real,    intent(in)    :: alpha, beta
        integer, intent(out)   :: exit_count
        integer :: i
        real    :: r

        !$omp single
        old_state  = state
        new_state  = old_state
        exit_count = 0
        if (old_state(L) == 1) then
            r = rng_real()
            if (r < beta) then
                new_state(L) = 0
                exit_count   = 1
            end if
        end if
        !$omp end single

        !$omp do
        do i = L-1, 1, -1
            if (old_state(i) == 1 .and. old_state(i+1) == 0) then
                new_state(i)   = 0
                new_state(i+1) = 1
            end if
        end do
        !$omp end do

        !$omp single
        if (new_state(1) == 0) then
            r = rng_real()
            if (r < alpha) then
                new_state(1) = 1
            end if
        end if
        state = new_state
        !$omp end single

    end subroutine tasep_step_omp

    function count_occupied(state, L) result(nocc)
        ! Count the number of occupied states
        integer, intent(in) :: L 
        integer, intent(in) :: state(L)
        integer             :: nocc 

        nocc = sum(state)
    end function count_occupied

    function compute_density(state, L) result(rho)
        ! compute the density = number of occupied states / L
        integer, intent(in) :: L 
        integer, intent(in) :: state(L)
        real                :: rho

        rho = real(sum(state))/real(L)
    end function compute_density

end module tasep_model 

    