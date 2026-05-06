module junction_mod
    !--------------------------------------------------------------------
    ! Junction evaluation: moves vehicles from holding cells (site L of
    ! each leg's inbound lane) into destination cells (site 1 of the
    ! relevant outbound lane), respecting the right-of-way hierarchy:
    !
    !   1. Physical block: never advance into an occupied destination cell.
    !   2. Yield to the right (R1): yield to the vehicle one step
    !      counter-clockwise (the CW-prev leg), when paths cross or when
    !      both vehicles target the same outbound.
    !   3. Right-turn yields oncoming (R2): a RIGHT-category vehicle yields
    !      to any LEFT or STRAIGHT vehicle from the opposite leg.
    !      (n_in == n_out == 4 only.)
    !   4. Opposite-leg conflict scan (R3): remaining delta=2 crossings —
    !        (S, L) opposite: straight has priority, left yields (UK rule)
    !        (R, R) opposite: mutual yields -> stochastic deadlock break
    !      (n_in == n_out == 4 only.)
    !   5. Parallel non-conflicting moves are allowed.
    !   6. Yield semantics: if you yield, you wait the FULL step regardless
    !      of whether your yield-target also moves this step.
    !   7. Symmetric standoff: stochastic tie-break selects one winner;
    !      remaining deadlocked candidates wait for the next timestep.
    !
    ! Move categories (derived from inbound and outbound leg indices):
    !   0 = UTURN    d = 0
    !   1 = LEFT     d = 1  (CW-next; UK easy turn)
    !   2 = STRAIGHT d = n/2
    !   3 = RIGHT    d = n-1 (CW-prev; cuts across)
    !  where d = modulo(out_idx - in_idx, n).
    !--------------------------------------------------------------------
    use vehicle_mod
    use road_network_mod
    implicit none
    private

    public :: evaluate_junctions, move_category

contains

    subroutine evaluate_junctions(net)
        type(road_network_t), intent(inout) :: net
        integer :: j
        do j = 1, size(net%junctions)
            call evaluate_one_junction(net, j)
        end do
    end subroutine evaluate_junctions

    !-----------------------------------------------------------------
    subroutine evaluate_one_junction(net, jid)
        type(road_network_t), intent(inout) :: net
        integer, intent(in) :: jid

        integer, parameter :: MAX_LEGS = 4

        integer :: N, k, m, src_road, src_lane, L_src
        integer :: dst_road(MAX_LEGS), dst_lane(MAX_LEGS), dst_idx(MAX_LEGS)
        integer :: holding(MAX_LEGS)
        logical :: phys_clear(MAX_LEGS), yields_to(MAX_LEGS, MAX_LEGS), approved(MAX_LEGS)
        real    :: r, cumsum

        N = net%junctions(jid)%n_in
        if (N > MAX_LEGS) then
            print *, 'evaluate_one_junction: N exceeds MAX_LEGS'
            stop 1
        end if

        holding    = V_EMPTY
        dst_road   = 0
        dst_lane   = 0
        dst_idx    = 0
        phys_clear = .false.
        yields_to  = .false.
        approved   = .false.

        ! 1. Read holding cells from the old snapshot.
        do k = 1, N
            src_road   = net%junctions(jid)%in_road(k)
            src_lane   = net%junctions(jid)%in_lane(k)
            L_src      = net%roads(src_road)%lane(src_lane)%length
            holding(k) = net%roads(src_road)%lane(src_lane)%old(L_src)%turning_intent
        end do

        ! 2. Sample destination from in_routes for each occupied leg.
        do k = 1, N
            if (holding(k) == 0) cycle
            call random_number(r)
            cumsum    = 0.0
            dst_idx(k) = net%junctions(jid)%n_out   ! fallback: last outbound
            do m = 1, net%junctions(jid)%n_out
                cumsum = cumsum + net%junctions(jid)%in_routes(k)%prob(m)
                if (r < cumsum) then
                    dst_idx(k) = m
                    exit
                end if
            end do
            dst_road(k) = net%junctions(jid)%out_road(dst_idx(k))
            dst_lane(k) = net%junctions(jid)%out_lane(dst_idx(k))
        end do

        ! 3. Physical-clear check (destination site 1 in the old snapshot).
        do k = 1, N
            if (holding(k) == V_EMPTY) cycle
            phys_clear(k) = (.not. net%roads(dst_road(k))%lane(dst_lane(k))%old(1)%has_car)
        end do

        ! 4. Build yield matrix.
        if (N == 4 .and. net%junctions(jid)%n_out == 4) then
            call build_yield_matrix_v2(net%junctions(jid), holding, dst_idx, yields_to)
        else
            call build_yield_matrix_asym(net%junctions(jid), holding, dst_idx, yields_to)
        end if

        ! 5. First approval pass.
        call approve_pass(holding, phys_clear, yields_to, approved)

        ! 6. Stochastic deadlock break for remaining symmetric standoffs.
        call resolve_deadlocks(holding, phys_clear, yields_to, approved)

        ! 7. Commit approved moves.
        do k = 1, N
            if (.not. approved(k)) cycle
            src_road = net%junctions(jid)%in_road(k)
            src_lane = net%junctions(jid)%in_lane(k)
            L_src    = net%roads(src_road)%lane(src_lane)%length
            net%roads(dst_road(k))%lane(dst_lane(k))%cells(1)%has_car = .true.
            net%roads(dst_road(k))%lane(dst_lane(k))%cells(1)%turning_intent = holding(k)
            net%roads(dst_road(k))%lane(dst_lane(k))%cells(1)%velocity = 0
            net%roads(src_road)%lane(src_lane)%cells(L_src)%has_car   = .false.
            net%roads(src_road)%lane(src_lane)%cells(L_src)%turning_intent = 0
            net%roads(src_road)%lane(src_lane)%cells(L_src)%velocity = 0
        end do
    end subroutine evaluate_one_junction

    !-----------------------------------------------------------------
    ! Move category: derived purely from inbound index, outbound index,
    ! and junction size.  No cell-level intent code needed.
    !-----------------------------------------------------------------
    pure function move_category(in_idx, out_idx, n) result(cat)
        integer, intent(in) :: in_idx, out_idx, n
        integer :: cat, d
        d = modulo(out_idx - in_idx, n)
        if      (d == 0)                    then; cat = 0   ! UTURN
        else if (d == 1)                    then; cat = 1   ! LEFT  (CW-next)
        else if (n >= 4 .and. d == n/2)     then; cat = 2   ! STRAIGHT
        else if (d == n - 1)                then; cat = 3   ! RIGHT (CW-prev)
        else                                    ; cat = -1
        end if
    end function move_category

    !-----------------------------------------------------------------
    ! Path-conflict predicate for symmetric (n_in == n_out) layouts.
    ! Returns .false. for n /= 4 (R2/R3 conditions gate on n==4 anyway).
    !-----------------------------------------------------------------
    pure function paths_conflict_sym(n, cat_i, cat_j, delta) result(c)
        integer, intent(in) :: n, cat_i, cat_j, delta
        logical :: c
        c = .false.
        if (n /= 4) return
        select case (delta)
        case (1, 3)   ! adjacent leg
            c = .not. (cat_i == 1 .and. cat_j == 1)   ! conflict unless both LEFT
        case (2)      ! opposite leg
            c = .not. (cat_i == 2 .and. cat_j == 2) &
                .and. .not. (cat_i == 1 .and. cat_j == 1)
        end select
    end function paths_conflict_sym

    !-----------------------------------------------------------------
    ! Yield matrix construction (category-based, Phase 3+).
    !-----------------------------------------------------------------
    subroutine build_yield_matrix_v2(jn, holding, dst_idx, yields_to)
        type(junction_t), intent(in)  :: jn
        integer,          intent(in)  :: holding(:), dst_idx(:)
        logical,          intent(out) :: yields_to(:,:)
        integer :: N, i, right_of_i, opp, cat_i, cat_opp, cat_right

        N = jn%n_in
        yields_to = .false.

        do i = 1, N
            if (holding(i) == V_EMPTY) cycle

            cat_i      = move_category(i, dst_idx(i), N)
            right_of_i = mod(i + N - 2, N) + 1

            ! R1 — yield to the right when paths conflict or same destination.
            if (holding(right_of_i) /= V_EMPTY) then
                cat_right = move_category(right_of_i, dst_idx(right_of_i), N)
                if (paths_conflict_sym(N, cat_i, cat_right, 3) .or. &
                    dst_idx(i) == dst_idx(right_of_i)) then
                    yields_to(i, right_of_i) = .true.
                end if
            end if

            ! R2 — right-turn yields to oncoming straight or left (4-way only).
            if (N == 4 .and. jn%n_out == 4 .and. cat_i == 3) then
                opp = mod(i + 1, N) + 1
                if (holding(opp) /= V_EMPTY) then
                    cat_opp = move_category(opp, dst_idx(opp), N)
                    if (cat_opp == 1 .or. cat_opp == 2) yields_to(i, opp) = .true.
                end if
            end if
        end do

        ! R3 — remaining delta=2 crossings not already resolved by R2 (4-way only).
        if (N == 4 .and. jn%n_out == 4) then
            do i = 1, N
                if (holding(i) == V_EMPTY) cycle
                opp = mod(i + 1, N) + 1
                if (holding(opp) == V_EMPTY) cycle
                cat_i   = move_category(i,   dst_idx(i),   N)
                cat_opp = move_category(opp, dst_idx(opp), N)
                if (.not. (paths_conflict_sym(N, cat_i, cat_opp, 2) .or. &
                           dst_idx(i) == dst_idx(opp))) cycle
                if (yields_to(i, opp) .or. yields_to(opp, i)) cycle   ! already resolved

                select case (cat_i * 4 + cat_opp)
                case (2*4 + 1)   ! i=STRAIGHT, opp=LEFT: LEFT yields to STRAIGHT
                    yields_to(opp, i) = .true.
                case (1*4 + 2)   ! i=LEFT, opp=STRAIGHT: LEFT yields to STRAIGHT
                    yields_to(i, opp) = .true.
                case default     ! (R,R) and others: mutual yields -> deadlock
                    yields_to(i,   opp) = .true.
                    yields_to(opp, i  ) = .true.
                end select
            end do
        end if
    end subroutine build_yield_matrix_v2

    !-----------------------------------------------------------------
    ! Yield matrix for asymmetric junctions (n_in /= 4 or n_out /= 4).
    ! Only R1 applies: yield to the CW-prev neighbour when paths cross
    ! (chord-crossing predicate) or when both target the same outbound.
    ! R2 and R3 are skipped.
    !-----------------------------------------------------------------
    subroutine build_yield_matrix_asym(jn, holding, dst_idx, yields_to)
        type(junction_t), intent(in)  :: jn
        integer,          intent(in)  :: holding(:), dst_idx(:)
        logical,          intent(out) :: yields_to(:,:)
        integer :: N, i, right_of_i, n_perim

        N       = jn%n_in
        n_perim = jn%n_in + jn%n_out
        yields_to = .false.

        do i = 1, N
            if (holding(i) == V_EMPTY) cycle
            right_of_i = mod(i + N - 2, N) + 1
            if (holding(right_of_i) == V_EMPTY) cycle
            ! R1: yield to right if chords cross or same destination.
            if (chords_cross(jn%in_perim(i),          jn%out_perim(dst_idx(i)), &
                             jn%in_perim(right_of_i), jn%out_perim(dst_idx(right_of_i)), &
                             n_perim) &
                .or. dst_idx(i) == dst_idx(right_of_i)) then
                yields_to(i, right_of_i) = .true.
            end if
        end do
    end subroutine build_yield_matrix_asym

    !-----------------------------------------------------------------
    ! Chord-crossing predicate on a cyclic perimeter of n_perim ports
    ! (0-indexed).  Two chords (pa1->pa2) and (pb1->pb2) cross iff
    ! exactly one of {pb1, pb2} lies in the open CW arc from pa1 to pa2.
    !-----------------------------------------------------------------
    pure function chords_cross(pa1, pa2, pb1, pb2, n_perim) result(c)
        integer, intent(in) :: pa1, pa2, pb1, pb2, n_perim
        logical :: c
        c = in_arc(pa1, pa2, pb1, n_perim) .neqv. in_arc(pa1, pa2, pb2, n_perim)
    end function chords_cross

    pure function in_arc(a, b, p, n) result(res)
        ! Is port p strictly inside the open CW arc from a to b?
        integer, intent(in) :: a, b, p, n
        logical :: res
        integer :: ab, ap
        if (a == b) then
            res = .false.
        else
            ab = modulo(b - a, n)
            ap = modulo(p - a, n)
            res = (ap > 0) .and. (ap < ab)
        end if
    end function in_arc

    !-----------------------------------------------------------------
    ! Approval pass.
    !
    ! Blocking rule: a candidate k is blocked if ANY of its yield targets
    ! has a vehicle in its holding cell (in the old snapshot), regardless
    ! of whether that target is itself approved this timestep.
    !
    ! This "full-step yield" semantics means: if you must yield, you
    ! wait the entire timestep even if your priority vehicle also moves.
    !-----------------------------------------------------------------
    pure subroutine approve_pass(holding, phys_clear, yields_to, approved)
        integer, intent(in)    :: holding(:)
        logical, intent(in)    :: phys_clear(:)
        logical, intent(in)    :: yields_to(:,:)
        logical, intent(inout) :: approved(:)
        integer :: i, N
        logical :: yields_live

        N = size(holding)
        do i = 1, N
            if (approved(i))           cycle
            if (holding(i) == V_EMPTY) cycle
            if (.not. phys_clear(i))   cycle

            ! Blocked if any yield target holds a vehicle (full-step yield).
            yields_live = any(yields_to(i, :) .and. (holding /= V_EMPTY))
            if (.not. yields_live) approved(i) = .true.
        end do
    end subroutine approve_pass

    !-----------------------------------------------------------------
    ! Deadlock breaker.
    !
    ! A candidate is "in a deadlock" if it is yield-blocked by a target
    ! that is itself non-approved and non-deadlocked (i.e., it cannot be
    ! unblocked without a tie-break).  Vehicles blocked only by an
    ! already-approved mover are simply waiting — no tie-break needed.
    !
    ! When a deadlock set is found:
    !   1. One member is chosen uniformly at random and approved.
    !   2. Every OTHER member of that set is marked `deadlocked` for
    !      this timestep and excluded from future iterations.
    !
    ! This "mark the losers" strategy ensures that a single tie-break
    ! resolves the entire standoff in one shot: the losers wait for the
    ! next network timestep and retry.  It prevents the cascade where
    ! successive iterations keep picking more winners from a ring.
    !-----------------------------------------------------------------
    subroutine resolve_deadlocks(holding, phys_clear, yields_to, approved)
        integer, intent(in)    :: holding(:)
        logical, intent(in)    :: phys_clear(:)
        logical, intent(in)    :: yields_to(:,:)
        logical, intent(inout) :: approved(:)
        integer :: N, k, idx, n_locked, pick
        integer, allocatable :: locked(:)
        logical :: deadlocked(size(holding))
        real    :: r
        logical :: in_deadlock

        N = size(holding)
        deadlocked = .false.
        allocate(locked(N))

        do
            n_locked = 0
            do k = 1, N
                if (approved(k))   cycle
                if (deadlocked(k)) cycle
                if (holding(k) == V_EMPTY) cycle
                if (.not. phys_clear(k))   cycle
                ! In a deadlock only if ALL blocking yield targets are non-approved
                ! and non-deadlocked.  A candidate that also yields to an already-
                ! approved vehicle must wait the full step; it is not a deadlock member.
                in_deadlock = any(yields_to(k, :) .and. (holding /= V_EMPTY) &
                                  .and. .not. approved .and. .not. deadlocked) &
                              .and. .not. any(yields_to(k, :) .and. (holding /= V_EMPTY) &
                                             .and. approved)
                if (in_deadlock) then
                    n_locked = n_locked + 1
                    locked(n_locked) = k
                end if
            end do
            if (n_locked == 0) exit

            ! Pick one winner stochastically.
            call random_number(r)
            pick = locked(min(int(r * n_locked) + 1, n_locked))
            approved(pick) = .true.

            ! Mark every other member of this deadlock set as blocked for
            ! this timestep.  They will retry at the next network step.
            do idx = 1, n_locked
                if (locked(idx) /= pick) deadlocked(locked(idx)) = .true.
            end do
        end do

        deallocate(locked)
    end subroutine resolve_deadlocks

end module junction_mod
