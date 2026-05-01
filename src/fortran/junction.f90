module junction_mod
    !--------------------------------------------------------------------
    ! Junction evaluation: moves vehicles from holding cells (site L of
    ! each leg's inbound lane) into destination cells (site 1 of the
    ! relevant outbound lane), respecting the right-of-way hierarchy:
    !
    !   1. Physical block: a vehicle never advances into an occupied cell.
    !   2. Yield to the right: yield to the holding-cell vehicle one
    !      step counter-clockwise (index i-1 wrap), but only if their
    !      paths actually cross.
    !   3. Right-turn yields oncoming: a vehicle turning right yields
    !      to any straight or left mover from the opposite leg.
    !   4. Parallel non-conflicting moves are allowed.
    !   5. Symmetric standoffs are broken by a single random pick.
    !
    ! All directional logic is modulo arithmetic on the clockwise leg
    ! index k (1..N).  The destination of intent X from leg k is:
    !
    !     STRAIGHT  ->  mod(k + N/2 - 1, N) + 1   (opposite road)
    !     LEFT      ->  mod(k,           N) + 1   (next clockwise)
    !     RIGHT     ->  mod(k + N - 2,   N) + 1   (previous clockwise)
    !
    ! No string compass labels are used inside the model.
    !--------------------------------------------------------------------
    use vehicle_mod
    use road_network_mod
    implicit none
    private

    public :: evaluate_junctions, compute_destination, paths_conflict

contains

    subroutine evaluate_junctions(net)
        type(road_network_t), intent(inout) :: net
        integer :: j
        do j = 1, size(net%junctions)
            call evaluate_one_junction(net, j)
        end do
    end subroutine evaluate_junctions

    !-----------------------------------------------------------------
    ! Per-junction evaluation.
    !-----------------------------------------------------------------
    subroutine evaluate_one_junction(net, jid)
        type(road_network_t), intent(inout) :: net
        integer, intent(in) :: jid

        integer, parameter :: MAX_LEGS = 4

        integer :: N, k
        integer :: src_road, src_lane, src_end, L_src
        integer :: dst_road(MAX_LEGS), dst_lane(MAX_LEGS), dst_end(MAX_LEGS)
        integer :: holding(MAX_LEGS)
        logical :: phys_clear(MAX_LEGS)
        logical :: yields_to(MAX_LEGS, MAX_LEGS)
        logical :: approved(MAX_LEGS)

        N = net%junctions(jid)%n_legs
        if (N > MAX_LEGS) then
            print *, 'evaluate_one_junction: N exceeds MAX_LEGS = ', MAX_LEGS
            stop 1
        end if

        holding   = V_EMPTY
        dst_road  = 0
        dst_lane  = 0
        dst_end   = 0
        phys_clear = .false.
        yields_to  = .false.
        approved   = .false.

        ! 1. Read holding cells (use the snapshot 'old' so junction sees
        !    a coherent pre-step view across all legs).
        do k = 1, N
            src_road = net%junctions(jid)%connected_road_ids(k)
            src_end  = net%junctions(jid)%end_at(k)
            src_lane = inbound_lane_idx(src_end)
            L_src    = net%roads(src_road)%lane(src_lane)%length
            holding(k) = net%roads(src_road)%lane(src_lane)%old(L_src)
        end do

        ! 2. Resolve destinations from indicators.
        do k = 1, N
            if (holding(k) == V_EMPTY) cycle
            call compute_destination(net%junctions(jid), k, holding(k), &
                                     dst_road(k), dst_lane(k), dst_end(k))
        end do

        ! 3. Physical-clear check on the destination cell.
        do k = 1, N
            if (holding(k) == V_EMPTY) cycle
            phys_clear(k) = (net%roads(dst_road(k))%lane(dst_lane(k))%old(1) == V_EMPTY)
        end do

        ! 4. Build yield matrix.
        call build_yield_matrix(net%junctions(jid), holding, yields_to)

        ! 5. First approval pass.
        call approve_pass(holding, phys_clear, yields_to, approved)

        ! 6. Resolve symmetric standoffs.
        call resolve_deadlocks(holding, phys_clear, yields_to, approved)

        ! 7. Commit approved moves.
        do k = 1, N
            if (.not. approved(k)) cycle
            src_road = net%junctions(jid)%connected_road_ids(k)
            src_end  = net%junctions(jid)%end_at(k)
            src_lane = inbound_lane_idx(src_end)
            L_src    = net%roads(src_road)%lane(src_lane)%length
            net%roads(dst_road(k))%lane(dst_lane(k))%cells(1) = holding(k)
            net%roads(src_road)%lane(src_lane)%cells(L_src)   = V_EMPTY
        end do
    end subroutine evaluate_one_junction

    !-----------------------------------------------------------------
    ! Destination resolution: pure modulo arithmetic on clockwise index.
    !-----------------------------------------------------------------
    pure subroutine compute_destination(jn, k, intent_code, dst_road, dst_lane, dst_end)
        type(junction_t), intent(in)  :: jn
        integer,          intent(in)  :: k, intent_code
        integer,          intent(out) :: dst_road, dst_lane, dst_end
        integer :: dst_k

        select case (intent_code)
        case (V_STRAIGHT)
            dst_k = mod(k + jn%n_legs/2 - 1, jn%n_legs) + 1
        case (V_LEFT)
            dst_k = mod(k,                  jn%n_legs) + 1
        case (V_RIGHT)
            dst_k = mod(k + jn%n_legs - 2,  jn%n_legs) + 1
        case default
            dst_k = k                                    ! invalid intent: self-loop, will be physical-blocked
        end select

        dst_road = jn%connected_road_ids(dst_k)
        dst_end  = jn%end_at(dst_k)
        dst_lane = outbound_lane_idx(dst_end)
    end subroutine compute_destination

    !-----------------------------------------------------------------
    ! Yield matrix: yields_to(i,j) = .true. means leg i must yield to
    ! leg j.  Combines R1 (yield-right, only if paths conflict) and R2
    ! (right-turn yields oncoming straight/left).
    !-----------------------------------------------------------------
    subroutine build_yield_matrix(jn, holding, yields_to)
        type(junction_t), intent(in)  :: jn
        integer,          intent(in)  :: holding(:)
        logical,          intent(out) :: yields_to(:,:)
        integer :: N, i, right_of_i, opposite_of_i

        N = jn%n_legs
        yields_to = .false.

        do i = 1, N
            if (holding(i) == V_EMPTY) cycle

            right_of_i = mod(i + N - 2, N) + 1                 ! clockwise previous

            ! R1: yield to the right, but only if paths actually conflict.
            if (holding(right_of_i) /= V_EMPTY .and. &
                paths_conflict(N, holding(i), holding(right_of_i), 3)) then
                yields_to(i, right_of_i) = .true.
            end if

            ! R2: right-turning vehicle yields to oncoming straight/left.
            if (N == 4 .and. holding(i) == V_RIGHT) then
                opposite_of_i = mod(i + 1, N) + 1              ! +N/2 wrap (only N=4)
                if (holding(opposite_of_i) == V_STRAIGHT .or. &
                    holding(opposite_of_i) == V_LEFT) then
                    yields_to(i, opposite_of_i) = .true.
                end if
            end if
        end do
    end subroutine build_yield_matrix

    !-----------------------------------------------------------------
    ! Path-crossing predicate for a 4-way junction.  delta = clockwise
    ! offset (j - i) mod N, so delta in {1, 2, 3}.
    ! (UK convention: left turn is the easy turn; right turn cuts across.)
    !-----------------------------------------------------------------
    pure function paths_conflict(N, ti, tj, delta) result(c)
        integer, intent(in) :: N, ti, tj, delta
        logical :: c

        c = .false.
        if (N /= 4) return                                     ! T-junction: separate logic later

        select case (delta)
        case (1)        ! j is one step clockwise from i (to i's left)
            c = .not. (ti == V_LEFT .and. tj == V_LEFT)
        case (2)        ! j is opposite i
            if      (ti == V_STRAIGHT .and. tj == V_STRAIGHT) then
                c = .false.                                    ! parallel through-traffic
            else if (ti == V_LEFT     .and. tj == V_LEFT)     then
                c = .false.                                    ! both swing UK-left, paths clear
            else
                c = .true.                                     ! everything else, incl. mutual right turn
            end if
        case (3)        ! j is one step counter-clockwise from i (to i's right)
            c = .not. (ti == V_LEFT .and. tj == V_LEFT)
        end select
    end function paths_conflict

    !-----------------------------------------------------------------
    ! Greedy approval pass: a candidate is approved if it has a vehicle,
    ! its destination is physically clear, and no live yield obligation
    ! points at another live candidate.
    !-----------------------------------------------------------------
    pure subroutine approve_pass(holding, phys_clear, yields_to, approved)
        integer, intent(in)    :: holding(:)
        logical, intent(in)    :: phys_clear(:)
        logical, intent(in)    :: yields_to(:,:)
        logical, intent(inout) :: approved(:)
        integer :: i, j, N
        logical :: yields_live

        N = size(holding)
        do i = 1, N
            if (approved(i))             cycle
            if (holding(i) == V_EMPTY)   cycle
            if (.not. phys_clear(i))     cycle

            yields_live = .false.
            do j = 1, N
                if (yields_to(i, j) .and. holding(j) /= V_EMPTY .and. .not. approved(j)) then
                    yields_live = .true.
                    exit
                end if
            end do
            if (.not. yields_live) approved(i) = .true.
        end do
    end subroutine approve_pass

    !-----------------------------------------------------------------
    ! Deadlock breaker: if any non-approved live candidate still has a
    ! yield obligation pointing at another non-approved live candidate,
    ! we are in a symmetric standoff.  Pick one stochastically, mark
    ! it approved, and clear its row/column from the yield matrix.
    !-----------------------------------------------------------------
    subroutine resolve_deadlocks(holding, phys_clear, yields_to, approved)
        integer, intent(in)    :: holding(:)
        logical, intent(in)    :: phys_clear(:)
        logical, intent(inout) :: yields_to(:,:), approved(:)
        integer :: N, k, n_locked, pick
        integer, allocatable :: locked(:)
        real    :: r
        logical :: still_yielding

        N = size(holding)
        allocate(locked(N))

        do
            n_locked = 0
            do k = 1, N
                if (approved(k))            cycle
                if (holding(k) == V_EMPTY)  cycle
                if (.not. phys_clear(k))    cycle
                still_yielding = any(yields_to(k, :) .and. holding /= V_EMPTY .and. .not. approved)
                if (still_yielding) then
                    n_locked = n_locked + 1
                    locked(n_locked) = k
                end if
            end do
            if (n_locked == 0) exit

            call random_number(r)
            pick = locked(min(int(r * n_locked) + 1, n_locked))
            approved(pick)     = .true.
            yields_to(pick, :) = .false.
            yields_to(:, pick) = .false.

            ! After the tie-break, re-run the regular approval pass: clearing
            ! pick's column may unblock vehicles that were yielding only to it.
            call approve_pass(holding, phys_clear, yields_to, approved)
        end do

        deallocate(locked)
    end subroutine resolve_deadlocks

end module junction_mod
