module junction_mod
    !--------------------------------------------------------------------
    ! Junction evaluation: moves vehicles from holding cells (site L of
    ! each leg's inbound lane) into destination cells (site 1 of the
    ! relevant outbound lane), respecting the right-of-way hierarchy:
    !
    !   1. Physical block: never advance into an occupied destination cell.
    !   2. Yield to the right (R1): yield to the vehicle one step
    !      counter-clockwise (index i-1 wrap), when paths cross.
    !   3. Right-turn yields oncoming (R2): a V_RIGHT vehicle yields to
    !      any V_STRAIGHT or V_LEFT vehicle from the opposite leg.
    !   4. Opposite-leg conflict scan (R3): for remaining delta=2 path
    !      crossings not resolved by R2 —
    !        (S, L) opposite: straight has priority, left yields (UK rule)
    !        (R, R) opposite: mutual yields -> stochastic deadlock break
    !   5. Parallel non-conflicting moves are allowed.
    !   6. Yield semantics: if you yield, you wait the FULL step regardless
    !      of whether your yield-target also moves this step.
    !   7. Symmetric standoff: stochastic tie-break selects one winner;
    !      remaining deadlocked candidates wait for the next timestep.
    !
    ! Destination formula (pure modulo on clockwise leg index k, 1..N):
    !   STRAIGHT  ->  mod(k + N/2 - 1, N) + 1
    !   LEFT      ->  mod(k,           N) + 1   (next clockwise; UK easy turn)
    !   RIGHT     ->  mod(k + N - 2,   N) + 1   (prev clockwise; cuts across)
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
    subroutine evaluate_one_junction(net, jid)
        type(road_network_t), intent(inout) :: net
        integer, intent(in) :: jid

        integer, parameter :: MAX_LEGS = 4

        integer :: N, k, src_road, src_lane, src_end, L_src
        integer :: dst_road(MAX_LEGS), dst_lane(MAX_LEGS), dst_end_arr(MAX_LEGS)
        integer :: holding(MAX_LEGS)
        logical :: phys_clear(MAX_LEGS), yields_to(MAX_LEGS, MAX_LEGS), approved(MAX_LEGS)

        N = net%junctions(jid)%n_legs
        if (N > MAX_LEGS) then
            print *, 'evaluate_one_junction: N exceeds MAX_LEGS'
            stop 1
        end if

        holding    = V_EMPTY
        dst_road   = 0
        dst_lane   = 0
        dst_end_arr = 0
        phys_clear = .false.
        yields_to  = .false.
        approved   = .false.

        ! 1. Read holding cells from the old snapshot.
        do k = 1, N
            src_road = net%junctions(jid)%connected_road_ids(k)
            src_end  = net%junctions(jid)%end_at(k)
            src_lane = inbound_lane_at_end(net%roads(src_road), src_end)
            L_src    = net%roads(src_road)%lane(src_lane)%length
            if (net%roads(src_road)%lane(src_lane)%old(L_src)%has_car) then
                holding(k) = net%roads(src_road)%lane(src_lane)%old(L_src)%turning_intent
            else
                holding(k) = V_EMPTY
            end if
        end do

        ! 2. Resolve destinations.
        do k = 1, N
            if (holding(k) == 0) cycle
            call compute_destination(net%junctions(jid), k, holding(k), &
                                     dst_road(k), dst_lane(k), dst_end_arr(k))
        end do

        ! 3. Physical-clear check (destination site 1 in the old snapshot).
        do k = 1, N
            if (holding(k) == V_EMPTY) cycle
            phys_clear(k) = (.not. net%roads(dst_road(k))%lane(dst_lane(k))%old(1)%has_car)
        end do

        ! 4. Build yield matrix (rules R1, R2, R3).
        call build_yield_matrix(net%junctions(jid), holding, yields_to)

        ! 5. First approval pass.
        call approve_pass(holding, phys_clear, yields_to, approved)

        ! 6. Stochastic deadlock break for remaining symmetric standoffs.
        call resolve_deadlocks(holding, phys_clear, yields_to, approved)

        ! 7. Commit approved moves.
        do k = 1, N
            if (.not. approved(k)) cycle
            src_road = net%junctions(jid)%connected_road_ids(k)
            src_end  = net%junctions(jid)%end_at(k)
            src_lane = inbound_lane_at_end(net%roads(src_road), src_end)
            L_src    = net%roads(src_road)%lane(src_lane)%length
            net%roads(dst_road(k))%lane(dst_lane(k))%cells(1)%has_car = .true.
            net%roads(dst_road(k))%lane(dst_lane(k))%cells(1)%turning_intent = holding(k)
            net%roads(dst_road(k))%lane(dst_lane(k))%cells(1)%velocity = 0
            net%roads(src_road)%lane(src_lane)%cells(L_src)%has_car   = .false.
            net%roads(src_road)%lane(src_lane)%cells(L_src)%turning_intent = V_EMPTY
            net%roads(src_road)%lane(src_lane)%cells(L_src)%velocity = 0
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
            dst_k = mod(k,                   jn%n_legs) + 1
        case (V_RIGHT)
            dst_k = mod(k + jn%n_legs - 2,   jn%n_legs) + 1
        case default
            dst_k = k                   ! invalid intent: self-loop, physical block will stop it
        end select

        dst_road = jn%connected_road_ids(dst_k)
        dst_end  = jn%end_at(dst_k)
        dst_lane = dst_end          ! outbound lane index == end index for 2-lane roads (Phase 1)
    end subroutine compute_destination

    !-----------------------------------------------------------------
    ! Yield matrix construction.
    !-----------------------------------------------------------------
    subroutine build_yield_matrix(jn, holding, yields_to)
        type(junction_t), intent(in)  :: jn
        integer,          intent(in)  :: holding(:)
        logical,          intent(out) :: yields_to(:,:)
        integer :: N, i, right_of_i, opp

        N = jn%n_legs
        yields_to = .false.

        do i = 1, N
            if (holding(i) == V_EMPTY) cycle

            right_of_i = mod(i + N - 2, N) + 1                 ! clockwise previous (delta=3 from i)

            ! R1 — yield to the right when paths conflict.
            if (holding(right_of_i) /= V_EMPTY .and. &
                paths_conflict(N, holding(i), holding(right_of_i), 3)) then
                yields_to(i, right_of_i) = .true.
            end if

            ! R2 — right-turn yields to oncoming straight or left.
            if (N == 4 .and. holding(i) == V_RIGHT) then
                opp = mod(i + 1, N) + 1                         ! opposite leg (delta=2, N=4)
                if (holding(opp) == V_STRAIGHT .or. holding(opp) == V_LEFT) then
                    yields_to(i, opp) = .true.
                end if
            end if
        end do

        ! R3 — scan remaining delta=2 path crossings not resolved by R2.
        ! For N=4 only; T-junctions (N=3) have no opposite leg.
        if (N == 4) then
            do i = 1, N
                if (holding(i) == V_EMPTY) cycle
                opp = mod(i + 1, N) + 1
                if (holding(opp) == V_EMPTY) cycle
                if (.not. paths_conflict(N, holding(i), holding(opp), 2)) cycle
                if (yields_to(i, opp) .or. yields_to(opp, i)) cycle   ! already resolved

                select case (holding(i) * 4 + holding(opp))         ! encode (ti, tj) pair
                case (V_STRAIGHT*4 + V_LEFT)
                    yields_to(opp, i) = .true.    ! UK: left yields to oncoming straight
                case (V_LEFT*4 + V_STRAIGHT)
                    yields_to(i, opp) = .true.    ! UK: left yields to oncoming straight
                case default
                    ! (R,R): mutual yields -> symmetric standoff, resolved by tie-break.
                    yields_to(i,   opp) = .true.
                    yields_to(opp, i  ) = .true.
                end select
            end do
        end if
    end subroutine build_yield_matrix

    !-----------------------------------------------------------------
    ! Path-crossing predicate for a 4-way junction.
    ! delta = clockwise offset mod N from leg i to leg j.
    ! UK convention: left turn is the easy swing; right turn cuts across.
    !-----------------------------------------------------------------
    pure function paths_conflict(N, ti, tj, delta) result(c)
        integer, intent(in) :: N, ti, tj, delta
        logical :: c

        c = .false.
        if (N /= 4) return

        select case (delta)
        case (1)        ! j is one step clockwise from i (to i's left)
            c = .not. (ti == V_LEFT .and. tj == V_LEFT)
        case (2)        ! j is opposite i
            if      (ti == V_STRAIGHT .and. tj == V_STRAIGHT) then
                c = .false.                        ! parallel through-traffic
            else if (ti == V_LEFT     .and. tj == V_LEFT)     then
                c = .false.                        ! both swing UK-left, no crossing
            else
                c = .true.                         ! all other opposite combinations cross
            end if
        case (3)        ! j is one step counter-clockwise from i (to i's right)
            c = .not. (ti == V_LEFT .and. tj == V_LEFT)
        end select
    end function paths_conflict

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
