module junction_mod
    ! Junction evaluation and right-of-way rules for road-network simulations.
    !
    ! This module evaluates vehicles waiting at junction holding cells and moves
    ! approved vehicles into their selected outbound roads. A holding cell is the
    ! final site of an inbound lane. A destination cell is the first site of an
    ! outbound lane.
    !
    ! The update is based on an old-state snapshot. This means all candidate
    ! movements are evaluated from the same previous network state before any
    ! approved moves are committed.
    !
    ! The right-of-way hierarchy is:
    !
    ! 1. A vehicle cannot move into an occupied destination cell.
    ! 2. A vehicle yields to the vehicle on its right when their paths conflict or
    !    when both vehicles target the same outbound road.
    ! 3. In a four-way junction, a right-turning vehicle yields to oncoming
    !    straight or left-turning vehicles.
    ! 4. Remaining opposite-leg conflicts are resolved using category-specific
    !    priority rules or stochastic deadlock breaking.
    ! 5. Non-conflicting moves may proceed in parallel.
    !
    ! Move categories are derived from the inbound leg index, outbound leg index
    ! and number of junction legs. Defining
    !
    ! .. math::
    !
    !    d = (i_{\mathrm{out}} - i_{\mathrm{in}}) \bmod n,
    !
    ! the categories are:
    !
    ! * ``0``: U-turn, ``d = 0``;
    ! * ``1``: left turn, ``d = 1``;
    ! * ``2``: straight, ``d = n/2`` for four-way junctions;
    ! * ``3``: right turn, ``d = n - 1``.
    !
    ! Yield semantics are full-step semantics: if a vehicle must yield to another
    ! candidate, it waits for the entire time step even if the priority vehicle
    ! also moves during that step.

    use vehicle_mod
    use road_network_mod

    implicit none
    private :: evaluate_one_junction, paths_conflict_sym, build_yield_matrix_v2, &
    build_yield_matrix_asym, chords_cross, in_arc, approve_pass, resolve_deadlocks

    public :: evaluate_junctions, move_category

contains

    subroutine evaluate_junctions(net)
        ! Evaluate all junctions in a road network.
        !
        ! This routine loops over every junction stored in ``net`` and applies the
        ! junction update. Each junction is evaluated using the network old-state
        ! arrays and then approved moves are committed to the live cell arrays.
        type(road_network_t), intent(inout) :: net
        ! Road network containing roads, lanes, cells and junction definitions.

        integer :: j
        ! Junction index.

        do j = 1, size(net%junctions)
            call evaluate_one_junction(net, j)
        end do
    end subroutine evaluate_junctions


    subroutine evaluate_one_junction(net, jid)
        ! Evaluate one junction and commit approved vehicle movements.
        !
        ! The routine performs the complete junction decision process for a single
        ! junction:
        !
        ! 1. read the occupied holding cells from the old network snapshot;
        ! 2. sample a destination outbound road for each occupied inbound leg;
        ! 3. reject candidates whose destination cell is physically occupied;
        ! 4. build a yield matrix representing right-of-way conflicts;
        ! 5. approve candidates that are not blocked by physical occupancy or
        !    live yield targets;
        ! 6. resolve symmetric deadlocks stochastically;
        ! 7. commit approved movements.
        !
        ! The yield matrix uses the convention:
        !
        ! .. math::
        !
        !    Y_{ij} =
        !    \begin{cases}
        !    1, & \text{candidate } i \text{ yields to candidate } j, \\
        !    0, & \text{otherwise}.
        !    \end{cases}
        type(road_network_t), intent(inout) :: net
        ! Road network containing the junction and road states.
        integer, intent(in) :: jid
        ! Index of the junction to evaluate.

        integer, parameter :: MAX_LEGS = 4
        ! Maximum number of inbound legs currently supported by this routine.

        integer :: N
        ! Number of inbound legs at the junction.
        integer :: k
        ! Candidate inbound-leg index.
        integer :: m
        ! Outbound-route index used during destination sampling.
        integer :: src_road
        ! Source road index for a candidate.
        integer :: src_lane
        ! Source lane index for a candidate.
        integer :: L_src
        ! Length of the source lane.
        integer :: dst_road(MAX_LEGS)
        ! Selected destination road for each candidate.
        integer :: dst_lane(MAX_LEGS)
        ! Selected destination lane for each candidate.
        integer :: dst_idx(MAX_LEGS)
        ! Selected outbound index within the junction definition.
        integer :: holding(MAX_LEGS)
        ! Holding-cell occupancy for each inbound leg.
        logical :: phys_clear(MAX_LEGS)
        ! Whether the destination cell is physically empty for each candidate.
        logical :: yields_to(MAX_LEGS, MAX_LEGS)
        ! Yield matrix. ``yields_to(i,j)`` is true if candidate ``i`` yields to ``j``.
        logical :: approved(MAX_LEGS)
        ! Approval flag for each candidate.
        real :: r
        ! Uniform random number used for route selection.
        real :: cumsum
        ! Cumulative probability used when sampling outbound routes.

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

        ! Read holding cells from the old snapshot.
        do k = 1, N
            src_road = net%junctions(jid)%in_road(k)
            src_lane = net%junctions(jid)%in_lane(k)
            L_src    = net%roads(src_road)%lane(src_lane)%length
            if (net%roads(src_road)%lane(src_lane)%old(L_src)%has_car) then
                holding(k) = V_OCCUPIED
            else
                holding(k) = V_EMPTY
            end if
        end do

        ! Sample a destination route for each occupied inbound leg.
        do k = 1, N
            if (holding(k) == 0) cycle

            call random_number(r)
            cumsum    = 0.0
            dst_idx(k) = net%junctions(jid)%n_out

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

        ! Check whether each destination site is physically clear.
        do k = 1, N
            if (holding(k) == V_EMPTY) cycle
            phys_clear(k) = (.not. net%roads(dst_road(k))%lane(dst_lane(k))%old(1)%has_car)
        end do

        ! Build the right-of-way yield matrix.
        if (N == 4 .and. net%junctions(jid)%n_out == 4) then
            call build_yield_matrix_v2(net%junctions(jid), holding, dst_idx, yields_to)
        else
            call build_yield_matrix_asym(net%junctions(jid), holding, dst_idx, yields_to)
        end if

        call approve_pass(holding, phys_clear, yields_to, approved)
        call resolve_deadlocks(holding, phys_clear, yields_to, approved)

        ! Commit approved moves to the live cell arrays.
        do k = 1, N
            if (.not. approved(k)) cycle

            src_road = net%junctions(jid)%in_road(k)
            src_lane = net%junctions(jid)%in_lane(k)
            L_src    = net%roads(src_road)%lane(src_lane)%length

            net%roads(dst_road(k))%lane(dst_lane(k))%cells(1)%has_car = .true.
            net%roads(dst_road(k))%lane(dst_lane(k))%cells(1)%velocity = 0

            net%roads(src_road)%lane(src_lane)%cells(L_src)%has_car   = .false.
            net%roads(src_road)%lane(src_lane)%cells(L_src)%velocity = 0
        end do
    end subroutine evaluate_one_junction


    pure function move_category(in_idx, out_idx, n) result(cat)
        ! Return the movement category for a junction candidate.
        !
        ! The category is derived only from the inbound leg index, outbound leg
        ! index and junction size. The relative displacement is:
        !
        ! .. math::
        !
        !    d = (i_{\mathrm{out}} - i_{\mathrm{in}}) \bmod n.
        !
        ! The returned category is:
        !
        ! * ``0`` for a U-turn;
        ! * ``1`` for a left turn;
        ! * ``2`` for straight ahead in a four-way junction;
        ! * ``3`` for a right turn;
        ! * ``-1`` if the movement does not match one of these categories.
        integer, intent(in) :: in_idx
        ! Inbound leg index.
        integer, intent(in) :: out_idx
        ! Outbound leg index.
        integer, intent(in) :: n
        ! Number of legs in the symmetric junction.
        integer :: cat
        ! Movement category.
        integer :: d
        ! Relative outbound displacement modulo ``n``.

        d = modulo(out_idx - in_idx, n)

        if      (d == 0)                then
            cat = 0
        else if (d == 1)                then
            cat = 1
        else if (n >= 4 .and. d == n/2) then
            cat = 2
        else if (d == n - 1)            then
            cat = 3
        else
            cat = -1
        end if
    end function move_category


    pure function paths_conflict_sym(n, cat_i, cat_j, delta) result(c)
        ! Return whether two movements conflict in a symmetric four-way junction.
        !
        ! This predicate is used by the category-based four-way junction rules.
        ! For adjacent legs, movements conflict unless both vehicles are taking
        ! the easy left turn. For opposite legs, straight-straight and left-left
        ! movements are treated as non-conflicting.
        integer, intent(in) :: n
        ! Number of junction legs.
        integer, intent(in) :: cat_i
        ! Movement category of the first candidate.
        integer, intent(in) :: cat_j
        ! Movement category of the second candidate.
        integer, intent(in) :: delta
        ! Relative leg separation between the candidates.
        logical :: c
        ! True if the two movement paths conflict.

        c = .false.
        if (n /= 4) return

        select case (delta)
        case (1, 3)
            c = .not. (cat_i == 1 .and. cat_j == 1)
        case (2)
            c = .not. (cat_i == 2 .and. cat_j == 2) &
                .and. .not. (cat_i == 1 .and. cat_j == 1)
        end select
    end function paths_conflict_sym


    subroutine build_yield_matrix_v2(jn, holding, dst_idx, yields_to)
        ! Build the category-based yield matrix for a four-way junction.
        !
        ! This routine implements the main right-of-way rules for symmetric
        ! four-way junctions.
        !
        ! Rule R1 applies yielding to the right. Candidate ``i`` yields to the
        ! clockwise-previous leg if their paths conflict or if both candidates
        ! target the same outbound.
        !
        ! Rule R2 makes right-turning vehicles yield to oncoming straight or
        ! left-turning vehicles.
        !
        ! Rule R3 handles remaining opposite-leg conflicts. Straight has priority
        ! over left in a straight-left conflict. Other unresolved conflicts are
        ! represented as mutual yields and are later passed to the stochastic
        ! deadlock resolver.
        type(junction_t), intent(in) :: jn
        ! Junction definition.
        integer, intent(in) :: holding(:)
        ! Holding-cell occupancy for each inbound leg.
        integer, intent(in) :: dst_idx(:)
        ! Selected outbound index for each inbound leg.
        logical, intent(out) :: yields_to(:,:)
        ! Output yield matrix.

        integer :: N
        ! Number of inbound legs.
        integer :: i
        ! Inbound candidate index.
        integer :: right_of_i
        ! Candidate immediately to the right of candidate ``i``.
        integer :: opp
        ! Candidate opposite candidate ``i``.
        integer :: cat_i
        ! Movement category of candidate ``i``.
        integer :: cat_opp
        ! Movement category of the opposite candidate.
        integer :: cat_right
        ! Movement category of the candidate to the right.

        N = jn%n_in
        yields_to = .false.

        do i = 1, N
            if (holding(i) == V_EMPTY) cycle

            cat_i      = move_category(i, dst_idx(i), N)
            right_of_i = mod(i + N - 2, N) + 1

            ! R1: yield to the right when paths conflict or the destination is shared.
            if (holding(right_of_i) /= V_EMPTY) then
                cat_right = move_category(right_of_i, dst_idx(right_of_i), N)
                if (paths_conflict_sym(N, cat_i, cat_right, 3) .or. &
                    dst_idx(i) == dst_idx(right_of_i)) then
                    yields_to(i, right_of_i) = .true.
                end if
            end if

            ! R2: right-turn yields to oncoming straight or left.
            if (N == 4 .and. jn%n_out == 4 .and. cat_i == 3) then
                opp = mod(i + 1, N) + 1
                if (holding(opp) /= V_EMPTY) then
                    cat_opp = move_category(opp, dst_idx(opp), N)
                    if (cat_opp == 1 .or. cat_opp == 2) yields_to(i, opp) = .true.
                end if
            end if
        end do

        ! R3: remaining opposite-leg conflicts.
        if (N == 4 .and. jn%n_out == 4) then
            do i = 1, N
                if (holding(i) == V_EMPTY) cycle

                opp = mod(i + 1, N) + 1
                if (holding(opp) == V_EMPTY) cycle

                cat_i   = move_category(i,   dst_idx(i),   N)
                cat_opp = move_category(opp, dst_idx(opp), N)

                if (.not. (paths_conflict_sym(N, cat_i, cat_opp, 2) .or. &
                           dst_idx(i) == dst_idx(opp))) cycle

                if (yields_to(i, opp) .or. yields_to(opp, i)) cycle

                select case (cat_i * 4 + cat_opp)
                case (2*4 + 1)
                    yields_to(opp, i) = .true.
                case (1*4 + 2)
                    yields_to(i, opp) = .true.
                case default
                    yields_to(i,   opp) = .true.
                    yields_to(opp, i  ) = .true.
                end select
            end do
        end if
    end subroutine build_yield_matrix_v2


    subroutine build_yield_matrix_asym(jn, holding, dst_idx, yields_to)
        ! Build the yield matrix for an asymmetric junction.
        !
        ! This routine is used when the junction is not a symmetric four-way
        ! junction. Only the yield-to-the-right rule is applied. Candidate ``i``
        ! yields to the clockwise-previous candidate if the corresponding path
        ! chords cross or if both candidates target the same outbound.
        type(junction_t), intent(in) :: jn
        ! Junction definition.
        integer, intent(in) :: holding(:)
        ! Holding-cell occupancy for each inbound leg.
        integer, intent(in) :: dst_idx(:)
        ! Selected outbound index for each inbound leg.
        logical, intent(out) :: yields_to(:,:)
        ! Output yield matrix.

        integer :: N
        ! Number of inbound legs.
        integer :: i
        ! Inbound candidate index.
        integer :: right_of_i
        ! Candidate immediately to the right of candidate ``i``.
        integer :: n_perim
        ! Number of perimeter ports used by the chord-crossing test.

        N       = jn%n_in
        n_perim = jn%n_in + jn%n_out
        yields_to = .false.

        do i = 1, N
            if (holding(i) == V_EMPTY) cycle

            right_of_i = mod(i + N - 2, N) + 1
            if (holding(right_of_i) == V_EMPTY) cycle

            if (chords_cross(jn%in_perim(i),          jn%out_perim(dst_idx(i)), &
                             jn%in_perim(right_of_i), jn%out_perim(dst_idx(right_of_i)), &
                             n_perim) &
                .or. dst_idx(i) == dst_idx(right_of_i)) then
                yields_to(i, right_of_i) = .true.
            end if
        end do
    end subroutine build_yield_matrix_asym


    pure function chords_cross(pa1, pa2, pb1, pb2, n_perim) result(c)
        ! Return whether two route chords cross on a cyclic junction perimeter.
        !
        ! Ports are treated as points on a cyclic perimeter. The two chords
        ! ``pa1 -> pa2`` and ``pb1 -> pb2`` cross if exactly one of ``pb1`` and
        ! ``pb2`` lies inside the open clockwise arc from ``pa1`` to ``pa2``:
        !
        ! .. math::
        !
        !    C =
        !    I(pa_1, pa_2, pb_1)
        !    \;\mathrm{xor}\;
        !    I(pa_1, pa_2, pb_2).
        integer, intent(in) :: pa1
        ! Start port of the first chord.
        integer, intent(in) :: pa2
        ! End port of the first chord.
        integer, intent(in) :: pb1
        ! Start port of the second chord.
        integer, intent(in) :: pb2
        ! End port of the second chord.
        integer, intent(in) :: n_perim
        ! Number of ports on the cyclic perimeter.
        logical :: c
        ! True if the two chords cross.

        c = in_arc(pa1, pa2, pb1, n_perim) .neqv. in_arc(pa1, pa2, pb2, n_perim)
    end function chords_cross


    pure function in_arc(a, b, p, n) result(res)
        ! Return whether a port lies in an open clockwise arc.
        !
        ! This helper tests whether port ``p`` lies strictly inside the open
        ! clockwise arc from port ``a`` to port ``b``. The endpoints are excluded.
        !
        ! .. math::
        !
        !    0 < (p - a) \bmod n < (b - a) \bmod n.
        integer, intent(in) :: a
        ! Start port of the open arc.
        integer, intent(in) :: b
        ! End port of the open arc.
        integer, intent(in) :: p
        ! Port to test.
        integer, intent(in) :: n
        ! Number of ports on the cyclic perimeter.
        logical :: res
        ! True if ``p`` lies strictly inside the open clockwise arc.
        integer :: ab
        ! Arc length from ``a`` to ``b`` modulo ``n``.
        integer :: ap
        ! Arc length from ``a`` to ``p`` modulo ``n``.

        if (a == b) then
            res = .false.
        else
            ab = modulo(b - a, n)
            ap = modulo(p - a, n)
            res = (ap > 0) .and. (ap < ab)
        end if
    end function in_arc


    pure subroutine approve_pass(holding, phys_clear, yields_to, approved)
        ! Approve candidates that are physically clear and not yield-blocked.
        !
        ! A candidate is approved if:
        !
        ! * it is occupied;
        ! * its destination cell is physically clear;
        ! * it is not required to yield to any occupied holding cell.
        !
        ! The full-step yield rule is:
        !
        ! .. math::
        !
        !    \mathrm{blocked}_i
        !    =
        !    \exists j:
        !    Y_{ij} = 1
        !    \land
        !    H_j = 1.
        !
        ! If ``blocked_i`` is true, candidate ``i`` waits for the entire time
        ! step, regardless of whether candidate ``j`` is also approved this step.
        integer, intent(in) :: holding(:)
        ! Holding-cell occupancy for each candidate.
        logical, intent(in) :: phys_clear(:)
        ! Physical-clear flag for each candidate destination.
        logical, intent(in) :: yields_to(:,:)
        ! Yield matrix.
        logical, intent(inout) :: approved(:)
        ! Approval flags updated in place.

        integer :: i
        ! Candidate index.
        integer :: N
        ! Number of candidates.
        logical :: yields_live
        ! True if the candidate yields to an occupied holding cell.

        N = size(holding)

        do i = 1, N
            if (approved(i))           cycle
            if (holding(i) == V_EMPTY) cycle
            if (.not. phys_clear(i))   cycle

            yields_live = any(yields_to(i, :) .and. (holding /= V_EMPTY))
            if (.not. yields_live) approved(i) = .true.
        end do
    end subroutine approve_pass


    subroutine resolve_deadlocks(holding, phys_clear, yields_to, approved)
        ! Resolve symmetric yield deadlocks stochastically.
        !
        ! A candidate is considered part of a deadlock if it is blocked only by
        ! other non-approved and non-deadlocked candidates. Vehicles blocked by an
        ! already-approved priority vehicle are simply waiting under the full-step
        ! yield rule and are not included in the deadlock set.
        !
        ! When a deadlock set is found, one member is selected uniformly at random
        ! and approved. All other members of that deadlock set are marked as
        ! deadlocked for the current junction step, so they wait until the next
        ! network update.
        !
        ! If a deadlock set contains ``n`` candidates, the probability that a
        ! given candidate is selected as the winner is:
        !
        ! .. math::
        !
        !    P(\mathrm{winner}) = \frac{1}{n}.
        integer, intent(in) :: holding(:)
        ! Holding-cell occupancy for each candidate.
        logical, intent(in) :: phys_clear(:)
        ! Physical-clear flag for each candidate destination.
        logical, intent(in) :: yields_to(:,:)
        ! Yield matrix.
        logical, intent(inout) :: approved(:)
        ! Approval flags updated in place.

        integer :: N
        ! Number of candidates.
        integer :: k
        ! Candidate index.
        integer :: idx
        ! Index into the temporary locked-candidate list.
        integer :: n_locked
        ! Number of candidates in the current deadlock set.
        integer :: pick
        ! Candidate selected as the stochastic deadlock winner.
        integer, allocatable :: locked(:)
        ! Temporary list of candidates in the current deadlock set.
        logical :: deadlocked(size(holding))
        ! Flags for candidates already marked as deadlock losers.
        real :: r
        ! Uniform random number used to choose a deadlock winner.
        logical :: in_deadlock
        ! True if the current candidate is in the current deadlock set.

        N = size(holding)
        deadlocked = .false.
        allocate(locked(N))

        do
            n_locked = 0

            do k = 1, N
                if (approved(k))             cycle
                if (deadlocked(k))           cycle
                if (holding(k) == V_EMPTY)   cycle
                if (.not. phys_clear(k))     cycle

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

            call random_number(r)
            pick = locked(min(int(r * n_locked) + 1, n_locked))
            approved(pick) = .true.

            do idx = 1, n_locked
                if (locked(idx) /= pick) deadlocked(locked(idx)) = .true.
            end do
        end do

        deallocate(locked)
    end subroutine resolve_deadlocks

end module junction_mod