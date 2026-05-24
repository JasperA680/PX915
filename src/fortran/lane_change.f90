module lane_change_mod
    !--------------------------------------------------------------------
    ! Lane-changing sub-step for the multi-lane Nagel-Schreckenberg lc_model.
    !
    ! Based on: Rickert, Nagel, Schreckenberg, Latour (1996),
    !           "Two lane traffic simulations using cellular automata",
    !           Physica A 231, 534-550.
    !
    ! Two lc_models are provided via the `lc_model` integer parameter:
    !   LC_SYMMETRIC  (0) - both lanes treated equally; vehicle changes
    !                       only when blocked (T1 required for any move).
    !   LC_ASYMMETRIC (1) - vehicles always prefer the rightmost lane
    !                       (lowest index); they may return right without
    !                       being blocked, but still need T1 to move left.
    !
    ! Lane change conditions for a vehicle at position i, velocity v:
    !   T1: gap_fwd(own lane, i) < v+1         [own lane is blocked]
    !   T2: gap_other_fwd(target lane, i) > v+1 [enough space ahead]
    !   T3: gap_other_back(target lane, i) > L_BACK [enough space behind]
    !   T4: rand() < p_change                   [stochastic acceptance]
    !
    ! Lane ordering convention: lane index 1 = rightmost, higher = more left.
    ! Only lanes with the same flow_direction within the same road can
    ! exchange vehicles.  Each vehicle may move at most one lane per step.
    ! Site L (the junction holding cell / open-outflow boundary) is excluded
    ! from lane-change participation.
    !
    ! All decisions read from the pre-snapshotted `old` arrays (strict
    ! parallel update).  Results are written to the `cells` arrays.
    ! The caller must call snapshot_network() before and after this step.
    !--------------------------------------------------------------------
    use road_network_mod
    implicit none
    private

    integer, parameter, public :: LC_DISABLED   = -1
    integer, parameter, public :: LC_SYMMETRIC  =  0
    integer, parameter, public :: LC_ASYMMETRIC =  1

    ! Backward look-ahead distance = V_MAX: a following car at maximum
    ! velocity needs at least L_BACK+1 clear cells to avoid a collision.

    public :: apply_lane_changes

contains

    !------------------------------------------------------------------
    ! Apply one lane-change sub-step to every road in the network.
    ! Reads from lane%old (caller must have snapshotted), writes cells.
    !------------------------------------------------------------------
    subroutine apply_lane_changes(net, lc_model, p_change, v_max)
        type(road_network_t), intent(inout) :: net
        integer, intent(in) :: lc_model, v_max
        real,    intent(in) :: p_change
        integer :: r
        do r = 1, size(net%roads)
            call road_lane_change(net%roads(r), lc_model, p_change, v_max)
        end do
    end subroutine apply_lane_changes

    !------------------------------------------------------------------
    ! Apply lane changes within one road.
    ! Iterates over all lanes and positions (excluding site L).
    ! Uses old state for all gap queries; writes to cells.
    !------------------------------------------------------------------
    subroutine road_lane_change(road, lc_model, p_change, v_max)
        type(road_t), intent(inout) :: road
        integer, intent(in) :: lc_model, v_max
        real,    intent(in) :: p_change

        integer :: n_lanes, k, L, i, v, tgt
        integer :: gf, go_fwd, go_back, L_BACK
        logical :: t1, t2r, t3r, t2l, t3l, want_right, want_left
        real    :: rnd

        L_BACK = v_max

        n_lanes = size(road%lane)
        if (n_lanes < 2) return

        do k = 1, n_lanes
            L = road%lane(k)%length
            do i = 1, L - 1   ! site L (holding cell / outflow) excluded
                if (.not. road%lane(k)%old(i)%has_car) cycle

                v   = road%lane(k)%old(i)%velocity
                tgt = 0

                ! T1: is the vehicle blocked in its own lane?
                gf = fwd_gap(road%lane(k)%old, i, L)
                t1 = (gf < v + 1)

                ! --- right neighbor: lane k-1 (toward rightmost) ---
                want_right = .false.
                if (k > 1 .and. same_dir(road, k, k-1)) then
                    go_fwd  = other_fwd_gap(road%lane(k-1)%old, i, L)
                    go_back = other_back_gap(road%lane(k-1)%old, i)
                    t2r = (go_fwd  > v + 1)                                 ! T2 - Is this vehicle going fast enough to cover the gap infront in the next step?
                    t3r = (go_back > L_BACK)                                ! T3 - Is the vehicle behind in the target lane going slow enough to allow a merge?
                    want_right = t2r .and. t3r
                end if

                ! --- left neighbor: lane k+1 (away from rightmost) ---
                want_left = .false.
                if (k < n_lanes .and. same_dir(road, k, k+1)) then
                    go_fwd  = other_fwd_gap(road%lane(k+1)%old, i, L)
                    go_back = other_back_gap(road%lane(k+1)%old, i)
                    t2l = (go_fwd  > v + 1)
                    t3l = (go_back > L_BACK)
                    want_left = t2l .and. t3l
                end if

                if (lc_model == LC_ASYMMETRIC) then
                    ! Rightward return does not require T1.
                    ! Leftward pass requires T1 (blocked on current lane).
                    if (want_right) then
                        call random_number(rnd)
                        if (rnd < p_change) tgt = k - 1
                    end if
                    if (tgt == 0 .and. want_left .and. t1) then
                        call random_number(rnd)
                        if (rnd < p_change) tgt = k + 1
                    end if
                else
                    ! LC_SYMMETRIC: T1 required for any lane change.
                    if (.not. t1) cycle
                    if (want_right .and. want_left) then
                        call random_number(rnd)
                        if (rnd < p_change) then
                            call random_number(rnd)
                            tgt = merge(k-1, k+1, rnd < 0.5)
                        end if
                    else if (want_right) then
                        call random_number(rnd)
                        if (rnd < p_change) tgt = k - 1
                    else if (want_left) then
                        call random_number(rnd)
                        if (rnd < p_change) tgt = k + 1
                    end if
                end if

                if (tgt /= 0) then
                    road%lane(k)%cells(i)%has_car    = .false.
                    road%lane(k)%cells(i)%velocity   = 0
                    road%lane(tgt)%cells(i)%has_car  = .true.
                    road%lane(tgt)%cells(i)%velocity = v
                end if
            end do
        end do
    end subroutine road_lane_change

    !------------------------------------------------------------------
    ! True iff lanes a and b within the same road share flow_direction.
    !------------------------------------------------------------------
    pure logical function same_dir(road, a, b)
        type(road_t), intent(in) :: road
        integer,      intent(in) :: a, b
        same_dir = (road%lane(a)%flow_direction == road%lane(b)%flow_direction)
    end function same_dir

    !------------------------------------------------------------------
    ! Forward gap in lane old state from position i toward site L.
    ! Returns (distance to next car - 1), or (L - i) if none.
    !------------------------------------------------------------------
    integer function fwd_gap(lane_old, i, L) result(gap)
        type(cell), intent(in) :: lane_old(:)
        integer,    intent(in) :: i, L
        integer :: j
        gap = L - i
        do j = i + 1, L
            if (lane_old(j)%has_car) then
                gap = j - i - 1
                return
            end if
        end do
    end function fwd_gap

    !------------------------------------------------------------------
    ! Forward gap in an adjacent lane old state at position i.
    ! Returns -1 if position i itself is occupied (cannot change lanes).
    !------------------------------------------------------------------
    integer function other_fwd_gap(other_old, i, L) result(gap)
        type(cell), intent(in) :: other_old(:)
        integer,    intent(in) :: i, L
        integer :: j
        if (other_old(i)%has_car) then
            gap = -1
            return
        end if
        gap = L - i
        do j = i + 1, L
            if (other_old(j)%has_car) then
                gap = j - i - 1
                return
            end if
        end do
    end function other_fwd_gap

    !------------------------------------------------------------------
    ! Backward gap in an adjacent lane old state at position i.
    ! Returns -1 if position i itself is occupied.
    ! Otherwise: (i - position_of_nearest_car_behind - 1), or (i-1) if none.
    !------------------------------------------------------------------
    integer function other_back_gap(other_old, i) result(gap)
        type(cell), intent(in) :: other_old(:)
        integer,    intent(in) :: i
        integer :: j
        if (other_old(i)%has_car) then
            gap = -1
            return
        end if
        gap = i - 1
        do j = i - 1, 1, -1
            if (other_old(j)%has_car) then
                gap = i - j - 1
                return
            end if
        end do
    end function other_back_gap

end module lane_change_mod
