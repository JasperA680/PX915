module lane_change_mod
    ! Lane-changing sub-step for the multi-lane Nagel-Schreckenberg model.
    !
    ! This module implements a lane-changing update for the multi-lane cellular
    ! automaton traffic model. The rules are based on Rickert, Nagel,
    ! Schreckenberg and Latour, "Two lane traffic simulations using cellular
    ! automata", Physica A 231, 534-550 (1996).
    !
    ! Two lane-changing models are supported through the ``lc_model`` parameter:
    !
    ! * ``LC_SYMMETRIC``: all lanes are treated equally. A vehicle changes lane
    !   only if it is blocked in its current lane.
    ! * ``LC_ASYMMETRIC``: vehicles prefer the rightmost lane, corresponding to
    !   the lowest lane index. A vehicle may return right without being blocked,
    !   but moving left still requires the current lane to be blocked.
    !
    ! For a vehicle at position ``i`` with velocity ``v``, the lane-change
    ! conditions are:
    !
    ! .. math::
    !
    !    T_1:\quad g_{\mathrm{own}}(i) < v + 1,
    !
    ! meaning that the vehicle is blocked in its current lane,
    !
    ! .. math::
    !
    !    T_2:\quad g_{\mathrm{target,fwd}}(i) > v + 1,
    !
    ! meaning that there is enough space ahead in the target lane,
    !
    ! .. math::
    !
    !    T_3:\quad g_{\mathrm{target,back}}(i) > L_{\mathrm{back}},
    !
    ! meaning that there is enough space behind in the target lane, and
    !
    ! .. math::
    !
    !    T_4:\quad r < p_{\mathrm{change}},
    !
    ! where ``r`` is a uniform random number.
    !
    ! Lane index ``1`` is treated as the rightmost lane, with larger lane indices
    ! representing lanes further to the left. Only lanes with the same
    ! ``flow_direction`` within the same road can exchange vehicles. Each vehicle
    ! may move at most one lane per step.
    !
    ! The final site ``L`` is excluded from lane changing because it acts as a
    ! junction holding cell or open outflow boundary.
    !
    ! The update is parallel: all decisions are read from the pre-snapshotted
    ! ``old`` arrays and all accepted lane changes are written to the live
    ! ``cells`` arrays. The caller should snapshot the network before and after
    ! this sub-step.

    use road_network_mod

    implicit none
    private

    integer, parameter, public :: LC_DISABLED = -1
    ! Disable the lane-changing sub-step.
    integer, parameter, public :: LC_SYMMETRIC = 0
    ! Symmetric lane-changing model. Lane changes require the current lane to be blocked.
    integer, parameter, public :: LC_ASYMMETRIC = 1
    ! Asymmetric lane-changing model with preference for the rightmost lane.

    public :: apply_lane_changes

contains

    subroutine apply_lane_changes(net, lc_model, p_change, v_max)
        ! Apply one lane-changing sub-step to every road in the network.
        !
        ! This routine loops over all roads in the network and applies the
        ! lane-changing rule within each road. It assumes that the caller has
        ! already snapshotted the network, so that ``lane%old`` contains the
        ! state used for all lane-change decisions.
        !
        ! The lane-change probability is applied through:
        !
        ! .. math::
        !
        !    r < p_{\mathrm{change}},
        !
        ! where ``r`` is a random number sampled uniformly from ``[0,1)``.
        type(road_network_t), intent(inout) :: net
        ! Road network to update in place.
        integer, intent(in) :: lc_model
        ! Lane-changing model selector: ``LC_DISABLED``, ``LC_SYMMETRIC`` or ``LC_ASYMMETRIC``.
        integer, intent(in) :: v_max
        ! Maximum vehicle speed used to set the backward safety distance.
        real, intent(in) :: p_change
        ! Probability of accepting a candidate lane change.

        integer :: r
        ! Road index.

        do r = 1, size(net%roads)
            call road_lane_change(net%roads(r), lc_model, p_change, v_max)
        end do
    end subroutine apply_lane_changes


    subroutine road_lane_change(road, lc_model, p_change, v_max)
        ! Apply lane changes within one road.
        !
        ! This internal routine scans all lanes and all eligible positions in a
        ! single road. The final site is excluded because it is treated as a
        ! holding cell or outflow boundary. Gap calculations read from
        ! ``lane%old``, while accepted moves are written into ``lane%cells``.
        !
        ! The backward safety distance is chosen as:
        !
        ! .. math::
        !
        !    L_{\mathrm{back}} = v_{\max}.
        !
        ! This means a following car travelling at the maximum speed needs at
        ! least ``L_BACK + 1`` clear cells to avoid a collision after the lane
        ! change.
        type(road_t), intent(inout) :: road
        ! Road whose lanes are updated.
        integer, intent(in) :: lc_model
        ! Lane-changing model selector.
        integer, intent(in) :: v_max
        ! Maximum vehicle speed.
        real, intent(in) :: p_change
        ! Probability of accepting a candidate lane change.

        integer :: n_lanes
        ! Number of lanes in this road.
        integer :: k
        ! Current lane index.
        integer :: L
        ! Length of the current lane.
        integer :: i
        ! Position index in the current lane.
        integer :: v
        ! Vehicle velocity in the old state.
        integer :: tgt
        ! Target lane index. Zero means no lane change.
        integer :: gf
        ! Forward gap in the current lane.
        integer :: go_fwd
        ! Forward gap in the candidate target lane.
        integer :: go_back
        ! Backward gap in the candidate target lane.
        integer :: L_BACK
        ! Backward safety distance.
        logical :: t1
        ! True if the vehicle is blocked in its current lane.
        logical :: t2r
        ! True if there is enough space ahead in the right target lane.
        logical :: t3r
        ! True if there is enough space behind in the right target lane.
        logical :: t2l
        ! True if there is enough space ahead in the left target lane.
        logical :: t3l
        ! True if there is enough space behind in the left target lane.
        logical :: want_right
        ! True if the right lane is a valid candidate target.
        logical :: want_left
        ! True if the left lane is a valid candidate target.
        real :: rnd
        ! Uniform random number used for stochastic acceptance.

        L_BACK = v_max

        n_lanes = size(road%lane)
        if (n_lanes < 2) return

        do k = 1, n_lanes
            L = road%lane(k)%length

            do i = 1, L - 1
                if (.not. road%lane(k)%old(i)%has_car) cycle

                v   = road%lane(k)%old(i)%velocity
                tgt = 0

                ! T1: current lane is blocked.
                gf = fwd_gap(road%lane(k)%old, i, L)
                t1 = (gf < v + 1)

                ! Candidate right lane: lane k-1, towards the rightmost lane.
                want_right = .false.
                if (k > 1 .and. same_dir(road, k, k-1)) then
                    go_fwd  = other_fwd_gap(road%lane(k-1)%old, i, L)
                    go_back = other_back_gap(road%lane(k-1)%old, i)

                    t2r = (go_fwd  > v + 1)
                    t3r = (go_back > L_BACK)

                    want_right = t2r .and. t3r
                end if

                ! Candidate left lane: lane k+1, away from the rightmost lane.
                want_left = .false.
                if (k < n_lanes .and. same_dir(road, k, k+1)) then
                    go_fwd  = other_fwd_gap(road%lane(k+1)%old, i, L)
                    go_back = other_back_gap(road%lane(k+1)%old, i)

                    t2l = (go_fwd  > v + 1)
                    t3l = (go_back > L_BACK)

                    want_left = t2l .and. t3l
                end if

                if (lc_model == LC_ASYMMETRIC) then
                    ! In the asymmetric model, returning right does not require T1.
                    if (want_right) then
                        call random_number(rnd)
                        if (rnd < p_change) tgt = k - 1
                    end if

                    ! Moving left still requires the vehicle to be blocked.
                    if (tgt == 0 .and. want_left .and. t1) then
                        call random_number(rnd)
                        if (rnd < p_change) tgt = k + 1
                    end if
                else
                    ! In the symmetric model, T1 is required for any lane change.
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


    pure logical function same_dir(road, a, b)
        ! Return whether two lanes in the same road have the same flow direction.
        !
        ! Lane changes are only permitted between lanes whose traffic moves in
        ! the same direction.
        type(road_t), intent(in) :: road
        ! Road containing the lanes.
        integer, intent(in) :: a
        ! First lane index.
        integer, intent(in) :: b
        ! Second lane index.

        same_dir = (road%lane(a)%flow_direction == road%lane(b)%flow_direction)
    end function same_dir


    integer function fwd_gap(lane_old, i, L) result(gap)
        ! Return the forward gap in the current lane.
        !
        ! The forward gap is the number of empty cells between position ``i`` and
        ! the next vehicle ahead. If no vehicle is ahead, the returned gap is the
        ! distance to the end of the road:
        !
        ! .. math::
        !
        !    g = L - i.
        !
        ! If the nearest vehicle ahead is at position ``j``, the returned gap is:
        !
        ! .. math::
        !
        !    g = j - i - 1.
        type(cell), intent(in) :: lane_old(:)
        ! Old-state cell array for the current lane.
        integer, intent(in) :: i
        ! Position of the vehicle whose forward gap is requested.
        integer, intent(in) :: L
        ! Lane length.
        integer :: j
        ! Search index.

        gap = L - i

        do j = i + 1, L
            if (lane_old(j)%has_car) then
                gap = j - i - 1
                return
            end if
        end do
    end function fwd_gap


    integer function other_fwd_gap(other_old, i, L) result(gap)
        ! Return the forward gap in a candidate target lane.
        !
        ! If the target position ``i`` is already occupied in the other lane, the
        ! function returns ``-1`` to indicate that a lane change is impossible.
        ! Otherwise it returns the number of empty cells ahead of position ``i``
        ! in the target lane.
        type(cell), intent(in) :: other_old(:)
        ! Old-state cell array for the candidate target lane.
        integer, intent(in) :: i
        ! Position at which the vehicle would enter the target lane.
        integer, intent(in) :: L
        ! Lane length.
        integer :: j
        ! Search index.

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


    integer function other_back_gap(other_old, i) result(gap)
        ! Return the backward gap in a candidate target lane.
        !
        ! If the target position ``i`` is already occupied in the other lane, the
        ! function returns ``-1``. Otherwise it returns the number of empty cells
        ! between position ``i`` and the nearest vehicle behind it.
        !
        ! If the nearest vehicle behind is at position ``j``, the returned gap is:
        !
        ! .. math::
        !
        !    g = i - j - 1.
        !
        ! If there is no vehicle behind, the returned gap is:
        !
        ! .. math::
        !
        !    g = i - 1.
        type(cell), intent(in) :: other_old(:)
        ! Old-state cell array for the candidate target lane.
        integer, intent(in) :: i
        ! Position at which the vehicle would enter the target lane.
        integer :: j
        ! Search index.

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