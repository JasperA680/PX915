program test_lane_change
    !--------------------------------------------------------------------
    ! Unit and integration tests for the lane-change sub-step.
    !
    ! Tests cover:
    !  1. test_no_change_not_blocked   - vehicle with free road does not change
    !  2. test_symmetric_change        - blocked vehicle switches lane
    !  3. t2_block_same_pos            - target cell occupied: T2 fails
    !  4. t3_block_close_follower      - close follower: T3 fails
    !  5. test_p_change_zero           - p_change=0 suppresses all changes
    !  6. test_asymmetric_return_right - non-rightmost returns right (no T1)
    !  7. test_asymmetric_left_blocked - blocked on left because path occupied
    !  8. test_asymmetric_no_t1_needed - confirms rightward requires no T1
    !  9. test_opposite_dir_no_change  - opposite-direction lanes: no exchange
    ! 10. test_mass_conserved_lc_only  - lane-change step preserves car count
    ! 11. test_mass_conserved_full_sim - multi-step sim with lane changes: mass
    ! 12. test_three_lane_symmetric    - three same-dir lanes, middle vehicle
    ! 13. test_three_lane_asymmetric   - three same-dir lanes, always-right rule
    !--------------------------------------------------------------------
    use road_network_mod
    use lane_change_mod
    use network_simulation_mod
    use network_init_mod
    implicit none

    call seed_rng(20260518)

    call test_no_change_not_blocked()
    call test_symmetric_change()
    call t2_block_same_pos()
    call t3_block_close_follower()
    call test_p_change_zero()
    call test_asymmetric_return_right()
    call test_asymmetric_left_blocked()
    call test_asymmetric_no_t1_needed()
    call test_opposite_dir_no_change()
    call test_mass_conserved_lc_only()
    call test_mass_conserved_full_sim()
    call test_three_lane_symmetric()
    call test_three_lane_asymmetric()

    print *, 'test_lane_change: ALL PASSED'

contains

    ! ----------------------------------------------------------------
    ! RNG helpers
    ! ----------------------------------------------------------------
    subroutine seed_rng(s)
        integer, intent(in) :: s
        integer :: n
        integer, allocatable :: seed(:)
        call random_seed(size = n)
        allocate(seed(n))
        seed = s
        call random_seed(put = seed)
        deallocate(seed)
    end subroutine seed_rng

    ! Force the RNG to always return 0.0 so p_change comparisons are
    ! predictable: set seed so next calls produce values near 0.
    ! Instead of hacking the RNG, for determinism tests we use p_change=0.

    ! ----------------------------------------------------------------
    ! Road-building helpers
    ! ----------------------------------------------------------------

    ! Build a network with a single road of n_lanes same-direction lanes,
    ! each of length L, all flow_direction = +1, no junctions.
    subroutine init_test_net(net, n_lanes, L)
        type(road_network_t), intent(out) :: net
        integer,              intent(in)  :: n_lanes, L
        integer :: k
        allocate(net%roads(1))
        allocate(net%junctions(0))
        net%roads(1)%id = 1
        net%roads(1)%end_junction = 0
        allocate(net%roads(1)%lane(n_lanes))
        do k = 1, n_lanes
            net%roads(1)%lane(k)%length         = L
            net%roads(1)%lane(k)%flow_direction = +1
            net%roads(1)%lane(k)%open_in        = .false.
            net%roads(1)%lane(k)%open_out       = .false.
            net%roads(1)%lane(k)%alpha          = 0.0
            net%roads(1)%lane(k)%beta           = 0.0
            allocate(net%roads(1)%lane(k)%cells(L))
            allocate(net%roads(1)%lane(k)%old(L))
            net%roads(1)%lane(k)%cells%has_car  = .false.
            net%roads(1)%lane(k)%cells%velocity = 0
            net%roads(1)%lane(k)%old%has_car    = .false.
            net%roads(1)%lane(k)%old%velocity   = 0
        end do
    end subroutine init_test_net

    subroutine free_test_net(net)
        type(road_network_t), intent(inout) :: net
        integer :: k
        do k = 1, size(net%roads(1)%lane)
            deallocate(net%roads(1)%lane(k)%cells)
            deallocate(net%roads(1)%lane(k)%old)
        end do
        deallocate(net%roads(1)%lane)
        deallocate(net%roads)
        deallocate(net%junctions)
    end subroutine free_test_net

    ! Place a vehicle and sync old state.
    subroutine place_car(net, lane, pos, vel)
        type(road_network_t), intent(inout) :: net
        integer,              intent(in)    :: lane, pos, vel
        net%roads(1)%lane(lane)%cells(pos)%has_car  = .true.
        net%roads(1)%lane(lane)%cells(pos)%velocity = vel
        net%roads(1)%lane(lane)%old(pos)%has_car    = .true.
        net%roads(1)%lane(lane)%old(pos)%velocity   = vel
    end subroutine place_car

    ! Count cars in lane k.
    integer function count_lane(net, k) result(n)
        type(road_network_t), intent(in) :: net
        integer,              intent(in) :: k
        n = count(net%roads(1)%lane(k)%cells%has_car)
    end function count_lane

    integer function total_cars(net) result(n)
        type(road_network_t), intent(in) :: net
        integer :: k
        n = 0
        do k = 1, size(net%roads(1)%lane)
            n = n + count(net%roads(1)%lane(k)%cells%has_car)
        end do
    end function total_cars

    ! ----------------------------------------------------------------
    ! Test 1: symmetric model – vehicle with clear road does not change
    ! ----------------------------------------------------------------
    subroutine test_no_change_not_blocked()
        type(road_network_t) :: net
        integer, parameter :: L = 20
        ! Vehicle at pos 8, velocity 3 – gap ahead = L-8 = 12, v+1 = 4
        ! gap(12) >= v+1(4) -> T1 false -> no change in symmetric model.
        call init_test_net(net, 2, L)
        call place_car(net, 1, 8, 3)
        call apply_lane_changes(net, LC_SYMMETRIC, 1.0)
        if (.not. net%roads(1)%lane(1)%cells(8)%has_car) then
            print *, 'FAIL test_no_change_not_blocked: vehicle left lane 1'
            stop 1
        end if
        if (net%roads(1)%lane(2)%cells(8)%has_car) then
            print *, 'FAIL test_no_change_not_blocked: ghost in lane 2'
            stop 1
        end if
        call free_test_net(net)
    end subroutine test_no_change_not_blocked

    ! ----------------------------------------------------------------
    ! Test 2: symmetric model – blocked vehicle switches lane when clear
    ! All four conditions (T1-T4) satisfied with p_change=1.
    !
    ! Setup: lane 1 has car at pos 8 (v=2) and blocker at pos 10.
    !        gap_fwd = 1 < v+1=3  -> T1 true
    !        Lane 2 is empty at pos 8: gap_other_fwd = L-8 = 12 > 3 -> T2 true
    !        No car behind pos 8 in lane 2: gap_other_back = 7 > L_BACK=5 -> T3 true
    !        p_change=1 -> T4 true  ->  vehicle moves to lane 2.
    ! ----------------------------------------------------------------
    subroutine test_symmetric_change()
        type(road_network_t) :: net
        integer, parameter :: L = 20
        call init_test_net(net, 2, L)
        call place_car(net, 1, 8, 2)    ! vehicle to change, v=2
        call place_car(net, 1, 10, 0)   ! blocker (gap_fwd = 1 < 3 = v+1)
        ! lane 2 is completely empty -> gap_other_fwd = 12, gap_other_back = 7
        call apply_lane_changes(net, LC_SYMMETRIC, 1.0)
        ! Vehicle must have moved from lane 1 pos 8 to lane 2 pos 8.
        if (net%roads(1)%lane(1)%cells(8)%has_car) then
            print *, 'FAIL test_symmetric_change: vehicle still in lane 1'
            stop 1
        end if
        if (.not. net%roads(1)%lane(2)%cells(8)%has_car) then
            print *, 'FAIL test_symmetric_change: vehicle not in lane 2'
            stop 1
        end if
        if (net%roads(1)%lane(2)%cells(8)%velocity /= 2) then
            print *, 'FAIL test_symmetric_change: velocity not preserved'
            stop 1
        end if
        ! Blocker must stay.
        if (.not. net%roads(1)%lane(1)%cells(10)%has_car) then
            print *, 'FAIL test_symmetric_change: blocker vanished'
            stop 1
        end if
        call free_test_net(net)
    end subroutine test_symmetric_change

    ! ----------------------------------------------------------------
    ! Test 3: T2 fails because the target cell (same position) is occupied
    ! ----------------------------------------------------------------
    subroutine t2_block_same_pos()
        type(road_network_t) :: net
        integer, parameter :: L = 20
        call init_test_net(net, 2, L)
        call place_car(net, 1, 8, 2)    ! vehicle
        call place_car(net, 1, 10, 0)   ! blocker: T1 true
        call place_car(net, 2, 8, 1)    ! occupies target cell: T2 fails
        call apply_lane_changes(net, LC_SYMMETRIC, 1.0)
        ! Both vehicles must remain exactly where they are.
        if (.not. net%roads(1)%lane(1)%cells(8)%has_car) then
            print *, 'FAIL t2_block_same_pos: vehicle moved despite T2 fail'
            stop 1
        end if
        if (.not. net%roads(1)%lane(2)%cells(8)%has_car) then
            print *, 'FAIL t2_block_same_pos: lane-2 car moved'
            stop 1
        end if
        call free_test_net(net)
    end subroutine t2_block_same_pos

    ! ----------------------------------------------------------------
    ! Test 4: T3 fails because a follower is too close in the target lane.
    ! gap_other_back = i - j - 1 = 8 - 4 - 1 = 3, which is NOT > L_BACK=5.
    ! ----------------------------------------------------------------
    subroutine t3_block_close_follower()
        type(road_network_t) :: net
        integer, parameter :: L = 20
        call init_test_net(net, 2, L)
        call place_car(net, 1, 8, 2)    ! vehicle (T1: gap=1 < 3)
        call place_car(net, 1, 10, 0)   ! blocker
        call place_car(net, 2, 4, 2)    ! close follower: gap_back = 3 <= 5
        call apply_lane_changes(net, LC_SYMMETRIC, 1.0)
        if (.not. net%roads(1)%lane(1)%cells(8)%has_car) then
            print *, 'FAIL t3_block_close_follower: vehicle moved despite T3 fail'
            stop 1
        end if
        call free_test_net(net)
    end subroutine t3_block_close_follower

    ! ----------------------------------------------------------------
    ! Test 5: p_change = 0 suppresses lane changes even when all T1-T3 met
    ! ----------------------------------------------------------------
    subroutine test_p_change_zero()
        type(road_network_t) :: net
        integer, parameter :: L = 20
        call init_test_net(net, 2, L)
        call place_car(net, 1, 8, 2)
        call place_car(net, 1, 10, 0)
        ! lane 2 is clear - all T1-T3 satisfied - but p_change=0
        call apply_lane_changes(net, LC_SYMMETRIC, 0.0)
        if (.not. net%roads(1)%lane(1)%cells(8)%has_car) then
            print *, 'FAIL test_p_change_zero: vehicle changed despite p_change=0'
            stop 1
        end if
        if (net%roads(1)%lane(2)%cells(8)%has_car) then
            print *, 'FAIL test_p_change_zero: ghost appeared in lane 2'
            stop 1
        end if
        call free_test_net(net)
    end subroutine test_p_change_zero

    ! ----------------------------------------------------------------
    ! Test 6: asymmetric – non-rightmost vehicle returns right even when
    ! NOT blocked in own lane (T1 not required for rightward return).
    !
    ! Vehicle in lane 2 (leftmost), no blocker ahead.
    ! gap_fwd_own = large -> T1 would be false in symmetric.
    ! Lane 1 (rightmost) is clear -> T2 & T3 satisfied.
    ! Asymmetric model: vehicle should return to lane 1.
    ! ----------------------------------------------------------------
    subroutine test_asymmetric_return_right()
        type(road_network_t) :: net
        integer, parameter :: L = 20
        call init_test_net(net, 2, L)
        ! Put vehicle in lane 2 at pos 8, velocity 2 (not blocked)
        call place_car(net, 2, 8, 2)
        ! Lane 1 is completely empty
        call apply_lane_changes(net, LC_ASYMMETRIC, 1.0)
        if (net%roads(1)%lane(2)%cells(8)%has_car) then
            print *, 'FAIL test_asymmetric_return_right: vehicle did not leave lane 2'
            stop 1
        end if
        if (.not. net%roads(1)%lane(1)%cells(8)%has_car) then
            print *, 'FAIL test_asymmetric_return_right: vehicle not in lane 1'
            stop 1
        end if
        call free_test_net(net)
    end subroutine test_asymmetric_return_right

    ! ----------------------------------------------------------------
    ! Test 7: asymmetric – rightward return blocked by car at same position.
    ! Vehicle in lane 2, lane 1 occupied at the same position -> T2 fails,
    ! own lane not blocked -> T1 false -> cannot go further left either.
    ! Vehicle must stay in lane 2.
    ! ----------------------------------------------------------------
    subroutine test_asymmetric_left_blocked()
        type(road_network_t) :: net
        integer, parameter :: L = 20
        call init_test_net(net, 2, L)
        call place_car(net, 2, 8, 2)    ! vehicle in left lane
        call place_car(net, 1, 8, 0)    ! occupies rightward target: T2 fails
        ! No blocker ahead in lane 2 -> T1 is false -> no leftward escape either
        call apply_lane_changes(net, LC_ASYMMETRIC, 1.0)
        if (.not. net%roads(1)%lane(2)%cells(8)%has_car) then
            print *, 'FAIL test_asymmetric_left_blocked: vehicle left lane 2'
            stop 1
        end if
        if (count_lane(net, 1) /= 1 .or. &
            .not. net%roads(1)%lane(1)%cells(8)%has_car) then
            print *, 'FAIL test_asymmetric_left_blocked: lane 1 changed unexpectedly'
            stop 1
        end if
        call free_test_net(net)
    end subroutine test_asymmetric_left_blocked

    ! ----------------------------------------------------------------
    ! Test 8: asymmetric – vehicle in rightmost lane (lane 1) does NOT
    ! move left unless explicitly blocked (T1 required for leftward moves).
    ! Vehicle at pos 8, v=2, large gap ahead -> T1 false -> no change.
    ! ----------------------------------------------------------------
    subroutine test_asymmetric_no_t1_needed()
        type(road_network_t) :: net
        integer, parameter :: L = 20
        call init_test_net(net, 2, L)
        call place_car(net, 1, 8, 2)    ! rightmost lane, unblocked
        ! lane 2 is clear
        call apply_lane_changes(net, LC_ASYMMETRIC, 1.0)
        ! Vehicle is in rightmost lane already; rightward move not possible (k=1).
        ! Leftward move requires T1, which is false. Vehicle stays.
        if (.not. net%roads(1)%lane(1)%cells(8)%has_car) then
            print *, 'FAIL test_asymmetric_no_t1_needed: unblocked vehicle moved left'
            stop 1
        end if
        if (net%roads(1)%lane(2)%cells(8)%has_car) then
            print *, 'FAIL test_asymmetric_no_t1_needed: ghost in lane 2'
            stop 1
        end if
        call free_test_net(net)
    end subroutine test_asymmetric_no_t1_needed

    ! ----------------------------------------------------------------
    ! Test 9: opposite-direction lanes must not exchange vehicles.
    ! The crossroad has lane 1 (outbound) and lane 2 (inbound) on each
    ! road; they have opposite flow_directions, so same_dir() is false.
    ! ----------------------------------------------------------------
    subroutine test_opposite_dir_no_change()
        type(road_network_t) :: net
        integer, parameter :: L = 12
        integer :: n_before_r1, n_before_r2, n_after_r1, n_after_r2

        call init_crossroad(net, L, 0.0, 0.0, 0.5, 0.0)

        ! Place vehicles in each lane on road 1.
        net%roads(1)%lane(1)%cells(4)%has_car = .true.
        net%roads(1)%lane(1)%cells(4)%velocity = 2
        net%roads(1)%lane(1)%old(4)%has_car = .true.
        net%roads(1)%lane(1)%old(4)%velocity = 2
        net%roads(1)%lane(2)%cells(4)%has_car = .true.
        net%roads(1)%lane(2)%cells(4)%velocity = 1
        net%roads(1)%lane(2)%old(4)%has_car = .true.
        net%roads(1)%lane(2)%old(4)%velocity = 1

        n_before_r1 = count(net%roads(1)%lane(1)%cells%has_car)
        n_before_r2 = count(net%roads(1)%lane(2)%cells%has_car)

        call apply_lane_changes(net, LC_SYMMETRIC, 1.0)

        n_after_r1 = count(net%roads(1)%lane(1)%cells%has_car)
        n_after_r2 = count(net%roads(1)%lane(2)%cells%has_car)

        if (n_before_r1 /= n_after_r1 .or. n_before_r2 /= n_after_r2) then
            print *, 'FAIL test_opposite_dir_no_change: vehicles crossed direction boundary'
            stop 1
        end if

        call free_network(net)
    end subroutine test_opposite_dir_no_change

    ! ----------------------------------------------------------------
    ! Test 10: the lane-change step alone must conserve vehicle count.
    ! ----------------------------------------------------------------
    subroutine test_mass_conserved_lc_only()
        type(road_network_t) :: net
        integer, parameter :: L = 30
        integer :: n_before, n_after, i

        call init_test_net(net, 3, L)
        ! Scatter vehicles across 3 lanes using a repeating pattern.
        do i = 2, L - 2, 3
            call place_car(net, 1, i,   1)
        end do
        do i = 3, L - 2, 4
            call place_car(net, 2, i,   2)
        end do
        do i = 5, L - 2, 5
            call place_car(net, 3, i,   0)
        end do

        n_before = total_cars(net)
        call apply_lane_changes(net, LC_ASYMMETRIC, 1.0)
        n_after  = total_cars(net)

        if (n_before /= n_after) then
            print *, 'FAIL test_mass_conserved_lc_only: count changed:', n_before, '->', n_after
            stop 1
        end if
        call free_test_net(net)
    end subroutine test_mass_conserved_lc_only

    ! ----------------------------------------------------------------
    ! Test 11: full simulation (network_step with lane changes) on a
    ! crossroad conserves mass at every step.
    ! ----------------------------------------------------------------
    subroutine test_mass_conserved_full_sim()
        type(road_network_t) :: net
        integer, parameter :: L = 15, N_STEPS = 200
        integer :: step, n_old, n_new, entries, exits, r, ln, Lk

        call init_crossroad(net, L, 0.5, 0.5, 0.3, 0.2)

        do step = 1, N_STEPS
            n_old = count_occupied_network(net)
            entries = 0
            exits   = 0

            call network_step(net, "NS", LC_SYMMETRIC, 1.0)

            n_new = count_occupied_network(net)
            do r = 1, size(net%roads)
                do ln = 1, size(net%roads(r)%lane)
                    if (net%roads(r)%lane(ln)%open_in) then
                        if (net%roads(r)%lane(ln)%cells(1)%has_car .and. &
                            .not. net%roads(r)%lane(ln)%old(1)%has_car) &
                            entries = entries + 1
                    end if
                    if (net%roads(r)%lane(ln)%open_out) then
                        Lk = net%roads(r)%lane(ln)%length
                        if (.not. net%roads(r)%lane(ln)%cells(Lk)%has_car .and. &
                            net%roads(r)%lane(ln)%old(Lk)%has_car) &
                            exits = exits + 1
                    end if
                end do
            end do

            if (n_new /= n_old + entries - exits) then
                print *, 'FAIL test_mass_conserved_full_sim: step', step, &
                         'n_old', n_old, 'entries', entries, 'exits', exits, 'n_new', n_new
                stop 1
            end if
        end do

        call free_network(net)
    end subroutine test_mass_conserved_full_sim

    ! ----------------------------------------------------------------
    ! Test 12: three same-direction lanes, symmetric model.
    ! Vehicle in middle lane (lane 2), blocked, checks both neighbors.
    ! Right (lane 1) is occupied at same position -> cannot go right.
    ! Left (lane 3) is clear -> should move left.
    ! ----------------------------------------------------------------
    subroutine test_three_lane_symmetric()
        type(road_network_t) :: net
        integer, parameter :: L = 20
        call init_test_net(net, 3, L)
        call place_car(net, 2, 8, 2)    ! vehicle in middle lane
        call place_car(net, 2, 10, 0)   ! blocker: gap=1 < 3=v+1 -> T1 true
        call place_car(net, 1, 8, 0)    ! right neighbor occupied -> T2 right fails
        ! lane 3 at pos 8 is empty, back is clear -> T2/T3 left pass
        call apply_lane_changes(net, LC_SYMMETRIC, 1.0)
        if (net%roads(1)%lane(2)%cells(8)%has_car) then
            print *, 'FAIL test_three_lane_symmetric: vehicle did not leave middle lane'
            stop 1
        end if
        if (.not. net%roads(1)%lane(3)%cells(8)%has_car) then
            print *, 'FAIL test_three_lane_symmetric: vehicle not in lane 3'
            stop 1
        end if
        call free_test_net(net)
    end subroutine test_three_lane_symmetric

    ! ----------------------------------------------------------------
    ! Test 13: three same-direction lanes, asymmetric model.
    ! Vehicle in lane 3 (leftmost), not blocked.
    ! Lane 2 is clear -> tries to move right to lane 2 (T1 not required).
    ! ----------------------------------------------------------------
    subroutine test_three_lane_asymmetric()
        type(road_network_t) :: net
        integer, parameter :: L = 20
        call init_test_net(net, 3, L)
        call place_car(net, 3, 8, 2)    ! vehicle in leftmost lane, unblocked
        ! lanes 1 and 2 are empty
        call apply_lane_changes(net, LC_ASYMMETRIC, 1.0)
        ! Asymmetric: vehicle always tries to go right. Lane 2 is one step closer.
        if (net%roads(1)%lane(3)%cells(8)%has_car) then
            print *, 'FAIL test_three_lane_asymmetric: vehicle did not leave lane 3'
            stop 1
        end if
        if (.not. net%roads(1)%lane(2)%cells(8)%has_car) then
            print *, 'FAIL test_three_lane_asymmetric: vehicle not in lane 2 (one step right)'
            stop 1
        end if
        ! Only ONE lane change per step: vehicle should be in lane 2, not lane 1.
        if (net%roads(1)%lane(1)%cells(8)%has_car) then
            print *, 'FAIL test_three_lane_asymmetric: vehicle jumped two lanes in one step'
            stop 1
        end if
        call free_test_net(net)
    end subroutine test_three_lane_asymmetric

end program test_lane_change
