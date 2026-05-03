program test_junction
    !--------------------------------------------------------------------
    ! Phase 3 validation suite for evaluate_junctions.
    !
    ! Each test_* subroutine constructs a fresh crossroad, plants
    ! vehicles in specific holding cells, runs evaluate_junctions, and
    ! checks the resulting state.  Tests use a fixed RNG seed so that
    ! stochastic tie-breaks are reproducible.
    !--------------------------------------------------------------------
    use vehicle_mod
    use road_network_mod
    use network_init_mod
    use junction_mod
    implicit none

    integer, parameter :: L = 10

    call seed_rng(42)

    call test_dest_math()
    call test_single_move()
    call test_physical_block()
    call test_yield_right()
    call test_right_turn_yields()
    call test_parallel_opposite_straight()
    call test_parallel_left_left()
    call test_mutual_right_deadlock()
    call test_4way_straight_deadlock()

    print *, 'test_junction: ALL OK'

contains

    !-----------------------------------------------------------------
    ! Helpers
    !-----------------------------------------------------------------
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

    subroutine fresh(net)
        type(road_network_t), intent(out) :: net
        call init_crossroad(net, L, 0.0, 0.0, 0.0, 0.0)        ! no inflow/outflow during these unit tests
    end subroutine fresh

    subroutine plant(net, leg, intent_code)
        ! Place a vehicle in the holding cell of clockwise leg `leg`.
        type(road_network_t), intent(inout) :: net
        integer, intent(in) :: leg, intent_code
        integer :: src_lane
        src_lane = inbound_lane_idx(net%junctions(1)%end_at(leg))
        net%roads(leg)%lane(src_lane)%cells(L) = intent_code
    end subroutine plant

    function holding_cell(net, leg) result(v)
        type(road_network_t), intent(in) :: net
        integer, intent(in) :: leg
        integer :: v, src_lane
        src_lane = inbound_lane_idx(net%junctions(1)%end_at(leg))
        v = net%roads(leg)%lane(src_lane)%cells(L)
    end function holding_cell

    function dest_cell(net, leg) result(v)
        ! Read the first cell of leg's outbound lane (where its arriving
        ! vehicle would be placed by the junction).
        type(road_network_t), intent(in) :: net
        integer, intent(in) :: leg
        integer :: v, dst_lane
        dst_lane = outbound_lane_idx(net%junctions(1)%end_at(leg))
        v = net%roads(leg)%lane(dst_lane)%cells(1)
    end function dest_cell

    subroutine assert(cond, msg)
        logical, intent(in) :: cond
        character(len=*), intent(in) :: msg
        if (.not. cond) then
            print *, 'ASSERT FAILED: ', trim(msg)
            stop 1
        end if
    end subroutine assert

    !-----------------------------------------------------------------
    ! 1. Destination math: hand-worked table.
    !-----------------------------------------------------------------
    subroutine test_dest_math()
        type(road_network_t) :: net
        integer :: dst_road, dst_lane, dst_end, k

        call fresh(net)

        ! From leg 1 (N): straight->3, right->4, left->2
        call compute_destination(net%junctions(1), 1, V_STRAIGHT, dst_road, dst_lane, dst_end)
        call assert(dst_road == 3, 'k=1 STRAIGHT')
        call compute_destination(net%junctions(1), 1, V_RIGHT,    dst_road, dst_lane, dst_end)
        call assert(dst_road == 4, 'k=1 RIGHT')
        call compute_destination(net%junctions(1), 1, V_LEFT,     dst_road, dst_lane, dst_end)
        call assert(dst_road == 2, 'k=1 LEFT')

        ! Spot-check the rest.
        call compute_destination(net%junctions(1), 2, V_STRAIGHT, dst_road, dst_lane, dst_end)
        call assert(dst_road == 4, 'k=2 STRAIGHT')
        call compute_destination(net%junctions(1), 2, V_RIGHT,    dst_road, dst_lane, dst_end)
        call assert(dst_road == 1, 'k=2 RIGHT')
        call compute_destination(net%junctions(1), 3, V_LEFT,     dst_road, dst_lane, dst_end)
        call assert(dst_road == 4, 'k=3 LEFT')
        call compute_destination(net%junctions(1), 4, V_STRAIGHT, dst_road, dst_lane, dst_end)
        call assert(dst_road == 2, 'k=4 STRAIGHT')

        ! All destinations should be outbound lane 1 (since end_at = 1 everywhere here).
        do k = 1, 4
            call compute_destination(net%junctions(1), k, V_STRAIGHT, dst_road, dst_lane, dst_end)
            call assert(dst_lane == 1, 'destination is outbound lane')
        end do

        call free_network(net)
        print *, '  test_dest_math OK'
    end subroutine test_dest_math

    !-----------------------------------------------------------------
    ! 2. Single vehicle moves with no contention.
    !-----------------------------------------------------------------
    subroutine test_single_move()
        type(road_network_t) :: net
        call fresh(net)
        call plant(net, 1, V_STRAIGHT)
        call snapshot_network(net)
        call evaluate_junctions(net)

        call assert(holding_cell(net, 1) == V_EMPTY,    'holding 1 not cleared')
        call assert(dest_cell(net, 3)    == V_STRAIGHT, 'road 3 site 1 not filled')

        call free_network(net)
        print *, '  test_single_move OK'
    end subroutine test_single_move

    !-----------------------------------------------------------------
    ! 3. Physical block: pre-fill destination -> no move.
    !-----------------------------------------------------------------
    subroutine test_physical_block()
        type(road_network_t) :: net
        call fresh(net)
        call plant(net, 1, V_STRAIGHT)
        net%roads(3)%lane(1)%cells(1) = V_STRAIGHT             ! pre-fill destination

        call snapshot_network(net)
        call evaluate_junctions(net)

        call assert(holding_cell(net, 1) == V_STRAIGHT, 'holding 1 should be unchanged')
        call assert(dest_cell(net, 3)    == V_STRAIGHT, 'destination should be unchanged')

        call free_network(net)
        print *, '  test_physical_block OK'
    end subroutine test_physical_block

    !-----------------------------------------------------------------
    ! 4. Yield to the right: legs 1 and 4 both straight.
    !    From leg 1's perspective leg 4 is to its right -> 1 yields to 4.
    !    They have different destinations (3 and 2) so 4 advances unhindered.
    !-----------------------------------------------------------------
    subroutine test_yield_right()
        type(road_network_t) :: net
        call fresh(net)
        call plant(net, 1, V_STRAIGHT)
        call plant(net, 4, V_STRAIGHT)
        call snapshot_network(net)
        call evaluate_junctions(net)

        call assert(holding_cell(net, 4) == V_EMPTY,    'leg 4 should advance')
        call assert(dest_cell(net, 2)    == V_STRAIGHT, 'leg 4 -> road 2')

        call assert(holding_cell(net, 1) == V_STRAIGHT, 'leg 1 must yield')
        call assert(dest_cell(net, 3)    == V_EMPTY,    'leg 1 destination unchanged')

        call free_network(net)
        print *, '  test_yield_right OK'
    end subroutine test_yield_right

    !-----------------------------------------------------------------
    ! 5. Right-turn vs oncoming straight.
    !    Leg 1 RIGHT, leg 3 STRAIGHT (oncoming).  Leg 1 must yield.
    !    Destinations: leg 1 -> road 4, leg 3 -> road 1 (different).
    !-----------------------------------------------------------------
    subroutine test_right_turn_yields()
        type(road_network_t) :: net
        call fresh(net)
        call plant(net, 1, V_RIGHT)
        call plant(net, 3, V_STRAIGHT)
        call snapshot_network(net)
        call evaluate_junctions(net)

        call assert(holding_cell(net, 3) == V_EMPTY,    'leg 3 should advance')
        call assert(dest_cell(net, 1)    == V_STRAIGHT, 'leg 3 -> road 1')

        call assert(holding_cell(net, 1) == V_RIGHT,    'leg 1 must yield to oncoming')

        call free_network(net)
        print *, '  test_right_turn_yields OK'
    end subroutine test_right_turn_yields

    !-----------------------------------------------------------------
    ! 6. Two opposite straights: parallel, both move.
    !-----------------------------------------------------------------
    subroutine test_parallel_opposite_straight()
        type(road_network_t) :: net
        call fresh(net)
        call plant(net, 1, V_STRAIGHT)
        call plant(net, 3, V_STRAIGHT)
        call snapshot_network(net)
        call evaluate_junctions(net)

        call assert(holding_cell(net, 1) == V_EMPTY, 'leg 1 should advance')
        call assert(holding_cell(net, 3) == V_EMPTY, 'leg 3 should advance')
        call assert(dest_cell(net, 3)    == V_STRAIGHT, 'leg 1 -> road 3')
        call assert(dest_cell(net, 1)    == V_STRAIGHT, 'leg 3 -> road 1')

        call free_network(net)
        print *, '  test_parallel_opposite_straight OK'
    end subroutine test_parallel_opposite_straight

    !-----------------------------------------------------------------
    ! 7. Two left turns from opposite legs: paths don't cross, both move.
    !    Leg 1 LEFT -> road 2; leg 3 LEFT -> road 4.
    !-----------------------------------------------------------------
    subroutine test_parallel_left_left()
        type(road_network_t) :: net
        call fresh(net)
        call plant(net, 1, V_LEFT)
        call plant(net, 3, V_LEFT)
        call snapshot_network(net)
        call evaluate_junctions(net)

        call assert(holding_cell(net, 1) == V_EMPTY, 'leg 1 should advance')
        call assert(holding_cell(net, 3) == V_EMPTY, 'leg 3 should advance')
        call assert(dest_cell(net, 2)    == V_LEFT,  'leg 1 -> road 2')
        call assert(dest_cell(net, 4)    == V_LEFT,  'leg 3 -> road 4')

        call free_network(net)
        print *, '  test_parallel_left_left OK'
    end subroutine test_parallel_left_left

    !-----------------------------------------------------------------
    ! 8. Mutual right turns from opposite legs: classic deadlock.
    !    Both want to cut across each other: exactly one wins via the
    !    stochastic tie-break.
    !-----------------------------------------------------------------
    subroutine test_mutual_right_deadlock()
        type(road_network_t) :: net
        integer :: n_moved
        call fresh(net)
        call plant(net, 1, V_RIGHT)
        call plant(net, 3, V_RIGHT)
        call snapshot_network(net)
        call evaluate_junctions(net)

        n_moved = 0
        if (holding_cell(net, 1) == V_EMPTY) n_moved = n_moved + 1
        if (holding_cell(net, 3) == V_EMPTY) n_moved = n_moved + 1
        call assert(n_moved == 1, 'mutual-right deadlock: exactly one should move')

        call free_network(net)
        print *, '  test_mutual_right_deadlock OK'
    end subroutine test_mutual_right_deadlock

    !-----------------------------------------------------------------
    ! 9. Four-way symmetric standoff: every leg STRAIGHT, every leg
    !    yields to its right.  Tie-break should release exactly one;
    !    the others remain blocked.
    !-----------------------------------------------------------------
    subroutine test_4way_straight_deadlock()
        type(road_network_t) :: net
        integer :: n_moved, k
        call fresh(net)
        do k = 1, 4
            call plant(net, k, V_STRAIGHT)
        end do
        call snapshot_network(net)
        call evaluate_junctions(net)

        n_moved = 0
        do k = 1, 4
            if (holding_cell(net, k) == V_EMPTY) n_moved = n_moved + 1
        end do
        call assert(n_moved == 1, '4-way standoff: exactly one should move')

        call free_network(net)
        print *, '  test_4way_straight_deadlock OK'
    end subroutine test_4way_straight_deadlock

end program test_junction
