program test_junction
    !--------------------------------------------------------------------
    ! Phase 3 validation suite for evaluate_junctions.
    !
    ! Each test_* subroutine constructs a fresh crossroad, plants
    ! V_OCCUPIED into specific holding cells, force-routes them to
    ! deterministic destinations via force_route (Kronecker-delta prob),
    ! runs evaluate_junctions, and checks the resulting state.
    ! Tests use a fixed RNG seed so that stochastic tie-breaks are
    ! reproducible.
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
        ! p_left=p_right=0 => all in_routes point to STRAIGHT (prob=1).
        call init_crossroad(net, L, 0.0, 0.0, 0.0, 0.0)
    end subroutine fresh

    subroutine plant(net, leg)
        ! Place V_OCCUPIED in the holding cell (site L of inbound lane 2) for leg.
        type(road_network_t), intent(inout) :: net
        integer, intent(in) :: leg
        net%roads(leg)%lane(2)%cells(L) = V_OCCUPIED
    end subroutine plant

    subroutine force_route(net, jid, leg, out_idx)
        ! Set in_routes(leg)%prob to a Kronecker delta on out_idx so
        ! destination sampling is deterministic for this leg.
        type(road_network_t), intent(inout) :: net
        integer, intent(in) :: jid, leg, out_idx
        net%junctions(jid)%in_routes(leg)%prob        = 0.0
        net%junctions(jid)%in_routes(leg)%prob(out_idx) = 1.0
    end subroutine force_route

    function holding_cell(net, leg) result(v)
        type(road_network_t), intent(in) :: net
        integer, intent(in) :: leg
        integer :: v
        v = net%roads(leg)%lane(2)%cells(L)
    end function holding_cell

    function dest_cell(net, leg) result(v)
        ! Read site 1 of leg's outbound lane (lane(1), flow_direction=+1).
        type(road_network_t), intent(in) :: net
        integer, intent(in) :: leg
        integer :: v
        v = net%roads(leg)%lane(1)%cells(1)
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
    ! 1. Destination math: exercise move_category directly.
    !    Categories: 1=LEFT, 2=STRAIGHT, 3=RIGHT, 0=UTURN.
    !    For n=4: d=modulo(out-in,4); d=1->LEFT, d=2->STRAIGHT, d=3->RIGHT.
    !-----------------------------------------------------------------
    subroutine test_dest_math()
        ! From leg 1 (N): LEFT->2, STRAIGHT->3, RIGHT->4
        call assert(move_category(1, 2, 4) == 1, 'k=1 out=2 LEFT')
        call assert(move_category(1, 3, 4) == 2, 'k=1 out=3 STRAIGHT')
        call assert(move_category(1, 4, 4) == 3, 'k=1 out=4 RIGHT')
        call assert(move_category(1, 1, 4) == 0, 'k=1 out=1 UTURN')
        ! Spot-check other legs.
        call assert(move_category(2, 4, 4) == 2, 'k=2 out=4 STRAIGHT')
        call assert(move_category(2, 1, 4) == 3, 'k=2 out=1 RIGHT')
        call assert(move_category(3, 4, 4) == 1, 'k=3 out=4 LEFT')
        call assert(move_category(4, 2, 4) == 2, 'k=4 out=2 STRAIGHT')
        print *, '  test_dest_math OK'
    end subroutine test_dest_math

    !-----------------------------------------------------------------
    ! 2. Single vehicle moves with no contention.
    !    Leg 1 -> out=3 (STRAIGHT); holding clears, road 3 site 1 filled.
    !-----------------------------------------------------------------
    subroutine test_single_move()
        type(road_network_t) :: net
        call fresh(net)
        call plant(net, 1)
        call force_route(net, 1, 1, 3)   ! leg 1 STRAIGHT -> outbound 3
        call snapshot_network(net)
        call evaluate_junctions(net)

        call assert(holding_cell(net, 1) == V_EMPTY,    'holding 1 not cleared')
        call assert(dest_cell(net, 3)    == V_OCCUPIED, 'road 3 site 1 not filled')

        call free_network(net)
        print *, '  test_single_move OK'
    end subroutine test_single_move

    !-----------------------------------------------------------------
    ! 3. Physical block: pre-fill destination -> no move.
    !-----------------------------------------------------------------
    subroutine test_physical_block()
        type(road_network_t) :: net
        call fresh(net)
        call plant(net, 1)
        call force_route(net, 1, 1, 3)
        net%roads(3)%lane(1)%cells(1) = V_OCCUPIED   ! pre-fill destination

        call snapshot_network(net)
        call evaluate_junctions(net)

        call assert(holding_cell(net, 1) == V_OCCUPIED, 'holding 1 should be unchanged')
        call assert(dest_cell(net, 3)    == V_OCCUPIED, 'destination should be unchanged')

        call free_network(net)
        print *, '  test_physical_block OK'
    end subroutine test_physical_block

    !-----------------------------------------------------------------
    ! 4. Yield to the right: legs 1 and 4 both STRAIGHT.
    !    Leg 4 is to leg 1's right -> leg 1 yields, leg 4 advances.
    !    Destinations differ (leg 1->3, leg 4->2) so no same-dest conflict.
    !-----------------------------------------------------------------
    subroutine test_yield_right()
        type(road_network_t) :: net
        call fresh(net)
        call plant(net, 1)
        call plant(net, 4)
        call force_route(net, 1, 1, 3)   ! leg 1 STRAIGHT -> outbound 3
        call force_route(net, 1, 4, 2)   ! leg 4 STRAIGHT -> outbound 2

        call snapshot_network(net)
        call evaluate_junctions(net)

        call assert(holding_cell(net, 4) == V_EMPTY,    'leg 4 should advance')
        call assert(dest_cell(net, 2)    == V_OCCUPIED, 'leg 4 -> road 2')
        call assert(holding_cell(net, 1) == V_OCCUPIED, 'leg 1 must yield')
        call assert(dest_cell(net, 3)    == V_EMPTY,    'leg 1 destination unchanged')

        call free_network(net)
        print *, '  test_yield_right OK'
    end subroutine test_yield_right

    !-----------------------------------------------------------------
    ! 5. Right-turn vs oncoming straight (R2).
    !    Leg 1 RIGHT (->out 4), leg 3 STRAIGHT (->out 1).
    !    Leg 1 yields to oncoming leg 3; leg 3 advances.
    !-----------------------------------------------------------------
    subroutine test_right_turn_yields()
        type(road_network_t) :: net
        call fresh(net)
        call plant(net, 1)
        call plant(net, 3)
        call force_route(net, 1, 1, 4)   ! leg 1 RIGHT  -> outbound 4
        call force_route(net, 1, 3, 1)   ! leg 3 STRAIGHT -> outbound 1

        call snapshot_network(net)
        call evaluate_junctions(net)

        call assert(holding_cell(net, 3) == V_EMPTY,    'leg 3 should advance')
        call assert(dest_cell(net, 1)    == V_OCCUPIED, 'leg 3 -> road 1')
        call assert(holding_cell(net, 1) == V_OCCUPIED, 'leg 1 must yield to oncoming')

        call free_network(net)
        print *, '  test_right_turn_yields OK'
    end subroutine test_right_turn_yields

    !-----------------------------------------------------------------
    ! 6. Two opposite straights: parallel, both move.
    !-----------------------------------------------------------------
    subroutine test_parallel_opposite_straight()
        type(road_network_t) :: net
        call fresh(net)
        call plant(net, 1)
        call plant(net, 3)
        call force_route(net, 1, 1, 3)   ! leg 1 STRAIGHT -> outbound 3
        call force_route(net, 1, 3, 1)   ! leg 3 STRAIGHT -> outbound 1

        call snapshot_network(net)
        call evaluate_junctions(net)

        call assert(holding_cell(net, 1) == V_EMPTY,    'leg 1 should advance')
        call assert(holding_cell(net, 3) == V_EMPTY,    'leg 3 should advance')
        call assert(dest_cell(net, 3)    == V_OCCUPIED, 'leg 1 -> road 3')
        call assert(dest_cell(net, 1)    == V_OCCUPIED, 'leg 3 -> road 1')

        call free_network(net)
        print *, '  test_parallel_opposite_straight OK'
    end subroutine test_parallel_opposite_straight

    !-----------------------------------------------------------------
    ! 7. Two left turns from opposite legs: paths don't cross, both move.
    !    Leg 1 LEFT->2; leg 3 LEFT->4.
    !-----------------------------------------------------------------
    subroutine test_parallel_left_left()
        type(road_network_t) :: net
        call fresh(net)
        call plant(net, 1)
        call plant(net, 3)
        call force_route(net, 1, 1, 2)   ! leg 1 LEFT -> outbound 2
        call force_route(net, 1, 3, 4)   ! leg 3 LEFT -> outbound 4

        call snapshot_network(net)
        call evaluate_junctions(net)

        call assert(holding_cell(net, 1) == V_EMPTY,    'leg 1 should advance')
        call assert(holding_cell(net, 3) == V_EMPTY,    'leg 3 should advance')
        call assert(dest_cell(net, 2)    == V_OCCUPIED, 'leg 1 -> road 2')
        call assert(dest_cell(net, 4)    == V_OCCUPIED, 'leg 3 -> road 4')

        call free_network(net)
        print *, '  test_parallel_left_left OK'
    end subroutine test_parallel_left_left

    !-----------------------------------------------------------------
    ! 8. Mutual right turns from opposite legs: classic deadlock.
    !    Leg 1 RIGHT->4, leg 3 RIGHT->2.  Exactly one wins the tie-break.
    !-----------------------------------------------------------------
    subroutine test_mutual_right_deadlock()
        type(road_network_t) :: net
        integer :: n_moved
        call fresh(net)
        call plant(net, 1)
        call plant(net, 3)
        call force_route(net, 1, 1, 4)   ! leg 1 RIGHT -> outbound 4
        call force_route(net, 1, 3, 2)   ! leg 3 RIGHT -> outbound 2

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
    !    yields to its right.  Tie-break releases exactly one.
    !-----------------------------------------------------------------
    subroutine test_4way_straight_deadlock()
        type(road_network_t) :: net
        integer :: n_moved, k
        call fresh(net)
        do k = 1, 4
            call plant(net, k)
        end do
        call force_route(net, 1, 1, 3)   ! leg 1 STRAIGHT -> outbound 3
        call force_route(net, 1, 2, 4)   ! leg 2 STRAIGHT -> outbound 4
        call force_route(net, 1, 3, 1)   ! leg 3 STRAIGHT -> outbound 1
        call force_route(net, 1, 4, 2)   ! leg 4 STRAIGHT -> outbound 2

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
