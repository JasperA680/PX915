program test_types
    !--------------------------------------------------------------------
    ! Phase 1 smoke test: build a road_network_t with 4 roads x 2 lanes
    ! x L cells, exercise the V_* codes and helper pure functions, then
    ! run a Monte Carlo check on sample_indicator's distribution.
    !--------------------------------------------------------------------
    use vehicle_mod
    use road_network_mod
    implicit none

    integer, parameter :: L = 100
    integer, parameter :: N_SAMPLES = 100000
    real,    parameter :: TOL = 0.01

    type(road_network_t) :: net
    integer :: r, ln, k, code
    integer :: counts(0:3)
    real    :: pL, pR, pS
    integer :: seed_size
    integer, allocatable :: seed(:)

    !------ Allocate a 4-road, 1-junction network ----------------------
    allocate(net%roads(4))
    allocate(net%junctions(1))

    net%junctions(1)%id     = 1
    net%junctions(1)%n_legs = 4
    allocate(net%junctions(1)%connected_road_ids(4))
    allocate(net%junctions(1)%end_at(4))
    net%junctions(1)%connected_road_ids = [1, 2, 3, 4]
    net%junctions(1)%end_at             = [1, 1, 1, 1]

    do r = 1, 4
        net%roads(r)%id           = r
        net%roads(r)%end_junction = [1, 0]
        allocate(net%roads(r)%lane(2))
        net%roads(r)%lane(1)%flow_direction = +1   ! end_1 -> end_2
        net%roads(r)%lane(2)%flow_direction = -1   ! end_2 -> end_1
        do ln = 1, 2
            net%roads(r)%lane(ln)%length = L
            allocate(net%roads(r)%lane(ln)%cells(L))
            allocate(net%roads(r)%lane(ln)%old(L))
            net%roads(r)%lane(ln)%cells%has_car = .false.
            net%roads(r)%lane(ln)%cells%turning_intent = V_EMPTY
            net%roads(r)%lane(ln)%cells%velocity = 0
            net%roads(r)%lane(ln)%old%has_car   = .false.
            net%roads(r)%lane(ln)%old%turning_intent = V_EMPTY
            net%roads(r)%lane(ln)%old%velocity = 0
        end do
    end do

    !------ is_occupied --------------------------------------------------
    call assert(.not. is_occupied(V_EMPTY),    'V_EMPTY occupied')
    call assert(      is_occupied(V_STRAIGHT), 'V_STRAIGHT not occupied')
    call assert(      is_occupied(V_LEFT),     'V_LEFT not occupied')
    call assert(      is_occupied(V_RIGHT),    'V_RIGHT not occupied')

    !------ Lane-at-end helpers (replaces old inbound/outbound_lane_idx) -
    ! For end 1: inbound lane has flow_direction=-1, which is lane 2.
    !            outbound lane has flow_direction=+1, which is lane 1.
    call assert(inbound_lane_at_end(net%roads(1), 1)  == 2, 'inbound_lane_at_end(road,1)')
    call assert(inbound_lane_at_end(net%roads(1), 2)  == 1, 'inbound_lane_at_end(road,2)')
    call assert(outbound_lane_at_end(net%roads(1), 1) == 1, 'outbound_lane_at_end(road,1)')
    call assert(outbound_lane_at_end(net%roads(1), 2) == 2, 'outbound_lane_at_end(road,2)')

    !------ count_occupied_network on an empty network ------------------
    call assert(count_occupied_network(net) == 0, 'empty network not zero')

    net%roads(1)%lane(2)%cells(L)%has_car = .true.
    net%roads(1)%lane(2)%cells(L)%turning_intent = V_STRAIGHT
    net%roads(1)%lane(2)%cells(L)%velocity = 0
    net%roads(3)%lane(2)%cells(L)%has_car = .true.
    net%roads(3)%lane(2)%cells(L)%turning_intent = V_LEFT
    net%roads(3)%lane(2)%cells(L)%velocity = 0
    call assert(count_occupied_network(net) == 2, 'two occupied not counted')

    !------ snapshot_network copies cells -> old ------------------------
    call snapshot_network(net)
    call assert(net%roads(1)%lane(2)%old(L)%has_car, 'snapshot road 1 occupied')
    call assert(net%roads(1)%lane(2)%old(L)%turning_intent == V_STRAIGHT, 'snapshot road 1')
    call assert(net%roads(3)%lane(2)%old(L)%has_car, 'snapshot road 3 occupied')
    call assert(net%roads(3)%lane(2)%old(L)%turning_intent == V_LEFT,     'snapshot road 3')

    !------ sample_indicator distribution -------------------------------
    call random_seed(size = seed_size)
    allocate(seed(seed_size))
    seed = 12345
    call random_seed(put = seed)

    counts = 0
    do k = 1, N_SAMPLES
        code = sample_indicator(0.25, 0.25)
        counts(code) = counts(code) + 1
    end do
    pL = real(counts(V_LEFT))     / real(N_SAMPLES)
    pR = real(counts(V_RIGHT))    / real(N_SAMPLES)
    pS = real(counts(V_STRAIGHT)) / real(N_SAMPLES)

    print '(a, f7.4, a, f7.4, a, f7.4)', &
        'P(left)=', pL, '  P(right)=', pR, '  P(straight)=', pS
    call assert(abs(pL - 0.25) < TOL, 'P(left) out of tolerance')
    call assert(abs(pR - 0.25) < TOL, 'P(right) out of tolerance')
    call assert(abs(pS - 0.50) < TOL, 'P(straight) out of tolerance')

    !------ Tidy up -----------------------------------------------------
    do r = 1, 4
        do ln = 1, size(net%roads(r)%lane)
            deallocate(net%roads(r)%lane(ln)%cells)
            deallocate(net%roads(r)%lane(ln)%old)
        end do
        deallocate(net%roads(r)%lane)
    end do
    deallocate(net%junctions(1)%connected_road_ids)
    deallocate(net%junctions(1)%end_at)
    deallocate(net%roads)
    deallocate(net%junctions)
    deallocate(seed)

    print *, 'test_types: OK'

contains

    subroutine assert(cond, msg)
        logical, intent(in) :: cond
        character(len=*), intent(in) :: msg
        if (.not. cond) then
            print *, 'ASSERT FAILED: ', trim(msg)
            stop 1
        end if
    end subroutine assert

end program test_types
