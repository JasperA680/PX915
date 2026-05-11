program test_types
    !--------------------------------------------------------------------
    ! Phase 3 smoke test: build a road_network_t with 4 roads x 2 lanes
    ! x L cells, exercise V_EMPTY/V_OCCUPIED, the lane-at-end helpers,
    ! count_occupied_network, and snapshot_network.
    !--------------------------------------------------------------------
    use vehicle_mod
    use road_network_mod
    implicit none

    integer, parameter :: L = 100

    type(road_network_t) :: net
    integer :: r, ln

    !------ Allocate a 4-road, 1-junction network ----------------------
    allocate(net%roads(4))
    allocate(net%junctions(1))

    net%junctions(1)%id   = 1
    net%junctions(1)%n_in = 4
    net%junctions(1)%n_out = 4
    allocate(net%junctions(1)%in_road(4),  net%junctions(1)%in_lane(4))
    allocate(net%junctions(1)%out_road(4), net%junctions(1)%out_lane(4))
    net%junctions(1)%in_road  = [1, 2, 3, 4]
    net%junctions(1)%in_lane  = [2, 2, 2, 2]
    net%junctions(1)%out_road = [1, 2, 3, 4]
    net%junctions(1)%out_lane = [1, 1, 1, 1]

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
            net%roads(r)%lane(ln)%cells%velocity = 0
            net%roads(r)%lane(ln)%old%has_car   = .false.
            net%roads(r)%lane(ln)%old%velocity = 0
        end do
    end do

    !------ is_occupied --------------------------------------------------
    call assert(.not. is_occupied(V_EMPTY),    'V_EMPTY should not be occupied')
    call assert(      is_occupied(V_OCCUPIED), 'V_OCCUPIED should be occupied')

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
    net%roads(3)%lane(2)%cells(L)%has_car = .true.
    call assert(count_occupied_network(net) == 2, 'two V_OCCUPIED cells not counted')

    !------ snapshot_network copies cells -> old ------------------------
    call snapshot_network(net)
    call assert(net%roads(1)%lane(2)%old(L)%has_car .eqv. .true., 'snapshot road 1')
    call assert(net%roads(3)%lane(2)%old(L)%has_car .eqv. .true., 'snapshot road 3')

    !------ Tidy up -----------------------------------------------------
    do r = 1, 4
        do ln = 1, size(net%roads(r)%lane)
            deallocate(net%roads(r)%lane(ln)%cells)
            deallocate(net%roads(r)%lane(ln)%old)
        end do
        deallocate(net%roads(r)%lane)
    end do
    deallocate(net%junctions(1)%in_road,  net%junctions(1)%in_lane)
    deallocate(net%junctions(1)%out_road, net%junctions(1)%out_lane)
    deallocate(net%roads)
    deallocate(net%junctions)

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
