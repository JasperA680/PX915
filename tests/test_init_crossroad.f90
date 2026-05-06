program test_init_crossroad
    !--------------------------------------------------------------------
    ! Phase 2 validation: build a crossroad via init_crossroad and verify
    ! every structural invariant the rest of the pipeline relies on.
    !--------------------------------------------------------------------
    use vehicle_mod
    use road_network_mod
    use network_init_mod
    implicit none

    integer, parameter :: L = 25
    type(road_network_t) :: net
    integer :: r, e

    call init_crossroad(net, L, 0.4, 0.5, 0.25, 0.25)

    !---- Container sizes ----------------------------------------------
    call assert(size(net%roads)     == 4, 'wrong number of roads')
    call assert(size(net%junctions) == 1, 'wrong number of junctions')

    !---- Junction wiring ----------------------------------------------
    call assert(net%junctions(1)%n_legs == 4, 'n_legs')
    call assert(all(net%junctions(1)%connected_road_ids == [1,2,3,4]), 'clockwise ids')
    call assert(all(net%junctions(1)%end_at             == [1,1,1,1]), 'end_at all 1')

    !---- Per-road wiring ----------------------------------------------
    do r = 1, 4
        call assert(net%roads(r)%id == r,            'road id mismatch')
        call assert(all(net%roads(r)%end_junction == [1, 0]), 'end_junction')

        call assert(size(net%roads(r)%lane) == 2,    'lane count not 2')

        call assert(net%roads(r)%lane(1)%length == L, 'lane 1 length')
        call assert(net%roads(r)%lane(2)%length == L, 'lane 2 length')

        call assert(net%roads(r)%lane(1)%flow_direction == +1, 'lane 1 flow_direction')
        call assert(net%roads(r)%lane(2)%flow_direction == -1, 'lane 2 flow_direction')

        call assert(all(.not. net%roads(r)%lane(1)%cells%has_car), 'lane 1 not empty')
        call assert(all(.not. net%roads(r)%lane(2)%cells%has_car), 'lane 2 not empty')

        ! Lane 2 is inbound to the junction (open inflow at site 1, far end).
        call assert(      net%roads(r)%lane(2)%open_in,  'lane 2 open_in')
        call assert(.not. net%roads(r)%lane(2)%open_out, 'lane 2 open_out')
        call assert(abs(net%roads(r)%lane(2)%alpha - 0.4) < 1e-6, 'alpha')

        ! Lane 1 is outbound from the junction (open outflow at site L, far end).
        call assert(      net%roads(r)%lane(1)%open_out, 'lane 1 open_out')
        call assert(.not. net%roads(r)%lane(1)%open_in,  'lane 1 open_in')
        call assert(abs(net%roads(r)%lane(1)%beta - 0.5) < 1e-6,  'beta')
    end do

    !---- lane_at_end helpers consistent with allocation order ---------
    do r = 1, 4
        e = net%junctions(1)%end_at(r)
        call assert(inbound_lane_at_end(net%roads(r), e)  == 2, 'inbound_lane_at_end not 2')
        call assert(outbound_lane_at_end(net%roads(r), e) == 1, 'outbound_lane_at_end not 1')
    end do

    !---- Print connectivity for human eyeballing -----------------------
    print *, 'Crossroad layout (clockwise):'
    print *, '  road 1  N  end_a -> junction'
    print *, '  road 2  E  end_a -> junction'
    print *, '  road 3  S  end_a -> junction'
    print *, '  road 4  W  end_a -> junction'
    print *, '  inbound lane = 2 (site L = junction holding cell)'
    print *, '  outbound lane = 1 (site 1 = junction destination)'

    call free_network(net)
    print *, 'test_init_crossroad: OK'

contains

    subroutine assert(cond, msg)
        logical, intent(in) :: cond
        character(len=*), intent(in) :: msg
        if (.not. cond) then
            print *, 'ASSERT FAILED: ', trim(msg)
            stop 1
        end if
    end subroutine assert

end program test_init_crossroad
