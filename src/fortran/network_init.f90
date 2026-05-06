module network_init_mod
    !--------------------------------------------------------------------
    ! Hard-coded initialisers for canonical road-network topologies.
    !
    ! init_crossroad: 4-arm crossroad, 1 junction, roads in CW order
    !   (1=N, 2=E, 3=S, 4=W), each connecting at end_1, open at end_2.
    !   Routing built from p_left / p_right: LEFT=mod(k,4)+1,
    !   STRAIGHT=mod(k+1,4)+1, RIGHT=mod(k+2,4)+1, UTURN=k (prob=0).
    !
    ! init_t_junction: T-junction (West arm bidirectional, South stem
    !   inbound-only, East arm outbound-only).
    !     Road 1 (west): 2 lanes, end_2 at junction.
    !     Road 2 (stem): 1 inbound lane, end_2 at junction.
    !     Road 3 (east): 1 outbound lane, end_1 at junction.
    !   n_in = 2, n_out = 2.  Perimeter ports (CW, 0-indexed):
    !     0=east_out, 1=stem_in, 2=west_in, 3=west_out.
    !--------------------------------------------------------------------
    use vehicle_mod
    use road_network_mod
    implicit none
    private

    public :: init_crossroad, init_t_junction, free_network

contains

    subroutine init_crossroad(net, lane_length, alpha, beta, p_left, p_right)
        type(road_network_t), intent(out) :: net
        integer, intent(in) :: lane_length
        real,    intent(in) :: alpha, beta, p_left, p_right
        integer :: i

        allocate(net%roads(4))
        allocate(net%junctions(1))

        net%junctions(1)%id   = 1
        net%junctions(1)%n_in = 4
        net%junctions(1)%n_out = 4
        allocate(net%junctions(1)%in_road(4),  net%junctions(1)%in_lane(4))
        allocate(net%junctions(1)%out_road(4), net%junctions(1)%out_lane(4))
        ! Roads in clockwise order (1=N, 2=E, 3=S, 4=W); all at end_1.
        net%junctions(1)%in_road  = [1, 2, 3, 4]
        net%junctions(1)%out_road = [1, 2, 3, 4]

        do i = 1, 4
            net%roads(i)%id           = i
            net%roads(i)%end_junction = [1, 0]

            allocate(net%roads(i)%lane(2))

            call init_lane(net%roads(i)%lane(1), lane_length, +1)  ! end_1->end_2, outbound
            call init_lane(net%roads(i)%lane(2), lane_length, -1)  ! end_2->end_1, inbound

            net%roads(i)%lane(2)%open_in  = .true.
            net%roads(i)%lane(2)%alpha    = alpha
            net%roads(i)%lane(1)%open_out = .true.
            net%roads(i)%lane(1)%beta     = beta

            net%junctions(1)%in_lane(i)  = inbound_lane_at_end(net%roads(i), 1)
            net%junctions(1)%out_lane(i) = outbound_lane_at_end(net%roads(i), 1)
        end do

        ! Per-inbound-leg route probability distributions.
        ! LEFT=mod(k,4)+1, STRAIGHT=mod(k+1,4)+1, RIGHT=mod(k+2,4)+1, UTURN=k (prob=0).
        allocate(net%junctions(1)%in_routes(4))
        do i = 1, 4
            allocate(net%junctions(1)%in_routes(i)%prob(4))
            net%junctions(1)%in_routes(i)%prob               = 0.0
            net%junctions(1)%in_routes(i)%prob(mod(i,   4)+1) = p_left
            net%junctions(1)%in_routes(i)%prob(mod(i+1, 4)+1) = 1.0 - p_left - p_right
            net%junctions(1)%in_routes(i)%prob(mod(i+2, 4)+1) = p_right
        end do
        ! in_perim / out_perim not allocated for the symmetric crossroad.
    end subroutine init_crossroad

    !-----------------------------------------------------------------
    ! T-junction: West arm (bidirectional), South stem (inbound only),
    ! East arm (outbound only).  n_in = 2, n_out = 2.
    !
    ! route_in_west(1:2) : probability vector for west-inbound vehicle;
    !   index 1 = west_out (U-turn back, set to 0 normally),
    !   index 2 = east_out.
    ! route_in_stem(1:2) : probability vector for stem-inbound vehicle;
    !   index 1 = west_out, index 2 = east_out.
    !
    ! Perimeter ports (CW, 0-indexed, total 4):
    !   0 = east_out, 1 = stem_in, 2 = west_in, 3 = west_out.
    !-----------------------------------------------------------------
    subroutine init_t_junction(net, lane_length, alpha, beta, route_in_west, route_in_stem)
        type(road_network_t), intent(out) :: net
        integer, intent(in) :: lane_length
        real,    intent(in) :: alpha, beta
        real,    intent(in) :: route_in_west(2), route_in_stem(2)

        allocate(net%roads(3))
        allocate(net%junctions(1))

        net%junctions(1)%id   = 1
        net%junctions(1)%n_in = 2
        net%junctions(1)%n_out = 2
        allocate(net%junctions(1)%in_road(2),  net%junctions(1)%in_lane(2))
        allocate(net%junctions(1)%out_road(2), net%junctions(1)%out_lane(2))
        allocate(net%junctions(1)%in_perim(2), net%junctions(1)%out_perim(2))
        ! in: [west_in=leg1, stem_in=leg2]; out: [west_out=outbound1, east_out=outbound2]
        net%junctions(1)%in_road  = [1, 2]
        net%junctions(1)%out_road = [1, 3]

        ! Road 1 (west arm): bidirectional, 2 lanes.
        ! end_junction = [0, 1]: end_2 is the junction end.
        net%roads(1)%id           = 1
        net%roads(1)%end_junction = [0, 1]
        allocate(net%roads(1)%lane(2))
        call init_lane(net%roads(1)%lane(1), lane_length, +1)  ! inbound: end_1->end_2 (flow +1)
        call init_lane(net%roads(1)%lane(2), lane_length, -1)  ! outbound: end_2->end_1 (flow -1)
        net%roads(1)%lane(1)%open_in  = .true.
        net%roads(1)%lane(1)%alpha    = alpha
        net%roads(1)%lane(2)%open_out = .true.
        net%roads(1)%lane(2)%beta     = beta

        ! Road 2 (south stem): inbound only, 1 lane.
        ! end_junction = [0, 1]: end_2 is the junction end.
        net%roads(2)%id           = 2
        net%roads(2)%end_junction = [0, 1]
        allocate(net%roads(2)%lane(1))
        call init_lane(net%roads(2)%lane(1), lane_length, +1)  ! inbound: end_1->end_2
        net%roads(2)%lane(1)%open_in  = .true.
        net%roads(2)%lane(1)%alpha    = alpha

        ! Road 3 (east arm): outbound only, 1 lane.
        ! end_junction = [1, 0]: end_1 is the junction end.
        net%roads(3)%id           = 3
        net%roads(3)%end_junction = [1, 0]
        allocate(net%roads(3)%lane(1))
        call init_lane(net%roads(3)%lane(1), lane_length, +1)  ! outbound: end_1->end_2
        net%roads(3)%lane(1)%open_out = .true.
        net%roads(3)%lane(1)%beta     = beta

        ! Junction flat lane lists.
        net%junctions(1)%in_lane(1)  = inbound_lane_at_end(net%roads(1), 2)   ! road1, end_2 in
        net%junctions(1)%in_lane(2)  = inbound_lane_at_end(net%roads(2), 2)   ! road2, end_2 in
        net%junctions(1)%out_lane(1) = outbound_lane_at_end(net%roads(1), 2)  ! road1, end_2 out
        net%junctions(1)%out_lane(2) = outbound_lane_at_end(net%roads(3), 1)  ! road3, end_1 out

        ! Perimeter port positions (0-indexed, CW, total 4):
        !   0 = east_out (out_idx=2), 1 = stem_in (in_idx=2),
        !   2 = west_in  (in_idx=1),  3 = west_out (out_idx=1).
        net%junctions(1)%in_perim  = [2, 1]   ! west_in=2, stem_in=1
        net%junctions(1)%out_perim = [3, 0]   ! west_out=3, east_out=0

        ! Route probability distributions.
        allocate(net%junctions(1)%in_routes(2))
        allocate(net%junctions(1)%in_routes(1)%prob(2))
        allocate(net%junctions(1)%in_routes(2)%prob(2))
        net%junctions(1)%in_routes(1)%prob = route_in_west
        net%junctions(1)%in_routes(2)%prob = route_in_stem
    end subroutine init_t_junction

    subroutine init_lane(ln, L, fd)
        type(lane_t), intent(out) :: ln
        integer, intent(in) :: L, fd
        ln%length         = L
        ln%flow_direction = fd
        allocate(ln%cells(L), ln%old(L))
        ln%cells%has_car = .false.
        ln%cells%turning_intent = V_EMPTY
        ln%cells%velocity = 0
        ln%old%has_car   = .false.
        ln%old%turning_intent = V_EMPTY
        ln%old%velocity = 0
    end subroutine init_lane

    subroutine free_network(net)
        type(road_network_t), intent(inout) :: net
        integer :: r, ln, k
        if (allocated(net%roads)) then
            do r = 1, size(net%roads)
                if (allocated(net%roads(r)%lane)) then
                    do ln = 1, size(net%roads(r)%lane)
                        if (allocated(net%roads(r)%lane(ln)%cells)) deallocate(net%roads(r)%lane(ln)%cells)
                        if (allocated(net%roads(r)%lane(ln)%old))   deallocate(net%roads(r)%lane(ln)%old)
                    end do
                    deallocate(net%roads(r)%lane)
                end if
            end do
            deallocate(net%roads)
        end if
        if (allocated(net%junctions)) then
            do r = 1, size(net%junctions)
                if (allocated(net%junctions(r)%in_road))   deallocate(net%junctions(r)%in_road)
                if (allocated(net%junctions(r)%in_lane))   deallocate(net%junctions(r)%in_lane)
                if (allocated(net%junctions(r)%out_road))  deallocate(net%junctions(r)%out_road)
                if (allocated(net%junctions(r)%out_lane))  deallocate(net%junctions(r)%out_lane)
                if (allocated(net%junctions(r)%in_perim))  deallocate(net%junctions(r)%in_perim)
                if (allocated(net%junctions(r)%out_perim)) deallocate(net%junctions(r)%out_perim)
                if (allocated(net%junctions(r)%in_routes)) then
                    do k = 1, size(net%junctions(r)%in_routes)
                        if (allocated(net%junctions(r)%in_routes(k)%prob)) &
                            deallocate(net%junctions(r)%in_routes(k)%prob)
                    end do
                    deallocate(net%junctions(r)%in_routes)
                end if
            end do
            deallocate(net%junctions)
        end if
    end subroutine free_network

end module network_init_mod
