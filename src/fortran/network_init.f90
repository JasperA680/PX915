module network_init_mod
    !--------------------------------------------------------------------
    ! Hard-coded initialisers for canonical road-network topologies.
    !
    ! init_crossroad builds a 4-arm crossroad with one central junction.
    ! Roads are placed in strict clockwise order about the junction
    ! (1=N, 2=E, 3=S, 4=W).  Every road connects to the junction at its
    ! end_a (end 1); end_b (end 2) is an open boundary with the supplied
    ! alpha (entry) and beta (exit) probabilities.
    !--------------------------------------------------------------------
    use vehicle_mod
    use road_network_mod
    implicit none
    private

    public :: init_crossroad, free_network

contains

    subroutine init_crossroad(net, lane_length, alpha, beta, p_left, p_right)
        type(road_network_t), intent(out) :: net
        integer, intent(in) :: lane_length
        real,    intent(in) :: alpha, beta, p_left, p_right
        integer :: i

        allocate(net%roads(4))
        allocate(net%junctions(1))

        net%junctions(1)%id     = 1
        net%junctions(1)%n_legs = 4
        allocate(net%junctions(1)%connected_road_ids(4))
        allocate(net%junctions(1)%end_at(4))
        net%junctions(1)%connected_road_ids = [1, 2, 3, 4]
        net%junctions(1)%end_at             = [1, 1, 1, 1]

        do i = 1, 4
            net%roads(i)%id           = i
            net%roads(i)%end_junction = [1, 0]

            allocate(net%roads(i)%lane(2))

            call init_lane(net%roads(i)%lane(1), lane_length, +1)  ! end_1->end_2, outbound from junction
            call init_lane(net%roads(i)%lane(2), lane_length, -1)  ! end_2->end_1, inbound to junction

            ! Inbound lane: site 1 sits at the open end_b -> open inflow with alpha.
            net%roads(i)%lane(2)%open_in  = .true.
            net%roads(i)%lane(2)%alpha    = alpha
            net%roads(i)%lane(2)%p_left   = p_left
            net%roads(i)%lane(2)%p_right  = p_right

            ! Outbound lane: site L sits at the open end_b -> open outflow with beta.
            net%roads(i)%lane(1)%open_out = .true.
            net%roads(i)%lane(1)%beta     = beta
        end do
    end subroutine init_crossroad

    subroutine init_lane(ln, L, fd)
        type(lane_t), intent(out) :: ln
        integer, intent(in) :: L, fd    ! fd: flow_direction (+1 or -1)
        ln%length         = L
        ln%flow_direction = fd
        allocate(ln%cells(L), ln%old(L))
        ln%cells%has_car = .false.
        ln%old%has_car   = .false.
    end subroutine init_lane

    subroutine free_network(net)
        type(road_network_t), intent(inout) :: net
        integer :: r, ln
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
                if (allocated(net%junctions(r)%connected_road_ids)) deallocate(net%junctions(r)%connected_road_ids)
                if (allocated(net%junctions(r)%end_at))             deallocate(net%junctions(r)%end_at)
            end do
            deallocate(net%junctions)
        end if
    end subroutine free_network

end module network_init_mod
