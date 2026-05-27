module network_init_mod
    ! Hard-coded initialisers for canonical road-network topologies.
    !
    ! This module provides helper routines for constructing standard road
    ! networks used in the cellular automaton simulations. The routines allocate
    ! and initialise ``road_network_t`` objects with roads, lanes, junctions,
    ! route probabilities and open-boundary parameters.
    !
    ! Two canonical topologies are provided:
    !
    ! * ``init_crossroad``: a symmetric four-arm crossroad with one junction and
    !   four bidirectional roads;
    ! * ``init_t_junction``: an asymmetric T-junction with a west arm, south stem
    !   and east arm.
    !
    ! The module also provides ``init_lane`` for allocating an individual lane
    ! and ``free_network`` for deallocating a network created by these
    ! initialisers.
    !
    ! Route probabilities should form a valid probability distribution for each
    ! inbound leg:
    !
    ! .. math::
    !
    !    \sum_{m=1}^{n_{\mathrm{out}}} p_{km} = 1,
    !
    ! where ``p_km`` is the probability that a vehicle arriving on inbound leg
    ! ``k`` chooses outbound leg ``m``.

    use vehicle_mod
    use road_network_mod

    implicit none
    private

    public :: init_crossroad, init_t_junction, free_network, init_lane
    public :: place_evenly

contains

    subroutine init_crossroad(net, lane_length, alpha, beta, p_left, p_right)
        ! Initialise a symmetric four-arm crossroad network.
        !
        ! The crossroad contains one central junction and four roads arranged in
        ! clockwise order:
        !
        ! * road 1: north arm;
        ! * road 2: east arm;
        ! * road 3: south arm;
        ! * road 4: west arm.
        !
        ! Each road connects to the junction at ``end_junction(1)`` and is open
        ! at ``end_junction(2)``. Each road has two lanes:
        !
        ! * lane 1 flows away from the junction and has an open outflow boundary;
        ! * lane 2 flows towards the junction and has an open inflow boundary.
        !
        ! The route probabilities are constructed from ``p_left`` and
        ! ``p_right``. For inbound leg ``k``, the probabilities are assigned as:
        !
        ! .. math::
        !
        !    P(\mathrm{left}) = p_{\mathrm{left}},
        !
        ! .. math::
        !
        !    P(\mathrm{right}) = p_{\mathrm{right}},
        !
        ! and
        !
        ! .. math::
        !
        !    P(\mathrm{straight})
        !    =
        !    1 - p_{\mathrm{left}} - p_{\mathrm{right}}.
        !
        ! U-turns are assigned probability zero.
        type(road_network_t), intent(out) :: net
        ! Allocated crossroad network returned by the routine.
        integer, intent(in) :: lane_length
        ! Number of cells in each lane.
        real, intent(in) :: alpha
        ! Entry probability for inbound open boundaries.
        real, intent(in) :: beta
        ! Exit probability for outbound open boundaries.
        real, intent(in) :: p_left
        ! Probability of choosing the left-turn route at the junction.
        real, intent(in) :: p_right
        ! Probability of choosing the right-turn route at the junction.

        integer :: i
        ! Road, lane or route index.

        allocate(net%roads(4))
        allocate(net%junctions(1))

        net%junctions(1)%id    = 1
        net%junctions(1)%n_in  = 4
        net%junctions(1)%n_out = 4

        allocate(net%junctions(1)%in_road(4),  net%junctions(1)%in_lane(4))
        allocate(net%junctions(1)%out_road(4), net%junctions(1)%out_lane(4))

        ! Roads are stored in clockwise order: 1=N, 2=E, 3=S, 4=W.
        net%junctions(1)%in_road  = [1, 2, 3, 4]
        net%junctions(1)%out_road = [1, 2, 3, 4]

        do i = 1, 4
            net%roads(i)%id           = i
            net%roads(i)%end_junction = [1, 0]

            allocate(net%roads(i)%lane(2))

            call init_lane(net%roads(i)%lane(1), lane_length, +1)
            call init_lane(net%roads(i)%lane(2), lane_length, -1)

            net%roads(i)%lane(2)%open_in  = .true.
            net%roads(i)%lane(2)%alpha    = alpha

            net%roads(i)%lane(1)%open_out = .true.
            net%roads(i)%lane(1)%beta     = beta

            net%junctions(1)%in_lane(i)  = inbound_lane_at_end(net%roads(i), 1)
            net%junctions(1)%out_lane(i) = outbound_lane_at_end(net%roads(i), 1)
        end do

        ! Build one outbound probability row for each inbound leg.
        allocate(net%junctions(1)%in_routes(4))

        do i = 1, 4
            allocate(net%junctions(1)%in_routes(i)%prob(4))

            net%junctions(1)%in_routes(i)%prob                 = 0.0
            net%junctions(1)%in_routes(i)%prob(mod(i,   4)+1) = p_left
            net%junctions(1)%in_routes(i)%prob(mod(i+1, 4)+1) = 1.0 - p_left - p_right
            net%junctions(1)%in_routes(i)%prob(mod(i+2, 4)+1) = p_right
        end do

        ! Perimeter-port arrays are not needed for the symmetric crossroad.
    end subroutine init_crossroad


    subroutine init_t_junction(net, lane_length, alpha, beta, route_in_west, route_in_stem)
        ! Initialise an asymmetric T-junction network.
        !
        ! The T-junction has one junction and three road arms:
        !
        ! * road 1: west arm, bidirectional with two lanes;
        ! * road 2: south stem, inbound only with one lane;
        ! * road 3: east arm, outbound only with one lane.
        !
        ! The junction has two inbound legs and two outbound legs:
        !
        ! * inbound leg 1: west-inbound lane;
        ! * inbound leg 2: stem-inbound lane;
        ! * outbound leg 1: west-outbound lane;
        ! * outbound leg 2: east-outbound lane.
        !
        ! The route arrays define the outbound probabilities for the two inbound
        ! legs:
        !
        ! .. math::
        !
        !    \mathbf{p}_{\mathrm{west}}
        !    =
        !    (P(\mathrm{west\_out}), P(\mathrm{east\_out})),
        !
        ! and
        !
        ! .. math::
        !
        !    \mathbf{p}_{\mathrm{stem}}
        !    =
        !    (P(\mathrm{west\_out}), P(\mathrm{east\_out})).
        !
        ! Perimeter-port positions are allocated because this is an asymmetric
        ! junction. The zero-indexed clockwise perimeter convention is:
        !
        ! * ``0``: east_out;
        ! * ``1``: stem_in;
        ! * ``2``: west_in;
        ! * ``3``: west_out.
        type(road_network_t), intent(out) :: net
        ! Allocated T-junction network returned by the routine.
        integer, intent(in) :: lane_length
        ! Number of cells in each lane.
        real, intent(in) :: alpha
        ! Entry probability for open inflow boundaries.
        real, intent(in) :: beta
        ! Exit probability for open outflow boundaries.
        real, intent(in) :: route_in_west(2)
        ! Route probabilities for vehicles entering from the west arm.
        real, intent(in) :: route_in_stem(2)
        ! Route probabilities for vehicles entering from the south stem.

        allocate(net%roads(3))
        allocate(net%junctions(1))

        net%junctions(1)%id    = 1
        net%junctions(1)%n_in  = 2
        net%junctions(1)%n_out = 2

        allocate(net%junctions(1)%in_road(2),  net%junctions(1)%in_lane(2))
        allocate(net%junctions(1)%out_road(2), net%junctions(1)%out_lane(2))
        allocate(net%junctions(1)%in_perim(2), net%junctions(1)%out_perim(2))

        ! Inbound legs: west_in, stem_in. Outbound legs: west_out, east_out.
        net%junctions(1)%in_road  = [1, 2]
        net%junctions(1)%out_road = [1, 3]

        ! Road 1: west arm, bidirectional.
        net%roads(1)%id           = 1
        net%roads(1)%end_junction = [0, 1]

        allocate(net%roads(1)%lane(2))

        call init_lane(net%roads(1)%lane(1), lane_length, +1)
        call init_lane(net%roads(1)%lane(2), lane_length, -1)

        net%roads(1)%lane(1)%open_in  = .true.
        net%roads(1)%lane(1)%alpha    = alpha

        net%roads(1)%lane(2)%open_out = .true.
        net%roads(1)%lane(2)%beta     = beta

        ! Road 2: south stem, inbound only.
        net%roads(2)%id           = 2
        net%roads(2)%end_junction = [0, 1]

        allocate(net%roads(2)%lane(1))

        call init_lane(net%roads(2)%lane(1), lane_length, +1)

        net%roads(2)%lane(1)%open_in = .true.
        net%roads(2)%lane(1)%alpha   = alpha

        ! Road 3: east arm, outbound only.
        net%roads(3)%id           = 3
        net%roads(3)%end_junction = [1, 0]

        allocate(net%roads(3)%lane(1))

        call init_lane(net%roads(3)%lane(1), lane_length, +1)

        net%roads(3)%lane(1)%open_out = .true.
        net%roads(3)%lane(1)%beta     = beta

        ! Junction flat lane lists.
        net%junctions(1)%in_lane(1)  = inbound_lane_at_end(net%roads(1), 2)
        net%junctions(1)%in_lane(2)  = inbound_lane_at_end(net%roads(2), 2)
        net%junctions(1)%out_lane(1) = outbound_lane_at_end(net%roads(1), 2)
        net%junctions(1)%out_lane(2) = outbound_lane_at_end(net%roads(3), 1)

        ! Perimeter port positions, zero-indexed and clockwise.
        net%junctions(1)%in_perim  = [2, 1]
        net%junctions(1)%out_perim = [3, 0]

        ! Route probability distributions.
        allocate(net%junctions(1)%in_routes(2))
        allocate(net%junctions(1)%in_routes(1)%prob(2))
        allocate(net%junctions(1)%in_routes(2)%prob(2))

        net%junctions(1)%in_routes(1)%prob = route_in_west
        net%junctions(1)%in_routes(2)%prob = route_in_stem
    end subroutine init_t_junction


    subroutine init_lane(ln, L, fd)
        ! Allocate and initialise a single lane.
        !
        ! This routine sets the lane length and flow direction, allocates both the
        ! live cell array and the old-state snapshot array, and initialises all
        ! cells to be empty with zero velocity.
        !
        ! The old-state array is used by parallel cellular automaton updates:
        !
        ! .. math::
        !
        !    \mathrm{old} \leftarrow \mathrm{cells}.
        !
        ! At initialisation both arrays are empty.
        type(lane_t), intent(out) :: ln
        ! Lane object to allocate and initialise.
        integer, intent(in) :: L
        ! Number of cells in the lane.
        integer, intent(in) :: fd
        ! Flow direction: ``+1`` or ``-1``.

        ln%length         = L
        ln%flow_direction = fd

        allocate(ln%cells(L), ln%old(L))

        ln%cells%has_car  = .false.
        ln%cells%velocity = 0

        ln%old%has_car    = .false.
        ln%old%velocity   = 0
    end subroutine init_lane


    subroutine place_evenly(ln, n_vehicles)
        ! Pre-populate a lane with ``n_vehicles`` cars at evenly-spaced sites
        ! and zero velocity.
        !
        ! Even spacing gives a deterministic starting density and avoids the
        ! initial transient that follows random or clumped placements. The
        ! placement is:
        !
        ! .. math::
        !
        !    i_k = 1 + \left\lfloor \frac{(k - 1)\, L}{n} \right\rfloor,
        !    \quad k = 1, \dots, n,
        !
        ! where ``L`` is the lane length and ``n`` is the (clamped) vehicle
        ! count. ``n_vehicles`` is clipped to ``[0, length]`` before placement,
        ! so passing 0 (or a negative value) is a no-op and passing more
        ! vehicles than cells fills the lane completely. Existing lane state is
        ! overwritten — any previous occupancy is cleared first.
        !
        ! Single source of truth used both by the periodic-ring
        ! fundamental-diagram measurement (``fundamental_diagram_mod``) and by
        ! the network builder when a lane specifies a non-zero ``n_vehicles``
        ! (typically a periodic NS ring set up via ``single_lane_periodic`` on
        ! the Python side).
        type(lane_t), intent(inout) :: ln
        ! Lane to populate. Must have its ``cells`` array already allocated
        ! (e.g. by ``init_lane``).
        integer, intent(in) :: n_vehicles
        ! Target number of vehicles, clipped to ``[0, ln%length]``.

        integer :: n, k, idx, L

        L = ln%length
        ln%cells%has_car  = .false.
        ln%cells%velocity = 0

        n = max(0, min(n_vehicles, L))
        if (n == 0) return

        do k = 1, n
            idx = 1 + ((k - 1) * L) / n
            if (idx > L) idx = L
            ln%cells(idx)%has_car  = .true.
            ln%cells(idx)%velocity = 0
        end do
    end subroutine place_evenly


    subroutine free_network(net)
        ! Deallocate all allocatable arrays owned by a road network.
        !
        ! This routine releases roads, lanes, cell arrays, junction arrays and
        ! route-probability arrays. It is intended for cleanup after a simulation
        ! or test case.
        !
        ! The routine checks each allocatable before deallocation, so it is safe
        ! to call on partially allocated networks created during setup or testing.
        type(road_network_t), intent(inout) :: net
        ! Network to deallocate in place.

        integer :: r
        ! Road or junction index.
        integer :: ln
        ! Lane index.
        integer :: k
        ! Route-row index.

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