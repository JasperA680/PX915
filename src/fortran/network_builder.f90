module network_builder_mod
    ! Build a road-network object from a flat network specification.
    !
    ! This module defines lightweight specification types that mirror the JSON
    ! configuration schema used by the road-network model. The main routine,
    ! ``build_network``, converts a ``network_spec_t`` into a fully allocated
    ! ``road_network_t`` that can be advanced by the network time-step routines.
    !
    ! The specification is intentionally flatter than the runtime network object.
    ! It contains:
    !
    ! * road specifications, including lane lists and boundary parameters;
    ! * junction specifications, including inbound legs, outbound legs and route
    !   probabilities;
    ! * optional perimeter-port indices for asymmetric junction conflict tests.
    !
    ! The route probabilities define a row-stochastic routing table for each
    ! junction:
    !
    ! .. math::
    !
    !    \sum_{m=1}^{n_{\mathrm{out}}} p_{km} = 1,
    !
    ! where ``p_{km}`` is the probability that a vehicle arriving from inbound
    ! leg ``k`` chooses outbound leg ``m``.
    !
    ! The builder reuses ``init_lane`` from ``network_init_mod`` so that lane
    ! allocation and initialisation are consistent with the rest of the network
    ! code.

    use road_network_mod
    use network_init_mod, only: init_lane

    implicit none
    private

    type, public :: lane_spec_t
        ! Specification for a single lane in a road.
        !
        ! This type stores the lane length, flow direction and open-boundary
        ! parameters read from the network configuration.

        integer :: length = 0
        ! Number of cells in the lane.
        integer :: flow_direction = 0
        ! Direction of traffic flow along the lane.
        real :: alpha = 0.0
        ! Entry probability for an open left/inflow boundary.
        real :: beta = 0.0
        ! Exit probability for an open right/outflow boundary.
        logical :: open_in = .false.
        ! Whether the lane has an open inflow boundary.
        logical :: open_out = .false.
        ! Whether the lane has an open outflow boundary.
    end type lane_spec_t


    type, public :: road_spec_t
        ! Specification for a road and its lanes.
        !
        ! A road connects two junctions and owns one or more lane specifications.
        ! The ``end_junction`` field stores the junction identifiers at the two
        ! ends of the road.

        integer :: id = 0
        ! Road identifier.
        integer :: end_junction(2) = 0
        ! Identifiers of the junctions connected by this road.
        type(lane_spec_t), allocatable :: lanes(:)
        ! Lane specifications belonging to this road.
    end type road_spec_t


    type, public :: route_row_t
        ! Route-probability row for one inbound junction leg.
        !
        ! The probability array has length ``n_out`` and stores the probabilities
        ! of choosing each outbound leg from a given inbound leg.

        real, allocatable :: prob(:)
        ! Route probabilities for one inbound leg, length ``n_out``.
    end type route_row_t


    type, public :: leg_spec_t
        ! Specification for one inbound or outbound junction leg.
        !
        ! A leg identifies a road-lane pair connected to a junction. The optional
        ! ``perim`` value gives the leg's port location on a cyclic perimeter for
        ! asymmetric junction conflict checks. A value of ``-1`` indicates that no
        ! perimeter index is supplied.

        integer :: road = 0
        ! Road index or identifier for this leg.
        integer :: lane = 0
        ! Lane index on the road.
        integer :: perim = -1
        ! Optional perimeter-port index. ``-1`` means no perimeter index.
    end type leg_spec_t


    type, public :: junc_spec_t
        ! Specification for one junction.
        !
        ! This type stores the flat inbound and outbound leg tables and the route
        ! probabilities associated with each inbound leg.

        integer :: id = 0
        ! Junction identifier.
        integer :: n_in = 0
        ! Number of inbound legs.
        integer :: n_out = 0
        ! Number of outbound legs.
        type(leg_spec_t), allocatable :: in_legs(:)
        ! Inbound leg specifications, length ``n_in``.
        type(leg_spec_t), allocatable :: out_legs(:)
        ! Outbound leg specifications, length ``n_out``.
        type(route_row_t), allocatable :: routes(:)
        ! Route-probability rows, one per inbound leg.
    end type junc_spec_t


    type, public :: network_spec_t
        ! Complete flat specification for a road network.
        !
        ! This type is the main input to ``build_network``. It contains all road
        ! and junction specifications required to construct the runtime
        ! ``road_network_t`` object.

        type(road_spec_t), allocatable :: roads(:)
        ! Road specifications.
        type(junc_spec_t), allocatable :: junctions(:)
        ! Junction specifications.
    end type network_spec_t

    public :: build_network, free_spec

contains

    subroutine build_network(spec, net)
        ! Build an allocated ``road_network_t`` from a ``network_spec_t``.
        !
        ! This routine allocates roads, lanes and junction arrays in ``net`` and
        ! copies the corresponding data from ``spec``. Lane storage is created by
        ! calling ``init_lane`` so that lane cells and old-state arrays are
        ! initialised consistently.
        !
        ! For each road, the builder:
        !
        ! * copies the road identifier and endpoint junctions;
        ! * allocates the lane array;
        ! * initialises each lane and copies its boundary parameters.
        !
        ! For each junction, the builder:
        !
        ! * copies inbound and outbound road-lane tables;
        ! * optionally allocates perimeter-port arrays if any port is specified;
        ! * allocates route-probability rows for all inbound legs.
        type(network_spec_t), intent(in) :: spec
        ! Flat input network specification.
        type(road_network_t), intent(out) :: net
        ! Allocated runtime road-network object ready for simulation.

        integer :: r
        ! Road index.
        integer :: l
        ! Lane index.
        integer :: j
        ! Junction index.
        integer :: k
        ! Leg or route index.
        integer :: has_perim_count
        ! Number of inbound or outbound legs with an explicit perimeter index.

        allocate(net%roads(size(spec%roads)))
        allocate(net%junctions(size(spec%junctions)))

        ! Build roads and lanes.
        do r = 1, size(spec%roads)
            net%roads(r)%id           = spec%roads(r)%id
            net%roads(r)%end_junction = spec%roads(r)%end_junction

            allocate(net%roads(r)%lane(size(spec%roads(r)%lanes)))

            do l = 1, size(spec%roads(r)%lanes)
                call init_lane(net%roads(r)%lane(l), &
                               spec%roads(r)%lanes(l)%length, &
                               spec%roads(r)%lanes(l)%flow_direction)

                net%roads(r)%lane(l)%alpha    = spec%roads(r)%lanes(l)%alpha
                net%roads(r)%lane(l)%beta     = spec%roads(r)%lanes(l)%beta
                net%roads(r)%lane(l)%open_in  = spec%roads(r)%lanes(l)%open_in
                net%roads(r)%lane(l)%open_out = spec%roads(r)%lanes(l)%open_out
            end do
        end do

        ! Build junctions.
        do j = 1, size(spec%junctions)
            net%junctions(j)%id    = spec%junctions(j)%id
            net%junctions(j)%n_in  = spec%junctions(j)%n_in
            net%junctions(j)%n_out = spec%junctions(j)%n_out

            allocate(net%junctions(j)%in_road(spec%junctions(j)%n_in))
            allocate(net%junctions(j)%in_lane(spec%junctions(j)%n_in))
            allocate(net%junctions(j)%out_road(spec%junctions(j)%n_out))
            allocate(net%junctions(j)%out_lane(spec%junctions(j)%n_out))

            do k = 1, spec%junctions(j)%n_in
                net%junctions(j)%in_road(k) = spec%junctions(j)%in_legs(k)%road
                net%junctions(j)%in_lane(k) = spec%junctions(j)%in_legs(k)%lane
            end do

            do k = 1, spec%junctions(j)%n_out
                net%junctions(j)%out_road(k) = spec%junctions(j)%out_legs(k)%road
                net%junctions(j)%out_lane(k) = spec%junctions(j)%out_legs(k)%lane
            end do

            ! Allocate perimeter-port arrays only when the specification uses them.
            has_perim_count = 0

            do k = 1, spec%junctions(j)%n_in
                if (spec%junctions(j)%in_legs(k)%perim >= 0) has_perim_count = has_perim_count + 1
            end do

            do k = 1, spec%junctions(j)%n_out
                if (spec%junctions(j)%out_legs(k)%perim >= 0) has_perim_count = has_perim_count + 1
            end do

            if (has_perim_count > 0) then
                allocate(net%junctions(j)%in_perim(spec%junctions(j)%n_in))
                allocate(net%junctions(j)%out_perim(spec%junctions(j)%n_out))

                do k = 1, spec%junctions(j)%n_in
                    net%junctions(j)%in_perim(k) = spec%junctions(j)%in_legs(k)%perim
                end do

                do k = 1, spec%junctions(j)%n_out
                    net%junctions(j)%out_perim(k) = spec%junctions(j)%out_legs(k)%perim
                end do
            end if

            ! Copy route-probability rows.
            allocate(net%junctions(j)%in_routes(spec%junctions(j)%n_in))

            do k = 1, spec%junctions(j)%n_in
                allocate(net%junctions(j)%in_routes(k)%prob(spec%junctions(j)%n_out))
                net%junctions(j)%in_routes(k)%prob = spec%junctions(j)%routes(k)%prob
            end do
        end do
    end subroutine build_network


    subroutine free_spec(spec)
        ! Deallocate all allocatable arrays owned by a network specification.
        !
        ! This routine releases nested arrays inside a ``network_spec_t``. It is
        ! intended for cleaning up specifications after they have been used to
        ! build a runtime network.
        !
        ! The runtime network object itself is not deallocated by this routine.
        type(network_spec_t), intent(inout) :: spec
        ! Network specification to deallocate in place.

        integer :: r
        ! Road index.
        integer :: j
        ! Junction index.
        integer :: k
        ! Route-row index.

        if (allocated(spec%roads)) then
            do r = 1, size(spec%roads)
                if (allocated(spec%roads(r)%lanes)) deallocate(spec%roads(r)%lanes)
            end do
            deallocate(spec%roads)
        end if

        if (allocated(spec%junctions)) then
            do j = 1, size(spec%junctions)
                if (allocated(spec%junctions(j)%in_legs))  deallocate(spec%junctions(j)%in_legs)
                if (allocated(spec%junctions(j)%out_legs)) deallocate(spec%junctions(j)%out_legs)

                if (allocated(spec%junctions(j)%routes)) then
                    do k = 1, size(spec%junctions(j)%routes)
                        if (allocated(spec%junctions(j)%routes(k)%prob)) &
                            deallocate(spec%junctions(j)%routes(k)%prob)
                    end do
                    deallocate(spec%junctions(j)%routes)
                end if
            end do
            deallocate(spec%junctions)
        end if
    end subroutine free_spec

end module network_builder_mod