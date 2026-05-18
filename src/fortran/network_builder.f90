module network_builder_mod
    !--------------------------------------------------------------------
    ! Generic network builder.  Takes a flat `network_spec_t` (mirroring
    ! the JSON config schema) and constructs a `road_network_t` ready
    ! for `network_step`.
    !
    ! The spec has three top-level arrays:
    !   road_spec_t  -- one entry per road; carries its lane list
    !   junc_spec_t  -- one entry per junction; carries its flat leg
    !                   tables and a route-prob row per inbound leg
    !
    ! `build_network` reuses init_lane (made public in network_init_mod).
    !--------------------------------------------------------------------
    use road_network_mod
    use network_init_mod, only: init_lane
    implicit none
    private

    type, public :: lane_spec_t
        integer :: length         = 0
        integer :: flow_direction = 0
        real    :: alpha          = 0.0
        real    :: beta           = 0.0
        logical :: open_in        = .false.
        logical :: open_out       = .false.
    end type lane_spec_t

    type, public :: road_spec_t
        integer :: id              = 0
        integer :: end_junction(2) = 0
        type(lane_spec_t), allocatable :: lanes(:)
    end type road_spec_t

    type, public :: route_row_t
        real, allocatable :: prob(:)        ! length n_out
    end type route_row_t

    type, public :: leg_spec_t
        integer :: road  = 0
        integer :: lane  = 0
        integer :: perim = -1              ! -1 = symmetric (no perim)
    end type leg_spec_t

    type, public :: junc_spec_t
        integer :: id    = 0
        integer :: n_in  = 0
        integer :: n_out = 0
        type(leg_spec_t),  allocatable :: in_legs(:)
        type(leg_spec_t),  allocatable :: out_legs(:)
        type(route_row_t), allocatable :: routes(:)    ! length n_in
    end type junc_spec_t

    type, public :: network_spec_t
        type(road_spec_t), allocatable :: roads(:)
        type(junc_spec_t), allocatable :: junctions(:)
    end type network_spec_t

    public :: build_network, free_spec

contains

    subroutine build_network(spec, net)
        type(network_spec_t), intent(in)  :: spec
        type(road_network_t), intent(out) :: net
        integer :: r, l, j, k, has_perim_count

        allocate(net%roads(size(spec%roads)))
        allocate(net%junctions(size(spec%junctions)))

        ! Roads + lanes.
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

        ! Junctions.
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

            ! Perimeter ports: allocate only if any leg has perim >= 0.
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

            ! Routes.
            allocate(net%junctions(j)%in_routes(spec%junctions(j)%n_in))
            do k = 1, spec%junctions(j)%n_in
                allocate(net%junctions(j)%in_routes(k)%prob(spec%junctions(j)%n_out))
                net%junctions(j)%in_routes(k)%prob = spec%junctions(j)%routes(k)%prob
            end do
        end do
    end subroutine build_network

    subroutine free_spec(spec)
        type(network_spec_t), intent(inout) :: spec
        integer :: r, j, k
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
