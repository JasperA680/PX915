module road_network_mod
    !--------------------------------------------------------------------
    ! Network topology and per-lane storage.
    !
    !   lane_t           a strictly unidirectional 1D TASEP lane
    !   road_t           a container of laterally-adjacent lanes
    !   leg_route_t      per-inbound-leg outbound probability distribution
    !   junction_t       a node connecting roads in clockwise order
    !   road_network_t   the master container of roads and junctions
    !
    ! Lane orientation:
    !   flow_direction = +1 : vehicles flow end_1 -> end_2 (site 1 at end_1,
    !                         site L at end_2)
    !   flow_direction = -1 : vehicles flow end_2 -> end_1 (site 1 at end_2,
    !                         site L at end_1)
    !
    ! In both cases site 1 is the entry end, site L is the holding cell
    ! toward whichever junction the lane flows into.
    !--------------------------------------------------------------------
    use vehicle_mod
    implicit none
    private

    type, public :: lane_t
        integer :: length         = 0
        integer :: flow_direction = 0     ! +1: end_1->end_2,  -1: end_2->end_1
        integer, allocatable :: cells(:)
        integer, allocatable :: old(:)    ! parallel-update snapshot
        real    :: alpha   = 0.0
        real    :: beta    = 0.0
        logical :: open_in  = .false.     ! site 1 is open inflow
        logical :: open_out = .false.     ! site L is open outflow
    end type lane_t

    type, public :: road_t
        integer :: id = 0
        type(lane_t), allocatable :: lane(:)   ! variable-count lane bundle
        integer :: end_junction(2) = 0         ! 0 = open boundary
    end type road_t

    type, public :: leg_route_t
        real, allocatable :: prob(:)           ! length n_out, sums to 1
    end type leg_route_t

    type, public :: junction_t
        integer :: id   = 0
        integer :: n_in = 0                    ! number of inbound lanes (clockwise-ordered)
        integer :: n_out = 0                   ! number of outbound lanes (clockwise-ordered)
        integer, allocatable :: in_road(:)     ! road id for inbound leg k   (length n_in)
        integer, allocatable :: in_lane(:)     ! lane index within road      (length n_in)
        integer, allocatable :: out_road(:)    ! road id for outbound leg k  (length n_out)
        integer, allocatable :: out_lane(:)    ! lane index within road      (length n_out)
        type(leg_route_t), allocatable :: in_routes(:)  ! length n_in
        ! Perimeter port positions (0-indexed, CW) for chord-crossing conflict predicate.
        ! Only allocated for asymmetric junctions (n_in /= 4 or n_out /= 4).
        integer, allocatable :: in_perim(:)    ! length n_in
        integer, allocatable :: out_perim(:)   ! length n_out
    end type junction_t

    type, public :: road_network_t
        type(road_t),     allocatable :: roads(:)
        type(junction_t), allocatable :: junctions(:)
    end type road_network_t

    public :: inbound_lane_at_end, outbound_lane_at_end, &
              snapshot_network, count_occupied_network

contains

    !--------------------------------------------------------------------
    ! Return the index of the lane that is inbound toward a junction at
    ! road end `junc_end` (1 or 2).  A lane is inbound to end e when its
    ! site L (holding cell) is at that end:
    !   flow_direction = -1  =>  site L at end_1  (inbound to end 1)
    !   flow_direction = +1  =>  site L at end_2  (inbound to end 2)
    ! Formula: inbound to end e when flow_direction == 2*e - 3.
    !--------------------------------------------------------------------
    pure function inbound_lane_at_end(road, junc_end) result(li)
        type(road_t), intent(in) :: road
        integer,      intent(in) :: junc_end
        integer :: li, k
        li = 0
        do k = 1, size(road%lane)
            if (road%lane(k)%flow_direction == 2*junc_end - 3) then
                li = k
                return
            end if
        end do
    end function inbound_lane_at_end

    !--------------------------------------------------------------------
    ! Return the index of the lane that is outbound from a junction at
    ! road end `junc_end`.  Outbound from end e when site 1 is at that end:
    !   flow_direction = +1  =>  site 1 at end_1  (outbound from end 1)
    !   flow_direction = -1  =>  site 1 at end_2  (outbound from end 2)
    ! Formula: outbound from end e when flow_direction == 3 - 2*e.
    !--------------------------------------------------------------------
    pure function outbound_lane_at_end(road, junc_end) result(li)
        type(road_t), intent(in) :: road
        integer,      intent(in) :: junc_end
        integer :: li, k
        li = 0
        do k = 1, size(road%lane)
            if (road%lane(k)%flow_direction == 3 - 2*junc_end) then
                li = k
                return
            end if
        end do
    end function outbound_lane_at_end

    subroutine snapshot_network(net)
        type(road_network_t), intent(inout) :: net
        integer :: r, l
        do r = 1, size(net%roads)
            do l = 1, size(net%roads(r)%lane)
                net%roads(r)%lane(l)%old = net%roads(r)%lane(l)%cells
            end do
        end do
    end subroutine snapshot_network

    function count_occupied_network(net) result(n)
        type(road_network_t), intent(in) :: net
        integer :: n, r, l, k
        n = 0
        do r = 1, size(net%roads)
            do l = 1, size(net%roads(r)%lane)
                do k = 1, net%roads(r)%lane(l)%length
                    if (is_occupied(net%roads(r)%lane(l)%cells(k))) n = n + 1
                end do
            end do
        end do
    end function count_occupied_network

end module road_network_mod
