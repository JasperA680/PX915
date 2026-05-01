module road_network_mod
    !--------------------------------------------------------------------
    ! Network topology and per-lane storage.
    !
    !   lane_t           a strictly unidirectional 1D TASEP lane
    !   road_t           a bidirectional container of two opposite lanes
    !   junction_t       a node connecting roads in clockwise order
    !   road_network_t   the master container of roads and junctions
    !
    ! Lane orientation convention (used everywhere downstream):
    !   lane(1) carries flow from road end_1 to end_2; vehicles travel
    !           in the +index direction so site 1 is at end_1, site L
    !           is at end_2.
    !   lane(2) carries flow from end_2 to end_1; site 1 at end_2,
    !           site L at end_1.
    !
    ! With this convention the inbound lane (toward a junction at end e)
    ! always has its holding cell at site L, and the outbound lane (away
    ! from that junction) always has its first cell at site 1.
    !--------------------------------------------------------------------
    use vehicle_mod
    implicit none
    private

    type, public :: lane_t
        integer :: length = 0
        integer, allocatable :: cells(:)
        integer, allocatable :: old(:)         ! parallel-update snapshot
        real    :: alpha   = 0.0
        real    :: beta    = 0.0
        real    :: p_left  = 0.0
        real    :: p_right = 0.0
        logical :: open_in  = .false.          ! site 1 is open inflow
        logical :: open_out = .false.          ! site L is open outflow
    end type lane_t

    type, public :: road_t
        integer      :: id = 0
        type(lane_t) :: lane(2)
        integer      :: end_junction(2) = 0    ! 0 = open boundary
    end type road_t

    type, public :: junction_t
        integer :: id     = 0
        integer :: n_legs = 0                  ! 3 (T) or 4 (cross)
        integer, allocatable :: connected_road_ids(:)   ! strict clockwise order
        integer, allocatable :: end_at(:)               ! 1 or 2: which end of each road we touch
    end type junction_t

    type, public :: road_network_t
        type(road_t),     allocatable :: roads(:)
        type(junction_t), allocatable :: junctions(:)
    end type road_network_t

    public :: inbound_lane_idx, outbound_lane_idx, snapshot_network, &
              count_occupied_network

contains

    pure function inbound_lane_idx(junction_end) result(li)
        integer, intent(in) :: junction_end       ! 1 or 2
        integer :: li
        li = 3 - junction_end                     ! end 1 -> lane 2, end 2 -> lane 1
    end function inbound_lane_idx

    pure function outbound_lane_idx(junction_end) result(li)
        integer, intent(in) :: junction_end
        integer :: li
        li = junction_end                         ! end 1 -> lane 1, end 2 -> lane 2
    end function outbound_lane_idx

    subroutine snapshot_network(net)
        type(road_network_t), intent(inout) :: net
        integer :: r, l
        do r = 1, size(net%roads)
            do l = 1, 2
                net%roads(r)%lane(l)%old = net%roads(r)%lane(l)%cells
            end do
        end do
    end subroutine snapshot_network

    function count_occupied_network(net) result(n)
        type(road_network_t), intent(in) :: net
        integer :: n, r, l, k
        n = 0
        do r = 1, size(net%roads)
            do l = 1, 2
                do k = 1, net%roads(r)%lane(l)%length
                    if (is_occupied(net%roads(r)%lane(l)%cells(k))) n = n + 1
                end do
            end do
        end do
    end function count_occupied_network

end module road_network_mod
