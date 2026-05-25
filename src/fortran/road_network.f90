module road_network_mod
    ! Network topology and per-lane storage for road-network simulations.
    !
    ! This module defines the core data structures used by the road-network
    ! cellular automaton models. A network is represented as a collection of
    ! roads and junctions. Each road contains one or more unidirectional lanes,
    ! and each junction connects inbound lanes to outbound lanes through route
    ! probabilities.
    !
    ! The main types are:
    !
    ! * ``cell``: one lattice cell, storing occupancy and vehicle velocity;
    ! * ``lane_t``: one strictly unidirectional one-dimensional lane;
    ! * ``road_t``: a bundle of laterally adjacent lanes;
    ! * ``leg_route_t``: route probabilities for one inbound junction leg;
    ! * ``junction_t``: a node connecting roads through inbound and outbound legs;
    ! * ``road_network_t``: the full network container.
    !
    ! Lane orientation is controlled by ``flow_direction``:
    !
    ! * ``flow_direction = +1`` means vehicles flow from ``end_junction(1)`` to
    !   ``end_junction(2)``;
    ! * ``flow_direction = -1`` means vehicles flow from ``end_junction(2)`` to
    !   ``end_junction(1)``.
    !
    ! In both cases, site ``1`` is the lane entry end and site ``L`` is the
    ! holding cell at the downstream end of the lane. The holding cell is the
    ! cell that interacts with a junction or open outflow boundary.
    !
    ! For a road end ``e`` equal to 1 or 2, a lane is inbound to that end when
    ! its holding cell is located there:
    !
    ! .. math::
    !
    !    \mathrm{flow\_direction} = 2e - 3.
    !
    ! A lane is outbound from road end ``e`` when its entry cell is located
    ! there:
    !
    ! .. math::
    !
    !    \mathrm{flow\_direction} = 3 - 2e.

    implicit none
    private

    type, public :: cell
        ! One lattice cell in a cellular automaton lane.
        !
        ! A cell stores whether it currently contains a vehicle and, for models
        ! with variable speed such as Nagel-Schreckenberg, the integer velocity
        ! of that vehicle.

        integer :: velocity = 0
        ! Vehicle velocity stored in this cell.
        logical :: has_car = .false.
        ! Whether this cell is occupied by a vehicle.
    end type cell


    type, public :: lane_t
        ! A strictly unidirectional one-dimensional traffic lane.
        !
        ! The lane stores its live cell state and a snapshot of the previous state.
        ! The snapshot is used for parallel updates, so that all movement decisions
        ! can be made from the same old configuration before the new configuration
        ! is committed.
        !
        ! The lane may also have open boundary conditions. Site ``1`` is the
        ! inflow end and site ``length`` is the downstream holding or outflow cell.

        integer :: length = 0
        ! Number of cells in the lane.
        integer :: flow_direction = 0
        ! Flow orientation. ``+1`` means road end 1 to end 2; ``-1`` means road end 2 to end 1.
        type(cell), allocatable :: cells(:)
        ! Live cell state.
        type(cell), allocatable :: old(:)
        ! Snapshot of the previous cell state used for parallel updates.
        real :: alpha = 0.0
        ! Open-boundary entry probability at site 1.
        real :: beta = 0.0
        ! Open-boundary exit probability at site ``length``.
        logical :: open_in = .false.
        ! Whether site 1 is an open inflow boundary.
        logical :: open_out = .false.
        ! Whether site ``length`` is an open outflow boundary.
    end type lane_t


    type, public :: road_t
        ! A road containing one or more laterally adjacent lanes.
        !
        ! Roads connect two junction endpoints. An endpoint value of zero denotes
        ! an open boundary rather than a junction.

        integer :: id = 0
        ! Road identifier.
        type(lane_t), allocatable :: lane(:)
        ! Lanes belonging to this road.
        integer :: end_junction(2) = 0
        ! Junction identifiers at the two road ends. ``0`` denotes an open boundary.
    end type road_t


    type, public :: leg_route_t
        ! Outbound route probabilities for one inbound junction leg.
        !
        ! The probability array has length ``n_out`` for the corresponding
        ! junction. It should satisfy:
        !
        ! .. math::
        !
        !    \sum_{m=1}^{n_{\mathrm{out}}} p_m = 1.

        real, allocatable :: prob(:)
        ! Probability of selecting each outbound leg.
    end type leg_route_t


    type, public :: junction_t
        ! A junction connecting inbound and outbound road-lane pairs.
        !
        ! A junction stores clockwise-ordered inbound and outbound legs. Each leg
        ! points to a road and lane. The ``in_routes`` array stores the outbound
        ! probability distribution for each inbound leg.
        !
        ! For asymmetric junctions, optional perimeter-port arrays can be allocated.
        ! These store zero-indexed clockwise port positions for chord-crossing
        ! conflict checks in the junction evaluator.

        integer :: id = 0
        ! Junction identifier.
        integer :: n_in = 0
        ! Number of inbound legs, ordered clockwise.
        integer :: n_out = 0
        ! Number of outbound legs, ordered clockwise.
        integer, allocatable :: in_road(:)
        ! Road identifier for each inbound leg, length ``n_in``.
        integer, allocatable :: in_lane(:)
        ! Lane index within the inbound road, length ``n_in``.
        integer, allocatable :: out_road(:)
        ! Road identifier for each outbound leg, length ``n_out``.
        integer, allocatable :: out_lane(:)
        ! Lane index within the outbound road, length ``n_out``.
        type(leg_route_t), allocatable :: in_routes(:)
        ! Route-probability row for each inbound leg, length ``n_in``.
        integer, allocatable :: in_perim(:)
        ! Optional zero-indexed clockwise perimeter port for each inbound leg.
        integer, allocatable :: out_perim(:)
        ! Optional zero-indexed clockwise perimeter port for each outbound leg.
    end type junction_t


    type, public :: road_network_t
        ! Master container for a road-network simulation.
        !
        ! The network contains all roads and all junctions. Simulation routines
        ! update the lane states and use the junction information to move vehicles
        ! between roads.

        type(road_t), allocatable :: roads(:)
        ! Roads in the network.
        type(junction_t), allocatable :: junctions(:)
        ! Junctions in the network.
    end type road_network_t

    public :: inbound_lane_at_end, outbound_lane_at_end
    public :: snapshot_network, count_occupied_network

contains

    pure function inbound_lane_at_end(road, junc_end) result(li)
        ! Return the lane index that is inbound to a given road end.
        !
        ! A lane is inbound to road end ``junc_end`` when its downstream holding
        ! cell, site ``L``, is located at that end.
        !
        ! The convention is:
        !
        ! * ``flow_direction = -1`` means site ``L`` is at road end 1;
        ! * ``flow_direction = +1`` means site ``L`` is at road end 2.
        !
        ! Therefore the matching condition is:
        !
        ! .. math::
        !
        !    \mathrm{flow\_direction} = 2e - 3,
        !
        ! where ``e`` is ``junc_end``. The function returns zero if no matching
        ! lane is found.
        type(road_t), intent(in) :: road
        ! Road to search.
        integer, intent(in) :: junc_end
        ! Road end to test, either 1 or 2.
        integer :: li
        ! Matching inbound lane index, or zero if no matching lane exists.
        integer :: k
        ! Lane loop index.

        li = 0

        do k = 1, size(road%lane)
            if (road%lane(k)%flow_direction == 2*junc_end - 3) then
                li = k
                return
            end if
        end do
    end function inbound_lane_at_end


    pure function outbound_lane_at_end(road, junc_end) result(li)
        ! Return the lane index that is outbound from a given road end.
        !
        ! A lane is outbound from road end ``junc_end`` when its entry cell,
        ! site ``1``, is located at that end.
        !
        ! The convention is:
        !
        ! * ``flow_direction = +1`` means site ``1`` is at road end 1;
        ! * ``flow_direction = -1`` means site ``1`` is at road end 2.
        !
        ! Therefore the matching condition is:
        !
        ! .. math::
        !
        !    \mathrm{flow\_direction} = 3 - 2e,
        !
        ! where ``e`` is ``junc_end``. The function returns zero if no matching
        ! lane is found.
        type(road_t), intent(in) :: road
        ! Road to search.
        integer, intent(in) :: junc_end
        ! Road end to test, either 1 or 2.
        integer :: li
        ! Matching outbound lane index, or zero if no matching lane exists.
        integer :: k
        ! Lane loop index.

        li = 0

        do k = 1, size(road%lane)
            if (road%lane(k)%flow_direction == 3 - 2*junc_end) then
                li = k
                return
            end if
        end do
    end function outbound_lane_at_end


    subroutine snapshot_network(net)
        ! Copy live cell states into old-state snapshots for the whole network.
        !
        ! Many cellular automaton updates require strict parallel semantics. This
        ! means decisions should be made from the old state, while accepted moves
        ! are written to the live ``cells`` arrays. This routine prepares the old
        ! state by copying:
        !
        ! .. math::
        !
        !    \mathrm{old} \leftarrow \mathrm{cells}
        !
        ! for every lane in every road.
        type(road_network_t), intent(inout) :: net
        ! Network whose lane snapshots are updated.

        integer :: r
        ! Road index.
        integer :: l
        ! Lane index.

        do r = 1, size(net%roads)
            do l = 1, size(net%roads(r)%lane)
                net%roads(r)%lane(l)%old = net%roads(r)%lane(l)%cells
            end do
        end do
    end subroutine snapshot_network


    function count_occupied_network(net) result(n)
        ! Count the total number of occupied cells in the whole network.
        !
        ! The returned value is:
        !
        ! .. math::
        !
        !    N =
        !    \sum_{\mathrm{roads}}
        !    \sum_{\mathrm{lanes}}
        !    \sum_{\mathrm{cells}}
        !    \tau_{r,l,k},
        !
        ! where ``tau`` is one for an occupied cell and zero for an empty cell.
        type(road_network_t), intent(in) :: net
        ! Network to count.
        integer :: n
        ! Total number of occupied cells.
        integer :: r
        ! Road index.
        integer :: l
        ! Lane index.
        integer :: k
        ! Cell index.

        n = 0

        do r = 1, size(net%roads)
            do l = 1, size(net%roads(r)%lane)
                do k = 1, net%roads(r)%lane(l)%length
                    if (net%roads(r)%lane(l)%cells(k)%has_car) then
                        n = n + 1
                    end if
                end do
            end do
        end do
    end function count_occupied_network

end module road_network_mod