program test_t_junction
    !--------------------------------------------------------------------
    ! Phase 4 validation suite for the T-junction topology.
    !
    ! Layout (looking down):
    !
    !    [open alpha/beta]   west arm   [junction]   east arm   [open beta]
    !    ←  road1 lane2 out  ←←←←←←←←←←←←←←←←←←←←←←←  road3 lane1 out  →
    !    →  road1 lane1 in   →→→→→→→→→→[jct]
    !                                      ↑
    !                             road2 lane1 in (stem)
    !                                      ↑
    !                                 [open alpha]
    !
    ! Roads: 1=west (bidirectional), 2=south-stem (inbound only),
    !        3=east (outbound only).
    ! Junction: n_in=2 [west_in=leg1, stem_in=leg2],
    !           n_out=2 [west_out=outbound1, east_out=outbound2].
    ! Perimeter ports (0-indexed CW): 0=east_out, 1=stem_in,
    !                                 2=west_in,  3=west_out.
    !
    ! Helper access:
    !   west-in  holding cell = road1 lane1 site L  (flow +1, site L at junction)
    !   stem-in  holding cell = road2 lane1 site L
    !   west-out destination  = road1 lane2 site 1  (flow -1, site 1 at junction)
    !   east-out destination  = road3 lane1 site 1  (flow +1, site 1 at junction)
    !--------------------------------------------------------------------
    use vehicle_mod
    use road_network_mod
    use network_init_mod
    use junction_mod
    use network_simulation_mod
    implicit none

    integer, parameter :: L       = 10
    integer, parameter :: N_STEPS = 2000
    real,    parameter :: ALPHA   = 0.4
    real,    parameter :: BETA    = 0.5

    call seed_rng(77)

    call test_t_init()
    call test_t_merge_yield()
    call test_t_diverge()
    call test_t_chord_no_conflict()
    call test_t_chord_conflict()
    call test_t_mass_conservation()

    print *, 'test_t_junction: ALL OK'

contains

    !-----------------------------------------------------------------
    subroutine seed_rng(s)
        integer, intent(in) :: s
        integer :: n
        integer, allocatable :: seed(:)
        call random_seed(size = n)
        allocate(seed(n))
        seed = s
        call random_seed(put = seed)
        deallocate(seed)
    end subroutine seed_rng

    subroutine fresh(net)
        type(road_network_t), intent(out) :: net
        ! west-in never U-turns; stem-in 50/50.  All boundaries closed.
        call init_t_junction(net, L, 0.0, 0.0, [0.0, 1.0], [0.5, 0.5])
    end subroutine fresh

    subroutine plant_west_in(net)
        type(road_network_t), intent(inout) :: net
        net%roads(1)%lane(1)%cells(L) = V_OCCUPIED
    end subroutine plant_west_in

    subroutine plant_stem_in(net)
        type(road_network_t), intent(inout) :: net
        net%roads(2)%lane(1)%cells(L) = V_OCCUPIED
    end subroutine plant_stem_in

    subroutine force_route(net, jid, leg, out_idx)
        type(road_network_t), intent(inout) :: net
        integer, intent(in) :: jid, leg, out_idx
        net%junctions(jid)%in_routes(leg)%prob         = 0.0
        net%junctions(jid)%in_routes(leg)%prob(out_idx) = 1.0
    end subroutine force_route

    function west_in_holding(net) result(v)
        type(road_network_t), intent(in) :: net
        integer :: v
        v = net%roads(1)%lane(1)%cells(L)
    end function west_in_holding

    function stem_in_holding(net) result(v)
        type(road_network_t), intent(in) :: net
        integer :: v
        v = net%roads(2)%lane(1)%cells(L)
    end function stem_in_holding

    function west_out_dest(net) result(v)
        ! road1 lane2 site 1 (outbound, flow -1, site 1 at junction end_2)
        type(road_network_t), intent(in) :: net
        integer :: v
        v = net%roads(1)%lane(2)%cells(1)
    end function west_out_dest

    function east_out_dest(net) result(v)
        ! road3 lane1 site 1 (outbound, flow +1, site 1 at junction end_1)
        type(road_network_t), intent(in) :: net
        integer :: v
        v = net%roads(3)%lane(1)%cells(1)
    end function east_out_dest

    subroutine assert(cond, msg)
        logical, intent(in) :: cond
        character(len=*), intent(in) :: msg
        if (.not. cond) then
            print *, 'ASSERT FAILED: ', trim(msg)
            stop 1
        end if
    end subroutine assert

    !-----------------------------------------------------------------
    ! 1. Structural initialisation checks.
    !-----------------------------------------------------------------
    subroutine test_t_init()
        type(road_network_t) :: net
        call fresh(net)

        call assert(size(net%roads)     == 3, 'T: wrong road count')
        call assert(size(net%junctions) == 1, 'T: wrong junction count')

        call assert(net%junctions(1)%n_in  == 2, 'T: n_in')
        call assert(net%junctions(1)%n_out == 2, 'T: n_out')
        call assert(all(net%junctions(1)%in_road  == [1, 2]), 'T: in_road')
        call assert(all(net%junctions(1)%out_road == [1, 3]), 'T: out_road')

        ! in_lane: road1 inbound at end_2 is lane1 (flow+1); road2 lane1.
        call assert(net%junctions(1)%in_lane(1) == 1, 'T: in_lane west')
        call assert(net%junctions(1)%in_lane(2) == 1, 'T: in_lane stem')
        ! out_lane: road1 outbound at end_2 is lane2 (flow-1); road3 lane1.
        call assert(net%junctions(1)%out_lane(1) == 2, 'T: out_lane west')
        call assert(net%junctions(1)%out_lane(2) == 1, 'T: out_lane east')

        call assert(allocated(net%junctions(1)%in_perim),  'T: in_perim allocated')
        call assert(allocated(net%junctions(1)%out_perim), 'T: out_perim allocated')
        call assert(size(net%junctions(1)%in_perim)  == 2, 'T: in_perim size')
        call assert(size(net%junctions(1)%out_perim) == 2, 'T: out_perim size')
        ! Perimeter: 0=east_out, 1=stem_in, 2=west_in, 3=west_out.
        call assert(net%junctions(1)%in_perim(1)  == 2, 'T: in_perim west_in=2')
        call assert(net%junctions(1)%in_perim(2)  == 1, 'T: in_perim stem_in=1')
        call assert(net%junctions(1)%out_perim(1) == 3, 'T: out_perim west_out=3')
        call assert(net%junctions(1)%out_perim(2) == 0, 'T: out_perim east_out=0')

        call assert(allocated(net%junctions(1)%in_routes),         'T: in_routes allocated')
        call assert(size(net%junctions(1)%in_routes) == 2,         'T: in_routes size')
        call assert(size(net%junctions(1)%in_routes(1)%prob) == 2, 'T: prob(1) size')
        call assert(size(net%junctions(1)%in_routes(2)%prob) == 2, 'T: prob(2) size')
        call assert(abs(sum(net%junctions(1)%in_routes(1)%prob) - 1.0) < 1e-6, 'T: prob(1) sums to 1')
        call assert(abs(sum(net%junctions(1)%in_routes(2)%prob) - 1.0) < 1e-6, 'T: prob(2) sums to 1')

        call free_network(net)
        print *, '  test_t_init OK'
    end subroutine test_t_init

    !-----------------------------------------------------------------
    ! 2. Merge conflict: both vehicles want east_out -> same destination
    !    -> mutual R1 yield -> exactly one advances.
    !-----------------------------------------------------------------
    subroutine test_t_merge_yield()
        type(road_network_t) :: net
        integer :: n_moved
        call fresh(net)
        call plant_west_in(net)
        call plant_stem_in(net)
        call force_route(net, 1, 1, 2)   ! west-in  -> east_out (out_idx=2)
        call force_route(net, 1, 2, 2)   ! stem-in  -> east_out (out_idx=2)

        call snapshot_network(net)
        call evaluate_junctions(net)

        n_moved = 0
        if (west_in_holding(net) == V_EMPTY) n_moved = n_moved + 1
        if (stem_in_holding(net) == V_EMPTY) n_moved = n_moved + 1
        call assert(n_moved == 1, 'T merge: exactly one should advance')

        call free_network(net)
        print *, '  test_t_merge_yield OK'
    end subroutine test_t_merge_yield

    !-----------------------------------------------------------------
    ! 3. Diverge: single stem vehicle routes to west_out.
    !    Holding clears; west_out destination gets V_OCCUPIED.
    !-----------------------------------------------------------------
    subroutine test_t_diverge()
        type(road_network_t) :: net
        call fresh(net)
        call plant_stem_in(net)
        call force_route(net, 1, 2, 1)   ! stem-in -> west_out (out_idx=1)

        call snapshot_network(net)
        call evaluate_junctions(net)

        call assert(stem_in_holding(net) == V_EMPTY,    'T diverge: holding not cleared')
        call assert(west_out_dest(net)   == V_OCCUPIED, 'T diverge: west_out not filled')
        call assert(east_out_dest(net)   == V_EMPTY,    'T diverge: east_out should be empty')

        call free_network(net)
        print *, '  test_t_diverge OK'
    end subroutine test_t_diverge

    !-----------------------------------------------------------------
    ! 4. Chord non-conflict: west-in->west_out and stem-in->east_out.
    !    Paths 2->3 and 1->0 on the perimeter do NOT cross.
    !    Both vehicles advance in the same step.
    !-----------------------------------------------------------------
    subroutine test_t_chord_no_conflict()
        type(road_network_t) :: net
        call fresh(net)
        call plant_west_in(net)
        call plant_stem_in(net)
        call force_route(net, 1, 1, 1)   ! west-in -> west_out (out_idx=1)
        call force_route(net, 1, 2, 2)   ! stem-in -> east_out (out_idx=2)

        call snapshot_network(net)
        call evaluate_junctions(net)

        call assert(west_in_holding(net) == V_EMPTY,    'T chord: west-in holding not cleared')
        call assert(stem_in_holding(net) == V_EMPTY,    'T chord: stem-in holding not cleared')
        call assert(west_out_dest(net)   == V_OCCUPIED, 'T chord: west_out not filled')
        call assert(east_out_dest(net)   == V_OCCUPIED, 'T chord: east_out not filled')

        call free_network(net)
        print *, '  test_t_chord_no_conflict OK'
    end subroutine test_t_chord_no_conflict

    !-----------------------------------------------------------------
    ! 5. Chord conflict: west-in->east_out (perimeter path 2->0) and
    !    stem-in->west_out (perimeter path 1->3).  These chords cross
    !    (port 3 is in the open CW arc from 2 to 0, but port 1 is not),
    !    so R1 fires for both legs -> mutual yield -> exactly one advances
    !    and the other waits in its holding cell.
    !-----------------------------------------------------------------
    subroutine test_t_chord_conflict()
        type(road_network_t) :: net
        integer :: n_moved
        call fresh(net)
        call plant_west_in(net)
        call plant_stem_in(net)
        call force_route(net, 1, 1, 2)   ! west-in -> east_out (perimeter: 2->0)
        call force_route(net, 1, 2, 1)   ! stem-in -> west_out (perimeter: 1->3)

        call snapshot_network(net)
        call evaluate_junctions(net)

        n_moved = 0
        if (west_in_holding(net) == V_EMPTY) n_moved = n_moved + 1
        if (stem_in_holding(net) == V_EMPTY) n_moved = n_moved + 1
        call assert(n_moved == 1, 'T chord conflict: exactly one should advance')
        ! Blocked vehicle must remain waiting in its holding cell.
        call assert(west_in_holding(net) == V_OCCUPIED .or. &
                    stem_in_holding(net) == V_OCCUPIED, &
                    'T chord conflict: blocked vehicle must wait in holding cell')

        call free_network(net)
        print *, '  test_t_chord_conflict OK'
    end subroutine test_t_chord_conflict

    !-----------------------------------------------------------------
    ! 6. Mass conservation soak: 2000 steps.
    !    At every step: count(t+1) == count(t) + entries - exits.
    !-----------------------------------------------------------------
    subroutine test_t_mass_conservation()
        type(road_network_t) :: net
        integer :: step, n_before, n_after, entries, exits
        integer :: seed_size
        integer, allocatable :: s(:)

        call random_seed(size = seed_size)
        allocate(s(seed_size))
        s = 12345
        call random_seed(put = s)
        deallocate(s)

        call init_t_junction(net, L, ALPHA, BETA, [0.0, 1.0], [0.5, 0.5])

        do step = 1, N_STEPS
            n_before = count_occupied_network(net)
            call network_step(net)

            entries = 0
            exits   = 0

            ! Entries: new vehicle at open-in site 1.
            if (net%roads(1)%lane(1)%cells(1) /= V_EMPTY .and. &
                net%roads(1)%lane(1)%old(1)   == V_EMPTY) entries = entries + 1
            if (net%roads(2)%lane(1)%cells(1) /= V_EMPTY .and. &
                net%roads(2)%lane(1)%old(1)   == V_EMPTY) entries = entries + 1

            ! Exits: vehicle cleared from open-out site L.
            if (net%roads(1)%lane(2)%cells(L) == V_EMPTY .and. &
                net%roads(1)%lane(2)%old(L)   /= V_EMPTY) exits = exits + 1
            if (net%roads(3)%lane(1)%cells(L) == V_EMPTY .and. &
                net%roads(3)%lane(1)%old(L)   /= V_EMPTY) exits = exits + 1

            n_after = count_occupied_network(net)
            if (n_after /= n_before + entries - exits) then
                print '(a,i6)', 'FAIL: T mass conservation violated at step ', step
                stop 1
            end if
        end do

        call free_network(net)
        print '(a,i6,a)', 'PASS: T mass conservation OK (all ', N_STEPS, ' steps)'
        print *, '  test_t_mass_conservation OK'
    end subroutine test_t_mass_conservation

end program test_t_junction
