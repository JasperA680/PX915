program test_ns_grid
    !--------------------------------------------------------------------
    ! Nagel-Schreckenberg validation on the 2D crossroad grid.
    ! Prints per-lane velocity snapshots for visual inspection.
    !--------------------------------------------------------------------
    use vehicle_mod
    use road_network_mod
    use network_init_mod
    use network_simulation_mod
    implicit none

    integer, parameter :: L = 12
    integer, parameter :: N_STEPS = 3
    integer, parameter :: V_MAX = 5

    type(road_network_t) :: net
    integer :: step

    call seed_rng(20260506)
    call init_crossroad(net, L, 0.0, 0.0, 0.0, 0.0)
    call close_boundaries(net)
    call seed_layout(net)

    print *, 'NS 2D grid (crossroad) test'
    call print_state(net, 0)

    do step = 1, N_STEPS
        call network_step(net)
        call check_step(net, step)
        call print_state(net, step)
    end do

    call free_network(net)
    print *, 'test_ns_grid: OK'

contains

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

    subroutine close_boundaries(net)
        type(road_network_t), intent(inout) :: net
        integer :: r, ln
        do r = 1, size(net%roads)
            do ln = 1, size(net%roads(r)%lane)
                net%roads(r)%lane(ln)%open_in  = .false.
                net%roads(r)%lane(ln)%open_out = .false.
                net%roads(r)%lane(ln)%alpha    = 0.0
                net%roads(r)%lane(ln)%beta     = 0.0
            end do
        end do
    end subroutine close_boundaries

    subroutine seed_layout(net)
        type(road_network_t), intent(inout) :: net
        call place_car(net%roads(1)%lane(1), 2, 1)
        call place_car(net%roads(1)%lane(1), 6, 0)
        call place_car(net%roads(2)%lane(2), 3, 2)
        call place_car(net%roads(3)%lane(1), 5, 1)
        call place_car(net%roads(4)%lane(2), 4, 0)
    end subroutine seed_layout

    subroutine place_car(ln, pos, vel)
        type(lane_t), intent(inout) :: ln
        integer, intent(in) :: pos, vel
        ln%cells(pos)%has_car = .true.
        ln%cells(pos)%turning_intent = V_STRAIGHT
        ln%cells(pos)%velocity = vel
    end subroutine place_car

    subroutine print_state(net, step)
        type(road_network_t), intent(in) :: net
        integer, intent(in) :: step
        character(len=L) :: out_lane, in_lane
        integer :: r

        print '(a,i3)', '--- Step ', step
        do r = 1, size(net%roads)
            call render_lane(net%roads(r)%lane(1)%cells, out_lane)
            call render_lane(net%roads(r)%lane(2)%cells, in_lane)
            print '(a,i1,a,1x,a)', ' Road ', r, ' outbound:', trim(out_lane)
            print '(a,i1,a,1x,a)', ' Road ', r, ' inbound :', trim(in_lane)
        end do
    end subroutine print_state

    subroutine render_lane(cells, line)
        type(cell), intent(in) :: cells(:)
        character(len=*), intent(out) :: line
        integer :: i, v
        do i = 1, size(cells)
            if (cells(i)%has_car) then
                v = min(cells(i)%velocity, 9)
                line(i:i) = achar(ichar('0') + v)
            else
                line(i:i) = '.'
            end if
        end do
    end subroutine render_lane

    subroutine check_step(net, step)
        type(road_network_t), intent(in) :: net
        integer, intent(in) :: step
        integer :: r, ln
        integer :: n_old, n_new

        n_old = count_occupied_old(net)
        n_new = count_occupied_network(net)
        if (n_old /= n_new) then
            print *, 'FAIL: mass not conserved at step ', step
            stop 1
        end if

        do r = 1, size(net%roads)
            do ln = 1, size(net%roads(r)%lane)
                call check_lane_ns(net%roads(r)%lane(ln)%old, &
                                   net%roads(r)%lane(ln)%cells, step, r, ln)
            end do
        end do
    end subroutine check_step

    integer function count_occupied_old(net) result(n)
        type(road_network_t), intent(in) :: net
        integer :: r, ln, k
        n = 0
        do r = 1, size(net%roads)
            do ln = 1, size(net%roads(r)%lane)
                do k = 1, net%roads(r)%lane(ln)%length
                    if (net%roads(r)%lane(ln)%old(k)%has_car) then
                        n = n + 1
                    end if
                end do
            end do
        end do
    end function count_occupied_old

    subroutine check_lane_ns(old_cells, new_cells, step, road_id, lane_id)
        type(cell), intent(in) :: old_cells(:), new_cells(:)
        integer, intent(in) :: step, road_id, lane_id
        integer :: n_old, n_new, i, idx, gap, v_calc, pos1, pos2, new_v, len
        integer, allocatable :: old_pos(:), old_v(:), new_pos(:)

        len = size(old_cells)
        n_old = count(old_cells%has_car)
        n_new = count(new_cells%has_car)
        if (n_old /= n_new) then
            print *, 'FAIL: lane count mismatch at step ', step, &
                     ' road ', road_id, ' lane ', lane_id
            stop 1
        end if
        if (n_old == 0) return

        allocate(old_pos(n_old), old_v(n_old), new_pos(n_new))
        idx = 0
        do i = 1, len
            if (old_cells(i)%has_car) then
                idx = idx + 1
                old_pos(idx) = i
                old_v(idx) = old_cells(i)%velocity
            end if
        end do
        idx = 0
        do i = 1, len
            if (new_cells(i)%has_car) then
                idx = idx + 1
                new_pos(idx) = i
            end if
        end do

        do i = 1, n_old
            if (new_pos(i) < old_pos(i)) then
                print *, 'FAIL: order inversion at step ', step, &
                         ' road ', road_id, ' lane ', lane_id
                stop 1
            end if
            if (i < n_old) then
                gap = old_pos(i+1) - old_pos(i) - 1
            else
                gap = len - old_pos(i)
            end if
            v_calc = min(old_v(i) + 1, V_MAX)
            v_calc = min(v_calc, gap)
            pos1 = old_pos(i) + v_calc
            pos2 = old_pos(i) + max(v_calc - 1, 0)
            if (.not. (new_pos(i) == pos1 .or. new_pos(i) == pos2)) then
                print *, 'FAIL: invalid move at step ', step, &
                         ' road ', road_id, ' lane ', lane_id
                stop 1
            end if
            new_v = new_cells(new_pos(i))%velocity
            if (new_v /= new_pos(i) - old_pos(i)) then
                print *, 'FAIL: velocity mismatch at step ', step, &
                         ' road ', road_id, ' lane ', lane_id
                stop 1
            end if
        end do

        deallocate(old_pos, old_v, new_pos)
    end subroutine check_lane_ns

end program test_ns_grid
