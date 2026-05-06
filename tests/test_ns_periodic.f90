program test_ns_periodic
    !--------------------------------------------------------------------
    ! Single-road Nagel-Schreckenberg test with periodic boundaries.
    ! Prints a ring snapshot for visual inspection.
    !--------------------------------------------------------------------
    use vehicle_mod
    use road_network_mod
    implicit none

    integer, parameter :: L = 18
    integer, parameter :: N_STEPS = 6
    integer, parameter :: V_MAX = 5
    real,    parameter :: P_SLOW = 0.2

    type(cell) :: state(L), old_state(L)
    integer :: step, n_expected

    call seed_rng(20260507)
    call init_ring(state)
    n_expected = count(state%has_car)

    print *, 'NS periodic ring test'
    call print_ring(state, 0)

    do step = 1, N_STEPS
        call periodic_step(state, old_state)
        call check_step(old_state, state, step, n_expected)
        call print_ring(state, step)
    end do

    print *, 'test_ns_periodic: OK'

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

    subroutine init_ring(state)
        type(cell), intent(out) :: state(:)
        state%has_car = .false.
        state%turning_intent = V_EMPTY
        state%velocity = 0
        call place_car(state, 3, 1)
        call place_car(state, 8, 2)
        call place_car(state, 14, 0)
    end subroutine init_ring

    subroutine place_car(state, pos, vel)
        type(cell), intent(inout) :: state(:)
        integer, intent(in) :: pos, vel
        state(pos)%has_car = .true.
        state(pos)%turning_intent = V_STRAIGHT
        state(pos)%velocity = vel
    end subroutine place_car

    subroutine periodic_step(state, old_state)
        type(cell), intent(inout) :: state(:)
        type(cell), intent(out) :: old_state(:)
        type(cell) :: new_state(size(state))
        integer :: i, len, gap, v, new_pos
        real :: rnd

        len = size(state)
        old_state = state
        new_state = old_state

        do i = 1, len
            if (old_state(i)%has_car) then
                new_state(i)%has_car = .false.
                new_state(i)%turning_intent = V_EMPTY
                new_state(i)%velocity = 0
            end if
        end do

        do i = 1, len
            if (.not. old_state(i)%has_car) cycle
            gap = gap_to_next(old_state, i)
            v = min(old_state(i)%velocity + 1, V_MAX)
            v = min(v, gap)
            call random_number(rnd)
            if (rnd < P_SLOW .and. v > 0) v = v - 1
            new_pos = 1 + mod(i - 1 + v, len)
            new_state(new_pos)%has_car = .true.
            new_state(new_pos)%turning_intent = old_state(i)%turning_intent
            new_state(new_pos)%velocity = v
        end do

        state = new_state
    end subroutine periodic_step

    integer function gap_to_next(old_state, pos) result(gap)
        type(cell), intent(in) :: old_state(:)
        integer, intent(in) :: pos
        integer :: len, offset, j
        len = size(old_state)
        gap = len - 1
        do offset = 1, len - 1
            j = 1 + mod(pos - 1 + offset, len)
            if (old_state(j)%has_car) then
                gap = offset - 1
                exit
            end if
        end do
    end function gap_to_next

    subroutine check_step(old_state, new_state, step, n_expected)
        type(cell), intent(in) :: old_state(:), new_state(:)
        integer, intent(in) :: step, n_expected
        integer :: len, i, gap, v_calc, pos1, pos2, new_pos, move_v
        logical, allocatable :: matched(:)

        len = size(old_state)
        if (count(old_state%has_car) /= n_expected) then
            print *, 'FAIL: old-state count mismatch at step ', step
            stop 1
        end if
        if (count(new_state%has_car) /= n_expected) then
            print *, 'FAIL: new-state count mismatch at step ', step
            stop 1
        end if

        allocate(matched(len))
        matched = .false.

        do i = 1, len
            if (.not. old_state(i)%has_car) cycle
            gap = gap_to_next(old_state, i)
            v_calc = min(old_state(i)%velocity + 1, V_MAX)
            v_calc = min(v_calc, gap)
            pos1 = 1 + mod(i - 1 + v_calc, len)
            pos2 = 1 + mod(i - 1 + max(v_calc - 1, 0), len)
            new_pos = 0

            if (new_state(pos1)%has_car .and. .not. matched(pos1)) then
                new_pos = pos1
            else if (new_state(pos2)%has_car .and. .not. matched(pos2)) then
                new_pos = pos2
            end if

            if (new_pos == 0) then
                print *, 'FAIL: periodic move mismatch at step ', step, ' pos ', i
                stop 1
            end if

            matched(new_pos) = .true.
            move_v = mod(new_pos - i + len, len)
            if (new_state(new_pos)%velocity /= move_v) then
                print *, 'FAIL: velocity mismatch at step ', step
                stop 1
            end if
        end do

        deallocate(matched)
    end subroutine check_step

    subroutine print_ring(state, step)
        type(cell), intent(in) :: state(:)
        integer, intent(in) :: step
        character(len=L) :: line
        call render_lane(state, line)
        print '(a,i3,1x,a)', 'Ring step', step, trim(line)
    end subroutine print_ring

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

end program test_ns_periodic
