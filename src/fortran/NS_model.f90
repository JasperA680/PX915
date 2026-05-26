module NS_model
    ! Update logic for one step in a Nagel-Schreckenberg type traffic model.
    !
    ! This module takes a road network and progresses all vehicles within it
    ! by one step under the rules of a Nagel-Schreckenberg model.
    !
    ! These rules are applied sequentially as follows:
    !
    !   1. Acceleration: The cars increase their velocity by 1, up to the speed limit, :math:`v_i^{(t+\frac{1}{3})} = \text{min}(v_i^t+1, v_{max})`
    !   2. Deceleration: The cars break if their velocity would cause a collision with the car in front, :math:`v_i^{(t+\frac{2}{3})} = \text{min}(v_i^{(t+\frac{1}{3})}, g_i^t)`
    !   3. Randomisation: With a random chance P, the cars decrease their velocity by 1, :math:`v_i^{(t+1)} = \text{max}(v_i^{(t+\frac{2}{3})}-1, 0)`
    ! The cars then update their positions, :math:`x^{(t+1)}_i = x_i^t + v_i^{(t+1)}`


    use vehicle_mod
    use road_network_mod

    implicit none
    private

    public :: NS_model_step


contains
    
    subroutine NS_model_step(net, v_max, p_slow)
        ! Per-lane Nagel-Schreckenberg step for the networked model.
        ! Takes a road network as input, and updates all car positions in-place for one step.
        !
        ! Inbound lane (junction at site L, open at site 1):
        !
        !   * Bulk hops use old state; site L is NOT explicitly exited here
        !     (the junction step handles that).
        !   * Open inflow at site 1 with alpha and stochastic indicator.
        !
        ! Outbound lane (junction at site 1, open at site L):
        !
        !   * Bulk hops use old state.
        !   * Open outflow at site L with beta.
        !   * Site 1 is NOT entered here (junction step places vehicles there).

        type(road_network_t), intent(inout) :: net
        integer, intent(in) :: v_max
        real, intent(in)    :: p_slow

        integer :: r, k, L, i_site, j, gap, v, new_pos
        real    :: rnd
        logical :: exit_now, junc_moved_L
        type(cell), allocatable :: updated(:)

        do r = 1, size(net%roads)
            do k = 1, size(net%roads(r)%lane)
                L = net%roads(r)%lane(k)%length
                allocate(updated(L))

                if (net%roads(r)%lane(k)%is_periodic) then
                    ! Periodic ring: site 1 follows site L with no junction or open
                    ! boundary interaction. Gap and position both wrap modulo L.
                    updated%has_car = .false.
                    updated%velocity = 0

                    do i_site = 1, L
                        if (.not. net%roads(r)%lane(k)%old(i_site)%has_car) cycle

                        gap = L - 1
                        do j = 1, L - 1
                            if (net%roads(r)%lane(k)%old(1 + mod(i_site - 1 + j, L))%has_car) then
                                gap = j - 1
                                exit
                            end if
                        end do

                        v = net%roads(r)%lane(k)%old(i_site)%velocity
                        v = min(v + 1, v_max)
                        v = min(v, gap)
                        if (v > 0) then
                            call random_number(rnd)
                            if (rnd < p_slow) v = v - 1
                        end if

                        new_pos = 1 + mod(i_site - 1 + v, L)
                        updated(new_pos)%has_car = .true.
                        updated(new_pos)%velocity = v
                    end do

                    net%roads(r)%lane(k)%cells = updated
                    deallocate(updated)
                    cycle
                end if

                updated = net%roads(r)%lane(k)%cells

                ! Detect whether the junction already moved the holding-cell vehicle out.
                ! If so, the lane step must not re-place it as a phantom.
                junc_moved_L = (net%roads(r)%lane(k)%old(L)%has_car .and. &
                                .not. net%roads(r)%lane(k)%cells(L)%has_car .and. &
                                .not. net%roads(r)%lane(k)%open_out)

                ! Clear positions that were occupied at the start of the step.
                do i_site = 1, L
                    if (net%roads(r)%lane(k)%old(i_site)%has_car) then
                        updated(i_site)%has_car = .false.
                        updated(i_site)%velocity = 0
                    end if
                end do

                ! Move cars using old state so the update is strictly parallel.
                do i_site = 1, L
                    if (.not. net%roads(r)%lane(k)%old(i_site)%has_car) cycle
                    ! Skip the holding cell when the junction already moved this vehicle.
                    if (i_site == L .and. junc_moved_L) cycle

                    exit_now = .false.
                    if (net%roads(r)%lane(k)%open_out .and. i_site == L) then
                        call random_number(rnd)
                        if (rnd < net%roads(r)%lane(k)%beta) exit_now = .true.
                    end if
                    if (exit_now) cycle

                    if (i_site == L) then
                        ! Boundary cell: exit handled via beta above, so no forward gap.
                        gap = 0
                    else
                        gap = L - i_site
                        do j = i_site + 1, L
                            if (net%roads(r)%lane(k)%old(j)%has_car) then
                                gap = j - i_site - 1
                                exit
                            end if
                        end do
                    end if

                    v = net%roads(r)%lane(k)%old(i_site)%velocity
                    v = min(v + 1, v_max)
                    v = min(v, gap)
                    call random_number(rnd)
                    if (rnd < p_slow .and. v > 0) v = v - 1

                    new_pos = min(i_site + v, L)
                    v = new_pos - i_site
                    updated(new_pos)%has_car = .true.
                    updated(new_pos)%velocity = v
                end do

                ! Open inflow at site 1 (only if site 1 was empty at step start).
                if (net%roads(r)%lane(k)%open_in .and. &
                    .not. net%roads(r)%lane(k)%old(1)%has_car .and. &
                    .not. updated(1)%has_car) then
                    call random_number(rnd)
                    if (rnd < net%roads(r)%lane(k)%alpha) then
                        updated(1)%has_car = .true.
                        updated(1)%velocity = 0
                    end if
                end if

                net%roads(r)%lane(k)%cells = updated
                deallocate(updated)
            end do
        end do
    end subroutine NS_model_step


end module NS_model