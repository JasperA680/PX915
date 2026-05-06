module network_simulation_mod
    !--------------------------------------------------------------------
    ! Top-level driver for one network timestep and a multi-step run.
    !
    ! Per-timestep order of operations (parallel update):
    !   1. snapshot_network  -- freeze all `old` arrays.
    !   2. evaluate_junctions -- move holding-cell vehicles to destinations.
    !   3. lane_internal_step -- bulk hops + open boundary inflow/outflow.
    !
    ! This ordering ensures coherent parallel-update semantics:
    !   * junction_step reads from `old` (pre-step state).
    !   * lane_internal_step also reads from `old` for hop decisions.
    !   * A vehicle placed at site 1 by the junction this step stays
    !     there until the next step (it was absent in `old`).
    !   * A vehicle that exits via the junction this step leaves a gap
    !     at site L; the vehicle behind it cannot advance into that gap
    !     this step (it sees site L occupied in `old`).
    !--------------------------------------------------------------------
    use vehicle_mod
    use road_network_mod
    use junction_mod
    implicit none
    private

    integer, parameter :: V_MAX  = 5
    real,    parameter :: P_SLOW = 0.2

    public :: network_step, run_network

contains

    subroutine network_step(net)
        type(road_network_t), intent(inout) :: net
        call snapshot_network(net)
        call evaluate_junctions(net)
        call lane_internal_step(net)
    end subroutine network_step

!-----------------------------------------------------------------
    ! Per-lane Nagel-Schreckenberg step for the networked model.
    !
    ! Inbound lane (junction at site L, open at site 1):
    !   * Bulk hops use old state; site L is NOT explicitly exited here
    !     (the junction step handles that).
    !   * Open inflow at site 1 with alpha and stochastic indicator.
    !
    ! Outbound lane (junction at site 1, open at site L):
    !   * Bulk hops use old state.
    !   * Open outflow at site L with beta.
    !   * Site 1 is NOT entered here (junction step places vehicles there).
    !-----------------------------------------------------------------
    subroutine lane_internal_step(net)
        type(road_network_t), intent(inout) :: net
        integer :: r, k, L, i_site, j, gap, v, new_pos
        real    :: rnd
        logical :: exit_now
        type(cell), allocatable :: updated(:)

        do r = 1, size(net%roads)
            do k = 1, size(net%roads(r)%lane)
                L = net%roads(r)%lane(k)%length
                allocate(updated(L))
                updated = net%roads(r)%lane(k)%cells

                ! Clear positions that were occupied at the start of the step.
                do i_site = 1, L
                    if (net%roads(r)%lane(k)%old(i_site)%has_car) then
                        updated(i_site)%has_car = .false.
                        updated(i_site)%turning_intent = V_EMPTY
                        updated(i_site)%velocity = 0
                    end if
                end do

                ! Move cars using old state so the update is strictly parallel.
                do i_site = 1, L
                    if (.not. net%roads(r)%lane(k)%old(i_site)%has_car) cycle

                    exit_now = .false.
                    if (net%roads(r)%lane(k)%open_out .and. i_site == L) then
                        call random_number(rnd)
                        if (rnd < net%roads(r)%lane(k)%beta) exit_now = .true.
                    end if
                    if (exit_now) cycle

                    if (i_site == L) then
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
                    v = min(v + 1, V_MAX)
                    v = min(v, gap)
                    call random_number(rnd)
                    if (rnd < P_SLOW .and. v > 0) v = v - 1

                    new_pos = i_site + v
                    updated(new_pos)%has_car = .true.
                    updated(new_pos)%turning_intent = net%roads(r)%lane(k)%old(i_site)%turning_intent
                    updated(new_pos)%velocity = v
                end do

                ! Open inflow at site 1 (only if site 1 was empty at step start).
                if (net%roads(r)%lane(k)%open_in .and. &
                    .not. net%roads(r)%lane(k)%old(1)%has_car) then
                    call random_number(rnd)
                    if (rnd < net%roads(r)%lane(k)%alpha) then
                        call random_number(rnd)
                        if (rnd < net%roads(r)%lane(k)%p_left) then      
                            updated(1)%turning_intent = V_LEFT
                        else if (rnd < net%roads(r)%lane(k)%p_left &
                            + net%roads(r)%lane(k)%p_right) then
                            updated(1)%turning_intent = V_RIGHT
                        else
                            updated(1)%turning_intent = V_STRAIGHT
                        end if
                        updated(1)%has_car = .true.
                        updated(1)%velocity = 0
                    end if
                end if

                net%roads(r)%lane(k)%cells = updated
                deallocate(updated)
            end do
        end do
    end subroutine lane_internal_step

    !-----------------------------------------------------------------
    ! Run the network for n_steps, returning per-road entry/exit counts.
    !-----------------------------------------------------------------
    subroutine run_network(net, n_steps, entries, exits)
        type(road_network_t), intent(inout) :: net
        integer,              intent(in)    :: n_steps
        integer,              intent(out)   :: entries(:), exits(:)
        integer :: step, r, k, Lk
        integer :: road_count

        road_count = size(net%roads)
        entries = 0
        exits   = 0

        do step = 1, n_steps
            call network_step(net)

            do r = 1, road_count
                do k = 1, size(net%roads(r)%lane)
                    ! Count entry: new vehicle appeared at open-in site 1.
                    if (net%roads(r)%lane(k)%open_in) then
                        if (net%roads(r)%lane(k)%cells(1)%has_car .and. &
                            .not. net%roads(r)%lane(k)%old(1)%has_car) then
                            entries(r) = entries(r) + 1
                        end if
                    end if
                    ! Count exit: vehicle disappeared from open-out site L.
                    if (net%roads(r)%lane(k)%open_out) then
                        Lk = net%roads(r)%lane(k)%length
                        if (.not. net%roads(r)%lane(k)%cells(Lk)%has_car &
                            .and. net%roads(r)%lane(k)%old(Lk)%has_car) then
                            exits(r) = exits(r) + 1
                        end if
                    end if
                end do
            end do
        end do
    end subroutine run_network

end module network_simulation_mod
