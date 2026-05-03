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

    public :: network_step, run_network

contains

    subroutine network_step(net)
        type(road_network_t), intent(inout) :: net
        call snapshot_network(net)
        call evaluate_junctions(net)
        call lane_internal_step(net)
    end subroutine network_step

!-----------------------------------------------------------------
    ! Per-lane TASEP step for the networked model.
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
        integer :: r, k, L, i_site
        real    :: rnd

        do r = 1, size(net%roads)
            do k = 1, 2
                L = net%roads(r)%lane(k)%length

                ! Open outflow at site L.
                if (net%roads(r)%lane(k)%open_out .and. &
                    net%roads(r)%lane(k)%old(L) /= V_EMPTY) then
                    call random_number(rnd)
                    if (rnd < net%roads(r)%lane(k)%beta) &
                        net%roads(r)%lane(k)%cells(L) = V_EMPTY
                end if

                ! Bulk hops: use old state so the update is strictly parallel.
                do i_site = L-1, 1, -1
                    if (net%roads(r)%lane(k)%old(i_site)   /= V_EMPTY .and. &
                        net%roads(r)%lane(k)%old(i_site+1) == V_EMPTY) then
                        net%roads(r)%lane(k)%cells(i_site+1) = net%roads(r)%lane(k)%old(i_site)
                        net%roads(r)%lane(k)%cells(i_site)   = V_EMPTY
                    end if
                end do

                ! Open inflow at site 1 (only if site 1 was empty at step start).
                if (net%roads(r)%lane(k)%open_in .and. &
                    net%roads(r)%lane(k)%old(1) == V_EMPTY) then
                    call random_number(rnd)
                    if (rnd < net%roads(r)%lane(k)%alpha) then
                        net%roads(r)%lane(k)%cells(1) = &
                            sample_indicator(net%roads(r)%lane(k)%p_left, &
                                             net%roads(r)%lane(k)%p_right)
                    end if
                end if
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
        integer :: step, r, L1
        integer :: road_count

        road_count = size(net%roads)
        entries = 0
        exits   = 0

        do step = 1, n_steps
            call network_step(net)

            do r = 1, road_count
                ! Inbound lane (lane 2): count new vehicles entering at site 1.
                if (net%roads(r)%lane(2)%open_in) then
                    if (net%roads(r)%lane(2)%cells(1) /= V_EMPTY .and. &
                        net%roads(r)%lane(2)%old(1)   == V_EMPTY) then
                        entries(r) = entries(r) + 1
                    end if
                end if
                ! Outbound lane (lane 1): count vehicles that exited at site L.
                if (net%roads(r)%lane(1)%open_out) then
                    L1 = net%roads(r)%lane(1)%length
                    if (net%roads(r)%lane(1)%cells(L1) == V_EMPTY .and. &
                        net%roads(r)%lane(1)%old(L1)   /= V_EMPTY) then
                        exits(r) = exits(r) + 1
                    end if
                end if
            end do
        end do
    end subroutine run_network

end module network_simulation_mod
