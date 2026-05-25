module tasep_model
    ! Network-wide open-boundary TASEP cellular automaton update module.
    !
    ! This module advances every lane in a ``road_network_t`` by one parallel
    ! TASEP time step. Each lane is a 1-D lattice of ``cell`` objects; each site
    ! is either empty or occupied by one vehicle. Vehicles move only toward the
    ! downstream end (site ``L``) and obey the exclusion rule: a vehicle may hop
    ! forward only if the next site is empty.
    !
    ! Boundary conditions on a lane are controlled by two flags:
    !
    ! - ``open_in``  -- site 1 accepts new vehicles with probability ``alpha``.
    ! - ``open_out`` -- site ``L`` releases vehicles with probability ``beta``.
    !
    ! When a flag is ``.false.``, the corresponding boundary is closed (the lane
    ! connects to a junction that handles cross-road transfers separately).
    !
    ! The lane density and time-averaged current are defined as:
    !
    ! .. math::
    !
    !    \rho = \frac{N}{L},
    !    \qquad
    !    J = \frac{1}{T}\sum_{t=1}^{T} n_{\mathrm{exit}}(t),
    !
    ! where ``N`` is the number of occupied sites, ``L`` is the lane length,
    ! and ``n_exit(t)`` is the number of vehicles that left via ``open_out``
    ! at step ``t``.

    use road_network_mod

    implicit none
    private

    public :: tasep_lane_step

contains
    subroutine tasep_lane_step(net)
        ! Advance every lane in the road network by one parallel TASEP step.
        !
        ! For each road and each lane within that road, this routine performs a
        ! parallel-update TASEP step. All movement decisions are based on the
        ! ``old`` snapshot that must have been set by ``snapshot_network`` before
        ! this routine is called.
        !
        ! The update for each lane proceeds in three stages:
        !
        ! 1. **Clear old positions.** Every site that was occupied in ``old`` is
        !    set to empty in the working copy ``updated``, ready to receive cars
        !    at their new positions.
        !
        ! 2. **Bulk movement.** Each car in ``old`` attempts to hop from site
        !    ``i`` to site ``i+1`` if site ``i+1`` was empty in ``old``. Cars
        !    that cannot hop remain at ``i``. If the lane has ``open_out = .true.``
        !    and the car is at site ``L``, it exits (is removed) with probability
        !    ``beta`` instead.
        !
        !    Cars at site ``L`` that were already moved out by the junction step
        !    (detected via ``junc_moved_L``) are skipped to avoid double-processing.
        !
        ! 3. **Open inflow.** If ``open_in = .true.`` and site 1 was empty in
        !    ``old`` *and* remains empty after the bulk update, a new vehicle
        !    enters at site 1 with probability ``alpha``.
        !
        ! The boundary rates satisfy:
        !
        ! .. math::
        !
        !    0 \leq \alpha \leq 1,
        !    \qquad
        !    0 \leq \beta \leq 1.
        !
        type(road_network_t), intent(inout) :: net
        integer :: r, k, L, i_site
        real    :: rnd
        logical :: exit_now, junc_moved_L
        type(cell), allocatable :: updated(:)

        do r = 1, size(net%roads)
            do k = 1, size(net%roads(r)%lane)
                L = net%roads(r)%lane(k)%length
                allocate(updated(L))
                updated = net%roads(r)%lane(k)%cells

                junc_moved_L = (net%roads(r)%lane(k)%old(L)%has_car .and. &
                                .not. net%roads(r)%lane(k)%cells(L)%has_car .and. &
                                .not. net%roads(r)%lane(k)%open_out)

                ! Clear positions occupied at step start.
                do i_site = 1, L
                    if (net%roads(r)%lane(k)%old(i_site)%has_car) then
                        updated(i_site)%has_car  = .false.
                        updated(i_site)%velocity = 0
                    end if
                end do

                ! TASEP movement from old state.
                do i_site = 1, L
                    if (.not. net%roads(r)%lane(k)%old(i_site)%has_car) cycle
                    if (i_site == L .and. junc_moved_L) cycle

                    exit_now = .false.
                    if (net%roads(r)%lane(k)%open_out .and. i_site == L) then
                        call random_number(rnd)
                        if (rnd < net%roads(r)%lane(k)%beta) exit_now = .true.
                    end if
                    if (exit_now) cycle

                    if (i_site < L .and. .not. net%roads(r)%lane(k)%old(i_site+1)%has_car) then
                        ! Hop right by one cell.
                        updated(i_site+1)%has_car  = .true.
                        updated(i_site+1)%velocity = 1
                    else
                        ! Blocked or holding cell: stay.
                        updated(i_site)%has_car  = .true.
                        updated(i_site)%velocity = 0
                    end if
                end do

                ! Open inflow at site 1 (only if empty at step start).
                if (net%roads(r)%lane(k)%open_in .and. &
                    .not. net%roads(r)%lane(k)%old(1)%has_car .and. &
                    .not. updated(1)%has_car) then
                    call random_number(rnd)
                    if (rnd < net%roads(r)%lane(k)%alpha) then
                        updated(1)%has_car  = .true.
                        updated(1)%velocity = 1
                    end if
                end if

                net%roads(r)%lane(k)%cells = updated
                deallocate(updated)
            end do
        end do
    end subroutine tasep_lane_step

end module tasep_model