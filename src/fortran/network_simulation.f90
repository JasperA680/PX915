module network_simulation_mod
    ! Top-level time-step driver for road-network cellular automaton simulations.
    !
    ! This module coordinates one full network update. It combines optional
    ! lateral lane changing, junction movement and longitudinal vehicle motion
    ! within each lane.
    !
    ! The network update is designed to preserve parallel-update semantics. Each
    ! sub-step reads from a snapshotted ``old`` state and writes to the live
    ! ``cells`` arrays. The per-time-step order is:
    !
    ! 1. ``snapshot_network`` freezes the current network state.
    ! 2. ``apply_lane_changes`` applies optional lateral lane changes.
    ! 3. ``snapshot_network`` is called again so that junction and longitudinal
    !    updates see the post-lane-change state.
    ! 4. ``evaluate_junctions`` moves approved holding-cell vehicles through
    !    junctions.
    ! 5. The selected longitudinal model updates all lanes, either using TASEP or
    !    Nagel-Schreckenberg dynamics.
    !
    ! Lane-changing is optional. If ``lc_model`` is omitted, the update defaults
    ! to ``LC_DISABLED`` for backward compatibility.
    !
    ! The available lane-change model constants are imported from
    ! ``lane_change_mod``:
    !
    ! * ``LC_DISABLED``: no lane changes;
    ! * ``LC_SYMMETRIC``: symmetric lane changing;
    ! * ``LC_ASYMMETRIC``: asymmetric lane changing with right-lane preference.
    !
    ! For the network TASEP update, the local occupation variable can be written
    ! as:
    !
    ! .. math::
    !
    !    \tau_i =
    !    \begin{cases}
    !    1, & \text{if site } i \text{ contains a vehicle}, \\
    !    0, & \text{otherwise}.
    !    \end{cases}
    !
    ! A vehicle moves from site ``i`` to site ``i+1`` if:
    !
    ! .. math::
    !
    !    \tau_i = 1
    !    \quad \text{and} \quad
    !    \tau_{i+1} = 0.
    !
    ! Open inflow and outflow boundaries are controlled by the lane parameters
    ! ``alpha`` and ``beta``.

    use vehicle_mod
    use road_network_mod
    use junction_mod
    use lane_change_mod
    use NS_model, only: NS_model_step

    implicit none
    private

    public :: network_step
    public :: LC_DISABLED, LC_SYMMETRIC, LC_ASYMMETRIC

contains

    subroutine network_step(net, model, lc_model, lc_p_change, v_max, p_slow)
        ! Advance a road network by one full time step.
        !
        ! This routine is the main network-level update. It applies the optional
        ! lane-changing sub-step, evaluates junction movements, and then performs
        ! the selected longitudinal lane update.
        !
        ! The update order is:
        !
        ! 1. snapshot the network before lane changing;
        ! 2. apply lane changes if ``lc_model`` is not ``LC_DISABLED``;
        ! 3. snapshot again so later sub-steps read the post-lane-change state;
        ! 4. evaluate junctions;
        ! 5. update lane interiors using either ``'TASEP'`` or ``'NS'``.
        !
        ! If ``model`` is not recognised, the routine currently falls back to the
        ! Nagel-Schreckenberg update.
        type(road_network_t), intent(inout) :: net
        ! Road network to update in place.
        character(len=16) :: model
        ! Longitudinal model selector. Supported values are ``'TASEP'`` and ``'NS'``.
        integer, optional, intent(in) :: lc_model
        ! Optional lane-change model selector.
        real, optional, intent(in) :: lc_p_change
        ! Optional probability of accepting a candidate lane change.
        integer, intent(in) :: v_max
        ! Maximum vehicle speed used by the Nagel-Schreckenberg and lane-change models.
        real, intent(in) :: p_slow
        ! Random slowing probability used by the Nagel-Schreckenberg model.

        integer :: model_int
        ! Lane-change model used for this update.
        real :: p_change_val
        ! Lane-change probability used for this update.

        model_int    = LC_DISABLED
        p_change_val = 1.0

        if (present(lc_model))    model_int    = lc_model
        if (present(lc_p_change)) p_change_val = lc_p_change

        call snapshot_network(net)

        if (model_int /= LC_DISABLED) then
            call apply_lane_changes(net, model_int, p_change_val, v_max)
            call snapshot_network(net)
        end if

        call evaluate_junctions(net)

        select case (trim(model))
        case ('TASEP')
            call tasep_lane_step(net)
        case ('NS')
            call NS_model_step(net, v_max, p_slow)
        case default
            call NS_model_step(net, v_max, p_slow)
        end select

    end subroutine network_step


    subroutine tasep_lane_step(net)
        ! Apply the networked TASEP update to every lane.
        !
        ! This routine performs the longitudinal TASEP part of the network update.
        ! It should be called after ``snapshot_network`` and after junction
        ! movements have been evaluated.
        !
        ! Each vehicle attempts to hop exactly one cell to the right. The move is
        ! accepted if the next cell was empty in the old state:
        !
        ! .. math::
        !
        !    \tau_i^n = 1,
        !    \quad
        !    \tau_{i+1}^n = 0
        !    \quad \Longrightarrow \quad
        !    \tau_{i+1}^{n+1} = 1.
        !
        ! For an open outflow lane, a vehicle at site ``L`` exits with probability
        ! ``beta``:
        !
        ! .. math::
        !
        !    P(\mathrm{exit}) = \beta.
        !
        ! For an open inflow lane, a vehicle enters site ``1`` with probability
        ! ``alpha`` if the site is empty:
        !
        ! .. math::
        !
        !    P(\mathrm{entry}) = \alpha.
        !
        ! The final site ``L`` may also act as a junction holding cell. If a
        ! junction has already moved a vehicle out of this holding cell during the
        ! current network step, the longitudinal TASEP update does not reinsert it.
        type(road_network_t), intent(inout) :: net
        ! Road network whose lanes are updated in place.

        integer :: r
        ! Road index.
        integer :: k
        ! Lane index within the road.
        integer :: L
        ! Length of the current lane.
        integer :: i_site
        ! Cell index within the current lane.
        real :: rnd
        ! Uniform random number used for stochastic open-boundary events.
        logical :: exit_now
        ! Whether a vehicle exits at the current open outflow boundary.
        logical :: junc_moved_L
        ! Whether the junction update already removed the holding-cell vehicle.
        type(cell), allocatable :: updated(:)
        ! Temporary next-state cell array for the current lane.

        do r = 1, size(net%roads)
            do k = 1, size(net%roads(r)%lane)
                L = net%roads(r)%lane(k)%length

                allocate(updated(L))
                updated = net%roads(r)%lane(k)%cells

                junc_moved_L = (net%roads(r)%lane(k)%old(L)%has_car .and. &
                                .not. net%roads(r)%lane(k)%cells(L)%has_car .and. &
                                .not. net%roads(r)%lane(k)%open_out)

                ! Clear positions occupied at the start of the step.
                do i_site = 1, L
                    if (net%roads(r)%lane(k)%old(i_site)%has_car) then
                        updated(i_site)%has_car  = .false.
                        updated(i_site)%velocity = 0
                    end if
                end do

                ! Apply one-cell TASEP movement from the old state.
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
                        ! Blocked vehicle or occupied holding cell remains in place.
                        updated(i_site)%has_car  = .true.
                        updated(i_site)%velocity = 0
                    end if
                end do

                ! Open inflow at site 1.
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

end module network_simulation_mod