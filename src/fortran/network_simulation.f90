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
    use tasep_model, only: tasep_lane_step

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

end module network_simulation_mod