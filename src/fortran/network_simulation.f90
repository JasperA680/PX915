module network_simulation_mod
    !--------------------------------------------------------------------
    ! Top-level driver for one network timestep and a multi-step run.
    !
    ! Per-timestep order of operations (parallel update):
    !   1. snapshot_network      -- freeze all `old` arrays.
    !   2. apply_lane_changes    -- lateral lane changes (optional).
    !      snapshot_network      -- re-freeze so junction/NS see post-LC state.
    !   3. evaluate_junctions    -- move holding-cell vehicles to destinations.
    !   4. lane_internal_step    -- bulk hops + open boundary inflow/outflow.
    !
    ! Lane changing (step 2) follows Rickert et al. (1996): it is a strict
    ! parallel sub-step where all decisions are read from the pre-LC `old`
    ! snapshot and results are written to `cells`.  Re-snapshotting before
    ! the junction/NS step ensures those steps see the post-LC configuration.
    !
    ! Lane-change model constants (exported from lane_change_mod):
    !   LC_DISABLED  (-1) -- no lane changes (default)
    !   LC_SYMMETRIC  (0) -- symmetric model (T1+T2+T3+T4 for any move)
    !   LC_ASYMMETRIC (1) -- asymmetric model (vehicles prefer rightmost lane)
    !
    ! network_step and run_network accept optional lc_model / lc_p_change
    ! arguments; omitting them disables lane changing for back-compatibility.
    !--------------------------------------------------------------------
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

    subroutine network_step(net, model, lc_model, lc_p_change)
        type(road_network_t), intent(inout) :: net
        integer, optional,    intent(in)    :: lc_model
        real,    optional,    intent(in)    :: lc_p_change
        character(len=16) :: model
        integer :: model_int
        real    :: p_change_val

        model_int    = LC_DISABLED
        p_change_val = 1.0
        if (present(lc_model))    model_int    = lc_model
        if (present(lc_p_change)) p_change_val = lc_p_change

        call snapshot_network(net)
        if (model_int /= LC_DISABLED) then
            call apply_lane_changes(net, model_int, p_change_val)
            call snapshot_network(net)   ! re-freeze with post-LC state
        end if
        call evaluate_junctions(net)

        select case (trim(model))
        case ('TASEP')
            call tasep_lane_step(net)
        case ('NS')
            call NS_model_step(net)
        case default                        ! Just default to NS model if input not recognised (could change to an error code or something?)
            call NS_model_step(net)
        end select

    end subroutine network_step

    !-----------------------------------------------------------------
    ! Per-lane TASEP step for the networked model.
    !
    ! Each vehicle hops exactly one cell to the right if the next cell
    ! was empty in the old (pre-step) state.  Open boundaries:
    !   - Inbound lanes: entry at site 1 with alpha.
    !   - Outbound lanes: exit at site L with beta.
    ! Junction holding-cell logic mirrors NS_model_step.
    !-----------------------------------------------------------------
    subroutine tasep_lane_step(net)
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

end module network_simulation_mod
