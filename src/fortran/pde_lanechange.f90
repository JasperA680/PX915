module pde_lanechange
  ! Conservative lane-change source terms for the multi-lane LWR PDE.
  !
  ! This module computes source terms that transfer traffic density between
  ! adjacent lanes. The source terms are conservative: for each spatial cell,
  ! the total gain over all lanes is equal to the total loss over all lanes.
  !
  ! Only adjacent lanes exchange density. For each neighbouring pair of lanes
  ! ``l`` and ``l+1``, the lane-change rate from lane ``l`` to lane ``l+1`` is
  ! modelled as:
  !
  ! .. math::
  !
  !    S_{l \rightarrow l+1, i}
  !    =
  !    k \rho_{l,i}
  !    \max\left(0, 1 - \frac{\rho_{l+1,i}}{\rho_{\max,l+1}}\right)
  !    \max\left(0, v_{l+1}(\rho_{l+1,i}) - v_l(\rho_{l,i})\right).
  !
  ! The reverse exchange is:
  !
  ! .. math::
  !
  !    S_{l+1 \rightarrow l, i}
  !    =
  !    k \rho_{l+1,i}
  !    \max\left(0, 1 - \frac{\rho_{l,i}}{\rho_{\max,l}}\right)
  !    \max\left(0, v_l(\rho_{l,i}) - v_{l+1}(\rho_{l+1,i})\right).
  !
  ! The velocity in each lane is a speed-limited Greenshields velocity:
  !
  ! .. math::
  !
  !    v_l(\rho_l)
  !    =
  !    \min\left[
  !      v_{\max,l}\left(1 - \frac{\rho_l}{\rho_{\max,l}}\right),
  !      v_{\mathrm{limit}}
  !    \right].
  !
  ! The source term is conservative because each pairwise exchange is added to
  ! one lane and subtracted from the other:
  !
  ! .. math::
  !
  !    \sum_l S_{l,i} = 0
  !
  ! for every spatial cell ``i``.
  !
  ! Stability note: large lane-change rates can make the source update unstable.
  ! A useful rule of thumb is ``k*dt < 1``. In a fully adaptive solver, this can
  ! be enforced by also bounding the time step by the source time scale.

  implicit none
  private

  public :: compute_lane_change_sources

contains

  subroutine compute_lane_change_sources(density, n_lanes, M, v_max_lanes, rho_max_lanes, &
                                          k, v_limit, source)
    ! Compute conservative lane-change source terms.
    !
    ! This routine evaluates pairwise density exchange between adjacent lanes.
    ! For each pair of neighbouring lanes, vehicles preferentially move towards
    ! the faster lane, provided that the receiving lane has available capacity.
    !
    ! At most one direction of exchange is active for a given lane pair and
    ! spatial cell, because the velocity-difference factors are asymmetric:
    !
    ! .. math::
    !
    !    \max(0, v_m - v_l)
    !
    ! and
    !
    ! .. math::
    !
    !    \max(0, v_l - v_m).
    !
    ! This means that if lane ``m`` is faster than lane ``l``, only the
    ! ``l -> m`` transfer is active. If lane ``l`` is faster than lane ``m``,
    ! only the ``m -> l`` transfer is active.
    !
    ! The output array ``source`` has the same lane and space indexing as
    ! ``density``. It can be added to the finite-volume update as:
    !
    ! .. math::
    !
    !    \rho_l^{n+1}
    !    =
    !    \rho_l^*
    !    +
    !    \Delta t\,S_l.
    integer, intent(in) :: n_lanes  ! Number of lanes.
    integer, intent(in) :: M        ! Number of spatial cells.
    real, intent(in) :: density(n_lanes, M)
    ! Lane-resolved traffic density, indexed as ``density(lane, i)``.
    real, intent(in) :: v_max_lanes(n_lanes)
    ! Maximum/free-flow velocity for each lane.
    real, intent(in) :: rho_max_lanes(n_lanes)
    ! Maximum/jam density for each lane.
    real, intent(in) :: k
    ! Lane-change rate coefficient.
    real, intent(in) :: v_limit
    ! Global speed cap. Set ``v_limit = v_max`` when no speed restriction is used.
    real, intent(out) :: source(n_lanes, M)
    ! Conservative lane-change source terms, indexed as ``source(lane, i)``.

    integer :: lane  ! Index of the lower lane in the adjacent lane pair.
    integer :: i     ! Spatial cell index.

    real :: rho_l    ! Density in lane ``lane`` at cell ``i``.
    real :: rho_m    ! Density in lane ``lane+1`` at cell ``i``.
    real :: v_l      ! Speed-limited velocity in lane ``lane``.
    real :: v_m      ! Speed-limited velocity in lane ``lane+1``.
    real :: s_lm     ! Source transfer rate from lane ``lane`` to lane ``lane+1``.
    real :: s_ml     ! Source transfer rate from lane ``lane+1`` to lane ``lane``.
    real :: gap_l    ! Available capacity fraction in lane ``lane``.
    real :: gap_m    ! Available capacity fraction in lane ``lane+1``.

    source = 0.0

    do lane = 1, n_lanes - 1
      do i = 1, M
        rho_l = density(lane,   i)
        rho_m = density(lane+1, i)

        ! Compute speed-limited Greenshields velocities.
        !
        ! In the free-flow speed-limited regime, both lanes may be capped at
        ! ``v_limit``. In that case the velocity difference is zero and there is
        ! no lane-change incentive from speed differences.
        v_l = min(v_max_lanes(lane)   * (1.0 - rho_l / rho_max_lanes(lane)),   v_limit)
        v_m = min(v_max_lanes(lane+1) * (1.0 - rho_m / rho_max_lanes(lane+1)), v_limit)

        ! Compute available capacity fractions in the receiving lanes.
        gap_l = max(0.0, 1.0 - rho_l / rho_max_lanes(lane))
        gap_m = max(0.0, 1.0 - rho_m / rho_max_lanes(lane+1))

        ! Lane ``lane`` to lane ``lane+1``: vehicles move if lane ``lane+1`` is
        ! faster and has spare capacity.
        s_lm = k * rho_l * gap_m * max(0.0, v_m - v_l)

        ! Lane ``lane+1`` to lane ``lane``: vehicles move if lane ``lane`` is
        ! faster and has spare capacity.
        s_ml = k * rho_m * gap_l * max(0.0, v_l - v_m)

        ! Conservative pairwise update. What one lane loses, the other gains.
        source(lane,   i) = source(lane,   i) + s_ml - s_lm
        source(lane+1, i) = source(lane+1, i) + s_lm - s_ml
      end do
    end do

  end subroutine compute_lane_change_sources

end module pde_lanechange