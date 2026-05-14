! pde_lanechange.f90 — Conservative lane-change source terms for the multi-lane LWR PDE.
!
! Model: asymmetric exchange between adjacent lanes only (lane l exchanges with l-1 and l+1).
! For each adjacent pair (l, l+1):
!
!   S_{l -> l+1, i} = k * rho_l * max(0, 1 - rho_{l+1}/rho_max(l+1))
!                       * max(0, v_{l+1}(rho_{l+1}) - v_l(rho_l))
!   S_{l+1 -> l, i} = k * rho_{l+1} * max(0, 1 - rho_l/rho_max(l))
!                       * max(0, v_l(rho_l) - v_{l+1}(rho_{l+1}))
!
! Velocities are capped at v_limit (the global speed limit):
!   v_eff = min(v_max*(1 - rho/rho_max), v_limit)
! In the speed-limited zone (rho < rho*) both lanes travel at v_limit, so
! Dv = 0 and no lane-changing occurs — a physically correct emergent property.
!
! At most one of S_{l->l+1} and S_{l+1->l} is nonzero per cell (asymmetric Dv factor).
!
! Conservation: sum(source(:, i)) == 0 for every cell i, guaranteed pairwise.
! This holds regardless of how velocities are computed (algebraic cancellation).
!
! Stability note: large k violates source-ODE stability. Rule of thumb: k * dt < 1.
! Consider bounding dt <= min(CFL*dx/max_v, 1/k) in the CFL check.
module pde_lanechange
  implicit none
  private
  public :: compute_lane_change_sources

contains

  subroutine compute_lane_change_sources(density, n_lanes, M, v_max_lanes, rho_max_lanes, &
                                          k, v_limit, source)
    integer, intent(in)  :: n_lanes, M
    real,    intent(in)  :: density(n_lanes, M)
    real,    intent(in)  :: v_max_lanes(n_lanes)
    real,    intent(in)  :: rho_max_lanes(n_lanes)
    real,    intent(in)  :: k
    real,    intent(in)  :: v_limit   ! global speed cap; v_limit = v_max means no restriction
    real,    intent(out) :: source(n_lanes, M)

    integer :: lane, i
    real    :: rho_l, rho_m, v_l, v_m, s_lm, s_ml, gap_l, gap_m

    source = 0.0

    do lane = 1, n_lanes - 1
      do i = 1, M
        rho_l = density(lane,   i)
        rho_m = density(lane+1, i)

        ! Speed-limited Greenshields velocity: min(v_max*(1-rho/rho_max), v_limit).
        ! Both lanes cap at v_limit in free-flow → Δv = 0 → no lane-change incentive.
        v_l = min(v_max_lanes(lane)   * (1.0 - rho_l / rho_max_lanes(lane)),   v_limit)
        v_m = min(v_max_lanes(lane+1) * (1.0 - rho_m / rho_max_lanes(lane+1)), v_limit)

        ! Headroom in the receiving lane (capacity fraction available)
        gap_l = max(0.0, 1.0 - rho_l / rho_max_lanes(lane))
        gap_m = max(0.0, 1.0 - rho_m / rho_max_lanes(lane+1))

        ! Lane l -> l+1: vehicles move right if lane l+1 is faster and has space
        s_lm = k * rho_l * gap_m * max(0.0, v_m - v_l)

        ! Lane l+1 -> l: vehicles move left if lane l is faster and has space
        s_ml = k * rho_m * gap_l * max(0.0, v_l - v_m)

        ! Conservative: what lane l loses is exactly what lane l+1 gains, and vice versa
        source(lane,   i) = source(lane,   i) + s_ml - s_lm
        source(lane+1, i) = source(lane+1, i) + s_lm - s_ml
      end do
    end do

  end subroutine compute_lane_change_sources

end module pde_lanechange
