module pde_flux
  implicit none
  private
  public :: v_of_rho, q_of_rho, dq_drho, rho_critical
  public :: q_newell, dq_drho_newell, rho_critical_newell
  public :: q_dispatch, dq_drho_dispatch, godunov_dispatch
  public :: lax_friedrichs_flux, godunov_flux

  real, parameter :: NEWELL_W = 0.5  ! congestion wave speed (default)

contains

  ! v_limit caps vehicle speed; v_limit = v_max means no restriction.
  elemental function v_of_rho(rho, v_max, rho_max, v_limit) result(v)
    real, intent(in) :: rho, v_max, rho_max, v_limit
    real :: v
    v = min(v_max * (1.0 - rho / rho_max), v_limit)
  end function v_of_rho

  ! Greenshields flow with speed limit:
  !   q(rho) = v_limit * rho              for rho < rho* = rho_max*(1 - v_limit/v_max)
  !   q(rho) = v_max * rho*(1-rho/rho_max) for rho >= rho*
  ! The diagram remains concave; q_max is still at rho_c = rho_max/2.
  elemental function q_of_rho(rho, v_max, rho_max, v_limit) result(q)
    real, intent(in) :: rho, v_max, rho_max, v_limit
    real :: q
    q = v_of_rho(rho, v_max, rho_max, v_limit) * rho
  end function q_of_rho

  ! Characteristic speed dq/dρ:
  !   v_limit                          for rho < rho* (speed-limited free-flow)
  !   v_max*(1 - 2*rho/rho_max)        for rho >= rho* (standard Greenshields)
  ! Monotonically decreasing, so the diagram is concave throughout.
  elemental function dq_drho(rho, v_max, rho_max, v_limit) result(dq)
    real, intent(in) :: rho, v_max, rho_max, v_limit
    real :: dq, rho_star
    rho_star = rho_max * (1.0 - v_limit / v_max)
    if (rho < rho_star) then
      dq = v_limit
    else
      dq = v_max * (1.0 - 2.0 * rho / rho_max)
    end if
  end function dq_drho

  ! Argmax of q for Greenshields (with or without speed limit): ρ_c = ρ_max / 2.
  elemental function rho_critical(rho_max) result(rc)
    real, intent(in) :: rho_max
    real :: rc
    rc = rho_max / 2.0
  end function rho_critical

  ! Lax-Friedrichs numerical flux — simple but diffusive, good for debugging.
  ! F_LF = [q(ρ_L) + q(ρ_R)] / 2  −  (Δx / 2Δt) * (ρ_R − ρ_L)
  function lax_friedrichs_flux(rho_L, rho_R, v_max, rho_max, v_limit, dx, dt) result(F)
    real, intent(in) :: rho_L, rho_R, v_max, rho_max, v_limit, dx, dt
    real :: F
    F = 0.5 * (q_of_rho(rho_L, v_max, rho_max, v_limit) &
              + q_of_rho(rho_R, v_max, rho_max, v_limit)) &
        - (dx / (2.0 * dt)) * (rho_R - rho_L)
  end function lax_friedrichs_flux

  ! Godunov numerical flux for a concave flux function.
  !
  ! With a speed limit the Greenshields diagram gains a kink at rho* but
  ! remains concave (slope is monotonically decreasing), so the same
  ! case structure applies.  The argmax is still rho_c = rho_max/2.
  !
  ! Closed-form cases (ρ_c = ρ_max/2):
  !   ρ_L ≤ ρ_R  — shock:
  !     both sub-critical:   F = q(ρ_L)
  !     both super-critical: F = q(ρ_R)
  !     ρ_L sub, ρ_R super:  F = min(q(ρ_L), q(ρ_R))
  !   ρ_L >  ρ_R  — rarefaction:
  !     both sub-critical:   F = q(ρ_L)
  !     both super-critical: F = q(ρ_R)
  !     sonic (ρ_R < ρ_c < ρ_L): F = q(ρ_c)
  !
  ! Reference: LeVeque, Finite Volume Methods for Hyperbolic Problems, §12.1
  function godunov_flux(rho_L, rho_R, v_max, rho_max, v_limit) result(F)
    real, intent(in) :: rho_L, rho_R, v_max, rho_max, v_limit
    real :: F, rc, qL, qR
    rc = rho_critical(rho_max)
    qL = q_of_rho(rho_L, v_max, rho_max, v_limit)
    qR = q_of_rho(rho_R, v_max, rho_max, v_limit)
    if (rho_L <= rho_R) then
      if (rho_L >= rc) then
        F = qR
      else if (rho_R <= rc) then
        F = qL
      else
        F = min(qL, qR)
      end if
    else
      if (rho_L <= rc) then
        F = qL
      else if (rho_R >= rc) then
        F = qR
      else
        F = q_of_rho(rc, v_max, rho_max, v_limit)
      end if
    end if
  end function godunov_flux

  ! ---------------------------------------------------------------
  ! Newell-Daganzo (triangular) fundamental diagram with speed limit
  !
  !   q(rho) = min( v_limit * rho,  w * (rho_max - rho) )
  !
  ! Speed limit replaces v_max in the free-flow branch, shifting the
  ! critical density:
  !   rho_c = w * rho_max / (v_limit + w)
  ! v_limit = v_max recovers the original Newell diagram exactly.
  ! ---------------------------------------------------------------

  elemental function q_newell(rho, rho_max, v_limit) result(q)
    real, intent(in) :: rho, rho_max, v_limit
    real :: q
    q = min(v_limit * rho, NEWELL_W * (rho_max - rho))
  end function q_newell

  ! Characteristic speed for Newell with speed limit (piecewise constant):
  !   +v_limit in free-flow,  -NEWELL_W in congested, 0 at the kink.
  elemental function dq_drho_newell(rho, rho_max, v_limit) result(dq)
    real, intent(in) :: rho, rho_max, v_limit
    real :: dq, rc
    rc = rho_critical_newell(rho_max, v_limit)
    if (rho < rc) then
      dq = v_limit
    else if (rho > rc) then
      dq = -NEWELL_W
    else
      dq = 0.0
    end if
  end function dq_drho_newell

  elemental function rho_critical_newell(rho_max, v_limit) result(rc)
    real, intent(in) :: rho_max, v_limit
    real :: rc
    rc = NEWELL_W * rho_max / (v_limit + NEWELL_W)
  end function rho_critical_newell

  ! Godunov flux for Newell-Daganzo with speed limit (concave, piecewise linear).
  function godunov_newell_flux(rho_L, rho_R, rho_max, v_limit) result(F)
    real, intent(in) :: rho_L, rho_R, rho_max, v_limit
    real :: F, rc, qL, qR
    rc = rho_critical_newell(rho_max, v_limit)
    qL = q_newell(rho_L, rho_max, v_limit)
    qR = q_newell(rho_R, rho_max, v_limit)
    if (rho_L <= rho_R) then
      if (rho_L >= rc) then
        F = qR
      else if (rho_R <= rc) then
        F = qL
      else
        F = min(qL, qR)
      end if
    else
      if (rho_L <= rc) then
        F = qL
      else if (rho_R >= rc) then
        F = qR
      else
        F = q_newell(rc, rho_max, v_limit)
      end if
    end if
  end function godunov_newell_flux

  ! ---------------------------------------------------------------
  ! Dispatchers — select the closure based on flux_type string
  ! ---------------------------------------------------------------
  ! Dispatchers use index() so that both 'newell' (Godunov) and 'newell_lf'
  ! (Lax-Friedrichs with Newell closure) route to the Newell q and dq/dρ.
  function q_dispatch(rho, v_max, rho_max, v_limit, flux_type) result(q)
    real,             intent(in) :: rho, v_max, rho_max, v_limit
    character(len=*), intent(in) :: flux_type
    real :: q
    if (index(trim(flux_type), 'newell') > 0) then
      q = q_newell(rho, rho_max, v_limit)
    else
      q = q_of_rho(rho, v_max, rho_max, v_limit)
    end if
  end function q_dispatch

  function dq_drho_dispatch(rho, v_max, rho_max, v_limit, flux_type) result(dq)
    real,             intent(in) :: rho, v_max, rho_max, v_limit
    character(len=*), intent(in) :: flux_type
    real :: dq
    if (index(trim(flux_type), 'newell') > 0) then
      dq = dq_drho_newell(rho, rho_max, v_limit)
    else
      dq = dq_drho(rho, v_max, rho_max, v_limit)
    end if
  end function dq_drho_dispatch

  ! godunov_dispatch is only reached when flux_type is 'godunov' or 'newell'
  ! (never 'newell_lf'), so the exact-match check is sufficient here.
  function godunov_dispatch(rho_L, rho_R, v_max, rho_max, v_limit, flux_type) result(F)
    real,             intent(in) :: rho_L, rho_R, v_max, rho_max, v_limit
    character(len=*), intent(in) :: flux_type
    real :: F
    if (trim(flux_type) == 'newell') then
      F = godunov_newell_flux(rho_L, rho_R, rho_max, v_limit)
    else
      F = godunov_flux(rho_L, rho_R, v_max, rho_max, v_limit)
    end if
  end function godunov_dispatch

end module pde_flux
