module pde_flux
  implicit none
  private
  public :: v_of_rho, q_of_rho, dq_drho, rho_critical
  public :: q_newell, dq_drho_newell, rho_critical_newell
  public :: q_dispatch, dq_drho_dispatch, godunov_dispatch
  public :: lax_friedrichs_flux, godunov_flux

  real, parameter :: NEWELL_W = 0.5  ! congestion wave speed (default)

contains

  elemental function v_of_rho(rho, v_max, rho_max) result(v)
    real, intent(in) :: rho, v_max, rho_max
    real :: v
    v = v_max * (1.0 - rho / rho_max)
  end function v_of_rho

  elemental function q_of_rho(rho, v_max, rho_max) result(q)
    real, intent(in) :: rho, v_max, rho_max
    real :: q
    q = v_max * rho * (1.0 - rho / rho_max)
  end function q_of_rho

  ! Characteristic speed dq/dρ = v_max*(1 - 2ρ/ρ_max).
  ! Positive in free-flow (ρ < ρ_c), zero at ρ_c, negative in congested (ρ > ρ_c).
  elemental function dq_drho(rho, v_max, rho_max) result(dq)
    real, intent(in) :: rho, v_max, rho_max
    real :: dq
    dq = v_max * (1.0 - 2.0 * rho / rho_max)
  end function dq_drho

  ! Argmax of q for the Greenshields diagram: ρ_c = ρ_max / 2.
  elemental function rho_critical(rho_max) result(rc)
    real, intent(in) :: rho_max
    real :: rc
    rc = rho_max / 2.0
  end function rho_critical

  ! Lax-Friedrichs numerical flux — simple but diffusive, good for debugging.
  ! F_LF = [q(ρ_L) + q(ρ_R)] / 2  −  (Δx / 2Δt) * (ρ_R − ρ_L)
  function lax_friedrichs_flux(rho_L, rho_R, v_max, rho_max, dx, dt) result(F)
    real, intent(in) :: rho_L, rho_R, v_max, rho_max, dx, dt
    real :: F
    F = 0.5 * (q_of_rho(rho_L, v_max, rho_max) + q_of_rho(rho_R, v_max, rho_max)) &
        - (dx / (2.0 * dt)) * (rho_R - rho_L)
  end function lax_friedrichs_flux

  ! Godunov numerical flux for a concave flux function.
  !
  ! F_G(ρ_L, ρ_R) = min_{[ρ_L,ρ_R]} q   if ρ_L ≤ ρ_R  (shock)
  !                 max_{[ρ_R,ρ_L]} q   if ρ_L >  ρ_R  (rarefaction)
  !
  ! Closed-form cases for Greenshields (ρ_c = ρ_max/2):
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
  function godunov_flux(rho_L, rho_R, v_max, rho_max) result(F)
    real, intent(in) :: rho_L, rho_R, v_max, rho_max
    real :: F, rc, qL, qR
    rc = rho_critical(rho_max)
    qL = q_of_rho(rho_L, v_max, rho_max)
    qR = q_of_rho(rho_R, v_max, rho_max)
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
        F = q_of_rho(rc, v_max, rho_max)
      end if
    end if
  end function godunov_flux

  ! ---------------------------------------------------------------
  ! Newell-Daganzo (triangular) fundamental diagram
  !
  !   q(rho) = min( v_f * rho,  w * (rho_max - rho) )
  !
  ! Two regimes meeting at the critical density:
  !   rho_c = w * rho_max / (v_f + w)
  ! Below rho_c: free-flow,    q grows linearly with slope v_f
  ! Above rho_c: congestion,   q falls linearly with slope -w
  !
  ! Note: v_f reuses the v_max parameter; w is hardcoded for now via
  ! NEWELL_W. To make w runtime-configurable, add it to pde_params_t
  ! and thread through the dispatchers.
  ! ---------------------------------------------------------------

  elemental function q_newell(rho, v_max, rho_max) result(q)
    real, intent(in) :: rho, v_max, rho_max
    real :: q
    q = min(v_max * rho, NEWELL_W * (rho_max - rho))
  end function q_newell

  ! Characteristic speed for Newell flux (piecewise constant):
  !   +v_max in free-flow,  -NEWELL_W in congested
  ! Returns 0 exactly at the kink to avoid CFL instability.
  elemental function dq_drho_newell(rho, v_max, rho_max) result(dq)
    real, intent(in) :: rho, v_max, rho_max
    real :: dq, rc
    rc = rho_critical_newell(v_max, rho_max)
    if (rho < rc) then
      dq = v_max
    else if (rho > rc) then
      dq = -NEWELL_W
    else
      dq = 0.0
    end if
  end function dq_drho_newell

  elemental function rho_critical_newell(v_max, rho_max) result(rc)
    real, intent(in) :: v_max, rho_max
    real :: rc
    rc = NEWELL_W * rho_max / (v_max + NEWELL_W)
  end function rho_critical_newell

  ! Godunov flux for Newell-Daganzo (concave, piecewise linear).
  ! Same case structure as the Greenshields version since q is concave.
  function godunov_newell_flux(rho_L, rho_R, v_max, rho_max) result(F)
    real, intent(in) :: rho_L, rho_R, v_max, rho_max
    real :: F, rc, qL, qR
    rc = rho_critical_newell(v_max, rho_max)
    qL = q_newell(rho_L, v_max, rho_max)
    qR = q_newell(rho_R, v_max, rho_max)
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
        F = q_newell(rc, v_max, rho_max)
      end if
    end if
  end function godunov_newell_flux

  ! ---------------------------------------------------------------
  ! Dispatchers — select the closure based on flux_type string
  ! ---------------------------------------------------------------
  function q_dispatch(rho, v_max, rho_max, flux_type) result(q)
    real,             intent(in) :: rho, v_max, rho_max
    character(len=*), intent(in) :: flux_type
    real :: q
    if (trim(flux_type) == 'newell') then
      q = q_newell(rho, v_max, rho_max)
    else
      q = q_of_rho(rho, v_max, rho_max)
    end if
  end function q_dispatch

  function dq_drho_dispatch(rho, v_max, rho_max, flux_type) result(dq)
    real,             intent(in) :: rho, v_max, rho_max
    character(len=*), intent(in) :: flux_type
    real :: dq
    if (trim(flux_type) == 'newell') then
      dq = dq_drho_newell(rho, v_max, rho_max)
    else
      dq = dq_drho(rho, v_max, rho_max)
    end if
  end function dq_drho_dispatch

  function godunov_dispatch(rho_L, rho_R, v_max, rho_max, flux_type) result(F)
    real,             intent(in) :: rho_L, rho_R, v_max, rho_max
    character(len=*), intent(in) :: flux_type
    real :: F
    if (trim(flux_type) == 'newell') then
      F = godunov_newell_flux(rho_L, rho_R, v_max, rho_max)
    else
      F = godunov_flux(rho_L, rho_R, v_max, rho_max)
    end if
  end function godunov_dispatch

end module pde_flux
