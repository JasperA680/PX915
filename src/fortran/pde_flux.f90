module pde_flux
  implicit none
  private
  public :: v_of_rho, q_of_rho, dq_drho, rho_critical
  public :: lax_friedrichs_flux, godunov_flux

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

end module pde_flux
