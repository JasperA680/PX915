module pde_flux
  ! Numerical flux functions and traffic flow closure relations for the PDE solver.
  !
  ! This module defines velocity, flux and characteristic speed functions used by
  ! the macroscopic traffic flow PDE solver. The default closure is the
  ! Greenshields fundamental diagram, where the traffic velocity decreases
  ! linearly with density and the flux is a concave quadratic function.
  !
  ! All Greenshields routines accept an optional speed-limit argument v_limit.
  ! Setting v_limit = v_max recovers the classical diagram. For v_limit < v_max
  ! the free-flow branch is capped, introducing a kink at rho* = rho_max*(1 -
  ! v_limit/v_max) while preserving concavity.
  !
  ! The module also provides Newell-Daganzo triangular flux routines (where
  ! v_limit replaces v_max in the free-flow branch), numerical flux functions,
  ! and dispatch routines for selecting a closure at runtime.

  implicit none
  private
  public :: v_of_rho, q_of_rho, dq_drho, rho_critical
  public :: q_newell, dq_drho_newell, rho_critical_newell
  public :: q_dispatch, dq_drho_dispatch, godunov_dispatch
  public :: lax_friedrichs_flux, godunov_flux

  real, parameter :: NEWELL_W = 0.5  ! congestion wave speed used by the Newell-Daganzo closure

contains

  elemental function v_of_rho(rho, v_max, rho_max, v_limit) result(v)
    ! Return the Greenshields velocity at a given density, capped by a speed limit.
    !
    ! The base Greenshields velocity decreases linearly with density:
    !
    !   v_gs(rho) = v_max * (1 - rho / rho_max).
    !
    ! The returned velocity is min(v_gs, v_limit), so free-flow speeds are
    ! capped at v_limit. Setting v_limit = v_max recovers the classical diagram.
    real, intent(in) :: rho      ! Traffic density.
    real, intent(in) :: v_max    ! Maximum/free-flow velocity (uncapped).
    real, intent(in) :: rho_max  ! Maximum/jam density.
    real, intent(in) :: v_limit  ! Speed limit; set equal to v_max for no restriction.
    real :: v                    ! Velocity corresponding to rho.

    v = min(v_max * (1.0 - rho / rho_max), v_limit)

  end function v_of_rho


  elemental function q_of_rho(rho, v_max, rho_max, v_limit) result(q)
    ! Return the Greenshields traffic flux at a given density.
    !
    ! The flux is density times the speed-limited velocity:
    !
    !   q(rho) = rho * v(rho, v_limit).
    !
    ! For v_limit = v_max this reduces to the classical quadratic:
    !
    !   q(rho) = v_max * rho * (1 - rho / rho_max).
    !
    ! The diagram remains concave in both cases; the maximum is still at
    ! rho_c = rho_max / 2.
    real, intent(in) :: rho      ! Traffic density.
    real, intent(in) :: v_max    ! Maximum/free-flow velocity (uncapped).
    real, intent(in) :: rho_max  ! Maximum/jam density.
    real, intent(in) :: v_limit  ! Speed limit; set equal to v_max for no restriction.
    real :: q                    ! Traffic flux corresponding to rho.

    q = v_of_rho(rho, v_max, rho_max, v_limit) * rho
  end function q_of_rho


  elemental function dq_drho(rho, v_max, rho_max, v_limit) result(dq)
    ! Return the characteristic speed for the Greenshields flux.
    !
    ! The characteristic speed is the derivative of the flux with respect to
    ! density. With a speed limit it is piecewise:
    !
    !   dq/drho = v_limit                        for rho < rho*
    !           = v_max * (1 - 2*rho/rho_max)    for rho >= rho*
    !
    ! where rho* = rho_max * (1 - v_limit/v_max) is the density at which the
    ! speed limit first bites. For v_limit = v_max, rho* = 0 and the expression
    ! reduces to the classical v_max * (1 - 2*rho/rho_max).
    real, intent(in) :: rho      ! Traffic density.
    real, intent(in) :: v_max    ! Maximum/free-flow velocity (uncapped).
    real, intent(in) :: rho_max  ! Maximum/jam density.
    real, intent(in) :: v_limit  ! Speed limit; set equal to v_max for no restriction.
    real :: dq                   ! Derivative dq/drho at rho.
    real :: rho_star             ! Density at which the speed limit activates.

    rho_star = rho_max * (1.0 - v_limit / v_max)
    if (rho < rho_star) then
      dq = v_limit
    else
      dq = v_max * (1.0 - 2.0 * rho / rho_max)
    end if
  end function dq_drho


  elemental function rho_critical(rho_max) result(rc)
    ! Return the critical density for the Greenshields flux.
    !
    ! The critical density is the density at which the Greenshields flux is
    ! maximised. For the quadratic flux (with or without speed limit), this is
    !
    !   rho_c = rho_max / 2.
    real, intent(in) :: rho_max  ! Maximum/jam density.
    real :: rc                   ! Critical density where q(rho) is maximised.

    rc = rho_max / 2.0
  end function rho_critical


  function lax_friedrichs_flux(rho_L, rho_R, v_max, rho_max, v_limit, dx, dt) result(F)
    ! Return the Lax-Friedrichs numerical flux for the Greenshields closure.
    !
    ! The Lax-Friedrichs flux averages the physical fluxes on the left and
    ! right states and adds numerical diffusion:
    !
    !   F_LF = 0.5 * (q(rho_L) + q(rho_R)) - (dx / (2*dt)) * (rho_R - rho_L)
    !
    ! This flux is robust and useful for debugging, but more diffusive than
    ! the Godunov flux. The speed-limited q is used when v_limit < v_max.
    real, intent(in) :: rho_L    ! Density on the left side of the cell interface.
    real, intent(in) :: rho_R    ! Density on the right side of the cell interface.
    real, intent(in) :: v_max    ! Maximum/free-flow velocity (uncapped).
    real, intent(in) :: rho_max  ! Maximum/jam density.
    real, intent(in) :: v_limit  ! Speed limit; set equal to v_max for no restriction.
    real, intent(in) :: dx       ! Spatial grid spacing.
    real, intent(in) :: dt       ! Time step size.
    real :: F                    ! Lax-Friedrichs numerical flux.

    F = 0.5 * (q_of_rho(rho_L, v_max, rho_max, v_limit) &
              + q_of_rho(rho_R, v_max, rho_max, v_limit)) &
        - (dx / (2.0 * dt)) * (rho_R - rho_L)
  end function lax_friedrichs_flux


  function godunov_flux(rho_L, rho_R, v_max, rho_max, v_limit) result(F)
    ! Return the Godunov numerical flux for the Greenshields closure.
    !
    ! The Godunov flux solves the local Riemann problem at a cell interface.
    ! The Greenshields diagram (with or without speed limit) is concave, so
    ! the argmax is always rho_c = rho_max / 2 and the closed-form cases are:
    !
    ! rho_L <= rho_R (shock):
    !   both sub-critical:   F = q(rho_L)
    !   both super-critical: F = q(rho_R)
    !   straddle rho_c:      F = min(q(rho_L), q(rho_R))
    !
    ! rho_L > rho_R (rarefaction):
    !   both sub-critical:   F = q(rho_L)
    !   both super-critical: F = q(rho_R)
    !   sonic (rho_R < rho_c < rho_L): F = q(rho_c)
    real, intent(in) :: rho_L    ! Density on the left side of the cell interface.
    real, intent(in) :: rho_R    ! Density on the right side of the cell interface.
    real, intent(in) :: v_max    ! Maximum/free-flow velocity (uncapped).
    real, intent(in) :: rho_max  ! Maximum/jam density.
    real, intent(in) :: v_limit  ! Speed limit; set equal to v_max for no restriction.
    real :: F                    ! Godunov numerical flux.
    real :: rc                   ! Critical density for the Greenshields flux.
    real :: qL                   ! Physical flux q(rho_L).
    real :: qR                   ! Physical flux q(rho_R).

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
  ! Newell-Daganzo triangular fundamental diagram
  ! ---------------------------------------------------------------

  elemental function q_newell(rho, rho_max, v_limit) result(q)
    ! Return the Newell-Daganzo triangular traffic flux.
    !
    ! The Newell-Daganzo closure uses a piecewise-linear fundamental diagram:
    !
    !   q(rho) = min(v_limit * rho, NEWELL_W * (rho_max - rho)).
    !
    ! The free-flow branch has slope v_limit (the speed limit). The congested
    ! branch has backward wave speed NEWELL_W. Setting v_limit = v_max
    ! recovers the classical Newell diagram.
    real, intent(in) :: rho      ! Traffic density.
    real, intent(in) :: rho_max  ! Maximum/jam density.
    real, intent(in) :: v_limit  ! Speed limit; slope of the free-flow branch.
    real :: q                    ! Newell-Daganzo traffic flux corresponding to rho.

    q = min(v_limit * rho, NEWELL_W * (rho_max - rho))
  end function q_newell


  elemental function dq_drho_newell(rho, rho_max, v_limit) result(dq)
    ! Return the characteristic speed for the Newell-Daganzo flux.
    !
    ! The derivative is piecewise constant: +v_limit in the free-flow branch
    ! and -NEWELL_W in the congested branch. Returns zero at the critical
    ! density where the triangular flux has a kink.
    real, intent(in) :: rho      ! Traffic density.
    real, intent(in) :: rho_max  ! Maximum/jam density.
    real, intent(in) :: v_limit  ! Speed limit; slope of the free-flow branch.
    real :: dq                   ! Derivative dq/drho at rho.
    real :: rc                   ! Critical density for the Newell-Daganzo flux.

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
    ! Return the critical density for the Newell-Daganzo flux.
    !
    ! The critical density is the density where the free-flow and congested
    ! branches of the triangular diagram meet:
    !
    !   rho_c = NEWELL_W * rho_max / (v_limit + NEWELL_W).
    real, intent(in) :: rho_max  ! Maximum/jam density.
    real, intent(in) :: v_limit  ! Speed limit; slope of the free-flow branch.
    real :: rc                   ! Critical density where the triangular flux is maximised.

    rc = NEWELL_W * rho_max / (v_limit + NEWELL_W)
  end function rho_critical_newell


  function godunov_newell_flux(rho_L, rho_R, rho_max, v_limit) result(F)
    ! Return the Godunov numerical flux for the Newell-Daganzo closure.
    !
    ! The Newell-Daganzo triangular fundamental diagram is concave, so this
    ! routine uses the same shock/rarefaction case structure as the
    ! Greenshields Godunov flux. The physical flux and critical density are
    ! evaluated using the Newell-Daganzo closure.
    real, intent(in) :: rho_L    ! Density on the left side of the cell interface.
    real, intent(in) :: rho_R    ! Density on the right side of the cell interface.
    real, intent(in) :: rho_max  ! Maximum/jam density.
    real, intent(in) :: v_limit  ! Speed limit; slope of the free-flow branch.
    real :: F                    ! Godunov numerical flux.
    real :: rc                   ! Critical density for the Newell-Daganzo flux.
    real :: qL                   ! Physical flux q(rho_L).
    real :: qR                   ! Physical flux q(rho_R).

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
  ! Dispatch routines
  ! ---------------------------------------------------------------

  function q_dispatch(rho, v_max, rho_max, v_limit, flux_type) result(q)
    ! Return the physical traffic flux selected by a flux-type string.
    !
    ! If flux_type contains 'newell' (matching 'newell' or 'newell_lf'), this
    ! routine evaluates the Newell-Daganzo triangular flux. Otherwise it falls
    ! back to the speed-limited Greenshields quadratic flux.
    real,             intent(in) :: rho       ! Traffic density.
    real,             intent(in) :: v_max     ! Maximum/free-flow velocity (uncapped).
    real,             intent(in) :: rho_max   ! Maximum/jam density.
    real,             intent(in) :: v_limit   ! Speed limit.
    character(len=*), intent(in) :: flux_type ! Closure name: 'newell', 'newell_lf', or Greenshields.
    real :: q                                 ! Selected physical traffic flux.

    if (index(trim(flux_type), 'newell') > 0) then
      q = q_newell(rho, rho_max, v_limit)
    else
      q = q_of_rho(rho, v_max, rho_max, v_limit)
    end if
  end function q_dispatch


  function dq_drho_dispatch(rho, v_max, rho_max, v_limit, flux_type) result(dq)
    ! Return the characteristic speed selected by a flux-type string.
    !
    ! If flux_type contains 'newell', this routine evaluates the derivative of
    ! the Newell-Daganzo triangular flux. Otherwise it falls back to the
    ! derivative of the speed-limited Greenshields quadratic flux.
    real,             intent(in) :: rho       ! Traffic density.
    real,             intent(in) :: v_max     ! Maximum/free-flow velocity (uncapped).
    real,             intent(in) :: rho_max   ! Maximum/jam density.
    real,             intent(in) :: v_limit   ! Speed limit.
    character(len=*), intent(in) :: flux_type ! Closure name: 'newell', 'newell_lf', or Greenshields.
    real :: dq                                ! Selected derivative dq/drho.

    if (index(trim(flux_type), 'newell') > 0) then
      dq = dq_drho_newell(rho, rho_max, v_limit)
    else
      dq = dq_drho(rho, v_max, rho_max, v_limit)
    end if
  end function dq_drho_dispatch


  function godunov_dispatch(rho_L, rho_R, v_max, rho_max, v_limit, flux_type) result(F)
    ! Return the Godunov numerical flux selected by a flux-type string.
    !
    ! If flux_type is 'newell', this routine evaluates the Godunov flux using
    ! the Newell-Daganzo triangular closure. For any other value it falls back
    ! to the speed-limited Greenshields Godunov flux.
    real,             intent(in) :: rho_L     ! Density on the left side of the cell interface.
    real,             intent(in) :: rho_R     ! Density on the right side of the cell interface.
    real,             intent(in) :: v_max     ! Maximum/free-flow velocity (uncapped).
    real,             intent(in) :: rho_max   ! Maximum/jam density.
    real,             intent(in) :: v_limit   ! Speed limit.
    character(len=*), intent(in) :: flux_type ! Closure name: 'newell' or Greenshields.
    real :: F                                 ! Selected Godunov numerical flux.

    if (trim(flux_type) == 'newell') then
      F = godunov_newell_flux(rho_L, rho_R, rho_max, v_limit)
    else
      F = godunov_flux(rho_L, rho_R, v_max, rho_max, v_limit)
    end if
  end function godunov_dispatch

end module pde_flux
