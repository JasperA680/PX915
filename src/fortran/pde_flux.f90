module pde_flux
  ! Numerical flux functions and traffic flow closure relations for the PDE solver.
  !
  ! This module defines velocity, flux and characteristic speed functions used by
  ! the macroscopic traffic flow PDE solver. The default closure is the 
  ! Greenshields fundamental diagram, where the traffic velocity decreases 
  ! linearly with density and the flux is a concave quadratic function.
  !
  ! The module also provides alternative Newell-Daganzo triangular flux routines,
  ! numerical flux functions, and dispatch routines for selecting a closure at 
  ! runtime.

  implicit none
  private
  public :: v_of_rho, q_of_rho, dq_drho, rho_critical
  public :: q_newell, dq_drho_newell, rho_critical_newell
  public :: q_dispatch, dq_drho_dispatch, godunov_dispatch
  public :: lax_friedrichs_flux, godunov_flux

  real, parameter :: NEWELL_W = 0.5  ! congestion wave speed used by the Newell-Daganzo closure

contains

  elemental function v_of_rho(rho, v_max, rho_max) result(v)
    ! Return the Greenshields velocity at a given density.
    !
    ! The Greenshields closure assumes that velocity decreases linearly with
    ! density,
    !
    !   v(rho) = v_max * (1 - rho / rho_max).
    !
    ! This gives free flow velocity v_max at zero density and zero velocity at
    ! the jam density rho_max.
    real, intent(in) :: rho      ! Traffic density
    real, intent(in) :: v_max    ! Maximum/free-flow velocity.
    real, intent(in) :: rho_max  ! Maximum/jam density.
    real :: v                    ! Velocity corresponding to rho.

    v = v_max * (1.0 - rho / rho_max)

  end function v_of_rho

  elemental function q_of_rho(rho, v_max, rho_max) result(q)
    ! Return the Greenshields traffic flux at a given density.
    ! 
    ! The Greenshields flux is defined as density times velocity:
    !
    !   q(rho) = rho * v(rho)
    !           = v_max * rho * (1 - rho / rho_max).
    !
    ! This quadratic fundamental diagram is zero at rho = 0 and rho = rho_max,
    ! with a maximum at the critical density rho_c = rho_max / 2.
    real, intent(in) :: rho      ! Traffic density.
    real, intent(in) :: v_max    ! Maximum/free-flow velocity.
    real, intent(in) :: rho_max  ! Maximum/jam density.
    real :: q                    ! Traffic flux corresponding to rho.

    q = v_max * rho * (1.0 - rho / rho_max)
  end function q_of_rho

  elemental function dq_drho(rho, v_max, rho_max) result(dq)
    ! Return the characteristic speed for the Greenshields flux.
    !
    ! For the scalar conservation law, the characteristic speed is the derivative
    ! of the flux with respect to density:
    !
    !   dq/drho = v_max * (1 - 2*rho/rho_max).
    !
    ! This is positive in free-flow traffic, zero at the critical density, and
    ! negative in congested traffic.
    real, intent(in) :: rho      ! Traffic density.
    real, intent(in) :: v_max    ! Maximum/free-flow velocity.
    real, intent(in) :: rho_max  ! Maximum/jam density.
    real :: dq                   ! Derivative dq/drho at rho.

    dq = v_max * (1.0 - 2.0 * rho / rho_max)
  end function dq_drho

  elemental function rho_critical(rho_max) result(rc)
    ! Return the critical density for the Greenshields flux.
    !
    ! The critical density is the density at which the Greenshields flux is
    ! maximised. For the quadratic flux, this is
    !
    !   rho_c = rho_max / 2.
    real, intent(in) :: rho_max  ! Maximum/jam density.
    real :: rc                   ! Critical density where q(rho) is maximised.

    rc = rho_max / 2.0
  end function rho_critical

  function lax_friedrichs_flux(rho_L, rho_R, v_max, rho_max, dx, dt) result(F)
    ! Return the Lax-Friedrichs numerical flux for the Greenshields closure.
    !
    ! The Lax-Friedrichs flux is a simple approximate numerical flux for the
    ! scalar conservation law. It averages the physical fluxes on the left and
    ! right states and adds numerical diffusion:
    !
    !   F_LF = 0.5 * (q(rho_L) + q(rho_R))
    !          - (dx / (2*dt)) * (rho_R - rho_L)
    !
    ! This flux is robust and useful for debugging, but it is more diffusive
    ! than the Godunov flux.
    real, intent(in) :: rho_L    ! Density on the left side of the cell interface.
    real, intent(in) :: rho_R    ! Density on the right side of the cell interface.
    real, intent(in) :: v_max    ! Maximum/free-flow velocity.
    real, intent(in) :: rho_max  ! Maximum/jam density.
    real, intent(in) :: dx       ! Spatial grid spacing.
    real, intent(in) :: dt       ! Time step size.
    real :: F                    ! Lax-Friedrichs numerical flux.

    F = 0.5 * (q_of_rho(rho_L, v_max, rho_max) + q_of_rho(rho_R, v_max, rho_max)) &
        - (dx / (2.0 * dt)) * (rho_R - rho_L)
  end function lax_friedrichs_flux


  function godunov_flux(rho_L, rho_R, v_max, rho_max) result(F)
    ! Return the Godunov numerical flux for the Greenshields closure.
    !
    ! The Godunov flux solves the local Riemann problem at a cell interface.
    ! For the Greenshields fundamental diagram, the flux is concave, so the
    ! closed-form cases can be written in terms of the critical density
    ! rho_c = rho_max / 2.
    !
    ! If rho_L <= rho_R, the Riemann problem corresponds to a shock. If both
    ! states are sub-critical, the interface flux is q(rho_L). If both states
    ! are super-critical, the interface flux is q(rho_R). If the two states
    ! straddle the critical density, the flux is min(q(rho_L), q(rho_R)).
    !
    ! If rho_L > rho_R, the Riemann problem corresponds to a rarefaction. If
    ! both states are sub-critical, the interface flux is q(rho_L). If both
    ! states are super-critical, the interface flux is q(rho_R). If the
    ! rarefaction crosses the critical density, the flux is q(rho_c).
    real, intent(in) :: rho_L    ! Density on the left side of the cell interface.
    real, intent(in) :: rho_R    ! Density on the right side of the cell interface.
    real, intent(in) :: v_max    ! Maximum/free-flow velocity.
    real, intent(in) :: rho_max  ! Maximum/jam density.
    real :: F                    ! Godunov numerical flux.
    real :: rc                   ! Critical density for the Greenshields flux.
    real :: qL                   ! Physical flux q(rho_L).
    real :: qR                   ! Physical flux q(rho_R).

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
  ! Newell-Daganzo triangular fundamental diagram
  ! ---------------------------------------------------------------

  elemental function q_newell(rho, v_max, rho_max) result(q)
    ! Return the Newell-Daganzo triangular traffic flux.
    !
    ! The Newell-Daganzo closure uses a piecewise-linear fundamental diagram:
    !
    !   q(rho) = min(v_max * rho, NEWELL_W * (rho_max - rho)).
    !
    ! The first branch describes free-flow traffic with slope v_max. The second
    ! branch describes congested traffic with backward wave speed NEWELL_W.
    real, intent(in) :: rho      ! Traffic density.
    real, intent(in) :: v_max    ! Free-flow velocity, used as the slope of the free-flow branch.
    real, intent(in) :: rho_max  ! Maximum/jam density.
    real :: q                    ! Newell-Daganzo traffic flux corresponding to rho.

    q = min(v_max * rho, NEWELL_W * (rho_max - rho))
  end function q_newell


  elemental function dq_drho_newell(rho, v_max, rho_max) result(dq)
    ! Return the characteristic speed for the Newell-Daganzo flux.
    !
    ! The derivative is piecewise constant. It is +v_max in the free-flow
    ! branch and -NEWELL_W in the congested branch. At the critical density,
    ! where the triangular flux has a kink, this routine returns zero.
    real, intent(in) :: rho      ! Traffic density.
    real, intent(in) :: v_max    ! Free-flow velocity.
    real, intent(in) :: rho_max  ! Maximum/jam density.
    real :: dq                   ! Derivative dq/drho at rho.
    real :: rc                   ! Critical density for the Newell-Daganzo flux.

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
    ! Return the critical density for the Newell-Daganzo flux.
    !
    ! The critical density is the density where the free-flow and congested
    ! branches of the triangular fundamental diagram meet:
    !
    !   rho_c = NEWELL_W * rho_max / (v_max + NEWELL_W).
    real, intent(in) :: v_max    ! Free-flow velocity.
    real, intent(in) :: rho_max  ! Maximum/jam density.
    real :: rc                   ! Critical density where the triangular flux is maximised.

    rc = NEWELL_W * rho_max / (v_max + NEWELL_W)
  end function rho_critical_newell


  function godunov_newell_flux(rho_L, rho_R, v_max, rho_max) result(F)
    ! Return the Godunov numerical flux for the Newell-Daganzo closure.
    !
    ! The Newell-Daganzo triangular fundamental diagram is concave, so this
    ! routine uses the same shock/rarefaction case structure as the
    ! Greenshields Godunov flux. The physical flux and critical density are
    ! evaluated using the Newell-Daganzo closure.
    real, intent(in) :: rho_L    ! Density on the left side of the cell interface.
    real, intent(in) :: rho_R    ! Density on the right side of the cell interface.
    real, intent(in) :: v_max    ! Free-flow velocity.
    real, intent(in) :: rho_max  ! Maximum/jam density.
    real :: F                    ! Godunov numerical flux.
    real :: rc                   ! Critical density for the Newell-Daganzo flux.
    real :: qL                   ! Physical flux q(rho_L).
    real :: qR                   ! Physical flux q(rho_R).

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
  ! Dispatch routines
  ! ---------------------------------------------------------------

  function q_dispatch(rho, v_max, rho_max, flux_type) result(q)
    ! Return the physical traffic flux selected by a flux-type string.
    !
    ! If flux_type is 'newell', this routine evaluates the Newell-Daganzo
    ! triangular flux. For any other value, it falls back to the Greenshields
    ! quadratic flux.
    real,             intent(in) :: rho       ! Traffic density.
    real,             intent(in) :: v_max     ! Maximum/free-flow velocity.
    real,             intent(in) :: rho_max   ! Maximum/jam density.
    character(*),     intent(in) :: flux_type ! Flux closure name: 'newell' or default Greenshields.
    real :: q                                 ! Selected physical traffic flux.

    if (trim(flux_type) == 'newell') then
      q = q_newell(rho, v_max, rho_max)
    else
      q = q_of_rho(rho, v_max, rho_max)
    end if
  end function q_dispatch


  function dq_drho_dispatch(rho, v_max, rho_max, flux_type) result(dq)
    ! Return the characteristic speed selected by a flux-type string.
    !
    ! If flux_type is 'newell', this routine evaluates the derivative of the
    ! Newell-Daganzo triangular flux. For any other value, it falls back to the
    ! derivative of the Greenshields quadratic flux.
    real,             intent(in) :: rho       ! Traffic density.
    real,             intent(in) :: v_max     ! Maximum/free-flow velocity.
    real,             intent(in) :: rho_max   ! Maximum/jam density.
    character(*), intent(in) :: flux_type     ! Flux closure name: 'newell' or default Greenshields.
    real :: dq                                ! Selected derivative dq/drho.

    if (trim(flux_type) == 'newell') then
      dq = dq_drho_newell(rho, v_max, rho_max)
    else
      dq = dq_drho(rho, v_max, rho_max)
    end if
  end function dq_drho_dispatch


  function godunov_dispatch(rho_L, rho_R, v_max, rho_max, flux_type) result(F)
    ! Return the Godunov numerical flux selected by a flux-type string.
    !
    ! If flux_type is 'newell', this routine evaluates the Godunov flux using
    ! the Newell-Daganzo triangular closure. For any other value, it falls back
    ! to the Greenshields Godunov flux.
    real,             intent(in) :: rho_L     ! Density on the left side of the cell interface.
    real,             intent(in) :: rho_R     ! Density on the right side of the cell interface.
    real,             intent(in) :: v_max     ! Maximum/free-flow velocity.
    real,             intent(in) :: rho_max   ! Maximum/jam density.
    character(*),     intent(in) :: flux_type ! Flux closure name: 'newell' or default Greenshields.
    real :: F                                 ! Selected Godunov numerical flux.

    if (trim(flux_type) == 'newell') then
      F = godunov_newell_flux(rho_L, rho_R, v_max, rho_max)
    else
      F = godunov_flux(rho_L, rho_R, v_max, rho_max)
    end if
  end function godunov_dispatch

end module pde_flux
