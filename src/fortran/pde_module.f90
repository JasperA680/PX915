! pde_module.f90 — LWR PDE solver, N-lane extension of the single-lane solver.
!
! Array convention (lane-major, Fortran column-major):
!   state%density(lane, i)              lane index first, spatial cell second
!   density_history(n_lanes, M, n_steps+1)   Fortran storage: (lane, x, time)
!
! NetCDF storage: dimensions [lane_dimid, x_dimid, time_dimid] (Fortran order).
! Python (netCDF4) reads in reversed order -> shape (time, x, lane).
! load_pde_netcdf in pde_runner.py transposes this to (time, lane, x).
!
! Backward compatibility: n_lanes=1, lane_change_rate=0 reproduces the
! previous single-lane solver bit-for-bit.
module pde_solver
  use pde_flux
  use pde_lanechange
  use netcdf
  implicit none
  private
  public :: pde_params_t, pde_state_t
  public :: pde_initialise, pde_step, pde_finalise
  public :: compute_dt, write_pde_netcdf
  public :: pde_setup_params
  real, parameter :: PI = 3.14159265358979

  type :: pde_params_t
    ! Scalar fields — kept as backward-compatible CLI aliases.
    ! pde_setup_params broadcasts them into the _lanes arrays.
    real    :: dx, dt, domain_length
    real    :: v_max, rho_max
    real    :: rho_left_bc, rho_right_bc, cfl_number
    integer :: M, n_steps, C_checkpoint
    character(len=16) :: ic_type   ! "constant", "riemann", "gaussian", "sine"
    character(len=16) :: bc_type   ! "open", "periodic", "sponge"
    character(len=16) :: flux_type ! "lf" or "godunov"
    logical :: use_adaptive_dt
    integer :: n_sponge
    real    :: sponge_damping

    ! Multi-lane extensions
    integer :: n_lanes                        ! number of lanes (default 1)
    real    :: lane_change_rate               ! k >= 0; 0 disables lane changing
    real, allocatable :: v_max_lanes(:)       ! (n_lanes)
    real, allocatable :: rho_max_lanes(:)     ! (n_lanes)
    real, allocatable :: rho_left_bc_lanes(:)  ! (n_lanes)
    real, allocatable :: rho_right_bc_lanes(:) ! (n_lanes)
  end type pde_params_t

  type :: pde_state_t
    real, allocatable :: density(:,:)  ! (n_lanes, M) physical cells
    real, allocatable :: rho_ext(:,:)  ! (n_lanes, 0:M+1) with ghost cells
    real, allocatable :: flux(:,:)     ! (n_lanes, 0:M) interface fluxes
    real    :: t_current
    integer :: step
    integer :: clip_count  ! cells clamped after lane-change source step
  end type pde_state_t

contains

  ! Broadcast scalar params into per-lane arrays.
  ! Call this after all scalar fields are set and before pde_initialise.
  subroutine pde_setup_params(params)
    type(pde_params_t), intent(inout) :: params
    integer :: n
    n = params%n_lanes

    if (.not. allocated(params%v_max_lanes)) then
      allocate(params%v_max_lanes(n))
      params%v_max_lanes = params%v_max
    end if
    if (.not. allocated(params%rho_max_lanes)) then
      allocate(params%rho_max_lanes(n))
      params%rho_max_lanes = params%rho_max
    end if
    if (.not. allocated(params%rho_left_bc_lanes)) then
      allocate(params%rho_left_bc_lanes(n))
      params%rho_left_bc_lanes = params%rho_left_bc
    end if
    if (.not. allocated(params%rho_right_bc_lanes)) then
      allocate(params%rho_right_bc_lanes(n))
      params%rho_right_bc_lanes = params%rho_right_bc
    end if
  end subroutine pde_setup_params


  subroutine pde_initialise(state, params)
    type(pde_state_t),  intent(out) :: state
    type(pde_params_t), intent(in)  :: params
    integer :: lane, i, M
    real    :: x, x_split, sigma

    M = params%M
    state%t_current = 0.0
    state%step      = 0
    state%clip_count = 0

    allocate(state%density(params%n_lanes, M))
    allocate(state%rho_ext(params%n_lanes, 0:M+1))
    allocate(state%flux(params%n_lanes, 0:M))

    do lane = 1, params%n_lanes
      select case (trim(params%ic_type))

      case ('constant')
        state%density(lane, :) = params%rho_left_bc_lanes(lane)

      case ('riemann')
        x_split = params%domain_length / 2.0
        do i = 1, M
          x = (real(i) - 0.5) * params%dx
          if (x < x_split) then
            state%density(lane, i) = params%rho_left_bc_lanes(lane)
          else
            state%density(lane, i) = params%rho_right_bc_lanes(lane)
          end if
        end do

      case ('gaussian')
        sigma = params%domain_length * 0.05
        do i = 1, M
          x = (real(i) - 0.5) * params%dx
          state%density(lane, i) = params%rho_max_lanes(lane) * 0.2 &
            + 0.6 * params%rho_max_lanes(lane) &
              * exp(-((x - params%domain_length*0.5)**2) / (2.0*sigma**2))
          state%density(lane, i) = min(state%density(lane, i), params%rho_max_lanes(lane))
        end do

      case ('sine')
        do i = 1, M
          x = (real(i) - 0.5) * params%dx
          state%density(lane, i) = params%rho_max_lanes(lane) * 0.5 &
            + 0.15 * params%rho_max_lanes(lane) &
              * sin(2.0 * PI * x / params%domain_length)
        end do

      ! Odd lanes start at rho_left_bc, even lanes at rho_right_bc.
      ! Uses the scalar aliases so rho_left_bc and rho_right_bc set
      ! different uniform densities in alternating lanes.
      case ('staggered')
        if (mod(lane, 2) == 1) then
          state%density(lane, :) = params%rho_left_bc
        else
          state%density(lane, :) = params%rho_right_bc
        end if

      case default
        state%density(lane, :) = params%rho_left_bc_lanes(lane)

      end select
    end do
  end subroutine pde_initialise


  ! Advance one step: per-lane longitudinal FV update, then conservative
  ! lateral lane-change source terms, then optional sponge damping.
  subroutine pde_step(state, params)
    type(pde_state_t),  intent(inout) :: state
    type(pde_params_t), intent(in)    :: params
    integer :: lane, i, M
    real    :: j_real, sigma_j, rho_new
    real, allocatable :: source(:,:)

    M = params%M

    ! --- Longitudinal finite-volume update per lane ---
    do lane = 1, params%n_lanes

      state%rho_ext(lane, 1:M) = state%density(lane, 1:M)

      if (trim(params%bc_type) == 'periodic') then
        state%rho_ext(lane, 0)   = state%density(lane, M)
        state%rho_ext(lane, M+1) = state%density(lane, 1)
      else
        state%rho_ext(lane, 0)   = params%rho_left_bc_lanes(lane)
        state%rho_ext(lane, M+1) = params%rho_right_bc_lanes(lane)
      end if

      ! Godunov (Greenshields or Newell via dispatch); LF falls back to Greenshields.
      if (trim(params%flux_type) == 'godunov' .or. trim(params%flux_type) == 'newell') then
        do i = 0, M
          state%flux(lane, i) = godunov_dispatch( &
            state%rho_ext(lane, i), state%rho_ext(lane, i+1), &
            params%v_max_lanes(lane), params%rho_max_lanes(lane), params%flux_type)
        end do
      else
        do i = 0, M
          state%flux(lane, i) = lax_friedrichs_flux( &
            state%rho_ext(lane, i), state%rho_ext(lane, i+1), &
            params%v_max_lanes(lane), params%rho_max_lanes(lane), &
            params%dx, params%dt)
        end do
      end if

      do i = 1, M
        state%density(lane, i) = state%density(lane, i) &
          - (params%dt / params%dx) * (state%flux(lane, i) - state%flux(lane, i-1))
      end do

    end do

    ! --- Conservative lane-change source terms ---
    ! Stability: large k can violate source-ODE stability; rule of thumb k*dt < 1.
    if (params%lane_change_rate > 0.0 .and. params%n_lanes > 1) then
      allocate(source(params%n_lanes, M))
      call compute_lane_change_sources( &
        state%density, params%n_lanes, M, &
        params%v_max_lanes, params%rho_max_lanes, &
        params%lane_change_rate, source)
      do lane = 1, params%n_lanes
        do i = 1, M
          state%density(lane, i) = state%density(lane, i) + params%dt * source(lane, i)
        end do
      end do
      ! Positivity / bounds clipping with warning counter
      do lane = 1, params%n_lanes
        do i = 1, M
          rho_new = max(0.0, min(params%rho_max_lanes(lane), state%density(lane, i)))
          if (abs(rho_new - state%density(lane, i)) > 1.0e-10) &
            state%clip_count = state%clip_count + 1
          state%density(lane, i) = rho_new
        end do
      end do
      deallocate(source)
    end if

    ! --- Sponge boundary layer (right side) per lane ---
    if (trim(params%bc_type) == 'sponge') then
      do lane = 1, params%n_lanes
        do i = 1, params%n_sponge
          j_real  = real(i) / real(params%n_sponge)
          sigma_j = 0.5 * (1.0 - cos(PI * j_real))
          state%density(lane, M - params%n_sponge + i) = &
            max(0.0, min(params%rho_max_lanes(lane), &
              state%density(lane, M - params%n_sponge + i) &
              * (1.0 - sigma_j * params%sponge_damping * params%dt)))
        end do
      end do
    end if

    state%t_current = state%t_current + params%dt
    state%step      = state%step + 1
  end subroutine pde_step


  subroutine pde_finalise(state)
    type(pde_state_t), intent(inout) :: state
    if (allocated(state%density)) deallocate(state%density)
    if (allocated(state%rho_ext)) deallocate(state%rho_ext)
    if (allocated(state%flux))    deallocate(state%flux)
  end subroutine pde_finalise


  ! CFL-based adaptive dt — takes the maximum wave speed over all lanes.
  ! For Greenshields the conservative bound per lane is dx*CFL/v_max(lane).
  subroutine compute_dt(state, params, dt_out)
    type(pde_state_t),  intent(in)  :: state
    type(pde_params_t), intent(in)  :: params
    real,               intent(out) :: dt_out
    real    :: max_speed, v_max_global, dt_conservative, lane_speed
    integer :: lane

    max_speed    = 0.0
    v_max_global = 0.0
    do lane = 1, params%n_lanes
      if (trim(params%flux_type) == 'newell') then
        lane_speed = maxval(abs(dq_drho_newell(state%density(lane,:), &
                                 params%v_max_lanes(lane), params%rho_max_lanes(lane))))
      else
        lane_speed = maxval(abs(dq_drho(state%density(lane,:), &
                                 params%v_max_lanes(lane), params%rho_max_lanes(lane))))
      end if
      if (lane_speed > max_speed) max_speed = lane_speed
      if (params%v_max_lanes(lane) > v_max_global) v_max_global = params%v_max_lanes(lane)
    end do

    dt_conservative = params%cfl_number * params%dx / v_max_global
    if (max_speed < 1.0e-10) then
      dt_out = dt_conservative
    else
      dt_out = min(params%cfl_number * params%dx / max_speed, dt_conservative)
    end if
  end subroutine compute_dt


  ! Write simulation output to NetCDF.
  !
  ! density_history(n_lanes, M, n_steps+1): Fortran lane-major storage.
  ! NetCDF dims [lane_dimid, x_dimid, time_dimid] -> Python reads (time, x, lane).
  ! pde_runner.load_pde_netcdf transposes that to (time, lane, x).
  !
  ! flow_history(n_lanes, n_steps+1): right-boundary flow per lane per step.
  ! NetCDF dims [lane_dimid, time_dimid] -> Python reads (time, lane).
  ! flow_total is the lane sum, stored as (time,) for backward compat.
  subroutine write_pde_netcdf(filename, params, density_history, flow_history)
    character(len=*),   intent(in) :: filename
    type(pde_params_t), intent(in) :: params
    real,               intent(in) :: density_history(:,:,:)  ! (n_lanes, M, n_steps+1)
    real,               intent(in) :: flow_history(:,:)       ! (n_lanes, n_steps+1)

    integer :: ncid, time_dimid, x_dimid, lane_dimid
    integer :: density_varid, flow_varid, flow_total_varid
    integer :: x_varid, time_varid, lane_varid
    integer :: i, n_times, M, n_lanes_out
    real, allocatable :: x_coord(:), time_coord(:), lane_coord(:), flow_total(:)

    n_lanes_out = size(density_history, 1)
    M           = size(density_history, 2)
    n_times     = size(density_history, 3)

    allocate(x_coord(M), time_coord(n_times), lane_coord(n_lanes_out), flow_total(n_times))

    do i = 1, M
      x_coord(i) = (real(i) - 0.5) * params%dx
    end do
    do i = 1, n_times
      time_coord(i) = real(i - 1) * params%dt
    end do
    do i = 1, n_lanes_out
      lane_coord(i) = real(i)
    end do
    do i = 1, n_times
      flow_total(i) = sum(flow_history(:, i))
    end do

    call nc_check(nf90_create(filename, NF90_CLOBBER, ncid))

    call nc_check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, time_dimid))
    call nc_check(nf90_def_dim(ncid, 'x',    M,              x_dimid))
    call nc_check(nf90_def_dim(ncid, 'lane', n_lanes_out,    lane_dimid))

    call nc_check(nf90_def_var(ncid, 'time',       NF90_FLOAT, [time_dimid],                      time_varid))
    call nc_check(nf90_def_var(ncid, 'x',          NF90_FLOAT, [x_dimid],                         x_varid))
    call nc_check(nf90_def_var(ncid, 'lane',       NF90_FLOAT, [lane_dimid],                      lane_varid))
    call nc_check(nf90_def_var(ncid, 'density',    NF90_FLOAT, [lane_dimid, x_dimid, time_dimid], density_varid))
    call nc_check(nf90_def_var(ncid, 'flow',       NF90_FLOAT, [lane_dimid, time_dimid],          flow_varid))
    call nc_check(nf90_def_var(ncid, 'flow_total', NF90_FLOAT, [time_dimid],                      flow_total_varid))

    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'model',            'LWR-Greenshields'))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'M',                params%M))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'n_lanes',          n_lanes_out))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'domain_length',    params%domain_length))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'dx',               params%dx))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'dt',               params%dt))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'v_max',            params%v_max))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'rho_max',          params%rho_max))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'rho_left_bc',      params%rho_left_bc))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'rho_right_bc',     params%rho_right_bc))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'ic_type',          trim(params%ic_type)))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'bc_type',          trim(params%bc_type)))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'flux_type',        trim(params%flux_type)))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'n_steps',          params%n_steps))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'n_sponge',         params%n_sponge))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'sponge_damping',   params%sponge_damping))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'lane_change_rate', params%lane_change_rate))
    ! Per-lane arrays stored as vector attributes so runs can be reproduced exactly
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'v_max_lanes',   params%v_max_lanes))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'rho_max_lanes', params%rho_max_lanes))

    call nc_check(nf90_enddef(ncid))

    call nc_check(nf90_put_var(ncid, x_varid,         x_coord))
    call nc_check(nf90_put_var(ncid, time_varid,       time_coord))
    call nc_check(nf90_put_var(ncid, lane_varid,       lane_coord))
    call nc_check(nf90_put_var(ncid, density_varid,    density_history))
    call nc_check(nf90_put_var(ncid, flow_varid,       flow_history))
    call nc_check(nf90_put_var(ncid, flow_total_varid, flow_total))

    call nc_check(nf90_close(ncid))
    deallocate(x_coord, time_coord, lane_coord, flow_total)
  end subroutine write_pde_netcdf


  subroutine nc_check(status)
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
      write(*, '(A,A)') 'NetCDF error: ', trim(nf90_strerror(status))
      stop 1
    end if
  end subroutine nc_check

end module pde_solver
