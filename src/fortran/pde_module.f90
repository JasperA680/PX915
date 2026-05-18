! pde_module.f90 -- LWR PDE solver with an N-lane extension.
!
! Array convention:
!   state%density(lane, i) stores density with lane index first and spatial cell second.
!   density_history(n_lanes, M, n_steps+1) stores lane, space and time in Fortran order.
!
! NetCDF convention:
!   The NetCDF dimensions are written as lane, x and time. Python netCDF4 reads this in
!   reversed order, giving shape (time, x, lane). The Python loading routines transpose
!   this to (time, lane, x).
!
! Backward compatibility:
!   Setting n_lanes = 1 and lane_change_rate = 0 gives the single-lane solver behaviour.

module pde_solver
  ! Finite volume PDE solver for single lane and multilane traffic flow.
  !
  ! This module advances a macroscopic traffic density model using a finite volume
  ! discretisation. Each lane is updated independently in the longitudinal direction,
  ! then optional conservative lane change source terms transfer density between
  ! neighbouring lanes.
  !
  ! The solver supports several initial conditions, boundary conditions, numerical 
  ! fluxes and NetCDF output. The single lane case is recovered by using one lane and
  ! disabling the lane change source term.

  use pde_flux
  use pde_lanechange
  use netcdf

  implicit none
  private
  public :: pde_params_t, pde_state_t
  public :: pde_initialise, pde_step, pde_finalise
  public :: compute_dt, write_pde_netcdf
  public :: pde_setup_params

  real, parameter :: PI = 3.14159265358979 ! Mathematical constant pi.

  type :: pde_params_t
    ! Parameters controlling the PDE solver.
    !
    ! This type stores both scalar single lane parameters and lane resolved arrays.
    ! The scalar fields are retained as backward compatible aliases. The routine
    ! pde_setup_params broadcasts scalar values into the lane resolved arrays when
    ! those arrays have not been provided explicitly.

    real    :: dx                ! Spatial grid spacing.
    real    :: dt                ! Time step size.
    real    :: domain_length     ! Physical length of the one dimensional road domain.
    real    :: v_max             ! Default maximum/free flow velocity.
    real    :: rho_max           ! Default maximum/jam density.
    real    :: v_limit           ! Speed cap; set equal to v_max for no speed limit.
    real    :: rho_left_bc       ! Default left boundary density.
    real    :: rho_right_bc      ! Default right boundary density.
    real    :: cfl_number        ! CFL number used for adaptive time step selection.
    
    integer :: M                 ! Number of physical spatial cells.
    integer :: n_steps           ! Number of time steps to run.
    integer :: C_checkpoint      ! Checkpoint/output interval used by the driver.

    character(len=16) :: ic_type   ! Initial condition: "constant", "riemann", "gaussian", "sine" or "staggered".
    character(len=16) :: bc_type   ! Boundary condition: "open", "periodic" or "sponge".
    character(len=16) :: flux_type ! Numerical flux or closure selector: "lf", "godunov" or "newell".

    logical :: use_adaptive_dt ! Whether the driver should update dt using compute_dt
    integer :: n_sponge        ! Number of cells in the right hand sponge layer.
    real    :: sponge_damping  ! Damping strength used by the sponge boundary layer.

    integer :: n_lanes                          ! Number of lanes in the PDE model.
    real    :: lane_change_rate                 ! Lane change rate: zero disables lane changing.
    real, allocatable :: v_max_lanes(:)         ! Lane resolved maximum/free flow velocities, shape (n_lanes).
    real, allocatable :: rho_max_lanes(:)       ! Lane resolved maximum/jam densities, shape (n_lanes).
    real, allocatable :: rho_left_bc_lanes(:)   ! Lane resolved left boundary densities, shape (n_lanes).
    real, allocatable :: rho_right_bc_lanes(:)  ! Lane resolved right boundary densities, shape (n_lanes).
  end type pde_params_t

  type :: pde_state_t
    ! Mutable state of a PDE simulation.
    !
    ! The density array contains only physical cells. The rho_ext array adds ome 
    ! ghost cell on each side and is rebuilt at every step before computing
    ! interface fluxes. The flux array stores numerical fluxes at cell interfaces.

    real, allocatable :: density(:,:)  ! Traffic density in physical cells, shape (n_lanes, M)
    real, allocatable :: rho_ext(:,:)  ! Density including ghost cells, shape (n_lanes, 0:M+1)
    real, allocatable :: flux(:,:)     ! Interface fluxes, shape (n_lanes, 0:M)
    real    :: t_current   ! Current simulation time.
    integer :: step        ! Current time step index
    integer :: clip_count  ! Number of density values clipped after lane change source updates.
  end type pde_state_t

contains

  subroutine pde_setup_params(params)
    ! Broadcast scalar solver parameters into lane resolved arrays.
    !
    ! This routine should be called after the scalar fields of params have been set
    ! and before pde_initialise is called. Any lane resolved arrays that are already
    ! allocated are left unchanged. Any missing arrays are allocated with length
    ! params%n_lanes and filled using the corresponding scalar default value.
    type(pde_params_t), intent(inout) :: params ! Solver parameters to complete in place.
    
    integer :: n ! number of lanes.

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

      ! Godunov (Greenshields or Newell via dispatch); LF (incl. newell_lf) falls through.
      if (trim(params%flux_type) == 'godunov' .or. trim(params%flux_type) == 'newell') then
        do i = 0, M
          state%flux(lane, i) = godunov_dispatch( &
            state%rho_ext(lane, i), state%rho_ext(lane, i+1), &
            params%v_max_lanes(lane), params%rho_max_lanes(lane), &
            params%v_limit, params%flux_type)
        end do
      else
        do i = 0, M
          state%flux(lane, i) = lax_friedrichs_flux( &
            state%rho_ext(lane, i), state%rho_ext(lane, i+1), &
            params%v_max_lanes(lane), params%rho_max_lanes(lane), &
            params%v_limit, params%dx, params%dt)
        end do
      end if

      do i = 1, M
        state%density(lane, i) = state%density(lane, i) &
          - (params%dt / params%dx) * (state%flux(lane, i) - state%flux(lane, i-1))
      end do

    end do

    ! --- Conservative lane-change source terms ---
    if (params%lane_change_rate > 0.0 .and. params%n_lanes > 1) then
      allocate(source(params%n_lanes, M))
      call compute_lane_change_sources( &
        state%density, params%n_lanes, M, &
        params%v_max_lanes, params%rho_max_lanes, &
        params%lane_change_rate, params%v_limit, source)
      do lane = 1, params%n_lanes
        do i = 1, M
          state%density(lane, i) = state%density(lane, i) + params%dt * source(lane, i)
        end do
      end do
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
  ! Uses index() for newell check so 'newell_lf' also routes to Newell dq/dρ.
  subroutine compute_dt(state, params, dt_out)
    type(pde_state_t),  intent(in)  :: state
    type(pde_params_t), intent(in)  :: params
    real,               intent(out) :: dt_out
    real    :: max_speed, v_max_global, dt_conservative, lane_speed
    integer :: lane

    max_speed    = 0.0
    v_max_global = 0.0
    do lane = 1, params%n_lanes
      if (index(trim(params%flux_type), 'newell') > 0) then
        lane_speed = maxval(abs(dq_drho_newell(state%density(lane,:), &
                                 params%rho_max_lanes(lane), params%v_limit)))
      else
        lane_speed = maxval(abs(dq_drho(state%density(lane,:), &
                                 params%v_max_lanes(lane), params%rho_max_lanes(lane), params%v_limit)))
      end if
      if (lane_speed > max_speed) max_speed = lane_speed
      if (params%v_max_lanes(lane) > v_max_global) v_max_global = params%v_max_lanes(lane)
    end do

    ! v_limit ≤ v_max, so min(..., v_limit) can only tighten or equal the bound.
    dt_conservative = params%cfl_number * params%dx / min(v_max_global, params%v_limit)
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
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'v_limit',          params%v_limit))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'rho_left_bc',      params%rho_left_bc))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'rho_right_bc',     params%rho_right_bc))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'ic_type',          trim(params%ic_type)))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'bc_type',          trim(params%bc_type)))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'flux_type',        trim(params%flux_type)))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'n_steps',          params%n_steps))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'n_sponge',         params%n_sponge))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'sponge_damping',   params%sponge_damping))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'lane_change_rate', params%lane_change_rate))
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
