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
    ! Allocate and initialise the PDE state.
    !
    ! The density field is initialised independently for each lane according to 
    ! params%ic_type. Supported initial conditions are constant, riemann, gaussian,
    ! sine and staggered. The ghost cel and flux arrays are allocated but are 
    ! populated later during pde_step.
    type(pde_state_t),  intent(out) :: state   ! Simulation state to allocate and initialise.
    type(pde_params_t), intent(in)  :: params  ! Solver parameters defining grid, lanes and initial condition.
    
    integer :: lane     ! Lane index.
    integer :: i        ! Spatial cell index.
    integer :: M        ! Number of physical spatial cells.
    real    :: x        ! Cell centre coordinate.
    real    :: x_split  ! Interface location for the Riemann initial condition.
    real    :: sigma    ! Width parameter for the Gaussian initial condition.

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


  subroutine pde_step(state, params)
    ! Advance the PDE state by one time step.
    !
    ! The update has three stages. First, each lane is advanced using a 
    ! longitudinal finite volume update. Second, conservative lane change source
    ! terms are applied when lane changing is enabled. Third, an optional sponge
    ! boundary layer damps densities near the right boundary.
    type(pde_state_t),  intent(inout) :: state    ! Simulation state to update in place.
    type(pde_params_t), intent(in)    :: params   ! Solver parameters for the update.

    integer :: lane     ! Lane index
    integer :: i        ! Spatial cell or interface index.
    integer :: M        ! Number of physical spatial cells.
    real    :: j_real   ! Normalised index inside the sponge layer.
    real    :: sigma_j  ! Sponge-layer damping profile value.
    real    :: rho_new  ! Clipped density after applying source terms.
    real, allocatable :: source(:,:) ! Lane change source terms, shape (n_lanes, M)

    M = params%M

    
    do lane = 1, params%n_lanes

      state%rho_ext(lane, 1:M) = state%density(lane, 1:M)

      if (trim(params%bc_type) == 'periodic') then
        state%rho_ext(lane, 0)   = state%density(lane, M)
        state%rho_ext(lane, M+1) = state%density(lane, 1)
      else
        state%rho_ext(lane, 0)   = params%rho_left_bc_lanes(lane)
        state%rho_ext(lane, M+1) = params%rho_right_bc_lanes(lane)
      end if

      
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
    ! Deallocate arrays owned by a PDE state.
    !
    ! This routine realeses the density, ghost cell and flux arrays if they are
    ! currently allocated. It is safe to call even if some arrays are already 
    ! deallocated.
    type(pde_state_t), intent(inout) :: state
    if (allocated(state%density)) deallocate(state%density)
    if (allocated(state%rho_ext)) deallocate(state%rho_ext)
    if (allocated(state%flux))    deallocate(state%flux)
  end subroutine pde_finalise

  subroutine compute_dt(state, params, dt_out)
    ! Compute an adaptive timestep from a CFL condition.
    !
    ! The time step is chosen from the largest characteristic speed over all lanes.
    ! A conservative upper bound based on the global maximum speed is also applied.
    ! This avoids overly large time steps when the local characteristic speed is 
    ! close to zero.
    type(pde_state_t),  intent(in)  :: state   ! Current simulation state.
    type(pde_params_t), intent(in)  :: params  ! Solver parameters and CFL time step bound.
    real,               intent(out) :: dt_out  ! Computed time step.

    real    :: max_speed          ! Maximum absolute characteristic speed over all lanes.
    real    :: v_max_global       ! Largest lane maximum velocity.
    real    :: dt_conservative    ! Conservative CFL time step bound.
    real    :: lane_speed         ! Maximum absolute characteristic speed in one lane.
    integer :: lane               ! Lane index.

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

    dt_conservative = params%cfl_number * params%dx / min(v_max_global, params%v_limit)

    if (max_speed < 1.0e-10) then
      dt_out = dt_conservative
    else
      dt_out = min(params%cfl_number * params%dx / max_speed, dt_conservative)
    end if
  end subroutine compute_dt


  subroutine write_pde_netcdf(filename, params, density_history, flow_history)
    ! Write PDE simulation output to a NetCDF file.
    !
    ! The density history is expected in Fortran oder with shape 
    ! (n_lanes, M, n_steps+1). The NetCDF variable is written with dimensions 
    ! (lane, x, time). Python netCDF4 reads this as (time, x, lane) and the 
    ! project Python loader transposes it to (time, lane, x).
    !
    ! The routine also writes lane-resolved flow, total flow, coordinate arrays
    ! and global metadata describing the simulaiton parameters.
    character(len=*),   intent(in) :: filename                ! Output NetCDF filenames.
    type(pde_params_t), intent(in) :: params                  ! Solver parameters to store as metadata.
    real,               intent(in) :: density_history(:,:,:)  ! Density history, shape (n_lanes, M, n_steps+1)
    real,               intent(in) :: flow_history(:,:)       ! Lane resolved flow history, shape (n_lanes, n_steps+1)

    integer :: ncid                          ! NetCDF file ID.
    integer :: time_dimid                    ! NetCDF dimension identifier for time.
    integer :: x_dimid                       ! NetCDF dimension identifier for x.
    integer :: lane_dimid                    ! NetCDF dimension identifier for lane.
    integer :: density_varid                 ! NetCDF variable identifier for density.
    integer :: flow_varid                    ! NetCDF variable identifier for lane resolved flow.
    integer :: flow_total_varid              ! NetCDF variable identifier for total flow.
    integer :: x_varid                       ! NetCDF variable identifier for x coordinate.
    integer :: time_varid                    ! NetCDF variable identifier for time coordinate.
    integer :: lane_varid                    ! NetCDF variable identifier for lane coordinate.
    integer :: i                             ! Loop index.
    integer :: n_times                        ! Number of time steps in the history arrays.
    integer :: M                             ! Number of physical spatial cells.
    integer :: n_lanes_out                   ! Number of lanes written to file.
    real, allocatable :: x_coord(:)          ! cell centre x coordinates, shape (M)
    real, allocatable :: time_coord(:)       ! output time coordinates, shape (n_times).
    real, allocatable :: lane_coord(:)       ! lane coordinates, shape (n_lanes_out).
    real, allocatable :: flow_total(:)       ! Total flow history, shape (n_times).

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
    ! Check a NetCDF status code and stop if error occurred.
    !
    ! This internal helper keeps NetCDF error handling concise. If status is not 
    ! nf90_noerr, the corresponding NetCDF error string is printed and execution
    ! stops.
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
      write(*, '(A,A)') 'NetCDF error: ', trim(nf90_strerror(status))
      stop 1
    end if
  end subroutine nc_check

end module pde_solver
