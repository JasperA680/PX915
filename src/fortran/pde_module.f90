module pde_solver
  use pde_flux
  use netcdf
  implicit none
  private
  public :: pde_params_t, pde_state_t
  public :: pde_initialise, pde_step, pde_finalise
  public :: compute_dt, write_pde_netcdf
  real, parameter :: PI = 3.14159265358979

  type :: pde_params_t
    real    :: dx, dt, domain_length, v_max, rho_max
    real    :: rho_left_bc, rho_right_bc, cfl_number
    integer :: M, n_steps, C_checkpoint
    character(len=16) :: ic_type   ! "constant", "riemann", "gaussian", "sine"
    character(len=16) :: bc_type   ! "open", "periodic", "sponge"
    character(len=16) :: flux_type ! "lf" (Lax-Friedrichs) or "godunov"
    logical :: use_adaptive_dt
    integer :: n_sponge          ! number of sponge cells at right boundary (default 10)
    real    :: sponge_damping    ! damping rate coefficient (default 5*v_max/domain_length)
  end type pde_params_t

  type :: pde_state_t
    real, allocatable :: density(:)  ! physical cells 1:M
    real, allocatable :: rho_ext(:)  ! ghost cells 0:M+1, scratch for pde_step
    real, allocatable :: flux(:)     ! interface fluxes 0:M, scratch for pde_step
    real    :: t_current
    integer :: step
  end type pde_state_t

contains

  subroutine pde_initialise(state, params)
    type(pde_state_t),  intent(out) :: state
    type(pde_params_t), intent(in)  :: params
    integer :: i, M
    real    :: x, x_split, sigma

    M = params%M
    state%t_current = 0.0
    state%step      = 0

    allocate(state%density(1:M))
    allocate(state%rho_ext(0:M+1))
    allocate(state%flux(0:M))

    select case (trim(params%ic_type))

    case ('constant')
      state%density = params%rho_left_bc

    case ('riemann')
      x_split = params%domain_length / 2.0
      do i = 1, M
        x = (real(i) - 0.5) * params%dx
        if (x < x_split) then
          state%density(i) = params%rho_left_bc
        else
          state%density(i) = params%rho_right_bc
        end if
      end do

    case ('gaussian')
      sigma = params%domain_length * 0.05
      do i = 1, M
        x = (real(i) - 0.5) * params%dx
        state%density(i) = params%rho_max * 0.2 &
          + 0.6 * params%rho_max &
            * exp(-((x - params%domain_length * 0.5)**2) / (2.0 * sigma**2))
        state%density(i) = min(state%density(i), params%rho_max)
      end do

    case ('sine')
      do i = 1, M
        x = (real(i) - 0.5) * params%dx
        state%density(i) = params%rho_max * 0.5 &
          + 0.15 * params%rho_max &
            * sin(2.0 * PI * x / params%domain_length)
      end do

    case default
      state%density = params%rho_left_bc

    end select
  end subroutine pde_initialise


  ! Advance one time step using the selected numerical flux.
  subroutine pde_step(state, params)
    type(pde_state_t),  intent(inout) :: state
    type(pde_params_t), intent(in)    :: params
    integer :: i, M
    real    :: j_real, sigma_j

    M = params%M

    ! Copy interior to extended array
    state%rho_ext(1:M) = state%density(1:M)

    ! Ghost cells
    if (trim(params%bc_type) == 'periodic') then
      state%rho_ext(0)   = state%density(M)
      state%rho_ext(M+1) = state%density(1)
    else
      ! Open boundaries: fix ghost cells to boundary densities
      state%rho_ext(0)   = params%rho_left_bc
      state%rho_ext(M+1) = params%rho_right_bc
    end if

    ! Compute interface fluxes  flux(i) = F_{i+1/2}
    if (trim(params%flux_type) == 'godunov') then
      do i = 0, M
        state%flux(i) = godunov_flux( &
          state%rho_ext(i), state%rho_ext(i+1), params%v_max, params%rho_max)
      end do
    else
      ! Default: Lax-Friedrichs
      do i = 0, M
        state%flux(i) = lax_friedrichs_flux( &
          state%rho_ext(i), state%rho_ext(i+1), &
          params%v_max, params%rho_max, params%dx, params%dt)
      end do
    end if

    ! Conservative update: ρ_i^{n+1} = ρ_i^n − (Δt/Δx)·[F_{i+1/2} − F_{i-1/2}]
    do i = 1, M
      state%density(i) = state%density(i) &
        - (params%dt / params%dx) * (state%flux(i) - state%flux(i-1))
    end do

    ! Sponge layer — applied after the conservative update so the interior
    ! remains fully conservative.
    !
    ! The buffer occupies the rightmost n_sponge physical cells:
    !   cell M-n_sponge+1  (interior edge, j=1) — weakest damping
    !   cell M             (right wall,    j=n_sponge) — strongest damping
    !
    ! Cosine taper: sigma(j) = 0.5 * (1 - cos(pi * j / n_sponge))
    !   -> sigma(1)        ~ 0   (nearly no damping at interior edge)
    !   -> sigma(n_sponge) = 1   (full damping at wall)
    !
    ! Density nudge towards zero:
    !   rho_i <- rho_i * (1 - sigma_j * sponge_damping * dt)
    ! clamped to [0, rho_max] to prevent overshoot.
    if (trim(params%bc_type) == 'sponge') then
      do i = 1, params%n_sponge
        j_real  = real(i) / real(params%n_sponge)
        sigma_j = 0.5 * (1.0 - cos(PI * j_real))
        state%density(M - params%n_sponge + i) = &
          max(0.0, min(params%rho_max, &
            state%density(M - params%n_sponge + i) &
            * (1.0 - sigma_j * params%sponge_damping * params%dt)))
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


  ! CFL-based adaptive time step.
  ! dt = cfl * dx / max_i |dq/dρ(ρ_i)|
  ! CFL-based adaptive time step, capped at the fully conservative bound.
  ! Near the sonic point (rho -> rho_c), dq/drho -> 0 and the naive adaptive
  ! dt -> inf, making LF a centred (unstable) scheme.  Capping at
  ! cfl*dx/v_max ensures the LF diffusion term always exceeds wave speed.
  subroutine compute_dt(state, params, dt_out)
    type(pde_state_t),  intent(in)  :: state
    type(pde_params_t), intent(in)  :: params
    real,               intent(out) :: dt_out
    real :: max_speed, dt_conservative
    max_speed      = maxval(abs(dq_drho(state%density, params%v_max, params%rho_max)))
    dt_conservative = params%cfl_number * params%dx / params%v_max
    if (max_speed < 1.0e-10) then
      dt_out = dt_conservative
    else
      dt_out = min(params%cfl_number * params%dx / max_speed, dt_conservative)
    end if
  end subroutine compute_dt


  subroutine write_pde_netcdf(filename, params, density_history, flow_history)
    character(len=*),   intent(in) :: filename
    type(pde_params_t), intent(in) :: params
    real,               intent(in) :: density_history(:,:)  ! (M, n_steps+1)
    real,               intent(in) :: flow_history(:)       ! (n_steps+1)

    integer :: ncid, time_dimid, x_dimid
    integer :: density_varid, flow_varid, x_varid, time_varid
    integer :: i, n_times, M
    real, allocatable :: x_coord(:), time_coord(:)

    M       = size(density_history, 1)
    n_times = size(density_history, 2)  ! n_steps + 1

    allocate(x_coord(M))
    allocate(time_coord(n_times))

    do i = 1, M
      x_coord(i) = (real(i) - 0.5) * params%dx
    end do
    do i = 1, n_times
      time_coord(i) = real(i - 1) * params%dt
    end do

    call nc_check(nf90_create(filename, NF90_CLOBBER, ncid))

    ! Dimensions
    call nc_check(nf90_def_dim(ncid, 'time', NF90_UNLIMITED, time_dimid))
    call nc_check(nf90_def_dim(ncid, 'x',    M,              x_dimid))

    ! Variables: Fortran (x, time) → Python/C reads as (time, x) after transpose
    call nc_check(nf90_def_var(ncid, 'time',    NF90_FLOAT, [time_dimid],             time_varid))
    call nc_check(nf90_def_var(ncid, 'x',       NF90_FLOAT, [x_dimid],               x_varid))
    call nc_check(nf90_def_var(ncid, 'density', NF90_FLOAT, [x_dimid, time_dimid],   density_varid))
    call nc_check(nf90_def_var(ncid, 'flow',    NF90_FLOAT, [time_dimid],            flow_varid))

    ! Global attributes (enough to reproduce the run exactly)
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'model',         'LWR-Greenshields'))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'M',             params%M))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'domain_length', params%domain_length))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'dx',            params%dx))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'dt',            params%dt))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'v_max',         params%v_max))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'rho_max',       params%rho_max))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'rho_left_bc',   params%rho_left_bc))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'rho_right_bc',  params%rho_right_bc))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'ic_type',       trim(params%ic_type)))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'bc_type',       trim(params%bc_type)))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'flux_type',     trim(params%flux_type)))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'n_steps',       params%n_steps))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'n_sponge',       params%n_sponge))
    call nc_check(nf90_put_att(ncid, NF90_GLOBAL, 'sponge_damping', params%sponge_damping))

    call nc_check(nf90_enddef(ncid))

    call nc_check(nf90_put_var(ncid, x_varid,       x_coord))
    call nc_check(nf90_put_var(ncid, time_varid,    time_coord))
    call nc_check(nf90_put_var(ncid, density_varid, density_history))
    call nc_check(nf90_put_var(ncid, flow_varid,    flow_history))

    call nc_check(nf90_close(ncid))
    deallocate(x_coord, time_coord)
  end subroutine write_pde_netcdf


  subroutine nc_check(status)
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
      write(*, '(A,A)') 'NetCDF error: ', trim(nf90_strerror(status))
      stop 1
    end if
  end subroutine nc_check

end module pde_solver
