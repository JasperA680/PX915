! PDE solver driver.
! Usage: pde_solver [M] [N_STEPS] [V_MAX] [RHO_MAX] [RHO_LEFT] [RHO_RIGHT] [IC_TYPE] [FLUX_TYPE] [OUTPUT]
! Defaults: 200 500 1.0 1.0 0.1 0.9 riemann lf data/output/pde_simulation.nc
program pde_driver
  use pde_solver
  use pde_flux, only: q_of_rho
  implicit none

  type(pde_state_t)  :: state
  type(pde_params_t) :: params

  real, allocatable :: density_history(:,:)  ! (M, n_steps+1)
  real, allocatable :: flow_history(:)       ! (n_steps+1)

  character(len=256) :: arg, output_file
  integer :: t

  ! ---------- defaults ----------
  params%M             = 200
  params%n_steps       = 500
  params%v_max         = 1.0
  params%rho_max       = 1.0
  params%rho_left_bc   = 0.1
  params%rho_right_bc  = 0.9
  params%ic_type       = 'riemann'
  params%bc_type       = 'open'
  params%flux_type     = 'lf'
  params%cfl_number    = 0.9
  params%domain_length = 1.0
  params%C_checkpoint  = 100
  params%use_adaptive_dt = .true.
  params%n_sponge      = 10
  params%sponge_damping  = 5.0 * params%v_max / params%domain_length
  output_file = 'data/output/pde_simulation.nc'

  ! ---------- command-line overrides ----------
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg); read(arg, *) params%M
  end if
  if (command_argument_count() >= 2) then
    call get_command_argument(2, arg); read(arg, *) params%n_steps
  end if
  if (command_argument_count() >= 3) then
    call get_command_argument(3, arg); read(arg, *) params%v_max
  end if
  if (command_argument_count() >= 4) then
    call get_command_argument(4, arg); read(arg, *) params%rho_max
  end if
  if (command_argument_count() >= 5) then
    call get_command_argument(5, arg); read(arg, *) params%rho_left_bc
  end if
  if (command_argument_count() >= 6) then
    call get_command_argument(6, arg); read(arg, *) params%rho_right_bc
  end if
  if (command_argument_count() >= 7) then
    call get_command_argument(7, arg); params%ic_type = trim(arg)
  end if
  if (command_argument_count() >= 8) then
    call get_command_argument(8, arg); params%flux_type = trim(arg)
  end if
  if (command_argument_count() >= 9) then
    call get_command_argument(9, arg); params%bc_type = trim(arg)
  end if
  if (command_argument_count() >= 10) then
    call get_command_argument(10, arg); output_file = trim(arg)
  end if

  ! ---------- derived parameters ----------
  params%dx = params%domain_length / real(params%M)
  ! Initial dt from CFL (may be updated each step if use_adaptive_dt)
  params%dt = params%cfl_number * params%dx / params%v_max

  write(*, '(A)')       '=== LWR PDE Solver ==================================='
  write(*, '(A,I0)')    'M           = ', params%M
  write(*, '(A,I0)')    'n_steps     = ', params%n_steps
  write(*, '(A,F6.3)')  'v_max       = ', params%v_max
  write(*, '(A,F6.3)')  'rho_max     = ', params%rho_max
  write(*, '(A,F6.3)')  'rho_left_bc = ', params%rho_left_bc
  write(*, '(A,F6.3)')  'rho_right_bc= ', params%rho_right_bc
  write(*, '(A,A)')     'ic_type     = ', trim(params%ic_type)
  write(*, '(A,A)')     'flux_type   = ', trim(params%flux_type)
  write(*, '(A,A)')     'output      = ', trim(output_file)
  write(*, '(A)')       '======================================================'

  ! ---------- allocate history ----------
  allocate(density_history(params%M, params%n_steps + 1))
  allocate(flow_history(params%n_steps + 1))

  ! ---------- initialise ----------
  call pde_initialise(state, params)

  ! Store initial condition (t=0)
  density_history(:, 1) = state%density
  flow_history(1) = q_of_rho(state%density(params%M), params%v_max, params%rho_max)

  ! ---------- time loop ----------
  do t = 1, params%n_steps
    if (params%use_adaptive_dt) call compute_dt(state, params, params%dt)
    call pde_step(state, params)
    density_history(:, t + 1) = state%density
    flow_history(t + 1) = q_of_rho(state%density(params%M), params%v_max, params%rho_max)

    if (mod(t, 100) == 0) then
      write(*, '(A,I0,A,F8.4,A,F8.4)') &
        'step ', t, '  t = ', state%t_current, &
        '  mean_rho = ', sum(state%density) / real(params%M)
    end if
  end do

  write(*, '(A,A)') 'Writing output to ', trim(output_file)
  call write_pde_netcdf(trim(output_file), params, density_history, flow_history)
  write(*, '(A)') 'Done.'

  call pde_finalise(state)
  deallocate(density_history, flow_history)

end program pde_driver
