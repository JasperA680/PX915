! PDE solver driver — multi-lane extension with speed limit.
!
! Usage (all positional, backward-compatible):
!   pde_solver [M] [N_STEPS] [V_MAX] [RHO_MAX] [RHO_LEFT] [RHO_RIGHT]
!              [IC_TYPE] [FLUX_TYPE] [BC_TYPE] [V_LIMIT] [OUTPUT]
!              [N_LANES] [LANE_CHANGE_RATE] [V_MAX_LANES] [RHO_MAX_LANES]
!
! Arg 10: V_LIMIT — speed cap ≤ V_MAX (default = V_MAX, no restriction)
! Arg 11: OUTPUT  — NetCDF output path
! Arg 12: N_LANES — number of lanes (default 1)
! Arg 13: LANE_CHANGE_RATE — k ≥ 0 (default 0)
! Arg 14: comma-separated v_max per lane, e.g. "1.0,1.5"  (optional)
! Arg 15: comma-separated rho_max per lane, e.g. "1.0,1.0" (optional)
program pde_driver
  use pde_solver
  use pde_flux, only: q_of_rho, q_dispatch
  implicit none

  type(pde_state_t)  :: state
  type(pde_params_t) :: params

  real, allocatable :: density_history(:,:,:)  ! (n_lanes, M, n_steps+1)
  real, allocatable :: flow_history(:,:)        ! (n_lanes, n_steps+1)

  character(len=256) :: arg, output_file
  integer :: t, lane, n_found

  ! ---------- defaults ----------
  params%M              = 200
  params%n_steps        = 500
  params%v_max          = 1.0
  params%rho_max        = 1.0
  params%rho_left_bc    = 0.1
  params%rho_right_bc   = 0.9
  params%ic_type        = 'riemann'
  params%bc_type        = 'open'
  params%flux_type      = 'lf'
  params%cfl_number     = 0.9
  params%domain_length  = 1.0
  params%C_checkpoint   = 100
  params%use_adaptive_dt = .true.
  params%n_sponge       = 20
  params%sponge_damping = 5.0 * params%v_max / params%domain_length
  params%v_limit        = params%v_max   ! default: no speed restriction
  params%n_lanes        = 1
  params%lane_change_rate = 0.0
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
    call get_command_argument(10, arg); read(arg, *) params%v_limit
  end if
  if (command_argument_count() >= 11) then
    call get_command_argument(11, arg); output_file = trim(arg)
  end if
  if (command_argument_count() >= 12) then
    call get_command_argument(12, arg); read(arg, *) params%n_lanes
  end if
  if (command_argument_count() >= 13) then
    call get_command_argument(13, arg); read(arg, *) params%lane_change_rate
  end if

  ! ---------- derived parameters ----------
  params%dx = params%domain_length / real(params%M)
  params%dt = params%cfl_number * params%dx / params%v_max

  ! Broadcast scalars into per-lane arrays (allocates the arrays)
  call pde_setup_params(params)

  ! ---------- per-lane overrides (args 14–15) ----------
  if (command_argument_count() >= 14) then
    call get_command_argument(14, arg)
    call parse_real_list(trim(arg), params%v_max_lanes, n_found)
    if (n_found /= params%n_lanes) then
      write(*, '(A,I0,A,I0)') 'ERROR: v_max_lanes list length (', n_found, &
        ') must equal n_lanes (', params%n_lanes, ')'
      stop 1
    end if
    params%dt = params%cfl_number * params%dx / maxval(params%v_max_lanes)
  end if

  if (command_argument_count() >= 15) then
    call get_command_argument(15, arg)
    call parse_real_list(trim(arg), params%rho_max_lanes, n_found)
    if (n_found /= params%n_lanes) then
      write(*, '(A,I0,A,I0)') 'ERROR: rho_max_lanes list length (', n_found, &
        ') must equal n_lanes (', params%n_lanes, ')'
      stop 1
    end if
  end if

  write(*, '(A)')       '=== LWR PDE Solver (multi-lane) ======================'
  write(*, '(A,I0)')    'M                = ', params%M
  write(*, '(A,I0)')    'n_steps          = ', params%n_steps
  write(*, '(A,I0)')    'n_lanes          = ', params%n_lanes
  write(*, '(A,F6.3)')  'v_max (lane 1)   = ', params%v_max_lanes(1)
  if (params%n_lanes > 1) then
    do lane = 2, params%n_lanes
      write(*, '(A,I0,A,F6.3)') 'v_max (lane ', lane, ')   = ', params%v_max_lanes(lane)
    end do
  end if
  write(*, '(A,F6.3)')  'v_limit          = ', params%v_limit
  write(*, '(A,F6.3)')  'rho_max          = ', params%rho_max
  write(*, '(A,F6.3)')  'rho_left_bc      = ', params%rho_left_bc
  write(*, '(A,F6.3)')  'rho_right_bc     = ', params%rho_right_bc
  write(*, '(A,F6.3)')  'lane_change_rate = ', params%lane_change_rate
  write(*, '(A,A)')     'ic_type          = ', trim(params%ic_type)
  write(*, '(A,A)')     'flux_type        = ', trim(params%flux_type)
  write(*, '(A,A)')     'bc_type          = ', trim(params%bc_type)
  write(*, '(A,A)')     'output           = ', trim(output_file)
  write(*, '(A)')       '======================================================'

  ! ---------- allocate history ----------
  allocate(density_history(params%n_lanes, params%M, params%n_steps + 1))
  allocate(flow_history(params%n_lanes, params%n_steps + 1))

  ! ---------- initialise ----------
  call pde_initialise(state, params)

  ! Store initial condition (t=0)
  density_history(:, :, 1) = state%density
  do lane = 1, params%n_lanes
    flow_history(lane, 1) = q_dispatch(state%density(lane, params%M), &
                              params%v_max_lanes(lane), params%rho_max_lanes(lane), &
                              params%v_limit, params%flux_type)
  end do

  ! ---------- time loop ----------
  do t = 1, params%n_steps
    if (params%use_adaptive_dt) call compute_dt(state, params, params%dt)
    call pde_step(state, params)
    density_history(:, :, t + 1) = state%density
    do lane = 1, params%n_lanes
      flow_history(lane, t + 1) = q_dispatch(state%density(lane, params%M), &
                                    params%v_max_lanes(lane), params%rho_max_lanes(lane), &
                                    params%v_limit, params%flux_type)
    end do

    if (mod(t, 100) == 0) then
      write(*, '(A,I0,A,F8.4,A,F8.4)') &
        'step ', t, '  t = ', state%t_current, &
        '  mean_rho = ', sum(state%density) / real(params%n_lanes * params%M)
    end if
  end do

  if (state%clip_count > 0) then
    write(*, '(A,I0,A)') 'WARNING: density clamped in ', state%clip_count, &
      ' cell-steps (consider reducing k or dt).'
  end if

  write(*, '(A,A)') 'Writing output to ', trim(output_file)
  call write_pde_netcdf(trim(output_file), params, density_history, flow_history)
  write(*, '(A)') 'Done.'

  call pde_finalise(state)
  deallocate(density_history, flow_history)

contains

  ! Parse a comma-separated string of reals into an array.
  subroutine parse_real_list(str, vals, n_found)
    character(len=*), intent(in)    :: str
    real,             intent(inout) :: vals(:)
    integer,          intent(out)   :: n_found

    integer :: pos, sep, slen
    character(len=64) :: token

    n_found = 0
    pos  = 1
    slen = len_trim(str)

    do while (pos <= slen)
      sep = index(str(pos:slen), ',')
      if (sep == 0) then
        token = str(pos:slen)
        pos   = slen + 1
      else
        token = str(pos : pos + sep - 2)
        pos   = pos + sep
      end if
      if (len_trim(token) > 0) then
        n_found = n_found + 1
        if (n_found <= size(vals)) read(token, *) vals(n_found)
      end if
    end do
  end subroutine parse_real_list

end program pde_driver
