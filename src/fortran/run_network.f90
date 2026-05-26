program run_network
    ! Command-line driver for running a road-network simulation.
    !
    ! This program reads a network configuration file, builds the corresponding
    ! ``road_network_t`` object, runs the requested cellular automaton simulation
    ! and writes the output to a NetCDF file.
    !
    ! Usage from the project root is:
    !
    ! .. code-block:: bash
    !
    !    ./build/run_network <config.nc> <out.nc>
    !
    ! The driver performs the following steps:
    !
    ! 1. read the simulation configuration using ``read_config``;
    ! 2. build a ``road_network_t`` using ``build_network``;
    ! 3. seed the Fortran random number generator using ``params%rng_seed``;
    ! 4. allocate occupancy and velocity history arrays;
    ! 5. advance the network for ``params%n_steps`` time steps;
    ! 6. record cell histories and road-level diagnostics;
    ! 7. write the results to NetCDF using ``write_network_netcdf``.
    !
    ! The main history arrays have shape:
    !
    ! .. math::
    !
    !    (\mathrm{max\_lane\_length},\ \mathrm{total\_lanes},\ n_{\mathrm{steps}}).
    !
    ! Progress is printed every ``max(1, n_steps/100)`` steps. This is useful for
    ! the Python GUI or launcher, which can parse the output stream and update a
    ! progress bar.

    use road_network_mod
    use network_init_mod,         only: free_network
    use network_simulation_mod,   only: network_step
    use network_builder_mod
    use nc_config_mod
    use network_io_mod

    implicit none

    character(len=1024) :: config_path
    ! Path to the input network configuration file.
    character(len=1024) :: out_path
    ! Path to the output NetCDF results file.
    integer :: argc
    ! Number of command-line arguments.

    type(network_spec_t) :: spec
    ! Flat network specification read from the configuration file.
    type(sim_params_t) :: params
    ! Simulation parameters read from the configuration file.
    type(road_network_t) :: net
    ! Runtime road-network object built from ``spec``.

    integer :: n_lanes
    ! Total number of lanes in the network after flattening.
    integer :: max_L
    ! Maximum lane length used as the cell dimension in history arrays.
    integer :: n_roads
    ! Number of roads in the network.
    integer :: n_steps
    ! Number of simulation steps.
    integer :: step
    ! Simulation time-step index.
    integer :: r
    ! Road index.
    integer :: k
    ! Lane index within a road.
    integer :: idx
    ! Flattened lane index.
    integer :: Lk
    ! Length of the current lane.
    integer :: n_occ
    ! Number of occupied cells on the current road.
    integer :: log_every
    ! Interval between progress messages.

    integer(kind=1), allocatable :: occupancy(:,:,:)
    ! Cell occupancy history, indexed as ``occupancy(cell, lane_index, time)``.
    integer(kind=1), allocatable :: velocity(:,:,:)
    ! Cell velocity history, indexed as ``velocity(cell, lane_index, time)``.
    real, allocatable :: road_density(:,:)
    ! Road-level density history, indexed as ``road_density(road, time)``.
    integer, allocatable :: road_entries(:,:)
    ! Road-level entry counts, indexed as ``road_entries(road, time)``.
    integer, allocatable :: road_exits(:,:)
    ! Road-level exit counts, indexed as ``road_exits(road, time)``.
    integer, allocatable :: lane_road_offset(:)
    ! Offset mapping from road index to the first flattened lane index minus one.
    integer, allocatable :: seed_arr(:)
    ! Seed array passed to the Fortran random number generator.
    integer :: seed_size
    ! Required size of the Fortran random seed array.

    integer, parameter :: lane_change_model = -1
    ! Legacy hard-coded lane-change selector. Currently not used because the
    ! runtime value is read from ``params%lc_model``.
    real, parameter :: lane_change_prob = 1
    ! Legacy hard-coded lane-change probability. Currently not used because the
    ! runtime value is read from ``params%lc_p_change``.

    argc = command_argument_count()

    if (argc < 2) then
        write(*,*) "usage: run_network <config.nc> <out.nc>"
        stop 1
    end if

    call get_command_argument(1, config_path)
    call get_command_argument(2, out_path)

    write(*,'(a,a)') "config: ", trim(config_path)
    write(*,'(a,a)') "output: ", trim(out_path)

    ! Read configuration file into a flat network specification and parameter set.
    call read_config(trim(config_path), spec, params)

    n_steps = params%n_steps

    ! Build the runtime road network.
    call build_network(spec, net)
    n_roads = size(net%roads)

    ! Compute total lane count and maximum lane length.
    n_lanes = 0
    max_L   = 0

    do r = 1, n_roads
        do k = 1, size(net%roads(r)%lane)
            n_lanes = n_lanes + 1
            max_L = max(max_L, net%roads(r)%lane(k)%length)
        end do
    end do

    if (params%max_lane_length > 0) max_L = max(max_L, params%max_lane_length)

    ! Build per-road offsets into the flattened lane index.
    allocate(lane_road_offset(n_roads))

    idx = 0
    do r = 1, n_roads
        lane_road_offset(r) = idx
        idx = idx + size(net%roads(r)%lane)
    end do

    ! Seed the Fortran random number generator deterministically.
    call random_seed(size=seed_size)
    allocate(seed_arr(seed_size))

    do k = 1, seed_size
        seed_arr(k) = params%rng_seed + 1013904223 * k
    end do

    call random_seed(put=seed_arr)

    ! Allocate history arrays.
    allocate(occupancy(max_L, n_lanes, n_steps), velocity(max_L, n_lanes, n_steps))
    occupancy = 0_1
    velocity  = 0_1

    allocate(road_density(n_roads, n_steps))
    allocate(road_entries(n_roads, n_steps), road_exits(n_roads, n_steps))

    road_density = 0.0
    road_entries = 0
    road_exits   = 0

    write(*,'(a,i0,a,i0,a,i0)') &
        "network: n_roads=", n_roads, " n_lanes=", n_lanes, " max_L=", max_L

    log_every = max(1, n_steps / 100)

    write(*,'(a,i0,a,f6.3)') &
        "lane change: model=", params%lc_model, " p_change=", params%lc_p_change

    ! Main simulation loop.
    do step = 1, n_steps
        call network_step(net, params%model, params%lc_model, params%lc_p_change, &
                          params%v_max, params%p_slow)

        ! Record state and road-level diagnostics.
        do r = 1, n_roads
            n_occ = 0

            do k = 1, size(net%roads(r)%lane)
                Lk = net%roads(r)%lane(k)%length
                idx = lane_road_offset(r) + k

                ! Occupancy and velocity for this lane. Cells beyond ``Lk`` remain
                ! zero because lanes are padded to ``max_L`` in the output arrays.
                where (net%roads(r)%lane(k)%cells(1:Lk)%has_car)
                    occupancy(1:Lk, idx, step) = 1_1
                end where

                velocity(1:Lk, idx, step) = int(net%roads(r)%lane(k)%cells(1:Lk)%velocity, kind=1)

                n_occ = n_occ + count(net%roads(r)%lane(k)%cells(1:Lk)%has_car)

                ! Count entries through open inflow boundaries.
                if (net%roads(r)%lane(k)%open_in) then
                    if (net%roads(r)%lane(k)%cells(1)%has_car .and. &
                        .not. net%roads(r)%lane(k)%old(1)%has_car) then
                        road_entries(r, step) = road_entries(r, step) + 1
                    end if
                end if

                ! Count exits through open outflow boundaries.
                if (net%roads(r)%lane(k)%open_out) then
                    if (.not. net%roads(r)%lane(k)%cells(Lk)%has_car .and. &
                        net%roads(r)%lane(k)%old(Lk)%has_car) then
                        road_exits(r, step) = road_exits(r, step) + 1
                    end if
                end if
            end do

            ! Per-road density is total occupied cells divided by total lane cells:
            !
            ! .. math::
            !
            !    \rho_r(t)
            !    =
            !    \frac{N_r(t)}
            !    {\sum_{\ell \in r} L_{\ell}}.
            road_density(r, step) = real(n_occ) / real(sum_lane_lengths(net%roads(r)%lane))
        end do

        if (mod(step, log_every) == 0 .or. step == n_steps) then
            write(*,'(a,i0,a,i0)') "step ", step, "/", n_steps
            flush(6)
        end if
    end do

    ! Write NetCDF output.
    call write_network_netcdf(trim(out_path), net, occupancy, velocity, &
                              road_density, road_entries, road_exits, &
                              n_steps, params%rng_seed, params%v_max, &
                              params%p_slow, params%dt)

    write(*,'(a,a)') "wrote: ", trim(out_path)

    call free_network(net)
    call free_spec(spec)

contains

    pure function sum_lane_lengths(lanes) result(s)
        ! Return the total number of cells across a road's lanes.
        !
        ! This helper is used to normalise road density. If a road has lanes with
        ! lengths ``L_1, L_2, ..., L_n``, the returned value is:
        !
        ! .. math::
        !
        !    S = \sum_{\ell=1}^{n} L_{\ell}.
        type(lane_t), intent(in) :: lanes(:)
        ! Lanes belonging to one road.
        integer :: s
        ! Sum of lane lengths.
        integer :: i
        ! Lane index.

        s = 0

        do i = 1, size(lanes)
            s = s + lanes(i)%length
        end do
    end function sum_lane_lengths

end program run_network