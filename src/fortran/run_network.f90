program run_network
    !--------------------------------------------------------------------
    ! Driver: read JSON config, build network, run simulation, write NC.
    !
    !   ./build/run_network <config.json> <out.nc>
    !
    ! 1. Parse JSON config (json_config_mod).
    ! 2. Build road_network_t from the spec (network_builder_mod).
    ! 3. Seed RNG deterministically from params%rng_seed.
    ! 4. Allocate occupancy/velocity history arrays sized to
    !    (max_lane_length, total_lanes, n_steps) bytes.
    ! 5. Loop network_step, recording state at the end of every step.
    !    Print 'step k/n' every max(1, n_steps/100) with flush(6) so
    !    the Python launcher can drive a progress bar.
    ! 6. Write NetCDF result (network_io_mod) including final_occupancy
    !    + rng_seed + raw config_json for future restart-readiness.
    !--------------------------------------------------------------------
    use road_network_mod
    use network_init_mod,         only: free_network
    use network_simulation_mod,   only: network_step
    use network_builder_mod
    use json_config_mod
    use network_io_mod
    implicit none

    character(len=1024) :: config_path, out_path
    integer :: argc
    type(network_spec_t)  :: spec
    type(sim_params_t)    :: params
    type(road_network_t)  :: net

    integer :: n_lanes, max_L, n_roads, n_steps
    integer :: step, r, k, idx, Lk, n_occ, log_every
    integer(kind=1), allocatable :: occupancy(:,:,:), velocity(:,:,:)
    real,            allocatable :: road_density(:,:)
    integer,         allocatable :: road_entries(:,:), road_exits(:,:)
    integer,         allocatable :: lane_road_offset(:)
    character(len=:), allocatable :: config_json
    integer, allocatable :: seed_arr(:)
    integer :: seed_size

    ! LANE CHANGING PARAMETERS - CURRENTLY HARDCODED
    integer, parameter   :: lane_change_model = -1           ! Defines what logic to use for lane changes: -1=DISABLED, 0=SYMMETRIC, 1=ASYMMETRIC
    real,    parameter   :: lane_change_prob = 1             ! Probability a car will switch lanes, given it meets other requirements

    argc = command_argument_count()
    if (argc < 2) then
        write(*,*) "usage: run_network <config.json> <out.nc>"
        stop 1
    end if
    call get_command_argument(1, config_path)
    call get_command_argument(2, out_path)

    write(*,'(a,a)') "config: ", trim(config_path)
    write(*,'(a,a)') "output: ", trim(out_path)

    ! Read config and also retain the raw JSON text for the NC attribute.
    call read_config(trim(config_path), spec, params)
    call load_text_file(trim(config_path), config_json)

    n_steps = params%n_steps

    ! Build network.
    call build_network(spec, net)
    n_roads = size(net%roads)

    ! Compute total lanes + max lane length (or use params override).
    n_lanes = 0
    max_L   = 0
    do r = 1, n_roads
        do k = 1, size(net%roads(r)%lane)
            n_lanes = n_lanes + 1
            max_L = max(max_L, net%roads(r)%lane(k)%length)
        end do
    end do
    if (params%max_lane_length > 0) max_L = max(max_L, params%max_lane_length)

    ! Per-road lane-flat-index offset (lane_road_offset(r) = first global lane index for road r, minus 1).
    allocate(lane_road_offset(n_roads))
    idx = 0
    do r = 1, n_roads
        lane_road_offset(r) = idx
        idx = idx + size(net%roads(r)%lane)
    end do

    ! Seed RNG deterministically.
    call random_seed(size=seed_size)
    allocate(seed_arr(seed_size))
    do k = 1, seed_size
        seed_arr(k) = params%rng_seed + 1013904223 * k
    end do
    call random_seed(put=seed_arr)

    ! Allocate history.
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

    ! Log the lane-change configuration so the GUI run log shows what's active.
    write(*,'(a,i0,a,f6.3)') &
        "lane change: model=", params%lc_model, " p_change=", params%lc_p_change

    ! ----- Simulation loop -----
    do step = 1, n_steps
        call network_step(net, params%model, params%lc_model, params%lc_p_change, params%v_max, params%p_slow)

        ! Record state and per-road counters.
        do r = 1, n_roads
            n_occ = 0
            do k = 1, size(net%roads(r)%lane)
                Lk = net%roads(r)%lane(k)%length
                idx = lane_road_offset(r) + k

                ! Occupancy + velocity for this lane (pad past Lk stays 0).
                where (net%roads(r)%lane(k)%cells(1:Lk)%has_car)
                    occupancy(1:Lk, idx, step) = 1_1
                end where
                velocity(1:Lk, idx, step) = int(net%roads(r)%lane(k)%cells(1:Lk)%velocity, kind=1)

                n_occ = n_occ + count(net%roads(r)%lane(k)%cells(1:Lk)%has_car)

                ! Entries: new car at site 1 of an open-in lane.
                if (net%roads(r)%lane(k)%open_in) then
                    if (net%roads(r)%lane(k)%cells(1)%has_car .and. &
                        .not. net%roads(r)%lane(k)%old(1)%has_car) then
                        road_entries(r, step) = road_entries(r, step) + 1
                    end if
                end if
                ! Exits: car disappeared from site L of an open-out lane.
                if (net%roads(r)%lane(k)%open_out) then
                    if (.not. net%roads(r)%lane(k)%cells(Lk)%has_car .and. &
                        net%roads(r)%lane(k)%old(Lk)%has_car) then
                        road_exits(r, step) = road_exits(r, step) + 1
                    end if
                end if
            end do
            ! Per-road density = total cars / total lane cells.
            road_density(r, step) = real(n_occ) / &
                real(sum_lane_lengths(net%roads(r)%lane))
        end do

        if (mod(step, log_every) == 0 .or. step == n_steps) then
            write(*,'(a,i0,a,i0)') "step ", step, "/", n_steps
            flush(6)
        end if
    end do

    ! ----- Write NetCDF -----
    call write_network_netcdf(trim(out_path), net, occupancy, velocity, &
                              road_density, road_entries, road_exits, &
                              n_steps, params%rng_seed, params%v_max, &
                              params%p_slow, params%dt, config_json)

    write(*,'(a,a)') "wrote: ", trim(out_path)

    call free_network(net)
    call free_spec(spec)

contains

    pure function sum_lane_lengths(lanes) result(s)
        type(lane_t), intent(in) :: lanes(:)
        integer :: s, i
        s = 0
        do i = 1, size(lanes)
            s = s + lanes(i)%length
        end do
    end function sum_lane_lengths

end program run_network
