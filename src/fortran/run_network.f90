program run_network
    !--------------------------------------------------------------------
    ! Driver: read JSON config, build network, run simulation, write NC.
    !
    !   ./build/run_network <config.json> <out.nc>
    !
    ! P1 STUB: parses config, builds network, prints a structural summary,
    ! frees and exits.  Simulation loop + NetCDF writer added in P2.
    !--------------------------------------------------------------------
    use road_network_mod
    use network_init_mod, only: free_network
    use network_builder_mod
    use json_config_mod
    implicit none

    character(len=1024) :: config_path, out_path
    integer :: argc, r, j
    type(network_spec_t) :: spec
    type(sim_params_t)   :: params
    type(road_network_t) :: net

    argc = command_argument_count()
    if (argc < 2) then
        write(*,*) "usage: run_network <config.json> <out.nc>"
        stop 1
    end if
    call get_command_argument(1, config_path)
    call get_command_argument(2, out_path)

    write(*,'(a,a)') "config: ", trim(config_path)
    write(*,'(a,a)') "output: ", trim(out_path)

    call read_config(trim(config_path), spec, params)

    write(*,'(a)') "--- params ---"
    write(*,'(a,i0)')   "  schema_version  = ", params%schema_version
    write(*,'(a,i0)')   "  n_steps         = ", params%n_steps
    write(*,'(a,f6.3)') "  dt              = ", params%dt
    write(*,'(a,i0)')   "  v_max           = ", params%v_max
    write(*,'(a,f6.3)') "  p_slow          = ", params%p_slow
    write(*,'(a,i0)')   "  rng_seed        = ", params%rng_seed
    write(*,'(a,i0)')   "  max_lane_length = ", params%max_lane_length

    call build_network(spec, net)

    write(*,'(a)') "--- network ---"
    write(*,'(a,i0)') "  n_roads     = ", size(net%roads)
    write(*,'(a,i0)') "  n_junctions = ", size(net%junctions)
    do r = 1, size(net%roads)
        write(*,'(a,i0,a,i0,a,i0,a,i0,a)') &
            "  road ", net%roads(r)%id, &
            " lanes=", size(net%roads(r)%lane), &
            " end_junction=[", net%roads(r)%end_junction(1), &
            ",", net%roads(r)%end_junction(2), "]"
    end do
    do j = 1, size(net%junctions)
        write(*,'(a,i0,a,i0,a,i0,a,l1)') &
            "  junc ", net%junctions(j)%id, &
            " n_in=", net%junctions(j)%n_in, &
            " n_out=", net%junctions(j)%n_out, &
            " has_perim=", allocated(net%junctions(j)%in_perim)
    end do

    write(*,'(a)') "P1 stub: build OK; simulation + NC writer come in P2."

    call free_network(net)
    call free_spec(spec)
end program run_network
