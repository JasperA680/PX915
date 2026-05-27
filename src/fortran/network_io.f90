module network_io_mod
    ! NetCDF output routines for the road-network simulator.
    !
    ! This module writes network simulation histories to a NetCDF file. The output
    ! file stores time-dependent cell occupancy and velocity, final-state arrays,
    ! lane metadata, road-level diagnostics and flattened junction tables.
    !
    ! The main history variables use the layout:
    !
    ! .. math::
    !
    !    \mathrm{occupancy}(c, l, t),
    !
    ! where ``c`` is the cell index, ``l`` is the flattened lane index and ``t``
    ! is the time index.
    !
    ! The output also includes enough information to reconstruct the final network
    ! state in a future restart workflow: final occupancy, final velocity, the
    ! random seed and the network metadata. A restart loader is not currently
    ! implemented.
    !
    ! The junction data are stored in flattened arrays with offset vectors. This
    ! avoids needing ragged NetCDF dimensions. For example, the inbound legs of
    ! junction ``j`` begin at:
    !
    ! .. math::
    !
    !    o_j = \mathrm{junction\_inleg\_offset}(j).
    !
    ! The number of inbound legs is stored separately as
    ! ``junction_n_in(j)``.

    use netcdf
    use road_network_mod

    implicit none
    private

    public :: write_network_netcdf, check

contains

    subroutine write_network_netcdf(filename, net, occupancy, velocity, &
                                     road_density, road_entries, road_exits, &
                                     n_steps, rng_seed, v_max, p_slow, dt)
        ! Write a complete road-network simulation history to a NetCDF file.
        !
        ! The file contains:
        !
        ! * cell-level occupancy and velocity histories;
        ! * final occupancy and velocity arrays;
        ! * lane lookup tables mapping flattened lane indices to roads and lanes
        !   (including the ``lane_is_periodic`` flag so closed-ring NS lanes can
        !   be identified after the fact);
        ! * road-level density, entry and exit histories;
        ! * flattened junction connectivity and route-probability tables;
        ! * global metadata describing the simulation.
        !
        ! The occupancy and velocity histories are expected to have shape
        ! ``(cell, lane_index, time)``. The road diagnostics are expected to have
        ! shape ``(road, time)``.
        !
        ! The flattened lane index is built by scanning roads in order and then
        ! scanning lanes within each road. The resulting lookup arrays make it
        ! possible to recover the original road and lane associated with each
        ! flattened lane index.
        character(len=*), intent(in) :: filename
        ! Name of the NetCDF file to create.
        type(road_network_t), intent(in) :: net
        ! Road network whose topology and final state are written.
        integer(kind=1), intent(in) :: occupancy(:,:,:)
        ! Occupancy history, indexed as ``occupancy(cell, lane_index, time)``.
        integer(kind=1), intent(in) :: velocity(:,:,:)
        ! Velocity history, indexed as ``velocity(cell, lane_index, time)``.
        real, intent(in) :: road_density(:,:)
        ! Road-level density history, indexed as ``road_density(road, time)``.
        integer, intent(in) :: road_entries(:,:)
        ! Number of entries into each road at each time step.
        integer, intent(in) :: road_exits(:,:)
        ! Number of exits from each road at each time step.
        integer, intent(in) :: n_steps
        ! Number of stored time steps.
        integer, intent(in) :: rng_seed
        ! Random seed associated with the simulation.
        integer, intent(in) :: v_max
        ! Maximum vehicle speed used by the simulation.
        real, intent(in) :: p_slow
        ! Random slowing probability used by the Nagel-Schreckenberg model.
        real, intent(in) :: dt
        ! Time step size stored as metadata.

        integer :: ncid
        ! NetCDF file identifier.
        integer :: dim_cell
        ! NetCDF dimension identifier for cells.
        integer :: dim_lane
        ! NetCDF dimension identifier for flattened lane index.
        integer :: dim_time
        ! NetCDF dimension identifier for time.
        integer :: dim_road
        ! NetCDF dimension identifier for roads.
        integer :: dim_junc
        ! NetCDF dimension identifier for junctions.
        integer :: dim_pair
        ! NetCDF dimension identifier of length two.
        integer :: dim_inleg
        ! NetCDF dimension identifier for flattened inbound junction legs.
        integer :: dim_outleg
        ! NetCDF dimension identifier for flattened outbound junction legs.
        integer :: dim_route
        ! NetCDF dimension identifier for flattened route probabilities.

        integer :: var_occ
        ! NetCDF variable identifier for occupancy history.
        integer :: var_vel
        ! NetCDF variable identifier for velocity history.
        integer :: var_focc
        ! NetCDF variable identifier for final occupancy.
        integer :: var_fvel
        ! NetCDF variable identifier for final velocity.
        integer :: var_l_roadid
        ! NetCDF variable identifier for lane road identifiers.
        integer :: var_l_within
        ! NetCDF variable identifier for lane indices within roads.
        integer :: var_l_len
        ! NetCDF variable identifier for lane lengths.
        integer :: var_l_fd
        ! NetCDF variable identifier for lane flow directions.
        integer :: var_l_alpha
        ! NetCDF variable identifier for lane alpha values.
        integer :: var_l_beta
        ! NetCDF variable identifier for lane beta values.
        integer :: var_l_oin
        ! NetCDF variable identifier for lane open-in flags.
        integer :: var_l_oout
        ! NetCDF variable identifier for lane open-out flags.
        integer :: var_l_period
        ! NetCDF variable identifier for lane is_periodic flags. (The
        ! companion ``lane_n_vehicles`` is an input-time field consumed by
        ! ``place_evenly`` at build time; it lives only in the config NetCDF.)
        integer :: var_r_endj
        ! NetCDF variable identifier for road endpoint junctions.
        integer :: var_r_dens
        ! NetCDF variable identifier for road density history.
        integer :: var_r_ent
        ! NetCDF variable identifier for road entry counts.
        integer :: var_r_ex
        ! NetCDF variable identifier for road exit counts.
        integer :: var_j_nin
        ! NetCDF variable identifier for number of inbound legs per junction.
        integer :: var_j_nout
        ! NetCDF variable identifier for number of outbound legs per junction.
        integer :: var_j_inoff
        ! NetCDF variable identifier for inbound-leg offsets.
        integer :: var_j_outoff
        ! NetCDF variable identifier for outbound-leg offsets.
        integer :: var_j_routeoff
        ! NetCDF variable identifier for route-probability offsets.
        integer :: var_j_inroad
        ! NetCDF variable identifier for flattened inbound road indices.
        integer :: var_j_inlane
        ! NetCDF variable identifier for flattened inbound lane indices.
        integer :: var_j_inperim
        ! NetCDF variable identifier for flattened inbound perimeter ports.
        integer :: var_j_outroad
        ! NetCDF variable identifier for flattened outbound road indices.
        integer :: var_j_outlane
        ! NetCDF variable identifier for flattened outbound lane indices.
        integer :: var_j_outperim
        ! NetCDF variable identifier for flattened outbound perimeter ports.
        integer :: var_j_prob
        ! NetCDF variable identifier for flattened route probabilities.

        integer :: n_lanes
        ! Number of flattened lanes in the output.
        integer :: max_L
        ! Maximum stored lane length, equal to the cell dimension.
        integer :: n_roads
        ! Number of roads in the network.
        integer :: n_juncs
        ! Number of junctions in the network.
        integer :: r
        ! Road index.
        integer :: j
        ! Junction index.
        integer :: k
        ! Lane, leg or route index.
        integer :: idx
        ! Flattened-array index.
        integer :: max_cell
        ! Number of cells in the first dimension of the history arrays.

        integer, allocatable :: lane_road_id(:)
        ! Road identifier for each flattened lane.
        integer, allocatable :: lane_within(:)
        ! Lane index within its road for each flattened lane.
        integer, allocatable :: lane_len(:)
        ! Length of each flattened lane.
        integer, allocatable :: lane_fd(:)
        ! Flow direction of each flattened lane.
        real, allocatable :: lane_alpha(:)
        ! Entry probability for each flattened lane.
        real, allocatable :: lane_beta(:)
        ! Exit probability for each flattened lane.
        integer(kind=1), allocatable :: lane_oin(:)
        ! Open-inflow flag for each flattened lane.
        integer(kind=1), allocatable :: lane_oout(:)
        ! Open-outflow flag for each flattened lane.
        integer(kind=1), allocatable :: lane_period(:)
        ! Periodic-boundary flag for each flattened lane.
        integer, allocatable :: road_endj(:,:)
        ! Endpoint junction identifiers for each road, shape ``(2, n_roads)``.

        integer, allocatable :: junc_nin(:)
        ! Number of inbound legs for each junction.
        integer, allocatable :: junc_nout(:)
        ! Number of outbound legs for each junction.
        integer, allocatable :: junc_inoff(:)
        ! Offset into flattened inbound-leg arrays for each junction.
        integer, allocatable :: junc_outoff(:)
        ! Offset into flattened outbound-leg arrays for each junction.
        integer, allocatable :: junc_routeoff(:)
        ! Offset into flattened route-probability array for each junction.
        integer, allocatable :: jin_road(:)
        ! Flattened inbound road identifiers.
        integer, allocatable :: jin_lane(:)
        ! Flattened inbound lane indices.
        integer, allocatable :: jin_perim(:)
        ! Flattened inbound perimeter-port indices.
        integer, allocatable :: jout_road(:)
        ! Flattened outbound road identifiers.
        integer, allocatable :: jout_lane(:)
        ! Flattened outbound lane indices.
        integer, allocatable :: jout_perim(:)
        ! Flattened outbound perimeter-port indices.
        real, allocatable :: junc_prob(:)
        ! Flattened junction route probabilities.

        integer :: total_inlegs
        ! Total number of inbound junction legs over all junctions.
        integer :: total_outlegs
        ! Total number of outbound junction legs over all junctions.
        integer :: total_route
        ! Total number of route-probability entries over all junctions.
        character(len=8) :: date_s
        ! Date string returned by ``date_and_time``.
        character(len=10) :: time_s
        ! Time string returned by ``date_and_time``.
        character(len=24) :: created_at
        ! ISO-style creation timestamp written as a global attribute.

        n_roads = size(net%roads)
        n_juncs = size(net%junctions)

        max_cell = size(occupancy, 1)
        n_lanes  = size(occupancy, 2)
        max_L    = max_cell

        ! Build per-lane lookup tables.
        allocate(lane_road_id(n_lanes), lane_within(n_lanes), lane_len(n_lanes), lane_fd(n_lanes))
        allocate(lane_alpha(n_lanes), lane_beta(n_lanes))
        allocate(lane_oin(n_lanes), lane_oout(n_lanes))
        allocate(lane_period(n_lanes))

        idx = 0
        do r = 1, n_roads
            do k = 1, size(net%roads(r)%lane)
                idx = idx + 1

                lane_road_id(idx) = net%roads(r)%id
                lane_within(idx)  = k
                lane_len(idx)     = net%roads(r)%lane(k)%length
                lane_fd(idx)      = net%roads(r)%lane(k)%flow_direction
                lane_alpha(idx)   = net%roads(r)%lane(k)%alpha
                lane_beta(idx)    = net%roads(r)%lane(k)%beta
                lane_oin(idx)     = merge(1_1, 0_1, net%roads(r)%lane(k)%open_in)
                lane_oout(idx)    = merge(1_1, 0_1, net%roads(r)%lane(k)%open_out)
                lane_period(idx)  = merge(1_1, 0_1, net%roads(r)%lane(k)%is_periodic)
            end do
        end do

        allocate(road_endj(2, n_roads))

        do r = 1, n_roads
            road_endj(:, r) = net%roads(r)%end_junction
        end do

        ! Flatten junction tables using offset arrays.
        allocate(junc_nin(n_juncs), junc_nout(n_juncs))
        allocate(junc_inoff(n_juncs), junc_outoff(n_juncs), junc_routeoff(n_juncs))

        total_inlegs  = 0
        total_outlegs = 0
        total_route   = 0

        do j = 1, n_juncs
            junc_nin(j)  = net%junctions(j)%n_in
            junc_nout(j) = net%junctions(j)%n_out

            junc_inoff(j)    = total_inlegs
            junc_outoff(j)   = total_outlegs
            junc_routeoff(j) = total_route

            total_inlegs  = total_inlegs  + net%junctions(j)%n_in
            total_outlegs = total_outlegs + net%junctions(j)%n_out
            total_route   = total_route   + net%junctions(j)%n_in * net%junctions(j)%n_out
        end do

        allocate(jin_road(max(1,total_inlegs)), jin_lane(max(1,total_inlegs)), jin_perim(max(1,total_inlegs)))
        allocate(jout_road(max(1,total_outlegs)), jout_lane(max(1,total_outlegs)), jout_perim(max(1,total_outlegs)))
        allocate(junc_prob(max(1,total_route)))

        jin_perim  = -1
        jout_perim = -1

        idx = 0
        do j = 1, n_juncs
            do k = 1, net%junctions(j)%n_in
                idx = idx + 1

                jin_road(idx) = net%junctions(j)%in_road(k)
                jin_lane(idx) = net%junctions(j)%in_lane(k)

                if (allocated(net%junctions(j)%in_perim)) jin_perim(idx) = net%junctions(j)%in_perim(k)
            end do
        end do

        idx = 0
        do j = 1, n_juncs
            do k = 1, net%junctions(j)%n_out
                idx = idx + 1

                jout_road(idx) = net%junctions(j)%out_road(k)
                jout_lane(idx) = net%junctions(j)%out_lane(k)

                if (allocated(net%junctions(j)%out_perim)) jout_perim(idx) = net%junctions(j)%out_perim(k)
            end do
        end do

        idx = 0
        do j = 1, n_juncs
            do k = 1, net%junctions(j)%n_in
                junc_prob(idx+1 : idx+net%junctions(j)%n_out) = net%junctions(j)%in_routes(k)%prob
                idx = idx + net%junctions(j)%n_out
            end do
        end do

        ! Construct an ISO-style creation timestamp.
        call date_and_time(date=date_s, time=time_s)

        created_at = date_s(1:4) // '-' // date_s(5:6) // '-' // date_s(7:8) // 'T' // &
                     time_s(1:2) // ':' // time_s(3:4) // ':' // time_s(5:6)

        ! Create the NetCDF file and define dimensions.
        call check( nf90_create(filename, NF90_CLOBBER, ncid) )

        call check( nf90_def_dim(ncid, 'cell',       max_L,         dim_cell) )
        call check( nf90_def_dim(ncid, 'lane_index', n_lanes,       dim_lane) )
        call check( nf90_def_dim(ncid, 'time',       n_steps,       dim_time) )
        call check( nf90_def_dim(ncid, 'road',       n_roads,       dim_road) )
        call check( nf90_def_dim(ncid, 'junction',   n_juncs,       dim_junc) )
        call check( nf90_def_dim(ncid, 'pair',       2,             dim_pair) )
        call check( nf90_def_dim(ncid, 'junction_inleg_total',  max(1,total_inlegs),  dim_inleg) )
        call check( nf90_def_dim(ncid, 'junction_outleg_total', max(1,total_outlegs), dim_outleg) )
        call check( nf90_def_dim(ncid, 'route_total',           max(1,total_route),   dim_route) )

        ! Write global metadata.
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'model',          'network_tasep_ns') )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'schema_version', 1) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'n_steps',        n_steps) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'dt',             dt) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'v_max',          v_max) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'p_slow',         p_slow) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'rng_seed',       rng_seed) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'created_at',     trim(created_at)) )

        ! Define history variables.
        call check( nf90_def_var(ncid, 'occupancy', NF90_BYTE, [dim_cell, dim_lane, dim_time], var_occ) )
        call check( nf90_def_var(ncid, 'velocity',  NF90_BYTE, [dim_cell, dim_lane, dim_time], var_vel) )
        call check( nf90_def_var(ncid, 'final_occupancy', NF90_BYTE, [dim_cell, dim_lane], var_focc) )
        call check( nf90_def_var(ncid, 'final_velocity',  NF90_BYTE, [dim_cell, dim_lane], var_fvel) )

        ! Define per-lane lookup variables.
        call check( nf90_def_var(ncid, 'lane_road_id',        NF90_INT,   [dim_lane], var_l_roadid) )
        call check( nf90_def_var(ncid, 'lane_within_road',    NF90_INT,   [dim_lane], var_l_within) )
        call check( nf90_def_var(ncid, 'lane_length',         NF90_INT,   [dim_lane], var_l_len) )
        call check( nf90_def_var(ncid, 'lane_flow_direction', NF90_INT,   [dim_lane], var_l_fd) )
        call check( nf90_def_var(ncid, 'lane_alpha',          NF90_FLOAT, [dim_lane], var_l_alpha) )
        call check( nf90_def_var(ncid, 'lane_beta',           NF90_FLOAT, [dim_lane], var_l_beta) )
        call check( nf90_def_var(ncid, 'lane_open_in',        NF90_BYTE,  [dim_lane], var_l_oin) )
        call check( nf90_def_var(ncid, 'lane_open_out',       NF90_BYTE,  [dim_lane], var_l_oout) )
        call check( nf90_def_var(ncid, 'lane_is_periodic',    NF90_BYTE,  [dim_lane], var_l_period) )

        ! Define per-road variables.
        call check( nf90_def_var(ncid, 'road_end_junction', NF90_INT,   [dim_pair, dim_road], var_r_endj) )
        call check( nf90_def_var(ncid, 'road_density',      NF90_FLOAT, [dim_road, dim_time], var_r_dens) )
        call check( nf90_def_var(ncid, 'road_entries',      NF90_INT,   [dim_road, dim_time], var_r_ent) )
        call check( nf90_def_var(ncid, 'road_exits',        NF90_INT,   [dim_road, dim_time], var_r_ex) )

        ! Define flattened junction variables.
        call check( nf90_def_var(ncid, 'junction_n_in',          NF90_INT, [dim_junc], var_j_nin) )
        call check( nf90_def_var(ncid, 'junction_n_out',         NF90_INT, [dim_junc], var_j_nout) )
        call check( nf90_def_var(ncid, 'junction_inleg_offset',  NF90_INT, [dim_junc], var_j_inoff) )
        call check( nf90_def_var(ncid, 'junction_outleg_offset', NF90_INT, [dim_junc], var_j_outoff) )
        call check( nf90_def_var(ncid, 'junction_route_offset',  NF90_INT, [dim_junc], var_j_routeoff) )
        call check( nf90_def_var(ncid, 'junction_in_road',       NF90_INT, [dim_inleg],  var_j_inroad) )
        call check( nf90_def_var(ncid, 'junction_in_lane',       NF90_INT, [dim_inleg],  var_j_inlane) )
        call check( nf90_def_var(ncid, 'junction_in_perim',      NF90_INT, [dim_inleg],  var_j_inperim) )
        call check( nf90_def_var(ncid, 'junction_out_road',      NF90_INT, [dim_outleg], var_j_outroad) )
        call check( nf90_def_var(ncid, 'junction_out_lane',      NF90_INT, [dim_outleg], var_j_outlane) )
        call check( nf90_def_var(ncid, 'junction_out_perim',     NF90_INT, [dim_outleg], var_j_outperim) )
        call check( nf90_def_var(ncid, 'junction_route_prob',    NF90_FLOAT, [dim_route], var_j_prob) )

        call check( nf90_enddef(ncid) )

        ! Write history and final-state variables.
        call check( nf90_put_var(ncid, var_occ,  occupancy) )
        call check( nf90_put_var(ncid, var_vel,  velocity) )
        call check( nf90_put_var(ncid, var_focc, occupancy(:,:,n_steps)) )
        call check( nf90_put_var(ncid, var_fvel, velocity(:,:,n_steps)) )

        ! Write per-lane lookup tables.
        call check( nf90_put_var(ncid, var_l_roadid, lane_road_id) )
        call check( nf90_put_var(ncid, var_l_within, lane_within) )
        call check( nf90_put_var(ncid, var_l_len,    lane_len) )
        call check( nf90_put_var(ncid, var_l_fd,     lane_fd) )
        call check( nf90_put_var(ncid, var_l_alpha,  lane_alpha) )
        call check( nf90_put_var(ncid, var_l_beta,   lane_beta) )
        call check( nf90_put_var(ncid, var_l_oin,    lane_oin) )
        call check( nf90_put_var(ncid, var_l_oout,   lane_oout) )
        call check( nf90_put_var(ncid, var_l_period, lane_period) )

        ! Write per-road diagnostics.
        call check( nf90_put_var(ncid, var_r_endj, road_endj) )
        call check( nf90_put_var(ncid, var_r_dens, road_density) )
        call check( nf90_put_var(ncid, var_r_ent,  road_entries) )
        call check( nf90_put_var(ncid, var_r_ex,   road_exits) )

        ! Write flattened junction tables.
        call check( nf90_put_var(ncid, var_j_nin,      junc_nin) )
        call check( nf90_put_var(ncid, var_j_nout,     junc_nout) )
        call check( nf90_put_var(ncid, var_j_inoff,    junc_inoff) )
        call check( nf90_put_var(ncid, var_j_outoff,   junc_outoff) )
        call check( nf90_put_var(ncid, var_j_routeoff, junc_routeoff) )
        call check( nf90_put_var(ncid, var_j_inroad,   jin_road) )
        call check( nf90_put_var(ncid, var_j_inlane,   jin_lane) )
        call check( nf90_put_var(ncid, var_j_inperim,  jin_perim) )
        call check( nf90_put_var(ncid, var_j_outroad,  jout_road) )
        call check( nf90_put_var(ncid, var_j_outlane,  jout_lane) )
        call check( nf90_put_var(ncid, var_j_outperim, jout_perim) )
        call check( nf90_put_var(ncid, var_j_prob,     junc_prob) )

        call check( nf90_close(ncid) )
    end subroutine write_network_netcdf


    subroutine check(status)
        ! Check the return status from a NetCDF operation.
        !
        ! If the NetCDF call failed, this routine prints the NetCDF error message
        ! and stops execution. This avoids silently producing incomplete or
        ! corrupted output files.
        integer, intent(in) :: status
        ! NetCDF status code returned by an ``nf90_*`` routine.

        if (status /= NF90_NOERR) then
            print *, 'NetCDF error: ', trim(nf90_strerror(status))
            stop 1
        end if
    end subroutine check

end module network_io_mod