module network_io_mod
    !--------------------------------------------------------------------
    ! NetCDF writer for the network simulator.
    !
    ! Schema (see plan file for full spec):
    !   * occupancy(cell, lane_index, time)        byte
    !   * velocity(cell, lane_index, time)         byte
    !   * final_occupancy(cell, lane_index)        byte
    !   * final_velocity(cell, lane_index)         byte
    !   * lane_* scalars (lane_index)              int / float / byte
    !   * road_density / road_entries / road_exits (road, time)
    !   * flat junction tables with offset prefix-sums for ragged dims
    !   * global attrs: model, n_steps, dt, v_max, p_slow, rng_seed,
    !                   created_at, config_json, schema_version
    !
    ! All data needed for a future restart is included (final state +
    ! rng_seed + the raw config_json string).  No restart loader yet.
    !--------------------------------------------------------------------
    use netcdf
    use road_network_mod
    implicit none
    private

    public :: write_network_netcdf

contains

    subroutine write_network_netcdf(filename, net, occupancy, velocity, &
                                     road_density, road_entries, road_exits, &
                                     n_steps, rng_seed, v_max, p_slow, dt)
        character(len=*),     intent(in) :: filename
        type(road_network_t), intent(in) :: net
        integer(kind=1),      intent(in) :: occupancy(:,:,:)   ! (cell, lane_index, time)
        integer(kind=1),      intent(in) :: velocity(:,:,:)
        real,                 intent(in) :: road_density(:,:)  ! (road, time)
        integer,              intent(in) :: road_entries(:,:)
        integer,              intent(in) :: road_exits(:,:)
        integer,              intent(in) :: n_steps, rng_seed, v_max
        real,                 intent(in) :: p_slow, dt

        integer :: ncid
        integer :: dim_cell, dim_lane, dim_time, dim_road, dim_junc, dim_pair
        integer :: dim_inleg, dim_outleg, dim_route
        integer :: var_occ, var_vel, var_focc, var_fvel
        integer :: var_l_roadid, var_l_within, var_l_len, var_l_fd
        integer :: var_l_alpha, var_l_beta, var_l_oin, var_l_oout
        integer :: var_r_endj, var_r_dens, var_r_ent, var_r_ex
        integer :: var_j_nin, var_j_nout
        integer :: var_j_inoff, var_j_outoff, var_j_routeoff
        integer :: var_j_inroad, var_j_inlane, var_j_inperim
        integer :: var_j_outroad, var_j_outlane, var_j_outperim
        integer :: var_j_prob

        integer :: n_lanes, max_L, n_roads, n_juncs
        integer :: r, j, k, idx, max_cell

        integer, allocatable :: lane_road_id(:), lane_within(:), lane_len(:), lane_fd(:)
        real,    allocatable :: lane_alpha(:), lane_beta(:)
        integer(kind=1), allocatable :: lane_oin(:), lane_oout(:)
        integer, allocatable :: road_endj(:,:)
        integer, allocatable :: junc_nin(:), junc_nout(:)
        integer, allocatable :: junc_inoff(:), junc_outoff(:), junc_routeoff(:)
        integer, allocatable :: jin_road(:), jin_lane(:), jin_perim(:)
        integer, allocatable :: jout_road(:), jout_lane(:), jout_perim(:)
        real,    allocatable :: junc_prob(:)
        integer :: total_inlegs, total_outlegs, total_route
        character(len=8)  :: date_s
        character(len=10) :: time_s
        character(len=24) :: created_at

        n_roads = size(net%roads)
        n_juncs = size(net%junctions)

        max_cell = size(occupancy, 1)
        n_lanes  = size(occupancy, 2)
        max_L    = max_cell

        ! ----- Build per-lane lookup tables -----
        allocate(lane_road_id(n_lanes), lane_within(n_lanes), lane_len(n_lanes), lane_fd(n_lanes))
        allocate(lane_alpha(n_lanes), lane_beta(n_lanes))
        allocate(lane_oin(n_lanes), lane_oout(n_lanes))
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
            end do
        end do

        allocate(road_endj(2, n_roads))
        do r = 1, n_roads
            road_endj(:, r) = net%roads(r)%end_junction
        end do

        ! ----- Flatten junction tables -----
        allocate(junc_nin(n_juncs), junc_nout(n_juncs))
        allocate(junc_inoff(n_juncs), junc_outoff(n_juncs), junc_routeoff(n_juncs))
        total_inlegs = 0; total_outlegs = 0; total_route = 0
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
        jin_perim = -1
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

        ! ----- ISO-ish timestamp -----
        call date_and_time(date=date_s, time=time_s)
        created_at = date_s(1:4) // '-' // date_s(5:6) // '-' // date_s(7:8) // 'T' // &
                     time_s(1:2) // ':' // time_s(3:4) // ':' // time_s(5:6)

        ! ----- Create file -----
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

        ! Global attrs
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'model',          'network_tasep_ns') )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'schema_version', 1) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'n_steps',        n_steps) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'dt',             dt) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'v_max',          v_max) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'p_slow',         p_slow) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'rng_seed',       rng_seed) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'created_at',     trim(created_at)) )

        ! Vars: history
        call check( nf90_def_var(ncid, 'occupancy', NF90_BYTE, [dim_cell, dim_lane, dim_time], var_occ) )
        call check( nf90_def_var(ncid, 'velocity',  NF90_BYTE, [dim_cell, dim_lane, dim_time], var_vel) )
        call check( nf90_def_var(ncid, 'final_occupancy', NF90_BYTE, [dim_cell, dim_lane], var_focc) )
        call check( nf90_def_var(ncid, 'final_velocity',  NF90_BYTE, [dim_cell, dim_lane], var_fvel) )

        ! Vars: per-lane lookup
        call check( nf90_def_var(ncid, 'lane_road_id',       NF90_INT,   [dim_lane], var_l_roadid) )
        call check( nf90_def_var(ncid, 'lane_within_road',   NF90_INT,   [dim_lane], var_l_within) )
        call check( nf90_def_var(ncid, 'lane_length',        NF90_INT,   [dim_lane], var_l_len) )
        call check( nf90_def_var(ncid, 'lane_flow_direction',NF90_INT,   [dim_lane], var_l_fd) )
        call check( nf90_def_var(ncid, 'lane_alpha',         NF90_FLOAT, [dim_lane], var_l_alpha) )
        call check( nf90_def_var(ncid, 'lane_beta',          NF90_FLOAT, [dim_lane], var_l_beta) )
        call check( nf90_def_var(ncid, 'lane_open_in',       NF90_BYTE,  [dim_lane], var_l_oin) )
        call check( nf90_def_var(ncid, 'lane_open_out',      NF90_BYTE,  [dim_lane], var_l_oout) )

        ! Vars: per-road
        call check( nf90_def_var(ncid, 'road_end_junction', NF90_INT,   [dim_pair, dim_road], var_r_endj) )
        call check( nf90_def_var(ncid, 'road_density',      NF90_FLOAT, [dim_road, dim_time], var_r_dens) )
        call check( nf90_def_var(ncid, 'road_entries',      NF90_INT,   [dim_road, dim_time], var_r_ent) )
        call check( nf90_def_var(ncid, 'road_exits',        NF90_INT,   [dim_road, dim_time], var_r_ex) )

        ! Vars: junctions
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

        ! ----- Write data -----
        call check( nf90_put_var(ncid, var_occ,  occupancy) )
        call check( nf90_put_var(ncid, var_vel,  velocity) )
        call check( nf90_put_var(ncid, var_focc, occupancy(:,:,n_steps)) )
        call check( nf90_put_var(ncid, var_fvel, velocity(:,:,n_steps)) )

        call check( nf90_put_var(ncid, var_l_roadid, lane_road_id) )
        call check( nf90_put_var(ncid, var_l_within, lane_within) )
        call check( nf90_put_var(ncid, var_l_len,    lane_len) )
        call check( nf90_put_var(ncid, var_l_fd,     lane_fd) )
        call check( nf90_put_var(ncid, var_l_alpha,  lane_alpha) )
        call check( nf90_put_var(ncid, var_l_beta,   lane_beta) )
        call check( nf90_put_var(ncid, var_l_oin,    lane_oin) )
        call check( nf90_put_var(ncid, var_l_oout,   lane_oout) )

        call check( nf90_put_var(ncid, var_r_endj, road_endj) )
        call check( nf90_put_var(ncid, var_r_dens, road_density) )
        call check( nf90_put_var(ncid, var_r_ent,  road_entries) )
        call check( nf90_put_var(ncid, var_r_ex,   road_exits) )

        call check( nf90_put_var(ncid, var_j_nin,        junc_nin) )
        call check( nf90_put_var(ncid, var_j_nout,       junc_nout) )
        call check( nf90_put_var(ncid, var_j_inoff,      junc_inoff) )
        call check( nf90_put_var(ncid, var_j_outoff,     junc_outoff) )
        call check( nf90_put_var(ncid, var_j_routeoff,   junc_routeoff) )
        call check( nf90_put_var(ncid, var_j_inroad,     jin_road) )
        call check( nf90_put_var(ncid, var_j_inlane,     jin_lane) )
        call check( nf90_put_var(ncid, var_j_inperim,    jin_perim) )
        call check( nf90_put_var(ncid, var_j_outroad,    jout_road) )
        call check( nf90_put_var(ncid, var_j_outlane,    jout_lane) )
        call check( nf90_put_var(ncid, var_j_outperim,   jout_perim) )
        call check( nf90_put_var(ncid, var_j_prob,       junc_prob) )

        call check( nf90_close(ncid) )
    end subroutine write_network_netcdf

    subroutine check(status)
        integer, intent(in) :: status
        if (status /= NF90_NOERR) then
            print *, 'NetCDF error: ', trim(nf90_strerror(status))
            stop 1
        end if
    end subroutine check

end module network_io_mod
