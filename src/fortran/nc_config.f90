module nc_config_mod
    !--------------------------------------------------------------------
    ! NetCDF config reader for the network simulator.
    !
    ! Reads `config.nc` (written by Python `io.write_config_netcdf`) into a
    ! `sim_params_t` + `network_spec_t`.  Variable names mirror the schema
    ! used by ``network_io_mod`` for the output writer — so the read/write
    ! conventions stay symmetric across the pipeline.
    !
    ! Public surface matches the old json_config_mod (``sim_params_t`` and
    ! ``read_config``) so callers don't change.
    !--------------------------------------------------------------------
    use netcdf
    use network_builder_mod
    implicit none
    private

    type, public :: sim_params_t
        integer :: schema_version  = 1
        integer :: n_steps         = 1000
        real    :: dt              = 1.0
        integer :: v_max           = 5
        real    :: p_slow          = 0.2
        integer :: rng_seed        = 42
        integer :: max_lane_length = 0      ! 0 => compute from lanes
        character(len=16) :: model = "NS"   ! "NS" or "TASEP"
        ! Lane-change behaviour:
        !   -1 LC_DISABLED, 0 LC_SYMMETRIC, 1 LC_ASYMMETRIC
        integer :: lc_model    = -1
        real    :: lc_p_change = 1.0
    end type sim_params_t

    public :: read_config

contains

    subroutine read_config(path, spec, params)
        character(len=*),     intent(in)  :: path
        type(network_spec_t), intent(out) :: spec
        type(sim_params_t),   intent(out) :: params

        integer :: ncid
        integer :: n_roads, n_juncs_dim, n_juncs, n_lanes
        integer :: n_inlegs, n_outlegs, n_routes
        integer :: dim_road, dim_junc, dim_lane, dim_in, dim_out, dim_route

        integer, allocatable :: road_id(:), road_end_junction(:,:)
        integer, allocatable :: lane_road_id(:), lane_within(:), lane_length(:)
        integer, allocatable :: lane_fd(:)
        real,    allocatable :: lane_alpha(:), lane_beta(:)
        integer(kind=1), allocatable :: lane_open_in(:), lane_open_out(:)

        integer, allocatable :: j_id(:), j_nin(:), j_nout(:)
        integer, allocatable :: j_inoff(:), j_outoff(:), j_routeoff(:)
        integer, allocatable :: jin_road(:), jin_lane(:), jin_perim(:)
        integer, allocatable :: jout_road(:), jout_lane(:), jout_perim(:)
        real,    allocatable :: j_prob(:)

        integer :: r, j, k, off

        call check( nf90_open(trim(path), NF90_NOWRITE, ncid) )

        ! ----- Scalar params from global attributes -----
        call get_att_int (ncid, 'schema_version',  params%schema_version)
        call get_att_int (ncid, 'n_steps',         params%n_steps)
        call get_att_real(ncid, 'dt',              params%dt)
        call get_att_int (ncid, 'v_max',           params%v_max)
        call get_att_real(ncid, 'p_slow',          params%p_slow)
        call get_att_int (ncid, 'rng_seed',        params%rng_seed)
        call get_att_int (ncid, 'max_lane_length', params%max_lane_length)
        call get_att_int (ncid, 'lc_model',        params%lc_model)
        call get_att_real(ncid, 'lc_p_change',     params%lc_p_change)
        call get_att_string(ncid, 'model',         params%model)

        ! ----- Dimension sizes -----
        call check( nf90_inq_dimid(ncid, 'road',                    dim_road) )
        call check( nf90_inquire_dimension(ncid, dim_road,    len=n_roads) )

        call check( nf90_inq_dimid(ncid, 'junction',                dim_junc) )
        call check( nf90_inquire_dimension(ncid, dim_junc,    len=n_juncs_dim) )

        call check( nf90_inq_dimid(ncid, 'lane_index',              dim_lane) )
        call check( nf90_inquire_dimension(ncid, dim_lane,    len=n_lanes) )

        call check( nf90_inq_dimid(ncid, 'junction_inleg_total',    dim_in) )
        call check( nf90_inquire_dimension(ncid, dim_in,      len=n_inlegs) )

        call check( nf90_inq_dimid(ncid, 'junction_outleg_total',   dim_out) )
        call check( nf90_inquire_dimension(ncid, dim_out,     len=n_outlegs) )

        call check( nf90_inq_dimid(ncid, 'route_total',             dim_route) )
        call check( nf90_inquire_dimension(ncid, dim_route,   len=n_routes) )

        ! "n_junctions" attribute disambiguates "real 0 junctions" from "1
        ! sentinel slot to satisfy NetCDF's non-zero dim requirement".
        call get_att_int(ncid, 'n_junctions', n_juncs)

        ! ----- Read per-road -----
        allocate(road_id(n_roads), road_end_junction(2, n_roads))
        call read_int_1d(ncid, 'road_id',           road_id)
        call read_int_2d(ncid, 'road_end_junction', road_end_junction)

        ! ----- Read per-lane -----
        allocate(lane_road_id(n_lanes), lane_within(n_lanes), lane_length(n_lanes), lane_fd(n_lanes))
        allocate(lane_alpha(n_lanes), lane_beta(n_lanes))
        allocate(lane_open_in(n_lanes), lane_open_out(n_lanes))
        call read_int_1d(ncid, 'lane_road_id',        lane_road_id)
        call read_int_1d(ncid, 'lane_within_road',    lane_within)
        call read_int_1d(ncid, 'lane_length',         lane_length)
        call read_int_1d(ncid, 'lane_flow_direction', lane_fd)
        call read_real_1d(ncid, 'lane_alpha',         lane_alpha)
        call read_real_1d(ncid, 'lane_beta',          lane_beta)
        call read_byte_1d(ncid, 'lane_open_in',       lane_open_in)
        call read_byte_1d(ncid, 'lane_open_out',      lane_open_out)

        ! ----- Read per-junction -----
        if (n_juncs > 0) then
            allocate(j_id(n_juncs_dim), j_nin(n_juncs_dim), j_nout(n_juncs_dim))
            allocate(j_inoff(n_juncs_dim), j_outoff(n_juncs_dim), j_routeoff(n_juncs_dim))
            call read_int_1d(ncid, 'junction_id',            j_id)
            call read_int_1d(ncid, 'junction_n_in',          j_nin)
            call read_int_1d(ncid, 'junction_n_out',         j_nout)
            call read_int_1d(ncid, 'junction_inleg_offset',  j_inoff)
            call read_int_1d(ncid, 'junction_outleg_offset', j_outoff)
            call read_int_1d(ncid, 'junction_route_offset',  j_routeoff)

            allocate(jin_road(n_inlegs), jin_lane(n_inlegs), jin_perim(n_inlegs))
            allocate(jout_road(n_outlegs), jout_lane(n_outlegs), jout_perim(n_outlegs))
            allocate(j_prob(n_routes))
            call read_int_1d(ncid, 'junction_in_road',     jin_road)
            call read_int_1d(ncid, 'junction_in_lane',     jin_lane)
            call read_int_1d(ncid, 'junction_in_perim',    jin_perim)
            call read_int_1d(ncid, 'junction_out_road',    jout_road)
            call read_int_1d(ncid, 'junction_out_lane',    jout_lane)
            call read_int_1d(ncid, 'junction_out_perim',   jout_perim)
            call read_real_1d(ncid, 'junction_route_prob', j_prob)
        end if

        call check( nf90_close(ncid) )

        ! ----- Populate spec -----
        allocate(spec%roads(n_roads))
        do r = 1, n_roads
            spec%roads(r)%id              = road_id(r)
            spec%roads(r)%end_junction(1) = road_end_junction(1, r)
            spec%roads(r)%end_junction(2) = road_end_junction(2, r)
            block
                integer :: n_lanes_r, kk
                n_lanes_r = count(lane_road_id == road_id(r))
                allocate(spec%roads(r)%lanes(n_lanes_r))
                kk = 0
                do k = 1, n_lanes
                    if (lane_road_id(k) /= road_id(r)) cycle
                    kk = kk + 1
                    spec%roads(r)%lanes(kk)%length         = lane_length(k)
                    spec%roads(r)%lanes(kk)%flow_direction = lane_fd(k)
                    spec%roads(r)%lanes(kk)%alpha          = lane_alpha(k)
                    spec%roads(r)%lanes(kk)%beta           = lane_beta(k)
                    spec%roads(r)%lanes(kk)%open_in        = (lane_open_in(k)  /= 0_1)
                    spec%roads(r)%lanes(kk)%open_out       = (lane_open_out(k) /= 0_1)
                end do
            end block
        end do

        allocate(spec%junctions(n_juncs))
        do j = 1, n_juncs
            spec%junctions(j)%id    = j_id(j)
            spec%junctions(j)%n_in  = j_nin(j)
            spec%junctions(j)%n_out = j_nout(j)
            allocate(spec%junctions(j)%in_legs(j_nin(j)))
            allocate(spec%junctions(j)%out_legs(j_nout(j)))
            allocate(spec%junctions(j)%routes(j_nin(j)))
            off = j_inoff(j)
            do k = 1, j_nin(j)
                spec%junctions(j)%in_legs(k)%road  = jin_road(off + k)
                spec%junctions(j)%in_legs(k)%lane  = jin_lane(off + k)
                spec%junctions(j)%in_legs(k)%perim = jin_perim(off + k)
            end do
            off = j_outoff(j)
            do k = 1, j_nout(j)
                spec%junctions(j)%out_legs(k)%road  = jout_road(off + k)
                spec%junctions(j)%out_legs(k)%lane  = jout_lane(off + k)
                spec%junctions(j)%out_legs(k)%perim = jout_perim(off + k)
            end do
            off = j_routeoff(j)
            do k = 1, j_nin(j)
                allocate(spec%junctions(j)%routes(k)%prob(j_nout(j)))
                spec%junctions(j)%routes(k)%prob = j_prob(off + 1 : off + j_nout(j))
                off = off + j_nout(j)
            end do
        end do
    end subroutine read_config

    ! ------------------------------------------------------------------
    ! Helpers
    ! ------------------------------------------------------------------

    subroutine get_att_int(ncid, name, dst)
        integer,          intent(in)    :: ncid
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: dst
        integer :: status
        status = nf90_get_att(ncid, NF90_GLOBAL, trim(name), dst)
        ! Missing optional attrs leave dst at its previous (default) value.
        if (status /= NF90_NOERR .and. status /= NF90_ENOTATT) &
            call check(status)
    end subroutine get_att_int

    subroutine get_att_real(ncid, name, dst)
        integer,          intent(in)    :: ncid
        character(len=*), intent(in)    :: name
        real,             intent(inout) :: dst
        integer :: status
        status = nf90_get_att(ncid, NF90_GLOBAL, trim(name), dst)
        if (status /= NF90_NOERR .and. status /= NF90_ENOTATT) &
            call check(status)
    end subroutine get_att_real

    subroutine get_att_string(ncid, name, dst)
        integer,          intent(in)    :: ncid
        character(len=*), intent(in)    :: name
        character(len=*), intent(inout) :: dst
        integer :: status
        character(len=64) :: buf
        buf = ''
        status = nf90_get_att(ncid, NF90_GLOBAL, trim(name), buf)
        if (status == NF90_NOERR) then
            dst = trim(buf)
        else if (status /= NF90_ENOTATT) then
            call check(status)
        end if
    end subroutine get_att_string

    subroutine read_int_1d(ncid, name, dst)
        integer,          intent(in)    :: ncid
        character(len=*), intent(in)    :: name
        integer,          intent(out)   :: dst(:)
        integer :: vid
        call check( nf90_inq_varid(ncid, trim(name), vid) )
        call check( nf90_get_var (ncid, vid, dst) )
    end subroutine read_int_1d

    subroutine read_int_2d(ncid, name, dst)
        integer,          intent(in)    :: ncid
        character(len=*), intent(in)    :: name
        integer,          intent(out)   :: dst(:,:)
        integer :: vid
        call check( nf90_inq_varid(ncid, trim(name), vid) )
        call check( nf90_get_var (ncid, vid, dst) )
    end subroutine read_int_2d

    subroutine read_real_1d(ncid, name, dst)
        integer,          intent(in)    :: ncid
        character(len=*), intent(in)    :: name
        real,             intent(out)   :: dst(:)
        integer :: vid
        call check( nf90_inq_varid(ncid, trim(name), vid) )
        call check( nf90_get_var (ncid, vid, dst) )
    end subroutine read_real_1d

    subroutine read_byte_1d(ncid, name, dst)
        integer,          intent(in)    :: ncid
        character(len=*), intent(in)    :: name
        integer(kind=1),  intent(out)   :: dst(:)
        integer :: vid
        call check( nf90_inq_varid(ncid, trim(name), vid) )
        call check( nf90_get_var (ncid, vid, dst) )
    end subroutine read_byte_1d

    subroutine check(status)
        integer, intent(in) :: status
        if (status /= NF90_NOERR) then
            print *, 'nc_config: NetCDF error: ', trim(nf90_strerror(status))
            stop 1
        end if
    end subroutine check

end module nc_config_mod
