module nc_config_mod
    ! NetCDF configuration reader for the road-network simulator.
    !
    ! Reads a ``config.nc`` file written by the Python frontend
    ! (``python.io.write_config_netcdf``) into a ``sim_params_t`` plus a
    ! ``network_spec_t`` ready for ``network_builder_mod::build_network`` to
    ! turn into a runtime ``road_network_t``. Variable and attribute names
    ! mirror the output-side schema used by ``network_io_mod``.
    !
    ! What's in a config NetCDF:
    !
    ! * **Global attributes** — the simulation-wide knobs: ``n_steps``,
    !   ``dt``, ``v_max``, ``p_slow``, ``rng_seed``, ``max_lane_length``,
    !   ``model``, ``lc_model``, ``lc_p_change``. Also ``n_junctions``
    !   (the real number of junctions; useful because the on-disk
    !   ``junction`` dimension is always at least one, so a network with
    !   no junctions still has a one-row sentinel in the tables below).
    ! * **Per-road**: ``road_id``, ``road_end_junction``.
    ! * **Per-lane**: ``lane_road_id``, ``lane_within_road``,
    !   ``lane_length``, ``lane_flow_direction``, ``lane_alpha``,
    !   ``lane_beta``, ``lane_open_in``, ``lane_open_out``, and the
    !   optional ``lane_is_periodic`` / ``lane_n_vehicles`` added when
    !   periodic-NS support was wired through the Python pipeline.
    ! * **Per-junction**: ``junction_id``, ``junction_n_in``,
    !   ``junction_n_out``, ``junction_inleg_offset``,
    !   ``junction_outleg_offset``, ``junction_route_offset``.
    ! * **Flattened junction tables**: ``junction_in_road``,
    !   ``junction_in_lane``, ``junction_in_perim``, ``junction_out_road``,
    !   ``junction_out_lane``, ``junction_out_perim``, ``junction_route_prob``.
    !
    ! The two optional per-lane fields are read via ``read_*_1d_opt``
    ! helpers that fall back to a caller-initialised default when the
    ! variable is absent, so NetCDFs written by older Python frontends
    ! keep loading without error.

    use netcdf
    use network_builder_mod

    implicit none
    private

    type, public :: sim_params_t
        ! Simulation-wide scalar parameters, populated from the config
        ! NetCDF's global attributes. Missing attributes leave the
        ! corresponding field at its default value below.

        integer :: schema_version = 1
        ! On-disk schema version of the config file.
        integer :: n_steps = 1000
        ! Number of simulation time steps to run.
        real :: dt = 1.0
        ! Time-step size (informational on the CA side; the update is unitless).
        integer :: v_max = 5
        ! Maximum vehicle speed for the Nagel-Schreckenberg model
        ! (ignored when ``model == 'TASEP'``).
        real :: p_slow = 0.2
        ! Random slowing probability for the Nagel-Schreckenberg model.
        integer :: rng_seed = 42
        ! Seed for the Fortran intrinsic RNG.
        integer :: max_lane_length = 0
        ! Pre-computed max lane length to size output history arrays;
        ! zero (default) tells the driver to compute it from the spec.
        character(len=16) :: model = "NS"
        ! Longitudinal model selector: ``"NS"`` or ``"TASEP"``.
        integer :: lc_model = -1
        ! Lane-change model: ``-1`` disabled, ``0`` symmetric, ``1`` asymmetric.
        real :: lc_p_change = 1.0
        ! Acceptance probability for a candidate lane change.
    end type sim_params_t

    public :: read_config

contains

    subroutine read_config(path, spec, params)
        ! Read a ``config.nc`` into a ``sim_params_t`` and a
        ! ``network_spec_t``. The optional lane fields ``lane_is_periodic``
        ! and ``lane_n_vehicles`` fall back to ``.false.`` / ``0`` when
        ! absent, so configs written by older Python frontends still load.
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
        integer(kind=1), allocatable :: lane_is_periodic(:)
        integer,         allocatable :: lane_n_vehicles(:)

        integer, allocatable :: j_id(:), j_nin(:), j_nout(:)
        integer, allocatable :: j_inoff(:), j_outoff(:), j_routeoff(:)
        integer, allocatable :: jin_road(:), jin_lane(:), jin_perim(:)
        integer, allocatable :: jout_road(:), jout_lane(:), jout_perim(:)
        real,    allocatable :: j_prob(:)

        integer :: r, j, k, off
        ! Lane-loop locals hoisted out of an inner BLOCK because the Sphinx
        ! crackfortran parser doesn't support F2008 BLOCK constructs.
        integer :: n_lanes_r, kk

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

        ! Optional lane fields — pre-initialised so legacy NetCDFs default
        ! to non-periodic, empty-initial-state.
        allocate(lane_is_periodic(n_lanes), lane_n_vehicles(n_lanes))
        lane_is_periodic = 0_1
        lane_n_vehicles  = 0
        call read_byte_1d_opt(ncid, 'lane_is_periodic', lane_is_periodic)
        call read_int_1d_opt (ncid, 'lane_n_vehicles',  lane_n_vehicles)

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
                spec%roads(r)%lanes(kk)%is_periodic    = (lane_is_periodic(k) /= 0_1)
                spec%roads(r)%lanes(kk)%n_vehicles     = lane_n_vehicles(k)
            end do
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
    ! Private NetCDF helpers — thin wrappers around the C library that add
    ! error reporting and (for the ``_opt`` variants) graceful handling of
    ! variables added in newer schema versions.
    ! ------------------------------------------------------------------

    subroutine get_att_int(ncid, name, dst)
        ! Read an optional integer global attribute; leave ``dst`` at its
        ! caller-initialised default if the attribute is absent.
        integer,          intent(in)    :: ncid
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: dst
        integer :: status
        status = nf90_get_att(ncid, NF90_GLOBAL, trim(name), dst)
        if (status /= NF90_NOERR .and. status /= NF90_ENOTATT) &
            call check(status)
    end subroutine get_att_int

    subroutine get_att_real(ncid, name, dst)
        ! Read an optional real global attribute; see ``get_att_int``.
        integer,          intent(in)    :: ncid
        character(len=*), intent(in)    :: name
        real,             intent(inout) :: dst
        integer :: status
        status = nf90_get_att(ncid, NF90_GLOBAL, trim(name), dst)
        if (status /= NF90_NOERR .and. status /= NF90_ENOTATT) &
            call check(status)
    end subroutine get_att_real

    subroutine get_att_string(ncid, name, dst)
        ! Read an optional string global attribute, trimmed flush-left.
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
        ! Read a required 1-D integer variable. Fatal error if absent.
        integer,          intent(in)    :: ncid
        character(len=*), intent(in)    :: name
        integer,          intent(out)   :: dst(:)
        integer :: vid
        call check( nf90_inq_varid(ncid, trim(name), vid) )
        call check( nf90_get_var (ncid, vid, dst) )
    end subroutine read_int_1d

    subroutine read_int_2d(ncid, name, dst)
        ! Read a required 2-D integer variable. NetCDF reverses dimension
        ! order between C and Fortran, so an on-disk ``(road, pair)``
        ! variable arrives as ``(pair, road)`` here.
        integer,          intent(in)    :: ncid
        character(len=*), intent(in)    :: name
        integer,          intent(out)   :: dst(:,:)
        integer :: vid
        call check( nf90_inq_varid(ncid, trim(name), vid) )
        call check( nf90_get_var (ncid, vid, dst) )
    end subroutine read_int_2d

    subroutine read_real_1d(ncid, name, dst)
        ! Read a required 1-D real variable. Fatal error if absent.
        integer,          intent(in)    :: ncid
        character(len=*), intent(in)    :: name
        real,             intent(out)   :: dst(:)
        integer :: vid
        call check( nf90_inq_varid(ncid, trim(name), vid) )
        call check( nf90_get_var (ncid, vid, dst) )
    end subroutine read_real_1d

    subroutine read_byte_1d(ncid, name, dst)
        ! Read a required 1-D ``integer(kind=1)`` variable (the boolean
        ! lane flags are stored as ``int8`` on disk). Fatal error if absent.
        integer,          intent(in)    :: ncid
        character(len=*), intent(in)    :: name
        integer(kind=1),  intent(out)   :: dst(:)
        integer :: vid
        call check( nf90_inq_varid(ncid, trim(name), vid) )
        call check( nf90_get_var (ncid, vid, dst) )
    end subroutine read_byte_1d

    subroutine read_byte_1d_opt(ncid, name, dst)
        ! Optional byte-array reader for fields added after schema v1:
        ! if the variable is absent ``dst`` is left untouched, so callers
        ! should pre-initialise it to the desired default before this call.
        integer,          intent(in)    :: ncid
        character(len=*), intent(in)    :: name
        integer(kind=1),  intent(inout) :: dst(:)
        integer :: vid, status
        status = nf90_inq_varid(ncid, trim(name), vid)
        if (status == NF90_ENOTVAR) return
        if (status /= NF90_NOERR) call check(status)
        call check( nf90_get_var(ncid, vid, dst) )
    end subroutine read_byte_1d_opt

    subroutine read_int_1d_opt(ncid, name, dst)
        ! Optional integer-array reader; see ``read_byte_1d_opt``.
        integer,          intent(in)    :: ncid
        character(len=*), intent(in)    :: name
        integer,          intent(inout) :: dst(:)
        integer :: vid, status
        status = nf90_inq_varid(ncid, trim(name), vid)
        if (status == NF90_ENOTVAR) return
        if (status /= NF90_NOERR) call check(status)
        call check( nf90_get_var(ncid, vid, dst) )
    end subroutine read_int_1d_opt

    subroutine check(status)
        ! Abort with a human-readable NetCDF error message on non-zero status.
        integer, intent(in) :: status
        if (status /= NF90_NOERR) then
            print *, 'nc_config: NetCDF error: ', trim(nf90_strerror(status))
            stop 1
        end if
    end subroutine check

end module nc_config_mod
