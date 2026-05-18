module json_config_mod
    !--------------------------------------------------------------------
    ! Minimal JSON parser tuned to the network simulator's config schema.
    !
    ! Python emits the JSON with sort_keys=True and ensure_ascii=True so
    ! we avoid Unicode escapes and unusual float formats.  The parser
    ! handles: object, array, string, number (int or real), true/false,
    ! null, and arbitrary nesting via skip_value.  It does NOT handle
    ! string escape sequences beyond \"; the schema has no string fields
    ! that need them (the "layout" sub-object is read into a single
    ! string blob via raw extraction and stored verbatim).
    !
    ! Public entry point: read_config(path, spec, params).
    !--------------------------------------------------------------------
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
        character(len=16) :: model = "NS"   ! "NS" (Nagel-Schreckenberg) or "TASEP" (skeleton)
    end type sim_params_t

    public :: read_config, load_text_file

contains

    !------------------------------------------------------------------
    ! Top-level entry: read JSON config file, populate spec and params.
    !------------------------------------------------------------------
    subroutine read_config(path, spec, params)
        character(len=*),     intent(in)  :: path
        type(network_spec_t), intent(out) :: spec
        type(sim_params_t),   intent(out) :: params

        character(len=:), allocatable :: buf
        integer :: pos

        call load_text_file(path, buf)
        pos = 1
        call parse_top(buf, pos, spec, params)
    end subroutine read_config

    !------------------------------------------------------------------
    ! Read a whole text file into an allocatable string.
    !------------------------------------------------------------------
    subroutine load_text_file(path, buf)
        character(len=*),              intent(in)  :: path
        character(len=:), allocatable, intent(out) :: buf
        integer :: u, ios, nbytes

        open(newunit=u, file=path, access='stream', form='unformatted', &
             status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(*,*) "ERROR: cannot open config file: ", trim(path)
            stop 1
        end if
        inquire(unit=u, size=nbytes)
        allocate(character(len=nbytes) :: buf)
        if (nbytes > 0) read(u) buf
        close(u)
    end subroutine load_text_file

    !------------------------------------------------------------------
    ! Top-level object dispatch.
    !------------------------------------------------------------------
    subroutine parse_top(buf, pos, spec, params)
        character(len=*),     intent(in)    :: buf
        integer,              intent(inout) :: pos
        type(network_spec_t), intent(out)   :: spec
        type(sim_params_t),   intent(out)   :: params

        character(len=:), allocatable :: key
        logical :: first

        call skip_ws(buf, pos)
        call expect_char(buf, pos, '{')
        first = .true.
        do
            call skip_ws(buf, pos)
            if (peek(buf, pos) == '}') then
                pos = pos + 1
                exit
            end if
            if (.not. first) then
                call expect_char(buf, pos, ',')
                call skip_ws(buf, pos)
            end if
            first = .false.

            call parse_string(buf, pos, key)
            call skip_ws(buf, pos)
            call expect_char(buf, pos, ':')
            call skip_ws(buf, pos)

            select case (key)
            case ("schema_version")
                params%schema_version = parse_int(buf, pos)
            case ("params")
                call parse_params(buf, pos, params)
            case ("roads")
                call parse_roads_array(buf, pos, spec)
            case ("junctions")
                call parse_junctions_array(buf, pos, spec)
            case ("layout")
                call skip_value(buf, pos)   ! ignored by Fortran
            case default
                call skip_value(buf, pos)
            end select
        end do
    end subroutine parse_top

    !------------------------------------------------------------------
    ! "params": { ... }
    !------------------------------------------------------------------
    subroutine parse_params(buf, pos, params)
        character(len=*),  intent(in)    :: buf
        integer,           intent(inout) :: pos
        type(sim_params_t), intent(inout) :: params

        character(len=:), allocatable :: key, sval
        logical :: first

        call expect_char(buf, pos, '{')
        first = .true.
        do
            call skip_ws(buf, pos)
            if (peek(buf, pos) == '}') then
                pos = pos + 1
                exit
            end if
            if (.not. first) then
                call expect_char(buf, pos, ',')
                call skip_ws(buf, pos)
            end if
            first = .false.

            call parse_string(buf, pos, key)
            call skip_ws(buf, pos)
            call expect_char(buf, pos, ':')
            call skip_ws(buf, pos)

            select case (key)
            case ("n_steps");         params%n_steps         = parse_int(buf, pos)
            case ("dt");              params%dt              = parse_real(buf, pos)
            case ("v_max");           params%v_max           = parse_int(buf, pos)
            case ("p_slow");          params%p_slow          = parse_real(buf, pos)
            case ("rng_seed");        params%rng_seed        = parse_int(buf, pos)
            case ("max_lane_length"); params%max_lane_length = parse_int(buf, pos)
            case ("model")
                call parse_string(buf, pos, sval)
                params%model = sval
            case default;             call skip_value(buf, pos)
            end select
        end do
    end subroutine parse_params

    !------------------------------------------------------------------
    ! "roads": [ {...}, {...}, ... ]
    !------------------------------------------------------------------
    subroutine parse_roads_array(buf, pos, spec)
        character(len=*),     intent(in)    :: buf
        integer,              intent(inout) :: pos
        type(network_spec_t), intent(inout) :: spec

        type(road_spec_t), allocatable :: tmp(:)
        integer :: n
        logical :: first

        n = 0
        allocate(tmp(0))

        call expect_char(buf, pos, '[')
        first = .true.
        do
            call skip_ws(buf, pos)
            if (peek(buf, pos) == ']') then
                pos = pos + 1
                exit
            end if
            if (.not. first) then
                call expect_char(buf, pos, ',')
                call skip_ws(buf, pos)
            end if
            first = .false.

            n = n + 1
            call grow_roads(tmp, n)
            call parse_road(buf, pos, tmp(n))
        end do
        spec%roads = tmp
    end subroutine parse_roads_array

    subroutine grow_roads(arr, new_size)
        type(road_spec_t), allocatable, intent(inout) :: arr(:)
        integer, intent(in) :: new_size
        type(road_spec_t), allocatable :: tmp(:)
        integer :: i, old_size
        old_size = size(arr)
        if (new_size <= old_size) return
        allocate(tmp(new_size))
        do i = 1, old_size
            tmp(i) = arr(i)
        end do
        call move_alloc(tmp, arr)
    end subroutine grow_roads

    subroutine parse_road(buf, pos, road)
        character(len=*),  intent(in)    :: buf
        integer,           intent(inout) :: pos
        type(road_spec_t), intent(out)   :: road

        character(len=:), allocatable :: key
        logical :: first

        call expect_char(buf, pos, '{')
        first = .true.
        do
            call skip_ws(buf, pos)
            if (peek(buf, pos) == '}') then
                pos = pos + 1
                exit
            end if
            if (.not. first) then
                call expect_char(buf, pos, ',')
                call skip_ws(buf, pos)
            end if
            first = .false.

            call parse_string(buf, pos, key)
            call skip_ws(buf, pos)
            call expect_char(buf, pos, ':')
            call skip_ws(buf, pos)

            select case (key)
            case ("id");
                road%id = parse_int(buf, pos)
            case ("end_junction")
                call parse_int_pair(buf, pos, road%end_junction)
            case ("lanes")
                call parse_lanes_array(buf, pos, road)
            case default
                call skip_value(buf, pos)
            end select
        end do
    end subroutine parse_road

    subroutine parse_int_pair(buf, pos, pair)
        character(len=*), intent(in)    :: buf
        integer,          intent(inout) :: pos
        integer,          intent(out)   :: pair(2)
        call expect_char(buf, pos, '[')
        call skip_ws(buf, pos)
        pair(1) = parse_int(buf, pos)
        call skip_ws(buf, pos)
        call expect_char(buf, pos, ',')
        call skip_ws(buf, pos)
        pair(2) = parse_int(buf, pos)
        call skip_ws(buf, pos)
        call expect_char(buf, pos, ']')
    end subroutine parse_int_pair

    subroutine parse_lanes_array(buf, pos, road)
        character(len=*),  intent(in)    :: buf
        integer,           intent(inout) :: pos
        type(road_spec_t), intent(inout) :: road

        type(lane_spec_t), allocatable :: tmp(:)
        integer :: n, i, old_size
        type(lane_spec_t), allocatable :: tmp2(:)
        logical :: first

        n = 0
        allocate(tmp(0))

        call expect_char(buf, pos, '[')
        first = .true.
        do
            call skip_ws(buf, pos)
            if (peek(buf, pos) == ']') then
                pos = pos + 1
                exit
            end if
            if (.not. first) then
                call expect_char(buf, pos, ',')
                call skip_ws(buf, pos)
            end if
            first = .false.

            n = n + 1
            old_size = size(tmp)
            allocate(tmp2(n))
            do i = 1, old_size
                tmp2(i) = tmp(i)
            end do
            call move_alloc(tmp2, tmp)
            call parse_lane(buf, pos, tmp(n))
        end do
        road%lanes = tmp
    end subroutine parse_lanes_array

    subroutine parse_lane(buf, pos, lane)
        character(len=*),  intent(in)    :: buf
        integer,           intent(inout) :: pos
        type(lane_spec_t), intent(out)   :: lane

        character(len=:), allocatable :: key
        logical :: first

        call expect_char(buf, pos, '{')
        first = .true.
        do
            call skip_ws(buf, pos)
            if (peek(buf, pos) == '}') then
                pos = pos + 1
                exit
            end if
            if (.not. first) then
                call expect_char(buf, pos, ',')
                call skip_ws(buf, pos)
            end if
            first = .false.

            call parse_string(buf, pos, key)
            call skip_ws(buf, pos)
            call expect_char(buf, pos, ':')
            call skip_ws(buf, pos)

            select case (key)
            case ("length");         lane%length         = parse_int(buf, pos)
            case ("flow_direction"); lane%flow_direction = parse_int(buf, pos)
            case ("alpha");          lane%alpha          = parse_real(buf, pos)
            case ("beta");           lane%beta           = parse_real(buf, pos)
            case ("open_in");        lane%open_in        = parse_bool(buf, pos)
            case ("open_out");       lane%open_out       = parse_bool(buf, pos)
            case default;            call skip_value(buf, pos)
            end select
        end do
    end subroutine parse_lane

    !------------------------------------------------------------------
    ! "junctions": [ {...}, ... ]
    !------------------------------------------------------------------
    subroutine parse_junctions_array(buf, pos, spec)
        character(len=*),     intent(in)    :: buf
        integer,              intent(inout) :: pos
        type(network_spec_t), intent(inout) :: spec

        type(junc_spec_t), allocatable :: tmp(:), tmp2(:)
        integer :: n, i, old_size
        logical :: first

        n = 0
        allocate(tmp(0))

        call expect_char(buf, pos, '[')
        first = .true.
        do
            call skip_ws(buf, pos)
            if (peek(buf, pos) == ']') then
                pos = pos + 1
                exit
            end if
            if (.not. first) then
                call expect_char(buf, pos, ',')
                call skip_ws(buf, pos)
            end if
            first = .false.

            n = n + 1
            old_size = size(tmp)
            allocate(tmp2(n))
            do i = 1, old_size
                tmp2(i) = tmp(i)
            end do
            call move_alloc(tmp2, tmp)
            call parse_junction(buf, pos, tmp(n))
        end do
        spec%junctions = tmp
    end subroutine parse_junctions_array

    subroutine parse_junction(buf, pos, junc)
        character(len=*),  intent(in)    :: buf
        integer,           intent(inout) :: pos
        type(junc_spec_t), intent(out)   :: junc

        character(len=:), allocatable :: key
        logical :: first

        call expect_char(buf, pos, '{')
        first = .true.
        do
            call skip_ws(buf, pos)
            if (peek(buf, pos) == '}') then
                pos = pos + 1
                exit
            end if
            if (.not. first) then
                call expect_char(buf, pos, ',')
                call skip_ws(buf, pos)
            end if
            first = .false.

            call parse_string(buf, pos, key)
            call skip_ws(buf, pos)
            call expect_char(buf, pos, ':')
            call skip_ws(buf, pos)

            select case (key)
            case ("id");        junc%id    = parse_int(buf, pos)
            case ("n_in");      junc%n_in  = parse_int(buf, pos)
            case ("n_out");     junc%n_out = parse_int(buf, pos)
            case ("in_legs");   call parse_legs_array(buf, pos, junc%in_legs)
            case ("out_legs");  call parse_legs_array(buf, pos, junc%out_legs)
            case ("routes");    call parse_routes_array(buf, pos, junc%routes)
            case default;       call skip_value(buf, pos)
            end select
        end do
    end subroutine parse_junction

    subroutine parse_legs_array(buf, pos, legs)
        character(len=*),               intent(in)    :: buf
        integer,                        intent(inout) :: pos
        type(leg_spec_t), allocatable,  intent(out)   :: legs(:)

        type(leg_spec_t), allocatable :: tmp(:), tmp2(:)
        integer :: n, i, old_size
        logical :: first

        n = 0
        allocate(tmp(0))

        call expect_char(buf, pos, '[')
        first = .true.
        do
            call skip_ws(buf, pos)
            if (peek(buf, pos) == ']') then
                pos = pos + 1
                exit
            end if
            if (.not. first) then
                call expect_char(buf, pos, ',')
                call skip_ws(buf, pos)
            end if
            first = .false.

            n = n + 1
            old_size = size(tmp)
            allocate(tmp2(n))
            do i = 1, old_size
                tmp2(i) = tmp(i)
            end do
            call move_alloc(tmp2, tmp)
            call parse_leg(buf, pos, tmp(n))
        end do
        legs = tmp
    end subroutine parse_legs_array

    subroutine parse_leg(buf, pos, leg)
        character(len=*), intent(in)    :: buf
        integer,          intent(inout) :: pos
        type(leg_spec_t), intent(out)   :: leg

        character(len=:), allocatable :: key
        logical :: first

        leg%perim = -1
        call expect_char(buf, pos, '{')
        first = .true.
        do
            call skip_ws(buf, pos)
            if (peek(buf, pos) == '}') then
                pos = pos + 1
                exit
            end if
            if (.not. first) then
                call expect_char(buf, pos, ',')
                call skip_ws(buf, pos)
            end if
            first = .false.

            call parse_string(buf, pos, key)
            call skip_ws(buf, pos)
            call expect_char(buf, pos, ':')
            call skip_ws(buf, pos)

            select case (key)
            case ("road");  leg%road  = parse_int(buf, pos)
            case ("lane");  leg%lane  = parse_int(buf, pos)
            case ("perim"); leg%perim = parse_int(buf, pos)
            case default;   call skip_value(buf, pos)
            end select
        end do
    end subroutine parse_leg

    subroutine parse_routes_array(buf, pos, routes)
        character(len=*),                intent(in)    :: buf
        integer,                         intent(inout) :: pos
        type(route_row_t), allocatable,  intent(out)   :: routes(:)

        type(route_row_t), allocatable :: tmp(:), tmp2(:)
        integer :: n, i, old_size, in_leg_idx
        real, allocatable :: prob(:)
        logical :: first

        n = 0
        allocate(tmp(0))

        call expect_char(buf, pos, '[')
        first = .true.
        do
            call skip_ws(buf, pos)
            if (peek(buf, pos) == ']') then
                pos = pos + 1
                exit
            end if
            if (.not. first) then
                call expect_char(buf, pos, ',')
                call skip_ws(buf, pos)
            end if
            first = .false.

            call parse_route_obj(buf, pos, in_leg_idx, prob)
            n = max(n, in_leg_idx)
            old_size = size(tmp)
            if (n > old_size) then
                allocate(tmp2(n))
                do i = 1, old_size
                    tmp2(i) = tmp(i)
                end do
                call move_alloc(tmp2, tmp)
            end if
            tmp(in_leg_idx)%prob = prob
        end do
        routes = tmp
    end subroutine parse_routes_array

    subroutine parse_route_obj(buf, pos, in_leg_idx, prob)
        character(len=*),               intent(in)    :: buf
        integer,                        intent(inout) :: pos
        integer,                        intent(out)   :: in_leg_idx
        real, allocatable,              intent(out)   :: prob(:)

        character(len=:), allocatable :: key
        logical :: first

        in_leg_idx = 0
        call expect_char(buf, pos, '{')
        first = .true.
        do
            call skip_ws(buf, pos)
            if (peek(buf, pos) == '}') then
                pos = pos + 1
                exit
            end if
            if (.not. first) then
                call expect_char(buf, pos, ',')
                call skip_ws(buf, pos)
            end if
            first = .false.

            call parse_string(buf, pos, key)
            call skip_ws(buf, pos)
            call expect_char(buf, pos, ':')
            call skip_ws(buf, pos)

            select case (key)
            case ("in_leg"); in_leg_idx = parse_int(buf, pos)
            case ("prob");   call parse_real_array(buf, pos, prob)
            case default;    call skip_value(buf, pos)
            end select
        end do
    end subroutine parse_route_obj

    subroutine parse_real_array(buf, pos, arr)
        character(len=*),  intent(in)    :: buf
        integer,           intent(inout) :: pos
        real, allocatable, intent(out)   :: arr(:)

        real, allocatable :: tmp(:), tmp2(:)
        integer :: n, i, old_size
        logical :: first

        n = 0
        allocate(tmp(0))

        call expect_char(buf, pos, '[')
        first = .true.
        do
            call skip_ws(buf, pos)
            if (peek(buf, pos) == ']') then
                pos = pos + 1
                exit
            end if
            if (.not. first) then
                call expect_char(buf, pos, ',')
                call skip_ws(buf, pos)
            end if
            first = .false.

            n = n + 1
            old_size = size(tmp)
            allocate(tmp2(n))
            do i = 1, old_size
                tmp2(i) = tmp(i)
            end do
            call move_alloc(tmp2, tmp)
            tmp(n) = parse_real(buf, pos)
        end do
        arr = tmp
    end subroutine parse_real_array

    !------------------------------------------------------------------
    ! Low-level lexer primitives.
    !------------------------------------------------------------------
    subroutine skip_ws(buf, pos)
        character(len=*), intent(in)    :: buf
        integer,          intent(inout) :: pos
        do while (pos <= len(buf))
            select case (buf(pos:pos))
            case (' ', achar(9), achar(10), achar(13))
                pos = pos + 1
            case default
                return
            end select
        end do
    end subroutine skip_ws

    pure function peek(buf, pos) result(c)
        character(len=*), intent(in) :: buf
        integer,          intent(in) :: pos
        character(len=1) :: c
        if (pos > len(buf)) then
            c = achar(0)
        else
            c = buf(pos:pos)
        end if
    end function peek

    subroutine expect_char(buf, pos, c)
        character(len=*), intent(in)    :: buf
        integer,          intent(inout) :: pos
        character(len=1), intent(in)    :: c
        call skip_ws(buf, pos)
        if (pos > len(buf) .or. buf(pos:pos) /= c) then
            write(*,'(a,a,a,i0,a,a)') "JSON parse error: expected '", c, &
                                       "' at pos ", pos, " got '", peek(buf, pos)
            stop 1
        end if
        pos = pos + 1
    end subroutine expect_char

    subroutine parse_string(buf, pos, s)
        character(len=*),              intent(in)    :: buf
        integer,                       intent(inout) :: pos
        character(len=:), allocatable, intent(out)   :: s
        integer :: start
        call skip_ws(buf, pos)
        if (peek(buf, pos) /= '"') then
            write(*,'(a,i0)') "JSON parse error: expected string at pos ", pos
            stop 1
        end if
        pos = pos + 1
        start = pos
        do while (pos <= len(buf) .and. buf(pos:pos) /= '"')
            ! Minimal: allow \" inside strings only.
            if (buf(pos:pos) == '\') then
                pos = pos + 2
            else
                pos = pos + 1
            end if
        end do
        s = buf(start:pos-1)
        pos = pos + 1
    end subroutine parse_string

    function parse_int(buf, pos) result(v)
        character(len=*), intent(in)    :: buf
        integer,          intent(inout) :: pos
        integer :: v, start, ios
        call skip_ws(buf, pos)
        start = pos
        if (peek(buf, pos) == '-' .or. peek(buf, pos) == '+') pos = pos + 1
        do while (pos <= len(buf))
            select case (buf(pos:pos))
            case ('0':'9'); pos = pos + 1
            case default;   exit
            end select
        end do
        read(buf(start:pos-1), *, iostat=ios) v
        if (ios /= 0) then
            write(*,'(a,a)') "JSON parse error: bad integer: ", buf(start:pos-1)
            stop 1
        end if
    end function parse_int

    function parse_real(buf, pos) result(v)
        character(len=*), intent(in)    :: buf
        integer,          intent(inout) :: pos
        real :: v
        integer :: start, ios
        call skip_ws(buf, pos)
        start = pos
        if (peek(buf, pos) == '-' .or. peek(buf, pos) == '+') pos = pos + 1
        do while (pos <= len(buf))
            select case (buf(pos:pos))
            case ('0':'9', '.', 'e', 'E', '+', '-'); pos = pos + 1
            case default; exit
            end select
        end do
        read(buf(start:pos-1), *, iostat=ios) v
        if (ios /= 0) then
            write(*,'(a,a)') "JSON parse error: bad number: ", buf(start:pos-1)
            stop 1
        end if
    end function parse_real

    function parse_bool(buf, pos) result(v)
        character(len=*), intent(in)    :: buf
        integer,          intent(inout) :: pos
        logical :: v
        call skip_ws(buf, pos)
        if (pos + 3 <= len(buf) .and. buf(pos:pos+3) == "true") then
            v = .true.
            pos = pos + 4
        else if (pos + 4 <= len(buf) .and. buf(pos:pos+4) == "false") then
            v = .false.
            pos = pos + 5
        else
            write(*,'(a,i0)') "JSON parse error: expected true/false at pos ", pos
            stop 1
        end if
    end function parse_bool

    !------------------------------------------------------------------
    ! Skip an arbitrary JSON value (object, array, primitive).
    !------------------------------------------------------------------
    recursive subroutine skip_value(buf, pos)
        character(len=*), intent(in)    :: buf
        integer,          intent(inout) :: pos
        character(len=1) :: c
        character(len=:), allocatable :: dummy
        integer :: depth
        logical :: in_string

        call skip_ws(buf, pos)
        c = peek(buf, pos)
        select case (c)
        case ('{', '[')
            depth = 0
            in_string = .false.
            do while (pos <= len(buf))
                c = buf(pos:pos)
                if (in_string) then
                    if (c == '\') then
                        pos = pos + 2
                        cycle
                    else if (c == '"') then
                        in_string = .false.
                    end if
                else
                    select case (c)
                    case ('"');             in_string = .true.
                    case ('{', '[');        depth = depth + 1
                    case ('}', ']')
                        depth = depth - 1
                        if (depth == 0) then
                            pos = pos + 1
                            return
                        end if
                    end select
                end if
                pos = pos + 1
            end do
        case ('"')
            call parse_string(buf, pos, dummy)
        case ('t', 'f')
            if (parse_bool(buf, pos)) then; end if
        case ('n')
            if (pos + 3 <= len(buf) .and. buf(pos:pos+3) == "null") then
                pos = pos + 4
            end if
        case default
            ! Treat as number (consume digits, sign, exponent, dot).
            do while (pos <= len(buf))
                select case (buf(pos:pos))
                case ('0':'9', '.', 'e', 'E', '+', '-'); pos = pos + 1
                case default; exit
                end select
            end do
        end select
    end subroutine skip_value

end module json_config_mod
