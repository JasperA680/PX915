module tasep_io
    use netcdf
    implicit none
    private

    public :: write_netcdf

contains

    subroutine write_netcdf(filename, L, n_steps, alpha, beta, &
                            history, density_history, current_history)
        character(len=*), intent(in) :: filename
        integer,          intent(in) :: L, n_steps
        real,             intent(in) :: alpha, beta
        integer,          intent(in) :: history(L, n_steps)
        real,             intent(in) :: density_history(n_steps)
        integer,          intent(in) :: current_history(n_steps)

        integer :: ncid, dim_site, dim_time
        integer :: varid_history, varid_density, varid_current

        call check( nf90_create(filename, NF90_CLOBBER, ncid) )

        ! Dimensions
        call check( nf90_def_dim(ncid, 'site',     L,       dim_site) )
        call check( nf90_def_dim(ncid, 'time',     n_steps, dim_time) )

        ! Global attributes
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'model',       '1D TASEP open boundary') )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'L',           L) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'n_steps',     n_steps) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'alpha',       alpha) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'beta',        beta) )

        ! Variables
        call check( nf90_def_var(ncid, 'history', NF90_INT, &
                                 [dim_site, dim_time], varid_history) )
        call check( nf90_put_att(ncid, varid_history, 'long_name', &
                                 'lattice occupancy (0=empty, 1=occupied)') )

        call check( nf90_def_var(ncid, 'density', NF90_FLOAT, &
                                 [dim_time], varid_density) )
        call check( nf90_put_att(ncid, varid_density, 'long_name', 'mean density rho = N/L') )

        call check( nf90_def_var(ncid, 'current', NF90_INT, &
                                 [dim_time], varid_current) )
        call check( nf90_put_att(ncid, varid_current, 'long_name', 'particle exits per timestep') )

        call check( nf90_enddef(ncid) )

        ! Write data
        call check( nf90_put_var(ncid, varid_history, history) )
        call check( nf90_put_var(ncid, varid_density, density_history) )
        call check( nf90_put_var(ncid, varid_current, current_history) )

        call check( nf90_close(ncid) )

    end subroutine write_netcdf

    subroutine check(status)
        integer, intent(in) :: status
        if (status /= NF90_NOERR) then
            print *, 'NetCDF error: ', trim(nf90_strerror(status))
            stop 1
        end if
    end subroutine check

end module tasep_io
