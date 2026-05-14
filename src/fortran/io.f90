module tasep_io
    use netcdf
    implicit none
    private

    public :: write_netcdf
    public :: write_fundamental_diagram_netcdf
    public :: write_benchmark_netcdf

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

    subroutine write_fundamental_diagram_netcdf(filename, L, n_points, n_burnin, n_measure, &
                                                 fixed_beta, fixed_alpha, &
                                                 alpha_param, alpha_rho, alpha_J, &
                                                 beta_param,  beta_rho,  beta_J)
        character(len=*), intent(in) :: filename
        integer,          intent(in) :: L, n_points, n_burnin, n_measure
        real,             intent(in) :: fixed_beta, fixed_alpha
        real,             intent(in) :: alpha_param(n_points), alpha_rho(n_points), alpha_J(n_points)
        real,             intent(in) :: beta_param(n_points),  beta_rho(n_points),  beta_J(n_points)

        integer :: ncid, dim_pts
        integer :: vid_ap, vid_ar, vid_aj, vid_bp, vid_br, vid_bj

        call check( nf90_create(filename, NF90_CLOBBER, ncid) )

        call check( nf90_def_dim(ncid, 'n_points', n_points, dim_pts) )

        call check( nf90_put_att(ncid, NF90_GLOBAL, 'model',       '1D TASEP fundamental diagram') )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'L',           L) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'n_points',    n_points) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'n_burnin',    n_burnin) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'n_measure',   n_measure) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'fixed_beta',  fixed_beta) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'fixed_alpha', fixed_alpha) )

        call check( nf90_def_var(ncid, 'alpha_param', NF90_FLOAT, [dim_pts], vid_ap) )
        call check( nf90_put_att(ncid, vid_ap, 'long_name', 'alpha values (beta sweep fixed at fixed_beta)') )

        call check( nf90_def_var(ncid, 'alpha_rho', NF90_FLOAT, [dim_pts], vid_ar) )
        call check( nf90_put_att(ncid, vid_ar, 'long_name', 'mean bulk density (alpha sweep)') )

        call check( nf90_def_var(ncid, 'alpha_J', NF90_FLOAT, [dim_pts], vid_aj) )
        call check( nf90_put_att(ncid, vid_aj, 'long_name', 'mean current (alpha sweep)') )

        call check( nf90_def_var(ncid, 'beta_param', NF90_FLOAT, [dim_pts], vid_bp) )
        call check( nf90_put_att(ncid, vid_bp, 'long_name', 'beta values (alpha sweep fixed at fixed_alpha)') )

        call check( nf90_def_var(ncid, 'beta_rho', NF90_FLOAT, [dim_pts], vid_br) )
        call check( nf90_put_att(ncid, vid_br, 'long_name', 'mean bulk density (beta sweep)') )

        call check( nf90_def_var(ncid, 'beta_J', NF90_FLOAT, [dim_pts], vid_bj) )
        call check( nf90_put_att(ncid, vid_bj, 'long_name', 'mean current (beta sweep)') )

        call check( nf90_enddef(ncid) )

        call check( nf90_put_var(ncid, vid_ap, alpha_param) )
        call check( nf90_put_var(ncid, vid_ar, alpha_rho) )
        call check( nf90_put_var(ncid, vid_aj, alpha_J) )
        call check( nf90_put_var(ncid, vid_bp, beta_param) )
        call check( nf90_put_var(ncid, vid_br, beta_rho) )
        call check( nf90_put_var(ncid, vid_bj, beta_J) )

        call check( nf90_close(ncid) )

    end subroutine write_fundamental_diagram_netcdf

    subroutine write_benchmark_netcdf(filename, n_L, n_N, n_rep, &
                                      alpha, beta, &
                                      L_values, N_values, &
                                      wallclock, steady_density, steady_current)
        character(len=*), intent(in) :: filename
        integer,          intent(in) :: n_L, n_N, n_rep
        real,             intent(in) :: alpha, beta
        integer,          intent(in) :: L_values(n_L)
        integer,          intent(in) :: N_values(n_N)
        real,             intent(in) :: wallclock(n_L, n_N, n_rep)
        real,             intent(in) :: steady_density(n_L)
        real,             intent(in) :: steady_current(n_L)

        integer :: ncid, dim_L, dim_N, dim_R
        integer :: vid_L, vid_N, vid_w, vid_rho, vid_J

        call check( nf90_create(filename, NF90_CLOBBER, ncid) )

        call check( nf90_def_dim(ncid, 'n_L',   n_L,   dim_L) )
        call check( nf90_def_dim(ncid, 'n_N',   n_N,   dim_N) )
        call check( nf90_def_dim(ncid, 'n_rep', n_rep, dim_R) )

        call check( nf90_put_att(ncid, NF90_GLOBAL, 'model',    '1D TASEP benchmark sweep') )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'alpha',    alpha) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'beta',     beta) )
        call check( nf90_put_att(ncid, NF90_GLOBAL, 'compiler', 'gfortran -O2') )

        call check( nf90_def_var(ncid, 'L_values', NF90_INT, [dim_L], vid_L) )
        call check( nf90_put_att(ncid, vid_L, 'long_name', 'lattice sizes swept') )

        call check( nf90_def_var(ncid, 'N_values', NF90_INT, [dim_N], vid_N) )
        call check( nf90_put_att(ncid, vid_N, 'long_name', 'n_steps values swept') )

        call check( nf90_def_var(ncid, 'wallclock', NF90_FLOAT, [dim_L, dim_N, dim_R], vid_w) )
        call check( nf90_put_att(ncid, vid_w, 'long_name', 'wallclock seconds per run') )
        call check( nf90_put_att(ncid, vid_w, 'units', 'seconds') )

        call check( nf90_def_var(ncid, 'steady_density', NF90_FLOAT, [dim_L], vid_rho) )
        call check( nf90_put_att(ncid, vid_rho, 'long_name', 'mean steady-state bulk density') )

        call check( nf90_def_var(ncid, 'steady_current', NF90_FLOAT, [dim_L], vid_J) )
        call check( nf90_put_att(ncid, vid_J, 'long_name', 'mean steady-state current') )

        call check( nf90_enddef(ncid) )

        call check( nf90_put_var(ncid, vid_L,   L_values) )
        call check( nf90_put_var(ncid, vid_N,   N_values) )
        call check( nf90_put_var(ncid, vid_w,   wallclock) )
        call check( nf90_put_var(ncid, vid_rho, steady_density) )
        call check( nf90_put_var(ncid, vid_J,   steady_current) )

        call check( nf90_close(ncid) )

    end subroutine write_benchmark_netcdf

    subroutine check(status)
        integer, intent(in) :: status
        if (status /= NF90_NOERR) then
            print *, 'NetCDF error: ', trim(nf90_strerror(status))
            stop 1
        end if
    end subroutine check

end module tasep_io
