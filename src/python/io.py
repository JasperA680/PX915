import netCDF4 as nc
import numpy as np


def load_netcdf(filename):
    """Load TASEP simulation output from a NetCDF file.

    Returns a dict with keys: L, n_steps, alpha, beta,
    history (L x n_steps int array), density (n_steps float array),
    current (n_steps int array).
    """
    with nc.Dataset(filename, 'r') as ds:
        data = {
            'L':        int(ds.getncattr('L')),
            'n_steps':  int(ds.getncattr('n_steps')),
            'alpha':    float(ds.getncattr('alpha')),
            'beta':     float(ds.getncattr('beta')),
            'history':  np.array(ds.variables['history'][:]).T,  # Fortran writes (site,time); Python reads (time,site) — transpose back
            'density':  np.array(ds.variables['density'][:]),
            'current':  np.array(ds.variables['current'][:]),
        }
    return data
