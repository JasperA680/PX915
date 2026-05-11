import json
from dataclasses import dataclass
from typing import Any, Dict, Optional

import netCDF4 as nc
import numpy as np


def load_netcdf(filename):
    """Load 1D TASEP simulation output from a NetCDF file.

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


@dataclass
class NetworkResult:
    """Loaded network simulation NetCDF.

    Array conventions follow the schema in network_io.f90.  Time-series
    arrays are time-major in Python: occupancy[t, lane, cell].
    """
    # History
    occupancy: np.ndarray      # uint8 (time, lane_index, cell)
    velocity:  np.ndarray      # int8  (time, lane_index, cell)
    final_occupancy: np.ndarray
    final_velocity:  np.ndarray

    # Per-lane lookup
    lane_road_id:       np.ndarray
    lane_within_road:   np.ndarray
    lane_length:        np.ndarray
    lane_flow_direction: np.ndarray
    lane_alpha:         np.ndarray
    lane_beta:          np.ndarray
    lane_open_in:       np.ndarray
    lane_open_out:      np.ndarray

    # Per-road
    road_end_junction:  np.ndarray   # (road, 2)
    road_density:       np.ndarray   # (time, road)
    road_entries:       np.ndarray
    road_exits:         np.ndarray

    # Junctions (flat with offsets)
    junction_n_in:        np.ndarray
    junction_n_out:       np.ndarray
    junction_inleg_offset:  np.ndarray
    junction_outleg_offset: np.ndarray
    junction_route_offset:  np.ndarray
    junction_in_road:    np.ndarray
    junction_in_lane:    np.ndarray
    junction_in_perim:   np.ndarray
    junction_out_road:   np.ndarray
    junction_out_lane:   np.ndarray
    junction_out_perim:  np.ndarray
    junction_route_prob: np.ndarray

    # Scalars
    n_steps: int
    v_max:   int
    p_slow:  float
    dt:      float
    rng_seed: int
    created_at: str
    schema_version: int
    config_json: str
    config: Dict[str, Any]    # parsed config_json (includes layout)


def load_network_netcdf(filename) -> NetworkResult:
    """Load a network simulation NetCDF file written by run_network.

    Time-series arrays come back time-major (time, lane_index, cell) etc.
    """
    with nc.Dataset(filename, 'r') as ds:
        config_json = str(ds.getncattr('config_json'))
        try:
            cfg = json.loads(config_json)
        except json.JSONDecodeError:
            cfg = {}

        v = ds.variables

        def arr(name):
            return np.array(v[name][:])

        result = NetworkResult(
            occupancy=arr('occupancy'),
            velocity=arr('velocity'),
            final_occupancy=arr('final_occupancy'),
            final_velocity=arr('final_velocity'),
            lane_road_id=arr('lane_road_id'),
            lane_within_road=arr('lane_within_road'),
            lane_length=arr('lane_length'),
            lane_flow_direction=arr('lane_flow_direction'),
            lane_alpha=arr('lane_alpha'),
            lane_beta=arr('lane_beta'),
            lane_open_in=arr('lane_open_in'),
            lane_open_out=arr('lane_open_out'),
            road_end_junction=arr('road_end_junction'),
            road_density=arr('road_density'),
            road_entries=arr('road_entries'),
            road_exits=arr('road_exits'),
            junction_n_in=arr('junction_n_in'),
            junction_n_out=arr('junction_n_out'),
            junction_inleg_offset=arr('junction_inleg_offset'),
            junction_outleg_offset=arr('junction_outleg_offset'),
            junction_route_offset=arr('junction_route_offset'),
            junction_in_road=arr('junction_in_road'),
            junction_in_lane=arr('junction_in_lane'),
            junction_in_perim=arr('junction_in_perim'),
            junction_out_road=arr('junction_out_road'),
            junction_out_lane=arr('junction_out_lane'),
            junction_out_perim=arr('junction_out_perim'),
            junction_route_prob=arr('junction_route_prob'),
            n_steps=int(ds.getncattr('n_steps')),
            v_max=int(ds.getncattr('v_max')),
            p_slow=float(ds.getncattr('p_slow')),
            dt=float(ds.getncattr('dt')),
            rng_seed=int(ds.getncattr('rng_seed')),
            created_at=str(ds.getncattr('created_at')),
            schema_version=int(ds.getncattr('schema_version')),
            config_json=config_json,
            config=cfg,
        )
    return result


def load_network_for_restart(filename) -> NetworkResult:
    """Load NC including final_occupancy + rng_seed for restart.

    Schema-ready entry point; the actual restart driver (--restart flag in
    the Fortran binary) is not yet implemented.  This stub exists so the
    GUI / scripts can wire to it once that's added.
    """
    return load_network_netcdf(filename)
