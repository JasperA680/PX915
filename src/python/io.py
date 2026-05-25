from dataclasses import dataclass
from pathlib import Path
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
    config: Dict[str, Any]    # reconstructed from the input config.nc (includes layout)


def load_network_netcdf(filename) -> NetworkResult:
    """Load a network simulation NetCDF file written by run_network.

    The companion ``config.nc`` is loaded from the same directory if present
    and used to populate ``result.config`` with the spec + layout (the layout
    coords aren't carried in result.nc itself).  When config.nc is missing
    (e.g. a legacy run) ``result.config`` falls back to an empty dict.
    """
    filename = Path(filename)
    with nc.Dataset(filename, 'r') as ds:
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
            config=_load_config_dict(filename.parent / "config.nc"),
        )
    return result


def _load_config_dict(config_path: Path) -> Dict[str, Any]:
    """Load a config.nc and reconstruct the layout/roads dict used by visualisation.

    Visualisation only reads ``layout`` (junction xy + road endpoints) and
    ``roads[*]['lanes'][*]['flow_direction']``, so that's all we materialise.
    """
    if not config_path.exists():
        return {}
    with nc.Dataset(config_path, 'r') as ds:
        v = ds.variables
        junction_ids = np.array(v['junction_id'][:]) if 'junction_id' in v else np.array([])
        junction_x = np.array(v['junction_x'][:]) if 'junction_x' in v else np.array([])
        junction_y = np.array(v['junction_y'][:]) if 'junction_y' in v else np.array([])
        road_ids   = np.array(v['road_id'][:])
        end1_x     = np.array(v['road_end1_x'][:])
        end1_y     = np.array(v['road_end1_y'][:])
        end2_x     = np.array(v['road_end2_x'][:])
        end2_y     = np.array(v['road_end2_y'][:])
        # Per-lane flow directions, indexed by lane_road_id.
        lane_road_id        = np.array(v['lane_road_id'][:])
        lane_flow_direction = np.array(v['lane_flow_direction'][:])

    roads = []
    for rid in road_ids:
        rid = int(rid)
        lane_mask = (lane_road_id == rid)
        fds = lane_flow_direction[lane_mask].tolist()
        roads.append({
            "id": rid,
            "lanes": [{"flow_direction": int(fd)} for fd in fds],
        })

    layout = {
        "junctions": [
            {"id": int(jid), "x": float(jx), "y": float(jy)}
            for jid, jx, jy in zip(junction_ids, junction_x, junction_y)
        ],
        "roads": [
            {"id": int(rid),
             "end_1": [float(x1), float(y1)],
             "end_2": [float(x2), float(y2)]}
            for rid, x1, y1, x2, y2 in zip(road_ids, end1_x, end1_y, end2_x, end2_y)
        ],
    }
    return {"layout": layout, "roads": roads}


def write_config_netcdf(spec, params, layout, path) -> Path:
    """Serialise a ``NetworkSpec`` + ``SimParams`` + ``LayoutSpec`` to a NetCDF
    config file the Fortran driver can ingest directly.

    Variable names match the schema used by ``network_io_mod`` in the Fortran
    output writer — so Fortran can read with identical naming.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    # Flatten per-lane and per-junction tables (matches network_io.f90 layout).
    n_roads = len(spec.roads)
    n_juncs = len(spec.junctions)
    road_ids = np.array([r.id for r in spec.roads], dtype=np.int32)
    road_end_junction = np.array([list(r.end_junction) for r in spec.roads], dtype=np.int32)

    lane_road_id, lane_within, lane_len, lane_fd = [], [], [], []
    lane_alpha, lane_beta, lane_open_in, lane_open_out = [], [], [], []
    for r in spec.roads:
        for k, ln in enumerate(r.lanes, start=1):
            lane_road_id.append(r.id)
            lane_within.append(k)
            lane_len.append(int(ln.length))
            lane_fd.append(int(ln.flow_direction))
            lane_alpha.append(float(ln.alpha))
            lane_beta.append(float(ln.beta))
            lane_open_in.append(1 if ln.open_in else 0)
            lane_open_out.append(1 if ln.open_out else 0)
    n_lanes = len(lane_len)

    j_nin, j_nout = [], []
    j_inoff, j_outoff, j_routeoff = [], [], []
    jin_road, jin_lane, jin_perim = [], [], []
    jout_road, jout_lane, jout_perim = [], [], []
    j_prob = []
    in_acc = 0; out_acc = 0; route_acc = 0
    for j in spec.junctions:
        j_nin.append(j.n_in); j_nout.append(j.n_out)
        j_inoff.append(in_acc); j_outoff.append(out_acc); j_routeoff.append(route_acc)
        for leg in j.in_legs:
            jin_road.append(leg.road); jin_lane.append(leg.lane); jin_perim.append(leg.perim)
        for leg in j.out_legs:
            jout_road.append(leg.road); jout_lane.append(leg.lane); jout_perim.append(leg.perim)
        for row in j.routes:
            j_prob.extend(float(v) for v in row)
        in_acc    += j.n_in
        out_acc   += j.n_out
        route_acc += j.n_in * j.n_out

    # Layout — pad to per-road / per-junction order matching road_ids / junction_ids.
    junction_ids = np.array([j.id for j in spec.junctions], dtype=np.int32)
    junction_x = np.array([layout.junctions.get(int(j.id), (0.0, 0.0))[0] for j in spec.junctions], dtype=np.float64)
    junction_y = np.array([layout.junctions.get(int(j.id), (0.0, 0.0))[1] for j in spec.junctions], dtype=np.float64)
    end1_x = np.zeros(n_roads, dtype=np.float64); end1_y = np.zeros(n_roads, dtype=np.float64)
    end2_x = np.zeros(n_roads, dtype=np.float64); end2_y = np.zeros(n_roads, dtype=np.float64)
    for i, r in enumerate(spec.roads):
        ep = layout.road_endpoints.get(int(r.id))
        if ep is not None:
            (x1, y1), (x2, y2) = ep
            end1_x[i] = x1; end1_y[i] = y1
            end2_x[i] = x2; end2_y[i] = y2

    with nc.Dataset(path, 'w', format='NETCDF4') as ds:
        # Dimensions.
        ds.createDimension('road',       n_roads)
        ds.createDimension('junction',   max(1, n_juncs))   # netCDF doesn't allow 0-sized
        ds.createDimension('pair',       2)
        ds.createDimension('lane_index', n_lanes)
        ds.createDimension('junction_inleg_total',  max(1, in_acc))
        ds.createDimension('junction_outleg_total', max(1, out_acc))
        ds.createDimension('route_total',           max(1, route_acc))

        # Scalar params as global attrs.
        ds.setncattr('schema_version',   1)
        ds.setncattr('n_steps',          int(params.n_steps))
        ds.setncattr('dt',               float(params.dt))
        ds.setncattr('v_max',            int(params.v_max))
        ds.setncattr('p_slow',           float(params.p_slow))
        ds.setncattr('rng_seed',         int(params.rng_seed))
        ds.setncattr('max_lane_length',  int(params.max_lane_length))
        ds.setncattr('model',            str(params.model))
        ds.setncattr('lc_model',         int(params.lc_model))
        ds.setncattr('lc_p_change',      float(params.lc_p_change))

        def _write(name, dim, arr, dtype):
            var = ds.createVariable(name, dtype, dim)
            var[:] = arr

        # Per-road.
        _write('road_id',           ('road',),         road_ids,            np.int32)
        # Fortran reads this into ``road_endj(2, n_roads)`` (col-major).  To
        # match that on disk we declare dims as (road, pair) in Python's
        # (row-major) view — NetCDF reverses dim order between conventions,
        # so Fortran sees (pair, road).
        _write('road_end_junction', ('road', 'pair'),  road_end_junction,   np.int32)
        _write('road_end1_x',       ('road',),         end1_x,              np.float64)
        _write('road_end1_y',       ('road',),         end1_y,              np.float64)
        _write('road_end2_x',       ('road',),         end2_x,              np.float64)
        _write('road_end2_y',       ('road',),         end2_y,              np.float64)

        # Per-lane.
        _write('lane_road_id',         ('lane_index',), lane_road_id, np.int32)
        _write('lane_within_road',     ('lane_index',), lane_within,  np.int32)
        _write('lane_length',          ('lane_index',), lane_len,     np.int32)
        _write('lane_flow_direction',  ('lane_index',), lane_fd,      np.int32)
        _write('lane_alpha',           ('lane_index',), lane_alpha,   np.float32)
        _write('lane_beta',            ('lane_index',), lane_beta,    np.float32)
        _write('lane_open_in',         ('lane_index',), lane_open_in, np.int8)
        _write('lane_open_out',        ('lane_index',), lane_open_out, np.int8)

        # Per-junction.  Even when n_juncs == 0 we still create dimensioned
        # vars of size 1 (with a sentinel) so the Fortran reader's
        # nf90_inq_dimid never NotFound's — easier than special-casing both
        # sides for empty networks.
        pad_int = lambda lst: lst if lst else [0]
        pad_real = lambda lst: lst if lst else [0.0]
        _write('junction_id',             ('junction',), pad_int(junction_ids.tolist()), np.int32)
        _write('junction_x',              ('junction',), pad_real(junction_x.tolist()),  np.float64)
        _write('junction_y',              ('junction',), pad_real(junction_y.tolist()),  np.float64)
        _write('junction_n_in',           ('junction',), pad_int(j_nin),                 np.int32)
        _write('junction_n_out',          ('junction',), pad_int(j_nout),                np.int32)
        _write('junction_inleg_offset',   ('junction',), pad_int(j_inoff),               np.int32)
        _write('junction_outleg_offset',  ('junction',), pad_int(j_outoff),              np.int32)
        _write('junction_route_offset',   ('junction',), pad_int(j_routeoff),            np.int32)
        _write('junction_in_road',        ('junction_inleg_total',),  pad_int(jin_road),   np.int32)
        _write('junction_in_lane',        ('junction_inleg_total',),  pad_int(jin_lane),   np.int32)
        _write('junction_in_perim',       ('junction_inleg_total',),  pad_int(jin_perim),  np.int32)
        _write('junction_out_road',       ('junction_outleg_total',), pad_int(jout_road),  np.int32)
        _write('junction_out_lane',       ('junction_outleg_total',), pad_int(jout_lane),  np.int32)
        _write('junction_out_perim',      ('junction_outleg_total',), pad_int(jout_perim), np.int32)
        _write('junction_route_prob',     ('route_total',),           pad_real(j_prob),    np.float32)

        # Sentinel telling Fortran which junction-related arrays are "empty"
        # (we always allocate >= 1 to satisfy NetCDF).
        ds.setncattr('n_junctions', int(n_juncs))

    return path


def load_network_for_restart(filename) -> NetworkResult:
    """Load NC including final_occupancy + rng_seed for restart.

    Schema-ready entry point; the actual restart driver (--restart flag in
    the Fortran binary) is not yet implemented.  This stub exists so the
    GUI / scripts can wire to it once that's added.
    """
    return load_network_netcdf(filename)
