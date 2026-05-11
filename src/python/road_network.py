"""Python-side network spec and preset builders.

This module mirrors the Fortran derived types (lane_spec_t, road_spec_t,
junc_spec_t, network_spec_t) so a network can be authored in Python,
validated, serialised to the JSON config schema, and shipped off to the
Fortran driver `./build/run_network`.

Four presets are exposed:
    single_lane(...)   - one-road 1D TASEP baseline
    crossroads(...)    - 4-way junction, mirrors init_crossroad
    roundabout(...)    - 4 T-junctions in a clockwise ring
    town(...)          - 2x2 grid of crossroads + 8 external arms

Each preset returns (NetworkSpec, LayoutSpec).  The layout is ignored by
Fortran but stored verbatim in the NC `config_json` attribute so the GUI
can re-render the diagram on load.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Dict, List, Tuple


# ---------------------------------------------------------------------------
# Dataclasses (mirror Fortran network_builder_mod)
# ---------------------------------------------------------------------------

@dataclass
class LaneSpec:
    length: int
    flow_direction: int       # +1 (end_1 -> end_2) or -1
    alpha: float = 0.0
    beta: float = 0.0
    open_in: bool = False
    open_out: bool = False


@dataclass
class RoadSpec:
    id: int
    end_junction: Tuple[int, int]   # 0 == open boundary
    lanes: List[LaneSpec]


@dataclass
class JunctionLegSpec:
    road: int
    lane: int                       # 1-indexed within road.lanes
    perim: int = -1                 # -1 = symmetric (no perim)


@dataclass
class JunctionSpec:
    id: int
    n_in: int
    n_out: int
    in_legs: List[JunctionLegSpec]
    out_legs: List[JunctionLegSpec]
    routes: List[List[float]]       # len n_in; each row len n_out, sums to 1


@dataclass
class LayoutSpec:
    junctions: Dict[int, Tuple[float, float]] = field(default_factory=dict)
    road_endpoints: Dict[int, Tuple[Tuple[float, float], Tuple[float, float]]] = field(default_factory=dict)


@dataclass
class NetworkSpec:
    roads: List[RoadSpec]
    junctions: List[JunctionSpec]


@dataclass
class SimParams:
    n_steps: int = 2000
    dt: float = 1.0
    v_max: int = 5
    p_slow: float = 0.2
    rng_seed: int = 42
    max_lane_length: int = 0        # 0 -> auto from lanes


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _lane_inbound_at(road: RoadSpec, junc_end: int) -> int:
    """Return 1-indexed lane index inbound to a junction at `junc_end` (1 or 2)."""
    target_fd = 2 * junc_end - 3
    for k, ln in enumerate(road.lanes, start=1):
        if ln.flow_direction == target_fd:
            return k
    raise ValueError(f"road {road.id}: no inbound lane at end {junc_end}")


def _lane_outbound_at(road: RoadSpec, junc_end: int) -> int:
    """Return 1-indexed lane index outbound from a junction at `junc_end`."""
    target_fd = 3 - 2 * junc_end
    for k, ln in enumerate(road.lanes, start=1):
        if ln.flow_direction == target_fd:
            return k
    raise ValueError(f"road {road.id}: no outbound lane at end {junc_end}")


def _bidir_lanes(L: int) -> List[LaneSpec]:
    return [LaneSpec(L, +1), LaneSpec(L, -1)]


# ---------------------------------------------------------------------------
# Validation + serialisation
# ---------------------------------------------------------------------------

def validate(spec: NetworkSpec) -> None:
    road_ids = [r.id for r in spec.roads]
    if len(set(road_ids)) != len(road_ids):
        raise ValueError("duplicate road ids")
    junc_ids = [j.id for j in spec.junctions]
    if len(set(junc_ids)) != len(junc_ids):
        raise ValueError("duplicate junction ids")

    road_by_id = {r.id: r for r in spec.roads}

    for r in spec.roads:
        if not (1 <= len(r.lanes) <= 2):
            raise ValueError(f"road {r.id}: lanes per road must be 1 or 2")
        if len(r.lanes) == 2 and r.lanes[0].flow_direction == r.lanes[1].flow_direction:
            raise ValueError(f"road {r.id}: bidirectional road must have distinct flow_directions")

    for j in spec.junctions:
        if len(j.in_legs) != j.n_in:
            raise ValueError(f"junction {j.id}: in_legs length != n_in")
        if len(j.out_legs) != j.n_out:
            raise ValueError(f"junction {j.id}: out_legs length != n_out")
        if len(j.routes) != j.n_in:
            raise ValueError(f"junction {j.id}: routes length != n_in")
        for k, row in enumerate(j.routes, start=1):
            if len(row) != j.n_out:
                raise ValueError(f"junction {j.id} leg {k}: route row length != n_out")
            s = sum(row)
            if abs(s - 1.0) > 1e-6:
                raise ValueError(f"junction {j.id} leg {k}: route probs sum to {s}, expected 1.0")

        for kind, legs in (("in", j.in_legs), ("out", j.out_legs)):
            for leg in legs:
                if leg.road not in road_by_id:
                    raise ValueError(f"junction {j.id} {kind}_leg: road {leg.road} not found")
                if not (1 <= leg.lane <= len(road_by_id[leg.road].lanes)):
                    raise ValueError(f"junction {j.id} {kind}_leg: lane {leg.lane} out of range for road {leg.road}")

        perims = [leg.perim for leg in j.in_legs] + [leg.perim for leg in j.out_legs]
        any_set = any(p >= 0 for p in perims)
        all_set = all(p >= 0 for p in perims)
        if any_set and not all_set:
            raise ValueError(f"junction {j.id}: perim values must be all-or-none")
        if all_set and len(set(perims)) != len(perims):
            raise ValueError(f"junction {j.id}: perim values must be distinct")


def spec_to_dict(spec: NetworkSpec, params: SimParams, layout: LayoutSpec) -> dict:
    """Serialise spec + params + layout to a JSON-ready dict."""
    return {
        "schema_version": 1,
        "params": asdict(params),
        "roads": [
            {
                "id": r.id,
                "end_junction": list(r.end_junction),
                "lanes": [asdict(ln) for ln in r.lanes],
            }
            for r in spec.roads
        ],
        "junctions": [
            {
                "id": j.id,
                "n_in": j.n_in,
                "n_out": j.n_out,
                "in_legs":  [asdict(leg) for leg in j.in_legs],
                "out_legs": [asdict(leg) for leg in j.out_legs],
                "routes":   [{"in_leg": k, "prob": list(row)}
                             for k, row in enumerate(j.routes, start=1)],
            }
            for j in spec.junctions
        ],
        "layout": {
            "junctions": [{"id": jid, "x": x, "y": y}
                          for jid, (x, y) in layout.junctions.items()],
            "roads":     [{"id": rid, "end_1": list(e1), "end_2": list(e2)}
                          for rid, (e1, e2) in layout.road_endpoints.items()],
        },
    }


# ---------------------------------------------------------------------------
# Presets
# ---------------------------------------------------------------------------

def single_lane(L: int = 50, alpha: float = 0.5, beta: float = 0.5) -> Tuple[NetworkSpec, LayoutSpec]:
    """One road, one lane, open boundaries.  1D-TASEP baseline."""
    road = RoadSpec(
        id=1, end_junction=(0, 0),
        lanes=[LaneSpec(L, +1, alpha=alpha, beta=beta, open_in=True, open_out=True)],
    )
    spec = NetworkSpec(roads=[road], junctions=[])
    layout = LayoutSpec(road_endpoints={1: ((0.0, 0.0), (1.0, 0.0))})
    return spec, layout


def crossroads(L: int = 20, alpha: float = 0.4, beta: float = 0.5,
               p_left: float = 0.25, p_right: float = 0.25) -> Tuple[NetworkSpec, LayoutSpec]:
    """4-way crossroad mirroring Fortran init_crossroad.

    Roads in CW order: 1=N, 2=E, 3=S, 4=W.  All end_1 connect to junction 1.
    Each road has 2 lanes [+1, -1].  Lane 2 (fd=-1) is inbound (open_in),
    lane 1 (fd=+1) is outbound (open_out).  This matches init_crossroad.
    """
    p_straight = 1.0 - p_left - p_right
    roads: List[RoadSpec] = []
    for i in range(1, 5):
        r = RoadSpec(id=i, end_junction=(1, 0), lanes=_bidir_lanes(L))
        r.lanes[1].open_in = True
        r.lanes[1].alpha = alpha
        r.lanes[0].open_out = True
        r.lanes[0].beta = beta
        roads.append(r)

    in_legs  = [JunctionLegSpec(i, _lane_inbound_at(roads[i - 1], 1))  for i in range(1, 5)]
    out_legs = [JunctionLegSpec(i, _lane_outbound_at(roads[i - 1], 1)) for i in range(1, 5)]

    routes: List[List[float]] = []
    for k in range(1, 5):
        row = [0.0] * 4
        row[k % 4]       = p_left
        row[(k + 1) % 4] = p_straight
        row[(k + 2) % 4] = p_right
        routes.append(row)

    junc = JunctionSpec(
        id=1, n_in=4, n_out=4,
        in_legs=in_legs, out_legs=out_legs, routes=routes,
    )

    layout = LayoutSpec(
        junctions={1: (0.0, 0.0)},
        road_endpoints={
            1: ((0.0, 0.0), (0.0,  1.0)),   # N
            2: ((0.0, 0.0), ( 1.0, 0.0)),   # E
            3: ((0.0, 0.0), (0.0, -1.0)),   # S
            4: ((0.0, 0.0), (-1.0, 0.0)),   # W
        },
    )
    return NetworkSpec(roads=roads, junctions=[junc]), layout


def roundabout(L_arm: int = 15, L_ring: int = 12,
               alpha: float = 0.4, beta: float = 0.5,
               p_exit: float = 0.25) -> Tuple[NetworkSpec, LayoutSpec]:
    """4 T-junctions in a CW ring, with 4 external arms (N, E, S, W).

    Roads:
        1..4 : ring R1..R4, single lane fd=+1, CW.  Rk: J_k -> J_{k%4+1}.
        5..8 : arms N, E, S, W (2 lanes each, open at external end).

    Junctions 1..4 each have 2 in-legs (ring_in, arm_in) and 2 out-legs
    (ring_out, arm_out).  Perim ports CW: [ring_in=0, arm_in=1, ring_out=2, arm_out=3].
    """
    p_continue = 1.0 - p_exit

    L = 1.0
    junc_xy = {
        1: (0.0,  L),
        2: ( L, 0.0),
        3: (0.0, -L),
        4: (-L, 0.0),
    }
    ext_xy = {
        5: (0.0,  2 * L),
        6: ( 2 * L, 0.0),
        7: (0.0, -2 * L),
        8: (-2 * L, 0.0),
    }
    ring_pairs = {1: (1, 2), 2: (2, 3), 3: (3, 4), 4: (4, 1)}
    arm_junc   = {5: 1, 6: 2, 7: 3, 8: 4}

    roads: List[RoadSpec] = []
    for rid, (a, b) in ring_pairs.items():
        roads.append(RoadSpec(id=rid, end_junction=(a, b), lanes=[LaneSpec(L_ring, +1)]))

    for rid, jid in arm_junc.items():
        lanes = _bidir_lanes(L_arm)
        lanes[0].open_in = True
        lanes[0].alpha = alpha
        lanes[1].open_out = True
        lanes[1].beta = beta
        roads.append(RoadSpec(id=rid, end_junction=(0, jid), lanes=lanes))

    road_by_id = {r.id: r for r in roads}

    junctions: List[JunctionSpec] = []
    for jid in (1, 2, 3, 4):
        ring_in_rid  = (jid - 2) % 4 + 1   # ring road feeding into J_jid
        ring_out_rid = jid                 # ring road leaving J_jid
        arm_rid      = jid + 4

        ring_in_road  = road_by_id[ring_in_rid]
        ring_out_road = road_by_id[ring_out_rid]
        arm_road      = road_by_id[arm_rid]

        ring_in_end  = 1 if ring_in_road.end_junction[0] == jid else 2
        ring_out_end = 1 if ring_out_road.end_junction[0] == jid else 2
        arm_end      = 2   # arms always have the junction at end_2 by construction

        in_legs = [
            JunctionLegSpec(ring_in_rid, _lane_inbound_at(ring_in_road, ring_in_end),  perim=0),
            JunctionLegSpec(arm_rid,     _lane_inbound_at(arm_road,    arm_end),       perim=1),
        ]
        out_legs = [
            JunctionLegSpec(ring_out_rid, _lane_outbound_at(ring_out_road, ring_out_end), perim=2),
            JunctionLegSpec(arm_rid,      _lane_outbound_at(arm_road,      arm_end),      perim=3),
        ]
        routes = [
            [p_continue, p_exit],
            [1.0,        0.0],
        ]
        junctions.append(JunctionSpec(
            id=jid, n_in=2, n_out=2,
            in_legs=in_legs, out_legs=out_legs, routes=routes,
        ))

    road_endpoints = {}
    for rid, (a, b) in ring_pairs.items():
        road_endpoints[rid] = (junc_xy[a], junc_xy[b])
    for rid, jid in arm_junc.items():
        road_endpoints[rid] = (ext_xy[rid], junc_xy[jid])
    layout = LayoutSpec(junctions=junc_xy, road_endpoints=road_endpoints)

    return NetworkSpec(roads=roads, junctions=junctions), layout


def town(L: int = 15, alpha: float = 0.35, beta: float = 0.5,
         p_left: float = 0.25, p_right: float = 0.25) -> Tuple[NetworkSpec, LayoutSpec]:
    """2x2 grid of crossroads with 8 external arms.

    Layout (unit spacing):
        J1 (0,1)  -- R1 -- J2 (1,1)
          |                  |
         R2                 R3
          |                  |
        J4 (0,0)  -- R4 -- J5 (1,0)

    Internal roads (2 lanes bidir): R1=J1-J2, R2=J1-J4, R3=J2-J5, R4=J4-J5.
    Each crossroad has 4 arms; the two non-internal arms are external
    (1 inbound + 1 outbound lane, alpha/beta on the open end).  8 externals total.
    """
    p_straight = 1.0 - p_left - p_right

    junc_xy = {1: (0.0, 1.0), 2: (1.0, 1.0), 4: (0.0, 0.0), 5: (1.0, 0.0)}

    internal = [
        (1, 1, 2),   # R1: J1 -> J2
        (2, 1, 4),   # R2: J1 -> J4
        (3, 2, 5),   # R3: J2 -> J5
        (4, 4, 5),   # R4: J4 -> J5
    ]
    externals = [
        # (rid, junction_id, dx, dy)
        (5,  1,  0.0,  1.0),   # NW-N
        (6,  1, -1.0,  0.0),   # NW-W
        (7,  2,  0.0,  1.0),   # NE-N
        (8,  2,  1.0,  0.0),   # NE-E
        (9,  4, -1.0,  0.0),   # SW-W
        (10, 4,  0.0, -1.0),   # SW-S
        (11, 5,  1.0,  0.0),   # SE-E
        (12, 5,  0.0, -1.0),   # SE-S
    ]

    roads: List[RoadSpec] = []
    road_endpoints: Dict[int, Tuple[Tuple[float, float], Tuple[float, float]]] = {}

    for rid, j1, j2 in internal:
        roads.append(RoadSpec(id=rid, end_junction=(j1, j2), lanes=_bidir_lanes(L)))
        road_endpoints[rid] = (junc_xy[j1], junc_xy[j2])

    for rid, jid, dx, dy in externals:
        jx, jy = junc_xy[jid]
        lanes = _bidir_lanes(L)
        lanes[0].open_in = True
        lanes[0].alpha = alpha
        lanes[1].open_out = True
        lanes[1].beta = beta
        roads.append(RoadSpec(id=rid, end_junction=(0, jid), lanes=lanes))
        road_endpoints[rid] = ((jx + dx, jy + dy), junc_xy[jid])

    road_by_id = {r.id: r for r in roads}

    # CW arms at each crossroad: [N, E, S, W].  Each entry: (road_id, junc_end_on_that_road).
    junction_arms = {
        1: [(5, 2), (1, 1), (2, 1), (6, 2)],
        2: [(7, 2), (8, 2), (3, 1), (1, 2)],
        4: [(2, 2), (4, 1), (10, 2), (9, 2)],
        5: [(3, 2), (11, 2), (12, 2), (4, 2)],
    }

    junctions: List[JunctionSpec] = []
    for jid, arms in junction_arms.items():
        in_legs: List[JunctionLegSpec] = []
        out_legs: List[JunctionLegSpec] = []
        for rid, junc_end in arms:
            r = road_by_id[rid]
            in_legs.append(JunctionLegSpec(rid, _lane_inbound_at(r, junc_end)))
            out_legs.append(JunctionLegSpec(rid, _lane_outbound_at(r, junc_end)))
        routes = []
        for k in range(1, 5):
            row = [0.0] * 4
            row[k % 4]       = p_left
            row[(k + 1) % 4] = p_straight
            row[(k + 2) % 4] = p_right
            routes.append(row)
        junctions.append(JunctionSpec(
            id=jid, n_in=4, n_out=4,
            in_legs=in_legs, out_legs=out_legs, routes=routes,
        ))

    layout = LayoutSpec(junctions=junc_xy, road_endpoints=road_endpoints)
    return NetworkSpec(roads=roads, junctions=junctions), layout


PRESETS = {
    "single_lane": single_lane,
    "crossroads":  crossroads,
    "roundabout":  roundabout,
    "town":        town,
}
