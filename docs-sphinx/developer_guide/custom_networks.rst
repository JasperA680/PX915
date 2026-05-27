Defining a custom road network
==============================

The CA / NS network simulator runs on any topology expressible through the
spec dataclasses in :mod:`road_network`. Seven preset builders cover the
common cases (single lane, two-lane, T-junction, crossroads, roundabout,
town, periodic ring), but anything outside that list — branches, bypasses,
asymmetric mergers, longer chains of junctions — has to be authored by
hand. This page documents the data model and the wiring rules so a
non-preset network can be assembled with confidence.

For a step-by-step worked example, see *Section 5: Custom Network —
Asymmetric Y-Junction* in the tutorial notebook
(``notebooks/tutorial.ipynb``).


The four-layer data model
-------------------------

A network is described by four nested dataclasses, all defined in
:mod:`road_network`:

* :class:`LaneSpec` — one directional chain of CA cells. Carries its
  length, flow direction (``+1`` or ``-1``), open-boundary inflow /
  outflow rates ``alpha`` / ``beta``, the ``open_in`` / ``open_out`` flags,
  and the optional ``is_periodic`` ring flag plus pre-placed
  ``n_vehicles`` count.
* :class:`RoadSpec` — a bundle of lanes between two endpoints. Each
  endpoint is either an *open boundary* (``end_junction = 0``) or
  *attached to a junction* by id. A bidirectional road is just a road
  with one lane of each flow direction.
* :class:`JunctionSpec` — a routing point. Carries lists of inbound and
  outbound :class:`JunctionLegSpec` legs (each leg references a specific
  road / lane index) and a routing-probability matrix.
* :class:`LayoutSpec` — purely cosmetic ``(x, y)`` coordinates for
  junctions and road endpoints. Ignored by the Fortran simulator but
  carried verbatim through the config NetCDF so the GUI and the
  ``plot_network_*`` helpers can re-render the geometry. An empty layout
  is allowed; the simulation will run, but visualisations will collapse
  every junction onto the origin.

Top-level :class:`NetworkSpec` holds ``roads`` and ``junctions``;
:class:`SimParams` carries the runtime parameters (``n_steps``, ``v_max``,
``p_slow``, ``rng_seed``, lane-change model, etc.). The full field-by-field
signatures are in the auto-generated API docs (see :doc:`../api/python`).


Wiring rules
------------

The presets in :mod:`road_network` encode several invariants silently. The
:func:`road_network.validate` function checks the easy ones at spec time;
the harder ones are runtime properties of the Fortran simulator and only
surface as crashes or nonsense flow. The rules below summarise both.

Lane flow direction
~~~~~~~~~~~~~~~~~~~

A lane's ``flow_direction`` determines which way vehicles travel along
the road:

* ``+1`` — vehicles travel ``end_1 → end_2``;
* ``-1`` — vehicles travel ``end_2 → end_1``.

A bidirectional road is built from two lanes with opposite
``flow_direction`` values. Several same-direction lanes on the same road
form a multilane carriageway and exercise the lane-change module when
``SimParams.lc_model >= 0``.

Road endpoints and junctions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``RoadSpec.end_junction = (a, b)`` records which junction (if any) sits at
each end of the road:

* ``0`` means *open boundary at this end*. Vehicles enter / exit the
  network at the open end at rates ``alpha`` / ``beta``; the lane must
  set ``open_in=True`` and / or ``open_out=True`` for those rates to be
  active.
* Any non-zero value is the ``id`` of the junction connected to that
  end. The same junction id may appear at one end, both ends, or neither.

A road with ``end_junction = (0, 0)`` and no junctions referencing it is
an isolated chain — useful for the single-lane and two-lane presets, and
for the closed periodic ring (``is_periodic=True``).

Junction legs and lane indices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each junction has a list of inbound legs (``in_legs``) and a list of
outbound legs (``out_legs``). A leg references a specific
``(road_id, lane_index)`` pair, where ``lane_index`` is the **1-based**
positional index into that road's ``lanes`` list (i.e. lane 1, lane 2,
…).

For a road whose end ``e`` is attached to a junction:

* an *inbound* leg of that junction references a lane whose
  ``flow_direction`` points **into** end ``e``;
* an *outbound* leg references a lane whose ``flow_direction`` points
  **out of** end ``e``.

The helpers ``_lane_inbound_at`` and ``_lane_outbound_at`` in
:mod:`road_network` are the canonical implementation of this rule and
the presets use them throughout — they're a useful reference when
writing a custom junction.

A multi-lane road feeding a junction contributes *one leg per lane* to
the appropriate side. A two-lane same-direction carriageway feeding into
a junction gives two inbound legs; the corresponding outflow road
contributes two outbound legs. ``n_in`` and ``n_out`` must equal the
lengths of ``in_legs`` and ``out_legs`` (``validate`` enforces this).

Routing matrix
~~~~~~~~~~~~~~

``JunctionSpec.routes`` is a list of length ``n_in``; each row has length
``n_out`` and gives the probability that a vehicle arriving on the
matching inbound leg leaves via each outbound leg. Each row must sum to
``1.0`` to within ``1e-6``; ``validate`` rejects rows that don't.

When the inbound and outbound legs of a junction reference the same
underlying road, the corresponding entry is a U-turn — by convention
these are set to ``0.0`` in every preset.

The matrix need not be square. A bifurcation (one inflow road, two
outflow roads) is ``n_in=1, n_out=2`` with a single-row routing matrix
``[[p_L, 1 - p_L]]`` — the Y-junction tutorial example. A merge is
``n_in=2, n_out=1``. Multi-lane same-direction roads contribute one leg
per lane, so the dimensions scale up accordingly.

Perimeter ports
~~~~~~~~~~~~~~~

Each :class:`JunctionLegSpec` carries an optional ``perim`` field. These
are clockwise port indices around the junction perimeter — they are used
by the Fortran asymmetric-yield code (``junction_mod::chords_cross``) to
decide which route chords cross when two candidates contend for the same
junction step.

The validation rules are:

* perim must be **all-or-none across a single junction**: if any leg
  sets ``perim >= 0``, *every* leg of that junction must.
* perim values within a junction must be **distinct integers** (one per
  port).

Symmetric four-way junctions (``n_in == 4 == n_out``) don't use perim:
the symmetric yield rules in ``build_yield_matrix_v2`` work from leg
indices alone. The :func:`road_network.crossroads` and
:func:`road_network.town` presets exploit this and leave ``perim`` at
its default of ``-1`` throughout.

**Any other junction shape (anything other than a symmetric 4-way) must
set perim on every leg.** The asymmetric yield routine
(``build_yield_matrix_asym``) dereferences ``in_perim(i)`` and
``out_perim(i)`` unconditionally; if those arrays were not allocated
(because no leg set perim) the run will segfault. Set perim ports
clockwise around the junction, starting from any convenient leg. The
exact starting position doesn't matter — only the cyclic order is
consulted.

The :func:`road_network.t_junction` and :func:`road_network.roundabout`
presets show the convention: each road takes two consecutive ports (in
then out) and roads are placed CW around the junction.

Periodic and pre-placed lanes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two optional :class:`LaneSpec` fields are consumed by the Fortran
builder rather than the runtime step:

* ``is_periodic = True`` — closes the lane into a ring of length
  ``length``. The periodic branch of the NS step (``NS_model_step`` in
  ``NS_model.f90``) ignores ``alpha`` / ``beta`` / ``open_in`` /
  ``open_out`` entirely. Mutually exclusive with open-boundary behaviour
  in practice. Used by the :func:`road_network.single_lane_periodic`
  preset and by the periodic-NS fundamental-diagram sweep.
* ``n_vehicles > 0`` — when set, ``network_init_mod::place_evenly``
  seeds that many cars at evenly-spaced sites on this lane before step
  1. Used to fix the global density on periodic rings; harmless on open
  lanes (just gives a non-empty initial state).


Spec field reference
--------------------

The tables below summarise each dataclass's fields. The full
auto-generated reference (including default values and docstrings) lives
in :doc:`../api/python`; this page is the quick lookup for hand-building.

LaneSpec
~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 22 16 62

   * - Field
     - Type
     - Meaning
   * - ``length``
     - ``int``
     - Number of CA cells in the lane.
   * - ``flow_direction``
     - ``int``
     - ``+1`` for end_1 → end_2, ``-1`` for end_2 → end_1.
   * - ``alpha``
     - ``float``
     - Open-boundary inflow rate. Ignored unless ``open_in`` is true.
   * - ``beta``
     - ``float``
     - Open-boundary outflow rate. Ignored unless ``open_out`` is true.
   * - ``open_in``
     - ``bool``
     - Open inflow at the upstream end of this lane.
   * - ``open_out``
     - ``bool``
     - Open outflow at the downstream end of this lane.
   * - ``is_periodic``
     - ``bool``
     - Treat the lane as a closed ring.
   * - ``n_vehicles``
     - ``int``
     - Number of pre-placed vehicles at step 0.

RoadSpec
~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 22 22 56

   * - Field
     - Type
     - Meaning
   * - ``id``
     - ``int``
     - Unique road identifier; referenced by every junction leg.
   * - ``end_junction``
     - ``tuple[int, int]``
     - ``(j_at_end_1, j_at_end_2)``; ``0`` = open boundary.
   * - ``lanes``
     - ``list[LaneSpec]``
     - One entry per lane on this road; at least one required.

JunctionLegSpec
~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 22 16 62

   * - Field
     - Type
     - Meaning
   * - ``road``
     - ``int``
     - Road id this leg attaches to.
   * - ``lane``
     - ``int``
     - 1-based lane index within that road.
   * - ``perim``
     - ``int``
     - Optional CW perimeter port index; ``-1`` = symmetric junction.

JunctionSpec
~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 22 28 50

   * - Field
     - Type
     - Meaning
   * - ``id``
     - ``int``
     - Unique junction identifier.
   * - ``n_in`` / ``n_out``
     - ``int``
     - Number of inbound / outbound legs.
   * - ``in_legs`` / ``out_legs``
     - ``list[JunctionLegSpec]``
     - Inbound and outbound leg specs.
   * - ``routes``
     - ``list[list[float]]``
     - ``n_in × n_out`` routing-probability matrix; each row sums to 1.

LayoutSpec
~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 28 36 36

   * - Field
     - Type
     - Meaning
   * - ``junctions``
     - ``dict[int, tuple[float, float]]``
     - ``{junction_id: (x, y)}``.
   * - ``road_endpoints``
     - ``dict[int, tuple[tuple[float, float], tuple[float, float]]]``
     - ``{road_id: ((x1, y1), (x2, y2))}``.


Validating, running, and inspecting a custom network
----------------------------------------------------

The flow is the same as for any preset. With a custom ``spec`` /
``layout`` in hand:

#. Call :func:`road_network.validate` on the spec. This checks duplicate
   road / junction ids, leg counts, route-row sums, and the perim
   all-or-none rule. ``validate`` raises ``ValueError`` for any
   violation — fix the spec rather than catching.
#. Build :class:`SimParams` with the desired ``n_steps``, ``v_max``,
   ``p_slow``, ``rng_seed`` and lane-change model (``lc_model = 0`` for
   symmetric lane changes on multi-lane carriageways; ``-1`` to
   disable).
#. Call :func:`python.run_simulation.run_simulation` with the spec,
   params, layout and an output directory. The runner writes
   ``config.nc`` (the spec + layout serialised by
   :func:`python.io.write_config_netcdf`), invokes
   ``./build/run_network`` as a subprocess, and returns the path to
   ``result.nc``.
#. Load the result with :func:`python.io.load_network_netcdf`. The
   returned :class:`NetworkResult` carries the full
   ``occupancy[t, lane_index, cell]`` history, ``velocity`` history,
   per-road density / entries / exits time series, and the lane-level
   metadata needed by the visualisation helpers.

The standard plotting helpers in :mod:`CA_visualisation` work on any
:class:`NetworkResult`:

* :func:`CA_visualisation.plot_network_layout` — the network diagram
  with optional per-lane occupancy heatmap (``occupancy_t=t`` for an
  instantaneous snapshot, or ``occupancy_t=(t0, t1)`` for a mean over a
  range);
* :func:`CA_visualisation.plot_network_density` — per-road density time
  series;
* :func:`CA_visualisation.plot_network_currents` — per-road exit-flow
  time series;
* :func:`CA_visualisation.plot_network_spacetime` — single-lane
  space-time diagram, keyed by ``road_id`` and ``lane_index``.

The heatmap colours and the space-time diagrams require the
:class:`LayoutSpec` to be populated; the density and currents plots do
not (they're indexed by road id alone).


Reusing the presets as templates
--------------------------------

The seven preset builders in :mod:`road_network` are the most
up-to-date worked examples in the codebase:

* :func:`road_network.single_lane` and
  :func:`road_network.single_lane_periodic` — one road, one lane, open
  or closed.
* :func:`road_network.two_lane` — one road, two same-direction lanes;
  the lane-change baseline.
* :func:`road_network.t_junction` — three roads at one junction, with
  perimeter ports set (good template for any asymmetric junction).
* :func:`road_network.crossroads` — four roads at one symmetric
  junction, no perim required.
* :func:`road_network.roundabout` — four T-junctions in a ring with
  four external arms.
* :func:`road_network.town` — 2×2 grid of crossroads, eight external
  arms.

For a topology not far off one of these, copying the preset source and
modifying it is usually faster than building from scratch. The
``crossroads`` source is the cleanest template for any symmetric multi-
arm junction; ``t_junction`` is the cleanest template when the junction
needs perim ports.
