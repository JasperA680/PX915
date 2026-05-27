PX915 Traffic Flow Modelling documentation
==========================================

This documentation describes the traffic-flow simulation code developed for
the PX915 group project.

The project contains two modelling approaches: a discrete cellular automaton
model (TASEP and Nagel–Schreckenberg) and a continuum PDE model (single-lane
and multilane LWR with conservative lane-changing source terms).  Both models
can be run from the command line or through a PyQt5 GUI, and both support
road-network topologies with junctions and probabilistic routing.

.. toctree::
   :maxdepth: 2
   :caption: User guide:

   user_guide/index

.. toctree::
   :maxdepth: 2
   :caption: Developer guide:

   developer_guide/index

.. toctree::
   :maxdepth: 2
   :caption: API reference:

   api/index