# Compiler
FC = gfortran

# NetCDF flags (Fortran wrapper + underlying C library).
# nf-config --flibs and nc-config --libs both append -lnetcdf, which double-links
# on this Mac.  Strip the trailing -lnetcdf from nf-config and let nc-config provide it.
NC_FFLAGS := $(shell nf-config --fflags)
NC_FLIBS   := $(shell nf-config --flibs | sed -E 's/-lnetcdf$$//') $(shell nc-config --libs)

# Compiler flags
FFLAGS = -Wall -O2 $(NC_FFLAGS)

# Directories
SRC_DIR   = src/fortran
BUILD_DIR = build


# Network-model sources
VEHICLE_SRC      = $(SRC_DIR)/vehicle.f90
NETWORK_SRC      = $(SRC_DIR)/road_network.f90
NETWORK_INIT_SRC = $(SRC_DIR)/network_init.f90
TASEP_SRC        = $(SRC_DIR)/tasep.f90
NS_MODEL_SRC     = $(SRC_DIR)/NS_model.f90
JUNCTION_SRC     = $(SRC_DIR)/junction.f90
LANE_CHANGE_SRC  = $(SRC_DIR)/lane_change.f90
NET_SIM_SRC      = $(SRC_DIR)/network_simulation.f90

# Network builder + NetCDF config reader + NetCDF result writer + driver
# (Python frontend writes config.nc via python.io.write_config_netcdf;
# Fortran reads it via nc_config_mod.read_config.)
NETWORK_BUILDER_SRC = $(SRC_DIR)/network_builder.f90
NC_CONFIG_SRC       = $(SRC_DIR)/nc_config.f90
NETWORK_IO_SRC      = $(SRC_DIR)/network_io.f90
RUN_NETWORK_SRC     = $(SRC_DIR)/run_network.f90

NETWORK_LIB_SRC = $(VEHICLE_SRC) $(NETWORK_SRC) $(NETWORK_INIT_SRC) \
                  $(LANE_CHANGE_SRC) $(JUNCTION_SRC) $(TASEP_SRC) $(NS_MODEL_SRC) \
                  $(NET_SIM_SRC)

# PDE source files (pde_flux and pde_lanechange must precede pde_module due to USE deps)
PDE_FLUX_SRC = $(SRC_DIR)/pde_flux.f90
PDE_LC_SRC   = $(SRC_DIR)/pde_lanechange.f90
PDE_MOD_SRC  = $(SRC_DIR)/pde_module.f90
PDE_DRV_SRC  = $(SRC_DIR)/pde_driver.f90

# Fundamental-diagram sweep.
# Single source file containing both ``fundamental_diagram_mod`` (measurement
# + NetCDF I/O, reusing the network primitives on a minimal 1-lane
# road_network_t) and the ``fd_sweep`` driver program. Shares the NetCDF
# ``check`` helper with network_io_mod.
FD_SRC          = $(SRC_DIR)/fundamental_diagram.f90

# Executables
PDE_EXE               = $(BUILD_DIR)/pde_solver
RUN_NETWORK_EXE       = $(BUILD_DIR)/run_network
FD_SWEEP_EXE          = $(BUILD_DIR)/fd_sweep

# Default simulation parameters (override: make run L=50 ALPHA=0.3 BETA=0.7)
L       ?= 10
N_STEPS ?= 100
ALPHA   ?= 0.5
BETA    ?= 0.5

# Default PDE parameters (override: make run-pde PDE_M=400 PDE_IC=riemann)
PDE_M      ?= 200
PDE_STEPS  ?= 500
PDE_VMAX   ?= 1.0
PDE_RHOMAX ?= 1.0
PDE_LEFT   ?= 0.1
PDE_RIGHT  ?= 0.9
PDE_IC     ?= riemann
PDE_FLUX   ?= lf
PDE_BC     ?= open
PDE_VLIMIT ?= 1.0
PDE_OUT    ?= data/output/pde_simulation.nc

# Default target: full setup + build all executables
all: setup run_network pde fd_sweep

# Python environment setup: install requirements and register a Jupyter
# kernel pointing at the active interpreter so `tutorial.ipynb` can find it.
# Run with the target venv activated (or PYTHON= pointing at it). pip is
# idempotent so re-running this on an up-to-date venv is a no-op.
PYTHON ?= python3

setup:
	$(PYTHON) -m pip install -r requirements.txt
	$(PYTHON) -m ipykernel install --user --name=px915 --display-name "Python (PX915)"

# Build PDE solver
$(PDE_EXE): $(PDE_FLUX_SRC) $(PDE_LC_SRC) $(PDE_MOD_SRC) $(PDE_DRV_SRC)
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) $^ -o $@ $(NC_FLIBS)

pde: $(PDE_EXE)

# Run targets

run-pde: $(PDE_EXE)
	./$(PDE_EXE) $(PDE_M) $(PDE_STEPS) $(PDE_VMAX) $(PDE_RHOMAX) $(PDE_LEFT) $(PDE_RIGHT) $(PDE_IC) $(PDE_FLUX) $(PDE_BC) $(PDE_VLIMIT) $(PDE_OUT)

# Python-frontend driver: reads NetCDF config, runs sim, writes NetCDF.
run_network: $(RUN_NETWORK_EXE)

$(RUN_NETWORK_EXE): $(NETWORK_LIB_SRC) $(NETWORK_BUILDER_SRC) $(NC_CONFIG_SRC) $(NETWORK_IO_SRC) $(RUN_NETWORK_SRC)
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) $^ -o $@ $(NC_FLIBS)

# Fundamental-diagram sweep driver (open-boundary 1D TASEP + periodic-ring NS).
# Drives the existing network TASEP/NS step routines on a 1-lane network, so
# it pulls in road_network, network_init, tasep, NS_model and network_io
# (for the shared NetCDF ``check`` helper). Junctions / lane-change /
# network_simulation are not needed.
#
# Built with -fopenmp: the two TASEP sweep loops (alpha and beta branches)
# and the NS density sweep are wrapped in !$omp parallel do. Each iteration
# is a self-contained steady-state measurement on its own road_network_t,
# and gfortran's random_number state is per-thread, so the only thing the
# parallel region needs is per-iteration RNG seeding (which the program
# already does via fundamental_diagram_mod::seed_iter_rng for bit-identical
# results across thread counts).
fd_sweep: $(FD_SWEEP_EXE)

$(FD_SWEEP_EXE): $(VEHICLE_SRC) $(NETWORK_SRC) $(NETWORK_INIT_SRC) \
                 $(TASEP_SRC) $(NS_MODEL_SRC) $(NETWORK_IO_SRC) \
                 $(FD_SRC)
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) -fopenmp -J$(BUILD_DIR) $^ -o $@ $(NC_FLIBS)

# Clean build files
clean:
	rm -rf $(BUILD_DIR)
	rm -f $(SRC_DIR)/*.mod *.mod
