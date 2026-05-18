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

# Source files (existing single-lattice TASEP)
MODEL_SRC = $(SRC_DIR)/tasep.f90
SIM_SRC   = $(SRC_DIR)/simulation.f90
IO_SRC    = $(SRC_DIR)/io.f90
TEST_SRC  = $(SRC_DIR)/test_simulation.f90
FD_SRC    = $(SRC_DIR)/fundamental_diagram.f90

# Network-model sources
VEHICLE_SRC      = $(SRC_DIR)/vehicle.f90
NETWORK_SRC      = $(SRC_DIR)/road_network.f90
NETWORK_INIT_SRC = $(SRC_DIR)/network_init.f90
NS_MODEL_SRC     = $(SRC_DIR)/NS_model.f90
JUNCTION_SRC     = $(SRC_DIR)/junction.f90
LANE_CHANGE_SRC  = $(SRC_DIR)/lane_change.f90
NET_SIM_SRC      = $(SRC_DIR)/network_simulation.f90

# Network builder + JSON config + NetCDF writer + driver (Python frontend integration)
NETWORK_BUILDER_SRC = $(SRC_DIR)/network_builder.f90
JSON_CONFIG_SRC     = $(SRC_DIR)/json_config.f90
NETWORK_IO_SRC      = $(SRC_DIR)/network_io.f90
RUN_NETWORK_SRC     = $(SRC_DIR)/run_network.f90

NETWORK_LIB_SRC = $(VEHICLE_SRC) $(NETWORK_SRC) $(NETWORK_INIT_SRC) \
                  $(LANE_CHANGE_SRC) $(JUNCTION_SRC) $(NET_SIM_SRC)

# PDE source files (pde_flux and pde_lanechange must precede pde_module due to USE deps)
PDE_FLUX_SRC = $(SRC_DIR)/pde_flux.f90
PDE_LC_SRC   = $(SRC_DIR)/pde_lanechange.f90
PDE_MOD_SRC  = $(SRC_DIR)/pde_module.f90
PDE_DRV_SRC  = $(SRC_DIR)/pde_driver.f90

# Executables
TEST_EXE              = $(BUILD_DIR)/test_simulation
FUND_EXE              = $(BUILD_DIR)/fundamental_diagram
PDE_EXE               = $(BUILD_DIR)/pde_solver
TEST_TYPES_EXE        = $(BUILD_DIR)/test_types
TEST_INIT_EXE         = $(BUILD_DIR)/test_init_crossroad
TEST_JUNCTION_EXE     = $(BUILD_DIR)/test_junction
TEST_NETWORK_RUN_EXE  = $(BUILD_DIR)/test_network_run
TEST_NS_GRID_EXE      = $(BUILD_DIR)/test_ns_grid
TEST_NS_PERIODIC_EXE  = $(BUILD_DIR)/test_ns_periodic
RUN_NETWORK_EXE       = $(BUILD_DIR)/run_network
TEST_LANE_CHANGE_EXE  = $(BUILD_DIR)/test_lane_change

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
PDE_OUT    ?= data/output/pde_simulation.nc

# Default target: build all executables
all: test test_types test_init_crossroad test_junction test_network_run test_ns_grid test_ns_periodic test_lane_change $(FUND_EXE) $(PDE_EXE)

# Build existing test program (needs NetCDF)
test: $(TEST_EXE)

$(TEST_EXE): $(NETWORK_LIB_SRC) $(MODEL_SRC) $(SIM_SRC) $(IO_SRC) $(TEST_SRC)
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) $^ -o $@ $(NC_FLIBS)

# Build fundamental diagram sweep
$(FUND_EXE): $(MODEL_SRC) $(SIM_SRC) $(IO_SRC) $(FD_SRC)
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) $^ -o $@ $(NC_FLIBS)

# Build PDE solver
$(PDE_EXE): $(PDE_FLUX_SRC) $(PDE_LC_SRC) $(PDE_MOD_SRC) $(PDE_DRV_SRC)
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) $^ -o $@ $(NC_FLIBS)

pde: $(PDE_EXE)

# Phase 1 smoke test
test_types: $(TEST_TYPES_EXE)

$(TEST_TYPES_EXE): $(VEHICLE_SRC) $(NETWORK_SRC) tests/test_types.f90
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) $^ -o $@

# Phase 2 init test
test_init_crossroad: $(TEST_INIT_EXE)

$(TEST_INIT_EXE): $(NETWORK_LIB_SRC) tests/test_init_crossroad.f90
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) $^ -o $@

# Run targets
run: test
	./$(TEST_EXE) $(L) $(N_STEPS) $(ALPHA) $(BETA)

run-fd: $(FUND_EXE)
	./$(FUND_EXE)

run-pde: $(PDE_EXE)
	./$(PDE_EXE) $(PDE_M) $(PDE_STEPS) $(PDE_VMAX) $(PDE_RHOMAX) $(PDE_LEFT) $(PDE_RIGHT) $(PDE_IC) $(PDE_FLUX) $(PDE_BC) $(PDE_OUT)

# Run the Phase 1 smoke test
run_test_types: test_types
	./$(TEST_TYPES_EXE)

# Phase 3 junction test
test_junction: $(TEST_JUNCTION_EXE)

$(TEST_JUNCTION_EXE): $(NETWORK_LIB_SRC) tests/test_junction.f90
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) $^ -o $@

# Run the Phase 2 init test
run_test_init_crossroad: test_init_crossroad
	./$(TEST_INIT_EXE)

# Run the Phase 3 junction tests
run_test_junction: test_junction
	./$(TEST_JUNCTION_EXE)

# Integration soak test
test_network_run: $(TEST_NETWORK_RUN_EXE)

$(TEST_NETWORK_RUN_EXE): $(NETWORK_LIB_SRC) tests/test_network_run.f90
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) $^ -o $@

run_test_network_run: test_network_run
	./$(TEST_NETWORK_RUN_EXE)

# NS 2D grid test
test_ns_grid: $(TEST_NS_GRID_EXE)

$(TEST_NS_GRID_EXE): $(NETWORK_LIB_SRC) tests/test_ns_grid.f90
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) $^ -o $@

run_test_ns_grid: test_ns_grid
	./$(TEST_NS_GRID_EXE)

# NS periodic boundary test
test_ns_periodic: $(TEST_NS_PERIODIC_EXE)

$(TEST_NS_PERIODIC_EXE): $(VEHICLE_SRC) $(NETWORK_SRC) tests/test_ns_periodic.f90
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) $^ -o $@

run_test_ns_periodic: test_ns_periodic
	./$(TEST_NS_PERIODIC_EXE)

# Python-frontend driver: reads JSON config, runs sim, writes NetCDF.
run_network: $(RUN_NETWORK_EXE)

$(RUN_NETWORK_EXE): $(NETWORK_LIB_SRC) $(NETWORK_BUILDER_SRC) $(JSON_CONFIG_SRC) $(NETWORK_IO_SRC) $(RUN_NETWORK_SRC)
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) $^ -o $@ $(NC_FLIBS)

# Lane-change unit/integration tests
test_lane_change: $(TEST_LANE_CHANGE_EXE)

$(TEST_LANE_CHANGE_EXE): $(NETWORK_LIB_SRC) tests/test_lane_change.f90
	mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) -J$(BUILD_DIR) $^ -o $@

run_test_lane_change: test_lane_change
	./$(TEST_LANE_CHANGE_EXE)

# Clean build files
clean:
	rm -rf $(BUILD_DIR)
	rm -f $(SRC_DIR)/*.mod *.mod
